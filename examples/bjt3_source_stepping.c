#include "dc_solver.h"
#include "source_stepping.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

/*
 * 电路三（Chua 四晶体管基准电路）的完整 MNA 模型。
 * 未知量 x[0] 到 x[13] 分别对应网表节点 x1 到 x14；
 * x[14] 到 x[17] 分别是 VCC1、VCC2、V1、V2 的支路电流。
 */
enum {
    BJT3_NODE_COUNT = 14,
    BJT3_VCC1_CURRENT = 14,
    BJT3_VCC2_CURRENT = 15,
    BJT3_V1_CURRENT = 16,
    BJT3_V2_CURRENT = 17,
    BJT3_DIMENSION = 18,
    BJT3_GROUND = -1
};

typedef struct {
    int positive_node;
    int negative_node;
    double resistance;
} Bjt3Resistor;

typedef struct {
    int positive_node;
    int negative_node;
    int branch_index;
    double voltage;
} Bjt3VoltageSource;

/* SPICE 晶体管端口顺序：collector, base, emitter。 */
typedef struct {
    int collector;
    int base;
    int emitter;
} Bjt3Transistor;

typedef struct {
    Bjt3Resistor resistors[14];
    Bjt3VoltageSource sources[4];
    Bjt3Transistor transistors[4];
    double is;
    double alpha_f;
    double alpha_r;
    double vt;
    double temperature_kelvin;
    double nominal_temperature_kelvin;
    double xti;
    double bandgap_ev;
} Bjt3Circuit;

typedef struct {
    double ic[4];
    double ib[4];
    double ie[4];
} Bjt3Values;

static int bjt3_node_index(int node)
{
    return node == BJT3_GROUND ? -1 : node - 1;
}

static double bjt3_node_voltage(const double *x, int node)
{
    const int index = bjt3_node_index(node);
    return index < 0 ? 0.0 : x[index];
}

/* 限制 exp() 的输入，防止不合适的 Newton 猜测导致溢出。 */
static double bjt3_safe_exp(double argument)
{
    if (argument > 40.0) {
        argument = 40.0;
    } else if (argument < -40.0) {
        argument = -40.0;
    }
    return exp(argument);
}

/* 与网表的 26 摄氏度 .op 条件一致的饱和电流温度修正。 */
static double bjt3_effective_is(const Bjt3Circuit *circuit)
{
    const double t = circuit->temperature_kelvin;
    const double tnom = circuit->nominal_temperature_kelvin;
    const double boltzmann_over_q = 8.617333262e-5;
    const double temperature_ratio = t / tnom;
    const double exponent = -circuit->bandgap_ev / boltzmann_over_q
                          * (1.0 / t - 1.0 / tnom);

    return circuit->is * pow(temperature_ratio, circuit->xti) * exp(exponent);
}

static void bjt3_add_residual(double *residual, int node, double value)
{
    const int index = bjt3_node_index(node);
    if (index >= 0) {
        residual[index] += value;
    }
}

static void bjt3_add_jacobian(
    double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS],
    int row_node,
    int column_node,
    double value)
{
    const int row = bjt3_node_index(row_node);
    const int column = bjt3_node_index(column_node);

    if (row >= 0 && column >= 0) {
        jacobian[row][column] += value;
    }
}

/* 将一个线性电阻的贡献写入残差和 Jacobian。 */
static void bjt3_stamp_resistor(
    const Bjt3Resistor *resistor,
    const double *x,
    double *residual,
    double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{
    const double conductance = 1.0 / resistor->resistance;
    const double voltage = bjt3_node_voltage(x, resistor->positive_node)
                         - bjt3_node_voltage(x, resistor->negative_node);

    bjt3_add_residual(residual, resistor->positive_node, conductance * voltage);
    bjt3_add_residual(residual, resistor->negative_node, -conductance * voltage);

    bjt3_add_jacobian(jacobian, resistor->positive_node,
                       resistor->positive_node, conductance);
    bjt3_add_jacobian(jacobian, resistor->positive_node,
                       resistor->negative_node, -conductance);
    bjt3_add_jacobian(jacobian, resistor->negative_node,
                       resistor->positive_node, -conductance);
    bjt3_add_jacobian(jacobian, resistor->negative_node,
                       resistor->negative_node, conductance);
}

/* 将 MNA 中独立电压源的贡献写入方程；source stepping 只缩放其约束右端。 */
static void bjt3_stamp_voltage_source(
    const Bjt3VoltageSource *source,
    const double *x,
    double lambda,
    double *residual,
    double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{
    const int positive = bjt3_node_index(source->positive_node);
    const int negative = bjt3_node_index(source->negative_node);
    const int branch = source->branch_index;
    const double current = x[branch];
    const double voltage = bjt3_node_voltage(x, source->positive_node)
                         - bjt3_node_voltage(x, source->negative_node);

    if (positive >= 0) {
        residual[positive] += current;
        jacobian[positive][branch] += 1.0;
        jacobian[branch][positive] += 1.0;
    }
    if (negative >= 0) {
        residual[negative] -= current;
        jacobian[negative][branch] -= 1.0;
        jacobian[branch][negative] -= 1.0;
    }
    residual[branch] += voltage - lambda * source->voltage;
}

/* 将一个 Ebers-Moll NPN 晶体管的电流和解析导数写入方程。 */
static void bjt3_stamp_transistor(
    const Bjt3Circuit *circuit,
    const Bjt3Transistor *transistor,
    const double *x,
    double *residual,
    double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS],
    Bjt3Values *values,
    int transistor_index)
{
    const int c = transistor->collector;
    const int b = transistor->base;
    const int e = transistor->emitter;
    const double vc = bjt3_node_voltage(x, c);
    const double vb = bjt3_node_voltage(x, b);
    const double ve = bjt3_node_voltage(x, e);
    const double ef = bjt3_safe_exp((vb - ve) / circuit->vt);
    const double er = bjt3_safe_exp((vb - vc) / circuit->vt);
    const double is = bjt3_effective_is(circuit);
    const double if_current = is * (ef - 1.0);
    const double ir_current = is * (er - 1.0);
    const double gf = is * ef / circuit->vt;
    const double gr = is * er / circuit->vt;
    /*
     * 这里的 IS、BF、BR 采用 Schmitt3_netlist 中 PSpice .model 的定义。
     * IS 是输运饱和电流，alpha_f=BF/(BF+1)，alpha_r=BR/(BR+1)。
     */
    const double ic = if_current - ir_current / circuit->alpha_r;
    const double ib = (1.0 / circuit->alpha_f - 1.0) * if_current
                    + (1.0 / circuit->alpha_r - 1.0) * ir_current;
    const double ie = if_current / circuit->alpha_f - ir_current;

    values->ic[transistor_index] = ic;
    values->ib[transistor_index] = ib;
    values->ie[transistor_index] = ie;

    /* KCL：集电极、基极电流流入晶体管；发射极电流从晶体管流出。 */
    bjt3_add_residual(residual, c, ic);
    bjt3_add_residual(residual, b, ib);
    bjt3_add_residual(residual, e, -ie);

    /* 集电极电流 IC = IF - IR / alpha_r 的导数。 */
    bjt3_add_jacobian(jacobian, c, c, gr / circuit->alpha_r);
    bjt3_add_jacobian(jacobian, c, b, gf - gr / circuit->alpha_r);
    bjt3_add_jacobian(jacobian, c, e, -gf);

    /* 基极电流 IB 的导数。 */
    bjt3_add_jacobian(jacobian, b, c,
                       -(1.0 / circuit->alpha_r - 1.0) * gr);
    bjt3_add_jacobian(jacobian, b, b,
                       (1.0 / circuit->alpha_f - 1.0) * gf
                     + (1.0 / circuit->alpha_r - 1.0) * gr);
    bjt3_add_jacobian(jacobian, b, e,
                       -(1.0 / circuit->alpha_f - 1.0) * gf);

    /* 发射极残差中为 -IE，IE = IF / alpha_f - IR。 */
    bjt3_add_jacobian(jacobian, e, c, -gr);
    bjt3_add_jacobian(jacobian, e, b, -gf / circuit->alpha_f + gr);
    bjt3_add_jacobian(jacobian, e, e, gf / circuit->alpha_f);
}

static void bjt3_build_system(
    const Bjt3Circuit *circuit,
    const double *x,
    double lambda,
    double *residual,
    double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS],
    Bjt3Values *values)
{
    memset(residual, 0, sizeof(double) * BJT3_DIMENSION);
    memset(jacobian, 0, sizeof(double) * DC_MAX_UNKNOWNS * DC_MAX_UNKNOWNS);
    memset(values, 0, sizeof(*values));

    for (int i = 0; i < 14; ++i) {
        bjt3_stamp_resistor(&circuit->resistors[i], x, residual, jacobian);
    }
    for (int i = 0; i < 4; ++i) {
        bjt3_stamp_voltage_source(&circuit->sources[i], x, lambda,
                                  residual, jacobian);
    }
    for (int i = 0; i < 4; ++i) {
        bjt3_stamp_transistor(circuit, &circuit->transistors[i], x,
                              residual, jacobian, values, i);
    }
}

static void bjt3_build_residual(
    const void *context,
    const double *x,
    double lambda,
    double *residual)
{
    double unused_jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS];
    Bjt3Values unused_values;

    bjt3_build_system(context, x, lambda, residual, unused_jacobian,
                      &unused_values);
}

static void bjt3_build_jacobian(
    const void *context,
    const double *x,
    double lambda,
    double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{
    double unused_residual[BJT3_DIMENSION];
    Bjt3Values unused_values;

    bjt3_build_system(context, x, lambda, unused_residual, jacobian,
                      &unused_values);
}

static Bjt3Circuit bjt3_default_circuit(void)
{
    return (Bjt3Circuit) {
        .resistors = {
            {2, 3, 10000.0}, {4, 13, 4000.0}, {5, 14, 4000.0},
            {1, BJT3_GROUND, 5000.0}, {4, 6, 30000.0},
            {8, BJT3_GROUND, 500.0}, {9, BJT3_GROUND, 500.0},
            {5, 7, 30000.0}, {6, BJT3_GROUND, 10100.0},
            {7, BJT3_GROUND, 10100.0}, {13, 10, 4000.0},
            {14, 12, 4000.0}, {10, 3, 30000.0}, {12, 11, 30000.0}
        },
        .sources = {
            {13, BJT3_GROUND, BJT3_VCC1_CURRENT, 12.0},
            {14, BJT3_GROUND, BJT3_VCC2_CURRENT, 12.0},
            {1, 2, BJT3_V1_CURRENT, 10.0},
            {11, 3, BJT3_V2_CURRENT, 2.0}
        },
        .transistors = {
            {4, 1, 8}, {10, 6, 8}, {5, 1, 9}, {12, 7, 9}
        },
        .is = 1.0e-9,
        .alpha_f = 100.0 / 101.0,
        .alpha_r = 0.5,
        /* 网表说明要求在 26 摄氏度运行，PSpice 的 TNOM 默认为 27 摄氏度。 */
        .vt = 8.617333262e-5 * (273.15 + 26.0),
        .temperature_kelvin = 273.15 + 26.0,
        .nominal_temperature_kelvin = 273.15 + 27.0,
        .xti = 3.0,
        .bandgap_ev = 1.11
    };
}

static void bjt3_make_problem(const Bjt3Circuit *circuit, DcProblem *problem)
{
    *problem = (DcProblem) {
        .dimension = BJT3_DIMENSION,
        .context = circuit,
        .build_residual = bjt3_build_residual,
        .build_jacobian = bjt3_build_jacobian
    };
}

static void bjt3_evaluate(
    const Bjt3Circuit *circuit,
    const double *x,
    Bjt3Values *values)
{
    double residual[BJT3_DIMENSION];
    double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS];

    bjt3_build_system(circuit, x, 1.0, residual, jacobian, values);
}

#ifndef BJT3_SOURCE_STEPPING_NO_MAIN
static void bjt3_print_step(
    double lambda,
    const double *x,
    int iterations,
    void *user_data)
{
    (void)user_data;
    printf(" %6.4f  %9.5f  %9.5f  %9.5f  %9.5f  %7d\n",
           lambda, x[0], x[3], x[5], x[7], iterations);
}

int main(void)
{
    const Bjt3Circuit circuit = bjt3_default_circuit();
    DcSolverOptions options = dc_solver_default_options();
    DcProblem problem;
    double x[DC_MAX_UNKNOWNS];
    Bjt3Values values;
    int total_iterations = 0;

    /* 多工作点电路使用较小初始步长，给延拓路径更多调整余量。 */
    options.initial_lambda_step = 0.02;
    options.maximum_lambda_step = 0.10;
    options.maximum_newton_iterations = 120;

    bjt3_make_problem(&circuit, &problem);
    puts(" lambda      X1(V)      X4(V)      X6(V)      X8(V)   Newton");

    if (!dc_source_stepping_solve(&problem, &options, x, bjt3_print_step,
                                  NULL, &total_iterations)) {
        fputs("Source stepping could not reach lambda = 1.\n", stderr);
        return 1;
    }

    bjt3_evaluate(&circuit, x, &values);
    printf("\nTotal Newton iterations across all successful lambda points: %d\n",
           total_iterations);
    puts("\nFinal DC operating point (lambda = 1):");
    for (int i = 0; i < BJT3_NODE_COUNT; ++i) {
        printf("X%-2d = % .9f V\n", i + 1, x[i]);
    }
    puts("\nVoltage-source branch currents:");
    printf("IVCC1 = % .9e A\nIVCC2 = % .9e A\n", x[BJT3_VCC1_CURRENT],
           x[BJT3_VCC2_CURRENT]);
    printf("IV1   = % .9e A\nIV2   = % .9e A\n", x[BJT3_V1_CURRENT],
           x[BJT3_V2_CURRENT]);
    puts("\nBJT terminal currents (positive into collector/base):");
    for (int i = 0; i < 4; ++i) {
        printf("T%d: IC = % .9e A, IB = % .9e A, IE = % .9e A\n", i + 1,
               values.ic[i], values.ib[i], values.ie[i]);
    }
    return 0;
}
#endif
