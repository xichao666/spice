#include "../dc_solver.h"
#include "../source_stepping.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

/* 五个未知节点电压：x1 发射极，x2/T1集电极，x3/T2集电极，x4/T2基极，x5/T1基极。 */
enum {
    BJT2_X1 = 0,
    BJT2_X2 = 1,
    BJT2_X3 = 2,
    BJT2_X4 = 3,
    BJT2_X5 = 4,
    BJT2_DIMENSION = 5
};

/* 第二个施密特触发器的元件参数和 Ebers-Moll 参数。 */
typedef struct {
    double r1;
    double r2;
    double r3;
    double r4;
    double rc1;
    double rc2;
    double re;
    double vcc;
    double is;
    double alpha_f;
    double alpha_r;
    double vt;
} Bjt2Circuit;

typedef struct {
    double ic1;
    double ie1;
    double ic2;
    double ie2;
} Bjt2Values;

typedef struct {
    double ic1;
    double ib1;
    double ie1;
    double ic2;
    double ib2;
    double ie2;
    double gf1;
    double gr1;
    double gf2;
    double gr2;
} Bjt2InternalValues;

/* 限制指数自变量，防止不合适的 Newton 猜测使 exp() 溢出。 */
static double bjt2_safe_exp(double x)
{
    if (x > 40.0) {
        x = 40.0;
    } else if (x < -40.0) {
        x = -40.0;
    }
    return exp(x);
}

/* 采用与第一个电路相同的 Ebers-Moll 电流符号约定。 */
static void bjt2_evaluate_internal(
    const Bjt2Circuit *c,
    const double *x,
    Bjt2InternalValues *b)
{
    const double vbe1 = x[BJT2_X5] - x[BJT2_X1];
    const double vbc1 = x[BJT2_X5] - x[BJT2_X2];
    const double vbe2 = x[BJT2_X4] - x[BJT2_X1];
    const double vbc2 = x[BJT2_X4] - x[BJT2_X3];
    const double ef1 = bjt2_safe_exp(vbe1 / c->vt);
    const double er1 = bjt2_safe_exp(vbc1 / c->vt);
    const double ef2 = bjt2_safe_exp(vbe2 / c->vt);
    const double er2 = bjt2_safe_exp(vbc2 / c->vt);
    const double if1 = c->is * (ef1 - 1.0);
    const double ir1 = c->is * (er1 - 1.0);
    const double if2 = c->is * (ef2 - 1.0);
    const double ir2 = c->is * (er2 - 1.0);

    b->ic1 = c->alpha_f * if1 - ir1;
    b->ib1 = (1.0 - c->alpha_f) * if1 + (1.0 - c->alpha_r) * ir1;
    b->ie1 = if1 - c->alpha_r * ir1;
    b->ic2 = c->alpha_f * if2 - ir2;
    b->ib2 = (1.0 - c->alpha_f) * if2 + (1.0 - c->alpha_r) * ir2;
    b->ie2 = if2 - c->alpha_r * ir2;
    b->gf1 = c->is * ef1 / c->vt;
    b->gr1 = c->is * er1 / c->vt;
    b->gf2 = c->is * ef2 / c->vt;
    b->gr2 = c->is * er2 / c->vt;
}

/*
 * 五个 KCL 残差：
 * x1: 公共发射极；x2: T1 集电极；x3: T2 集电极；
 * x4: T2 基极；x5: T1 基极。
 * VCC 在 Source Stepping 中缩放为 lambda * VCC。
 */
static void bjt2_build_residual(
    const void *context,
    const double *x,
    double lambda,
    double *residual)
{
    const Bjt2Circuit *c = context;
    const double vcc = lambda * c->vcc;
    Bjt2InternalValues b;

    bjt2_evaluate_internal(c, x, &b);
    residual[BJT2_X1] = x[BJT2_X1] / c->re - b.ie1 - b.ie2;
    residual[BJT2_X2] = (x[BJT2_X2] - x[BJT2_X4]) / c->r1
                       + (x[BJT2_X2] - vcc) / c->rc1 + b.ic1;
    residual[BJT2_X3] = (x[BJT2_X3] - vcc) / c->rc2 + b.ic2;
    residual[BJT2_X4] = (x[BJT2_X4] - x[BJT2_X2]) / c->r1
                       + x[BJT2_X4] / c->r4 + b.ib2;
    residual[BJT2_X5] = (x[BJT2_X5] - vcc) / c->r2
                       + x[BJT2_X5] / c->r3 + b.ib1;
}

/* 构造上述五个 KCL 方程对 x1 到 x5 的解析 Jacobian。 */
static void bjt2_build_jacobian(
    const void *context,
    const double *x,
    double lambda,
    double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{
    const Bjt2Circuit *c = context;
    Bjt2InternalValues b;

    (void)lambda; /* 本电路的 Jacobian 不直接包含独立源的幅值。 */
    bjt2_evaluate_internal(c, x, &b);
    memset(jacobian, 0, sizeof(double) * DC_MAX_UNKNOWNS * DC_MAX_UNKNOWNS);

    jacobian[BJT2_X1][BJT2_X1] = 1.0 / c->re + b.gf1 + b.gf2;
    jacobian[BJT2_X1][BJT2_X2] = -c->alpha_r * b.gr1;
    jacobian[BJT2_X1][BJT2_X3] = -c->alpha_r * b.gr2;
    jacobian[BJT2_X1][BJT2_X4] = -b.gf2 + c->alpha_r * b.gr2;
    jacobian[BJT2_X1][BJT2_X5] = -b.gf1 + c->alpha_r * b.gr1;

    jacobian[BJT2_X2][BJT2_X1] = -c->alpha_f * b.gf1;
    jacobian[BJT2_X2][BJT2_X2] = 1.0 / c->r1 + 1.0 / c->rc1 + b.gr1;
    jacobian[BJT2_X2][BJT2_X4] = -1.0 / c->r1;
    jacobian[BJT2_X2][BJT2_X5] = c->alpha_f * b.gf1 - b.gr1;

    jacobian[BJT2_X3][BJT2_X1] = -c->alpha_f * b.gf2;
    jacobian[BJT2_X3][BJT2_X3] = 1.0 / c->rc2 + b.gr2;
    jacobian[BJT2_X3][BJT2_X4] = c->alpha_f * b.gf2 - b.gr2;

    jacobian[BJT2_X4][BJT2_X2] = -1.0 / c->r1;
    jacobian[BJT2_X4][BJT2_X1] = -(1.0 - c->alpha_f) * b.gf2;
    jacobian[BJT2_X4][BJT2_X3] = -(1.0 - c->alpha_r) * b.gr2;
    jacobian[BJT2_X4][BJT2_X4] = 1.0 / c->r1 + 1.0 / c->r4
                               + (1.0 - c->alpha_f) * b.gf2
                               + (1.0 - c->alpha_r) * b.gr2;

    jacobian[BJT2_X5][BJT2_X1] = -(1.0 - c->alpha_f) * b.gf1;
    jacobian[BJT2_X5][BJT2_X2] = -(1.0 - c->alpha_r) * b.gr1;
    jacobian[BJT2_X5][BJT2_X5] = 1.0 / c->r2 + 1.0 / c->r3
                               + (1.0 - c->alpha_f) * b.gf1
                               + (1.0 - c->alpha_r) * b.gr1;
}

/* 图中给定的第二个施密特触发器参数。 */
static Bjt2Circuit bjt2_default_circuit(void)
{
    return (Bjt2Circuit) {
        .r1 = 10000.0,
        .r2 = 5000.0,
        .r3 = 1250.0,
        .r4 = 1.0e6,
        .rc1 = 1500.0,
        .rc2 = 1000.0,
        .re = 100.0,
        .vcc = 10.0,
        .is = 1.0e-16,
        .alpha_f = 0.99,
        .alpha_r = 0.5,
        .vt = 1.0 / 38.78
    };
}

/* 将 BJT2 电路注册为通用 DC 问题。 */
static void bjt2_make_problem(const Bjt2Circuit *circuit, DcProblem *problem)
{
    *problem = (DcProblem) {
        .dimension = BJT2_DIMENSION,
        .context = circuit,
        .build_residual = bjt2_build_residual,
        .build_jacobian = bjt2_build_jacobian
    };
}

/* 根据最终节点电压计算并输出两个晶体管的电流。 */
static void bjt2_evaluate(
    const Bjt2Circuit *circuit,
    const double *x,
    Bjt2Values *values)
{
    Bjt2InternalValues internal;

    bjt2_evaluate_internal(circuit, x, &internal);
    values->ic1 = internal.ic1;
    values->ie1 = internal.ie1;
    values->ic2 = internal.ic2;
    values->ie2 = internal.ie2;
}

#ifndef BJT2_SOURCE_STEPPING_NO_MAIN
/* 打印每一个成功的延拓点。 */
static void bjt2_print_step(
    double lambda,
    const double *x,
    int iterations,
    void *user_data)
{
    (void)user_data;
    printf(" %6.4f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %7d\n",
           lambda, x[BJT2_X1], x[BJT2_X2], x[BJT2_X3],
           x[BJT2_X4], x[BJT2_X5], iterations);
}

int main(void)
{
    const Bjt2Circuit circuit = bjt2_default_circuit();
    const DcSolverOptions options = dc_solver_default_options();
    DcProblem problem;
    double x[DC_MAX_UNKNOWNS];
    Bjt2Values values;
    int total_newton_iterations = 0;

    bjt2_make_problem(&circuit, &problem);
    puts(" lambda       X1(V)       X2(V)       X3(V)       X4(V)       X5(V)   Newton");

    if (!dc_source_stepping_solve(
            &problem, &options, x, bjt2_print_step, NULL,
            &total_newton_iterations)) {
        fputs("Source stepping could not reach lambda = 1.\n", stderr);
        return 1;
    }

    bjt2_evaluate(&circuit, x, &values);
    printf("\nTotal Newton iterations across all successful lambda points: %d\n",
           total_newton_iterations);
    puts("\nFinal DC operating point (lambda = 1):");
    printf("X1 = %.9f V\nX2 = %.9f V\nX3 = %.9f V\nX4 = %.9f V\nX5 = %.9f V\n",
           x[BJT2_X1], x[BJT2_X2], x[BJT2_X3], x[BJT2_X4], x[BJT2_X5]);
    printf("IC1 = %.9e A, IE1 = %.9e A\n", values.ic1, values.ie1);
    printf("IC2 = %.9e A, IE2 = %.9e A\n", values.ic2, values.ie2);
    return 0;
}
#endif
