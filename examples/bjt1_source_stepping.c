#include "dc_solver.h"
#include "source_stepping.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

/* 四个未知节点电压在解向量 x 中的位置。 */
enum {
    BJT1_V1 = 0,
    BJT1_V2 = 1,
    BJT1_V3 = 2,
    BJT1_V4 = 3,
    BJT1_DIMENSION = 4
};

/* 当前 BJT1 施密特触发器的元件参数与晶体管模型参数。 */
typedef struct {
    double rc1;
    double rc2;
    double r3;
    double re;
    double vcc;
    double vin;
    double is1;
    double is2;
    double alpha_f;
    double alpha_r;
    double vt;
} Bjt1Circuit;

/* 最终工作点中需要输出的晶体管电流。 */
typedef struct {
    double ic1;
    double ie1;
    double ic2;
    double ie2;
} Bjt1Values;

/* 内部变量：同时保存电流和 Jacobian 需要的微分电导。 */
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
} Bjt1InternalValues;

/* 限制指数自变量，避免 Newton 的异常猜测使 exp() 溢出。 */
static double safe_exp(double x)
{
    if (x > 40.0) {
        x = 40.0;
    } else if (x < -40.0) {
        x = -40.0;
    }
    return exp(x);
}

/* 依据 Ebers-Moll 模型计算两个 BJT 的电流与微分电导。 */
static void bjt1_evaluate_internal(
    const Bjt1Circuit *c,
    const double *x,
    double lambda,
    Bjt1InternalValues *b)
{
    const double vbe1 = lambda * c->vin - x[BJT1_V2];
    const double vbc1 = lambda * c->vin - x[BJT1_V1];
    const double vbe2 = x[BJT1_V4] - x[BJT1_V2];
    const double vbc2 = x[BJT1_V4] - x[BJT1_V3];
    const double ef1 = safe_exp(vbe1 / c->vt);
    const double er1 = safe_exp(vbc1 / c->vt);
    const double ef2 = safe_exp(vbe2 / c->vt);
    const double er2 = safe_exp(vbc2 / c->vt);
    const double if1 = c->is1 * (ef1 - 1.0);
    const double ir1 = c->is1 * (er1 - 1.0);
    const double if2 = c->is2 * (ef2 - 1.0);
    const double ir2 = c->is2 * (er2 - 1.0);

    b->ic1 = c->alpha_f * if1 - ir1;
    b->ib1 = (1.0 - c->alpha_f) * if1 + (1.0 - c->alpha_r) * ir1;
    b->ie1 = if1 - c->alpha_r * ir1;
    b->ic2 = c->alpha_f * if2 - ir2;
    b->ib2 = (1.0 - c->alpha_f) * if2 + (1.0 - c->alpha_r) * ir2;
    b->ie2 = if2 - c->alpha_r * ir2;
    b->gf1 = c->is1 * ef1 / c->vt;
    b->gr1 = c->is1 * er1 / c->vt;
    b->gf2 = c->is2 * ef2 / c->vt;
    b->gr2 = c->is2 * er2 / c->vt;
}

/* 构造四个节点的 KCL 残差，采用“从节点流出的电流代数和”。 */
static void bjt1_build_residual(
    const void *context,
    const double *x,
    double lambda,
    double *residual)
{
    const Bjt1Circuit *c = context;
    Bjt1InternalValues b;

    bjt1_evaluate_internal(c, x, lambda, &b);
    residual[BJT1_V1] = (x[BJT1_V1] - lambda * c->vcc) / c->rc1
                       + (x[BJT1_V1] - x[BJT1_V4]) / c->r3
                       + b.ic1;
    residual[BJT1_V2] = x[BJT1_V2] / c->re - b.ie1 - b.ie2;
    residual[BJT1_V3] = (x[BJT1_V3] - lambda * c->vcc) / c->rc2 + b.ic2;
    residual[BJT1_V4] = (x[BJT1_V4] - x[BJT1_V1]) / c->r3 + b.ib2;
}

/* 构造残差对四个节点电压的解析偏导数矩阵。 */
static void bjt1_build_jacobian(
    const void *context,
    const double *x,
    double lambda,
    double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{
    const Bjt1Circuit *c = context;
    Bjt1InternalValues b;

    bjt1_evaluate_internal(c, x, lambda, &b);
    memset(jacobian, 0, sizeof(double) * DC_MAX_UNKNOWNS * DC_MAX_UNKNOWNS);

    jacobian[BJT1_V1][BJT1_V1] = 1.0 / c->rc1 + 1.0 / c->r3 + b.gr1;
    jacobian[BJT1_V1][BJT1_V2] = -c->alpha_f * b.gf1;
    jacobian[BJT1_V1][BJT1_V4] = -1.0 / c->r3;

    jacobian[BJT1_V2][BJT1_V1] = -c->alpha_r * b.gr1;
    jacobian[BJT1_V2][BJT1_V2] = 1.0 / c->re + b.gf1 + b.gf2;
    jacobian[BJT1_V2][BJT1_V3] = -c->alpha_r * b.gr2;
    jacobian[BJT1_V2][BJT1_V4] = -b.gf2 + c->alpha_r * b.gr2;

    jacobian[BJT1_V3][BJT1_V2] = -c->alpha_f * b.gf2;
    jacobian[BJT1_V3][BJT1_V3] = 1.0 / c->rc2 + b.gr2;
    jacobian[BJT1_V3][BJT1_V4] = c->alpha_f * b.gf2 - b.gr2;

    jacobian[BJT1_V4][BJT1_V1] = -1.0 / c->r3;
    jacobian[BJT1_V4][BJT1_V2] = -(1.0 - c->alpha_f) * b.gf2;
    jacobian[BJT1_V4][BJT1_V3] = -(1.0 - c->alpha_r) * b.gr2;
    jacobian[BJT1_V4][BJT1_V4] = 1.0 / c->r3
                               + (1.0 - c->alpha_f) * b.gf2
                               + (1.0 - c->alpha_r) * b.gr2;
}

/* 返回题目给定的 BJT1 电路参数。 */
static Bjt1Circuit bjt1_default_circuit(void)
{
    return (Bjt1Circuit) {
        .rc1 = 2000.0,
        .rc2 = 1000.0,
        .r3 = 10000.0,
        .re = 100.0,
        .vcc = 10.0,
        .vin = 1.50,
        .is1 = 1.0e-16,
        .is2 = 1.0e-16,
        .alpha_f = 100.0 / 101.0,
        .alpha_r = 1.0 / 2.0,
        .vt = 1.0 / 38.78
    };
}

/* 将 BJT1 电路的残差/Jacobian 注册为通用 DcProblem。 */
static void bjt1_make_problem(const Bjt1Circuit *circuit, DcProblem *problem)
{
    *problem = (DcProblem) {
        .dimension = BJT1_DIMENSION,
        .context = circuit,
        .build_residual = bjt1_build_residual,
        .build_jacobian = bjt1_build_jacobian
    };
}

/* 根据收敛后的节点电压计算需要显示的晶体管电流。 */
static void bjt1_evaluate(
    const Bjt1Circuit *circuit,
    const double *x,
    double lambda,
    Bjt1Values *values)
{
    Bjt1InternalValues internal;

    bjt1_evaluate_internal(circuit, x, lambda, &internal);
    values->ic1 = internal.ic1;
    values->ie1 = internal.ie1;
    values->ic2 = internal.ic2;
    values->ie2 = internal.ie2;
}

#ifndef BJT1_SOURCE_STEPPING_NO_MAIN
/* 每个成功的 Source Stepping 点调用一次，用于打印延拓轨迹。 */
static void print_step(
    double lambda,
    const double *x,
    int iterations,
    void *user_data)
{
    (void)user_data;
    printf(" %6.4f  %10.6f  %10.6f  %10.6f  %10.6f  %7d\n",
           lambda,
           x[BJT1_V1], x[BJT1_V2], x[BJT1_V3], x[BJT1_V4],
           iterations);
}

int main(void)
{
    const Bjt1Circuit circuit = bjt1_default_circuit();
    const DcSolverOptions options = dc_solver_default_options();
    DcProblem problem;
    double x[DC_MAX_UNKNOWNS];
    Bjt1Values values;
    int total_newton_iterations = 0;

    bjt1_make_problem(&circuit, &problem);
    puts(" lambda       V1(V)       V2(V)       V3(V)       V4(V)   Newton");

    if (!dc_source_stepping_solve(
            &problem, &options, x, print_step, NULL,
            &total_newton_iterations)) {
        fputs("Source stepping could not reach lambda = 1.\n", stderr);
        return 1;
    }

    bjt1_evaluate(&circuit, x, 1.0, &values);
    printf("\nTotal Newton iterations across all successful lambda points: %d\n",
           total_newton_iterations);
    puts("\nFinal DC operating point (lambda = 1):");
    printf("V1 = %.9f V\nV2 = %.9f V\nV3 = %.9f V\nV4 = %.9f V\n",
           x[BJT1_V1], x[BJT1_V2], x[BJT1_V3], x[BJT1_V4]);
    printf("IC1 = %.9e A, IE1 = %.9e A\n", values.ic1, values.ie1);
    printf("IC2 = %.9e A, IE2 = %.9e A\n", values.ic2, values.ie2);
    return 0;
}
#endif
