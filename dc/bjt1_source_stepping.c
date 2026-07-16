#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

/* 节点连接：
 *   V1: T1 collector / RC1 / R3
 *   V2: common-emitter / Re
 *   V3: T2 collector / RC2
 *   V4: T2 base / R3
 *
 * Source stepping 改变两个激励
 *   VCC(lambda) = lambda * VCC
 *   VIN(lambda) = lambda * VIN.
 *
 * BJT采用 Ebers-Moll 模型，同时考虑 VBE 与 VBC，能够描述饱和区。
 */

enum { V1 = 0, V2 = 1, V3 = 2, V4 = 3, N = 4 };/* 常量 */

/* 电路参数结构体 */
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
} Circuit;

/* BJT值结构体*/
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
} BjtValues;

/* 防止 Newton 迭代过程中 exp() 溢出 */
static double safe_exp(double x)
{
    
    if (x > 40.0) {
        x = 40.0;
    } else if (x < -40.0) {
        x = -40.0;
    }

    return exp(x);
}

/* Calculate Ebers-Moll currents.  IC and IB enter the collector and base;
 * IE leaves the emitter and therefore flows into the common-emitter node. */
static void evaluate_bjts(const Circuit *c, const double x[N],
                          double lambda, BjtValues *b)
{
    const double vbe1 = lambda * c->vin - x[V2];
    const double vbe2 = x[V4] - x[V2];
    const double vbc1 = lambda * c->vin - x[V1];
    const double vbc2 = x[V4] - x[V3];
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

/* 构建残差向量,每个节点方向为流出的电流和 */
static void build_residual(const Circuit *c, const double x[N],
                           double lambda, double f[N])
{
    BjtValues b;
    evaluate_bjts(c, x, lambda, &b);

    f[V1] = (x[V1] - lambda * c->vcc) / c->rc1
          + (x[V1] - x[V4]) / c->r3
          + b.ic1;

    f[V2] = x[V2] / c->re - b.ie1 - b.ie2;

    f[V3] = (x[V3] - lambda * c->vcc) / c->rc2 + b.ic2;

    f[V4] = (x[V4] - x[V1]) / c->r3 + b.ib2;
}

/* 构建Jacobian矩阵 */
static void build_jacobian(const Circuit *c, const double x[N],
                           double lambda, double j[N][N])
{
    BjtValues b;
    evaluate_bjts(c, x, lambda, &b);
    memset(j, 0, sizeof(double) * N * N);

    j[V1][V1] = 1.0 / c->rc1 + 1.0 / c->r3 + b.gr1;
    j[V1][V2] = -c->alpha_f * b.gf1;
    j[V1][V4] = -1.0 / c->r3;

    j[V2][V1] = -c->alpha_r * b.gr1;
    j[V2][V2] = 1.0 / c->re + b.gf1 + b.gf2;
    j[V2][V3] = -c->alpha_r * b.gr2;
    j[V2][V4] = -b.gf2 + c->alpha_r * b.gr2;

    j[V3][V2] = -c->alpha_f * b.gf2;
    j[V3][V3] = 1.0 / c->rc2 + b.gr2;
    j[V3][V4] = c->alpha_f * b.gf2 - b.gr2;

    j[V4][V1] = -1.0 / c->r3;
    j[V4][V2] = -(1.0 - c->alpha_f) * b.gf2;
    j[V4][V3] = -(1.0 - c->alpha_r) * b.gr2;
    j[V4][V4] = 1.0 / c->r3
              + (1.0 - c->alpha_f) * b.gf2
              + (1.0 - c->alpha_r) * b.gr2;
}

/*返回max vi*/
static double infinity_norm(const double v[N])
{
    double maximum = 0.0;
    for (int i = 0; i < N; ++i) {
        if (!isfinite(v[i])) {
            return NAN;
        }
        const double value = fabs(v[i]);
        if (value > maximum) {
            maximum = value;
        }
    }
    return maximum;
}

/* 高斯消元法求矩阵方程 Ax = b, A为Jacobian矩阵，b为-F，x为Delta x */
static bool solve_4x4(double a_in[N][N], const double b_in[N],
                      double solution[N])
{
    double a[N][N];
    double b[N];
    memcpy(a, a_in, sizeof(a));
    memcpy(b, b_in, sizeof(b));

    for (int column = 0; column < N; ++column) {
        int pivot = column;
        for (int row = column + 1; row < N; ++row) {
            if (fabs(a[row][column]) > fabs(a[pivot][column])) {
                pivot = row;
            }
        }
        if (!isfinite(a[pivot][column]) || fabs(a[pivot][column]) < DBL_EPSILON) {
            return false;
        }
        if (pivot != column) {
            for (int k = column; k < N; ++k) {
                const double temporary = a[column][k];
                a[column][k] = a[pivot][k];
                a[pivot][k] = temporary;
            }
            const double temporary = b[column];
            b[column] = b[pivot];
            b[pivot] = temporary;
        }/* 当最佳主元不在当前行时，交换矩阵两行 */
        for (int row = column + 1; row < N; ++row) {
            const double factor = a[row][column] / a[column][column];
            for (int k = column; k < N; ++k) {
                a[row][k] -= factor * a[column][k];
            }
            b[row] -= factor * b[column];
        }
    }
    for (int row = N - 1; row >= 0; --row) {
        double sum = b[row];
        for (int column = row + 1; column < N; ++column) {
            sum -= a[row][column] * solution[column];
        }
        solution[row] = sum / a[row][row];
        if (!isfinite(solution[row])) {
            return false;
        }
    }
    return true;
}

/* 固定 lambda，使用 Newton-Raphson 求工作点。*/
static bool newton_solve(const Circuit *c, double lambda, double x[N],
                         int *iteration_count)
{
    const int maximum_iterations = 80;
    const double residual_tolerance = 1.0e-11;
    const double voltage_tolerance = 1.0e-10;
    const double maximum_component_step = 0.20;

    for (int iteration = 0; iteration < maximum_iterations; ++iteration) {
        double f[N], j[N][N], rhs[N], delta[N];
        build_residual(c, x, lambda, f);
        if (!isfinite(infinity_norm(f))) {
            return false;
        }
        if (infinity_norm(f) < residual_tolerance) {
            *iteration_count = iteration;
            return true;
        }

        build_jacobian(c, x, lambda, j);
        for (int i = 0; i < N; ++i) {
            rhs[i] = -f[i];
        }
        if (!solve_4x4(j, rhs, delta)) {
            return false;
        }

        /* 限制最大的单个节点电压变化 */
        double scale = 1.0;
        for (int i = 0; i < N; ++i) {
            if (fabs(delta[i]) > maximum_component_step) {
                const double candidate = maximum_component_step / fabs(delta[i]);
                if (candidate < scale) {
                    scale = candidate;
                }
            }
        }

        const double old_norm = infinity_norm(f);
        bool accepted = false;
        double actual_step[N] = {0.0, 0.0, 0.0, 0.0};
        for (double alpha = scale; alpha >= scale / 1024.0; alpha *= 0.5) {
            double candidate[N], candidate_f[N];
            for (int i = 0; i < N; ++i) {
                candidate[i] = x[i] + alpha * delta[i];
            }
            build_residual(c, candidate, lambda, candidate_f);
            if (isfinite(infinity_norm(candidate_f)) &&
                infinity_norm(candidate_f) < old_norm) {
                for (int i = 0; i < N; ++i) {
                    actual_step[i] = candidate[i] - x[i];
                    x[i] = candidate[i];
                }
                accepted = true;
                break;
            }
        }
        if (!accepted) {
            return false;
        }
        if (infinity_norm(actual_step) < voltage_tolerance) {
            build_residual(c, x, lambda, f);
            if (infinity_norm(f) < residual_tolerance) {
                *iteration_count = iteration + 1;
                return true;
            }
        }
    }
    return false;
}

static bool source_stepping_dc(const Circuit *c, double x[N])
{
    const double minimum_step = 1.0e-6;
    const double maximum_step = 0.25;
    double lambda = 0.0;
    double lambda_step = 0.05;
    int iterations = 0;

    memset(x, 0, sizeof(double) * N);
    if (!newton_solve(c, lambda, x, &iterations)) {
        return false;
    }

    puts(" lambda       V1(V)       V2(V)       V3(V)       V4(V)   Newton");
    printf(" %6.4f  %10.6f  %10.6f  %10.6f  %10.6f  %7d\n",
           lambda, x[V1], x[V2], x[V3], x[V4], iterations);

    while (lambda < 1.0 - 1.0e-15) {
        const double target = fmin(lambda + lambda_step, 1.0);
        double trial[N];
        memcpy(trial, x, sizeof(trial));

        if (newton_solve(c, target, trial, &iterations)) {
            lambda = target;
            memcpy(x, trial, sizeof(trial));
            printf(" %6.4f  %10.6f  %10.6f  %10.6f  %10.6f  %7d\n",
                   lambda, x[V1], x[V2], x[V3], x[V4], iterations);
            if (iterations <= 5) {
                lambda_step = fmin(lambda_step * 1.5, maximum_step);
            } else if (iterations >= 15) {
                lambda_step *= 0.5;
            }/* Newton 成功 */
        } else {
            lambda_step *= 0.5;
            fprintf(stderr, "lambda=%.8f failed; retrying with step %.3e\n",
                    target, lambda_step);
            if (lambda_step < minimum_step) {
                return false;
            }/* Newton 失败，尝试减小步长 */
        }
    }
    return true;
}

int main(void)
{
    /* 电路参数 */
    const Circuit circuit = {
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
    double x[N];
    BjtValues b;

    if (!source_stepping_dc(&circuit, x)) {
        fputs("Source stepping could not reach lambda = 1.\n", stderr);
        return 1;
    }

    evaluate_bjts(&circuit, x, 1.0, &b);
    puts("\nFinal DC operating point (lambda = 1):");
    printf("V1 = %.9f V\nV2 = %.9f V\nV3 = %.9f V\nV4 = %.9f V\n",
           x[V1], x[V2], x[V3], x[V4]);
    printf("IC1 = %.9e A, IE1 = %.9e A\n", b.ic1, b.ie1);
    printf("IC2 = %.9e A, IE2 = %.9e A\n", b.ic2, b.ie2);
    return 0;
}
