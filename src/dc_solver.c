#include "dc_solver.h"

#include <float.h>
#include <math.h>
#include <string.h>

/* 返回向量的无穷范数 max(|v_i|)，同时检查 NaN/无穷大。 */
static double infinity_norm(const double *vector, int dimension)
{
    double maximum = 0.0;

    for (int i = 0; i < dimension; ++i) {
        if (!isfinite(vector[i])) {
            return NAN;
        }
        if (fabs(vector[i]) > maximum) {
            maximum = fabs(vector[i]);
        }
    }
    return maximum;
}

/*
 * 带部分主元选取的稠密高斯消元法，求解 A * solution = b。
 * 输入矩阵会先复制到局部数组，避免破坏调用者保存的 Jacobian。
 */
static bool solve_dense_system(
    int dimension,
    double input_a[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS],
    const double *input_b,
    double *solution)
{
    double a[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS];
    double b[DC_MAX_UNKNOWNS];

    memcpy(a, input_a, sizeof(a));
    memcpy(b, input_b, sizeof(b));

    /* 前向消元：依次把主元下方元素消为 0，得到上三角矩阵。 */
    for (int column = 0; column < dimension; ++column) {
        int pivot = column;

        /* 在当前列选取绝对值最大的主元，提高数值稳定性。 */
        for (int row = column + 1; row < dimension; ++row) {
            if (fabs(a[row][column]) > fabs(a[pivot][column])) {
                pivot = row;
            }
        }

        if (!isfinite(a[pivot][column]) ||
            fabs(a[pivot][column]) < DBL_EPSILON) {
            return false;
        }

        /* 若最佳主元不在当前行，则同时交换 A 和 b 的对应行。 */
        if (pivot != column) {
            for (int k = column; k < dimension; ++k) {
                const double temporary = a[column][k];
                a[column][k] = a[pivot][k];
                a[pivot][k] = temporary;
            }
            const double temporary = b[column];
            b[column] = b[pivot];
            b[pivot] = temporary;
        }

        /* 用当前主元行消去其下方各行。 */
        for (int row = column + 1; row < dimension; ++row) {
            const double factor = a[row][column] / a[column][column];

            for (int k = column; k < dimension; ++k) {
                a[row][k] -= factor * a[column][k];
            }
            b[row] -= factor * b[column];
        }
    }

    /* 回代：从最后一个未知量开始得到完整的 Newton 修正量。 */
    for (int row = dimension - 1; row >= 0; --row) {
        double sum = b[row];

        for (int column = row + 1; column < dimension; ++column) {
            sum -= a[row][column] * solution[column];
        }
        solution[row] = sum / a[row][row];

        if (!isfinite(solution[row])) {
            return false;
        }
    }
    return true;
}

/* 默认数值参数可由调用者复制后按不同电路的收敛特性调整。 */
DcSolverOptions dc_solver_default_options(void)
{
    return (DcSolverOptions) {
        .maximum_newton_iterations = 80,
        .residual_tolerance = 1.0e-11,
        .voltage_tolerance = 1.0e-10,
        .maximum_component_step = 0.20,
        .initial_lambda_step = 0.05,
        .minimum_lambda_step = 1.0e-6,
        .maximum_lambda_step = 0.25
    };
}

bool dc_newton_solve(
    const DcProblem *problem,
    const DcSolverOptions *options,
    double lambda,
    double *x,
    int *iteration_count)
{
    const int n = problem->dimension; /* 当前电路未知量个数。 */

    if (n <= 0 || n > DC_MAX_UNKNOWNS ||
        problem->build_residual == NULL || problem->build_jacobian == NULL) {
        return false;
    }

    /* 固定 lambda 后的 Newton-Raphson 内层循环。 */
    for (int iteration = 0;
         iteration < options->maximum_newton_iterations;
         ++iteration) {
        double residual[DC_MAX_UNKNOWNS];
        double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS];
        double rhs[DC_MAX_UNKNOWNS];
        double delta[DC_MAX_UNKNOWNS] = {0.0};
        double actual_step[DC_MAX_UNKNOWNS] = {0.0};

        /* 1. 在当前猜测 x 上计算非线性方程残差 F(x,lambda)。 */
        problem->build_residual(problem->context, x, lambda, residual);
        const double old_norm = infinity_norm(residual, n);

        if (!isfinite(old_norm)) {
            return false;
        }
        if (old_norm < options->residual_tolerance) {
            *iteration_count = iteration;
            return true;
        }

        /* 2. 构造 Jacobian，并建立 J * delta = -F。 */
        problem->build_jacobian(problem->context, x, lambda, jacobian);
        for (int i = 0; i < n; ++i) {
            rhs[i] = -residual[i];
        }

        if (!solve_dense_system(n, jacobian, rhs, delta)) {
            return false;
        }

        /* 3. 先限制任何单个未知量的最大变化量。 */
        double scale = 1.0;
        for (int i = 0; i < n; ++i) {
            if (fabs(delta[i]) > options->maximum_component_step) {
                const double candidate =
                    options->maximum_component_step / fabs(delta[i]);
                if (candidate < scale) {
                    scale = candidate;
                }
            }
        }

        /*
         * 4. 回溯线搜索：若完整 Newton 步不能降低残差，则逐次减半。
         *    只有新点的残差更小，才接受本次更新。
         */
        bool accepted = false;
        for (double alpha = scale;
             alpha >= scale / 1024.0;
             alpha *= 0.5) {
            double candidate[DC_MAX_UNKNOWNS];
            double candidate_residual[DC_MAX_UNKNOWNS];

            for (int i = 0; i < n; ++i) {
                candidate[i] = x[i] + alpha * delta[i];
            }

            problem->build_residual(
                problem->context, candidate, lambda, candidate_residual);

            if (isfinite(infinity_norm(candidate_residual, n)) &&
                infinity_norm(candidate_residual, n) < old_norm) {
                for (int i = 0; i < n; ++i) {
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

        /* 5. 变量变化很小时，再检查残差，避免假收敛。 */
        if (infinity_norm(actual_step, n) < options->voltage_tolerance) {
            problem->build_residual(problem->context, x, lambda, residual);
            if (infinity_norm(residual, n) < options->residual_tolerance) {
                *iteration_count = iteration + 1;
                return true;
            }
        }
    }
    return false;
}
