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
        .maximum_junction_voltage_step = 0.10,
        .initial_lambda_step = 0.25,
        .minimum_lambda_step = 1.0e-6,
        .maximum_lambda_step = 0.5,
        .fast_newton_iteration_threshold = 5,
        .slow_newton_iteration_threshold = 15,
        .lambda_step_growth_factor = 1.5,
        .lambda_step_shrink_factor = 0.5,
        .source_step_policy = DC_SOURCE_STEP_LEGACY
    };
}

/* 固定 lambda 执行 Newton-Raphson，并统计线搜索回退情况。 */
bool dc_newton_solve_with_report(
    const DcProblem *problem,
    const DcSolverOptions *options,
    double lambda,
    double *x,
    DcNewtonReport *report)
{
    const int n = problem->dimension; /* 当前电路未知量个数。 */

    if (report == NULL) {
        return false;
    }
    *report = (DcNewtonReport) { 0, 0, NAN };

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
            report->iterations = iteration + 1;
            report->final_residual_norm = old_norm;
            return false;
        }
        if (old_norm < options->residual_tolerance) {
            report->iterations = iteration;
            report->final_residual_norm = old_norm;
            return true;
        }

        /* 2. 构造 Jacobian，并建立 J * delta = -F。 */
        problem->build_jacobian(problem->context, x, lambda, jacobian);
        for (int i = 0; i < n; ++i) {
            rhs[i] = -residual[i];
        }

        if (!solve_dense_system(n, jacobian, rhs, delta)) {
            report->iterations = iteration + 1;
            report->final_residual_norm = old_norm;
            return false;
        }

        /*
         * 3. First try the ordinary damped Newton direction.  If that cannot
         * reduce the residual, retry once with the circuit-specific junction
         * limiter.  This keeps normal solves fast while retaining a recovery
         * path for difficult exponential-device steps.
         */
        bool accepted = false;
        int reductions_this_iteration = 0;
        const int limiter_attempts =
            problem->limit_newton_step != NULL &&
            options->maximum_junction_voltage_step > 0.0 ? 2 : 1;

        for (int attempt = 0;
             attempt < limiter_attempts && !accepted;
             ++attempt) {
            if (attempt == 1) {
                problem->limit_newton_step(
                    problem->context, x, delta,
                    options->maximum_junction_voltage_step);
            }

            double scale = 1.0;
            for (int i = 0; i < n; ++i) {
                if (fabs(delta[i]) > options->maximum_component_step) {
                    const double candidate =
                        options->maximum_component_step / fabs(delta[i]);
                    if (candidate < scale) scale = candidate;
                }
            }

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
                ++reductions_this_iteration;
            }
        }
        report->line_search_reductions += reductions_this_iteration;

        if (!accepted) {
            report->iterations = iteration + 1;
            report->final_residual_norm = old_norm;
            return false;
        }

        /* 5. 变量变化很小时，再检查残差，避免假收敛。 */
        if (infinity_norm(actual_step, n) < options->voltage_tolerance) {
            problem->build_residual(problem->context, x, lambda, residual);
            const double final_norm = infinity_norm(residual, n);
            if (final_norm < options->residual_tolerance) {
                report->iterations = iteration + 1;
                report->final_residual_norm = final_norm;
                return true;
            }
        }
    }
    double final_residual[DC_MAX_UNKNOWNS];
    problem->build_residual(problem->context, x, lambda, final_residual);
    report->iterations = options->maximum_newton_iterations;
    report->final_residual_norm = infinity_norm(final_residual, n);
    return false;
}

/* 保留原有简化接口，供已有示例程序继续使用。 */
bool dc_newton_solve(
    const DcProblem *problem,
    const DcSolverOptions *options,
    double lambda,
    double *x,
    int *iteration_count)
{
    DcNewtonReport report;
    const bool converged =
        dc_newton_solve_with_report(problem, options, lambda, x, &report);

    if (iteration_count != NULL) {
        *iteration_count = report.iterations;
    }
    return converged;
}
