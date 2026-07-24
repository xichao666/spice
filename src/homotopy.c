#include "homotopy.h"

#include <float.h>
#include <math.h>
#include <string.h>

/* 状态 y=[lambda;x] 比原问题多一个分量。 */
#define HOMOTOPY_MAX_STATE (DC_MAX_UNKNOWNS + 1)

static double infinity_norm(const double *values, int count)
{
    double maximum = 0.0;
    for (int i = 0; i < count; ++i) {
        if (!isfinite(values[i])) return NAN;
        if (fabs(values[i]) > maximum) maximum = fabs(values[i]);
    }
    return maximum;
}

static double dot_product(const double *left, const double *right, int count)
{
    double sum = 0.0;
    for (int i = 0; i < count; ++i) sum += left[i] * right[i];
    return sum;
}

static double two_norm(const double *values, int count)
{
    return sqrt(dot_product(values, values, count));
}

/* 带部分主元的稠密求解；输入矩阵由调用者提供的局部副本。 */
static bool solve_dense(int n, double matrix[HOMOTOPY_MAX_STATE]
                                      [HOMOTOPY_MAX_STATE],
                        const double *input_rhs, double *solution)
{
    double rhs[HOMOTOPY_MAX_STATE];
    memcpy(rhs, input_rhs, sizeof(double) * (size_t)n);

    for (int column = 0; column < n; ++column) {
        int pivot = column;
        for (int row = column + 1; row < n; ++row) {
            if (fabs(matrix[row][column]) > fabs(matrix[pivot][column])) {
                pivot = row;
            }
        }
        if (!isfinite(matrix[pivot][column]) ||
            fabs(matrix[pivot][column]) < DBL_EPSILON) return false;
        if (pivot != column) {
            for (int k = column; k < n; ++k) {
                const double temporary = matrix[column][k];
                matrix[column][k] = matrix[pivot][k];
                matrix[pivot][k] = temporary;
            }
            const double temporary = rhs[column];
            rhs[column] = rhs[pivot];
            rhs[pivot] = temporary;
        }
        for (int row = column + 1; row < n; ++row) {
            const double factor = matrix[row][column] / matrix[column][column];
            for (int k = column; k < n; ++k) {
                matrix[row][k] -= factor * matrix[column][k];
            }
            rhs[row] -= factor * rhs[column];
        }
    }
    for (int row = n - 1; row >= 0; --row) {
        double sum = rhs[row];
        for (int column = row + 1; column < n; ++column) {
            sum -= matrix[row][column] * solution[column];
        }
        solution[row] = sum / matrix[row][row];
        if (!isfinite(solution[row])) return false;
    }
    return true;
}

static double scale_at(const DcHomotopyOptions *options, int index)
{
    return options->diagonal_scale == NULL ? options->uniform_scale :
                                              options->diagonal_scale[index];
}

/*
 * 计算根页面公式的 H、Hx 和 Hlambda。DcProblem 的 lambda=1 是原电路，
 * 因而同伦 lambda 不会改变器件或独立源 stamp 的原有定义。
 */
static void evaluate_homotopy(
    const DcProblem *problem, const DcHomotopyOptions *options,
    const double *x, double lambda, double *h, double *h_lambda,
    double h_x[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{
    const int n = problem->dimension;
    double f[DC_MAX_UNKNOWNS];
    double jf[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS];

    problem->build_residual(problem->context, x, 1.0, f);
    problem->build_jacobian(problem->context, x, 1.0, jf);
    for (int row = 0; row < n; ++row) {
        const double g = scale_at(options, row);
        h[row] = (1.0 - lambda) * g *
                     (x[row] - options->starting_point[row]) + lambda * f[row];
        h_lambda[row] = f[row] - g * (x[row] - options->starting_point[row]);
        for (int column = 0; column < n; ++column) {
            h_x[row][column] = lambda * jf[row][column];
        }
        h_x[row][row] += (1.0 - lambda) * g;
    }
}

/*
 * 对 H(y(s))=0 求导：DH * y'=0，DH=[Hlambda Hx]。
 * 额外加入一个方向约束，得到唯一切线；随后归一化并令其与旧切线同向。
 */
static bool evaluate_tangent(const DcProblem *problem,
                             const DcHomotopyOptions *options,
                             const double *state, const double *old_tangent,
                             double *tangent)
{
    const int n = problem->dimension;
    const int m = n + 1;
    double h[DC_MAX_UNKNOWNS];
    double h_lambda[DC_MAX_UNKNOWNS];
    double h_x[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS];
    double matrix[HOMOTOPY_MAX_STATE][HOMOTOPY_MAX_STATE] = {{0.0}};
    double rhs[HOMOTOPY_MAX_STATE] = {0.0};

    evaluate_homotopy(problem, options, state + 1, state[0], h, h_lambda, h_x);
    for (int row = 0; row < n; ++row) {
        matrix[row][0] = h_lambda[row];
        for (int column = 0; column < n; ++column) {
            matrix[row][column + 1] = h_x[row][column];
        }
    }
    if (old_tangent == NULL) {
        /* 起点的初始方向选择 lambda 递增，后续允许它自然折返。 */
        matrix[n][0] = 1.0;
        rhs[n] = 1.0;
    } else {
        for (int column = 0; column < m; ++column) {
            matrix[n][column] = old_tangent[column];
        }
        rhs[n] = 1.0;
    }
    if (!solve_dense(m, matrix, rhs, tangent)) return false;

    const double norm = two_norm(tangent, m);
    if (!isfinite(norm) || norm < DBL_EPSILON) return false;
    for (int i = 0; i < m; ++i) tangent[i] /= norm;
    if (old_tangent != NULL && dot_product(tangent, old_tangent, m) < 0.0) {
        for (int i = 0; i < m; ++i) tangent[i] = -tangent[i];
    }
    return true;
}

/* 将 lambda=1 附近的插值点交给原始 Newton，并按变量无穷范数去重。 */
static void record_crossing(
    const DcProblem *problem, const DcSolverOptions *newton_options,
    const double *left, const double *right, double solutions[DC_MAX_UNKNOWNS]
                                                    [DC_MAX_UNKNOWNS],
    int *count, DcHomotopySolutionCallback callback, void *user_data)
{
    const int n = problem->dimension;
    const double denominator = right[0] - left[0];
    const double ratio = fabs(denominator) < DBL_EPSILON ? 0.0 :
                         (1.0 - left[0]) / denominator;
    double guess[DC_MAX_UNKNOWNS];
    double residual[DC_MAX_UNKNOWNS];
    DcNewtonReport report;

    for (int i = 0; i < n; ++i) {
        guess[i] = left[i + 1] + ratio * (right[i + 1] - left[i + 1]);
    }
    if (!dc_newton_solve_with_report(problem, newton_options, 1.0, guess,
                                     &report)) return;
    problem->build_residual(problem->context, guess, 1.0, residual);
    const double residual_norm = infinity_norm(residual, n);
    if (!isfinite(residual_norm)) return;

    for (int solution = 0; solution < *count; ++solution) {
        double distance = 0.0;
        for (int i = 0; i < n; ++i) {
            distance = fmax(distance, fabs(guess[i] - solutions[solution][i]));
        }
        if (distance < 1.0e-6) return;
    }
    if (*count < DC_MAX_UNKNOWNS) {
        memcpy(solutions[*count], guess, sizeof(double) * (size_t)n);
        if (callback != NULL) callback(*count, guess, residual_norm, user_data);
        ++*count;
    }
}

DcHomotopyOptions dc_homotopy_default_options(const double *starting_point)
{
    return (DcHomotopyOptions) {
        .starting_point = starting_point,
        .diagonal_scale = NULL,
        .uniform_scale = 1.0e-3,
        .initial_arc_step = 0.01,
        .minimum_arc_step = 1.0e-8,
        .maximum_arc_step = 0.10,
        .local_error_tolerance = 1.0e-6,
        .crossing_tolerance = 1.0e-6,
        .maximum_path_steps = 10000
    };
}

bool dc_fixed_point_homotopy_solve(
    const DcProblem *problem, const DcSolverOptions *newton_options,
    const DcHomotopyOptions *options, DcHomotopyPathCallback path_callback,
    DcHomotopySolutionCallback solution_callback, void *user_data,
    int *solution_count)
{
    const int n = problem == NULL ? 0 : problem->dimension;
    const int m = n + 1;
    double state[HOMOTOPY_MAX_STATE] = {0.0};
    double tangent[HOMOTOPY_MAX_STATE] = {0.0};
    double solutions[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS] = {{0.0}};
    double arc_step;
    double arc_length = 0.0;
    int count = 0;

    if (solution_count != NULL) *solution_count = 0;
    if (n <= 0 || n > DC_MAX_UNKNOWNS || newton_options == NULL ||
        options == NULL || options->starting_point == NULL ||
        options->uniform_scale <= 0.0 || options->initial_arc_step <= 0.0 ||
        options->minimum_arc_step <= 0.0 ||
        options->maximum_arc_step < options->minimum_arc_step ||
        options->local_error_tolerance <= 0.0 || options->maximum_path_steps <= 0) {
        return false;
    }
    for (int i = 0; i < n; ++i) {
        if (!isfinite(scale_at(options, i)) || scale_at(options, i) <= 0.0) {
            return false;
        }
    }

    /* lambda=0 时 H=G(x-a)，故 [0;a] 无需额外 Newton 即为精确起点。 */
    memcpy(state + 1, options->starting_point, sizeof(double) * (size_t)n);
    if (!evaluate_tangent(problem, options, state, NULL, tangent)) return false;
    arc_step = options->initial_arc_step;

    for (int step = 1; step <= options->maximum_path_steps; ++step) {
        double predicted[HOMOTOPY_MAX_STATE];
        double corrected[HOMOTOPY_MAX_STATE];
        double predicted_tangent[HOMOTOPY_MAX_STATE];
        double local_error_vector[HOMOTOPY_MAX_STATE];
        double local_error;
        bool accepted = false;

        /* Euler predictor + Adams-Moulton(梯形) corrector，与参考启动段同类。 */
        while (arc_step >= options->minimum_arc_step) {
            for (int i = 0; i < m; ++i) predicted[i] = state[i] + arc_step * tangent[i];
            if (!evaluate_tangent(problem, options, predicted, tangent,
                                  predicted_tangent)) {
                arc_step *= 0.25;
                continue;
            }
            for (int i = 0; i < m; ++i) {
                corrected[i] = state[i] + 0.5 * arc_step *
                               (tangent[i] + predicted_tangent[i]);
                local_error_vector[i] = corrected[i] - predicted[i];
            }
            local_error = infinity_norm(local_error_vector, m);
            const double tolerance = options->local_error_tolerance *
                                     fmax(two_norm(state, m), 1.0);
            if (isfinite(local_error) && local_error <= tolerance) {
                accepted = true;
                break;
            }
            arc_step *= fmax(0.1, 0.8 * sqrt(tolerance / fmax(local_error, DBL_MIN)));
        }
        if (!accepted) return false;

        /* 穿过 lambda=1 后，以原始 DC 方程复核，而非只相信积分点。 */
        if ((state[0] - 1.0) * (corrected[0] - 1.0) <= 0.0 ||
            fabs(corrected[0] - 1.0) <= options->crossing_tolerance) {
            record_crossing(problem, newton_options, state, corrected, solutions,
                            &count, solution_callback, user_data);
        }

        arc_length += arc_step;
        memcpy(state, corrected, sizeof(double) * (size_t)m);
        if (!evaluate_tangent(problem, options, state, tangent, tangent)) return false;
        if (path_callback != NULL) {
            path_callback(step, arc_length, state[0], arc_step, local_error,
                          state + 1, user_data);
        }

        /* 一阶 predictor-corrector 的局部误差阶为 O(h^2)，据此调节下一步。 */
        const double tolerance = options->local_error_tolerance *
                                 fmax(two_norm(state, m), 1.0);
        const double factor = local_error < DBL_MIN ? 2.0 :
            fmin(2.0, fmax(0.5, 0.8 * sqrt(tolerance / local_error)));
        arc_step = fmin(options->maximum_arc_step, arc_step * factor);
    }
    if (solution_count != NULL) *solution_count = count;
    return true;
}
