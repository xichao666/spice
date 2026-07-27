#include "source_stepping.h"

#include <math.h>
#include <string.h>

/* 返回向量中各分量绝对值的最大值。 */
static double infinity_norm(const double *values, int count)
{
    double maximum = 0.0;
    for (int i = 0; i < count; ++i) {
        const double value = fabs(values[i]);
        if (value > maximum) maximum = value;
    }
    return maximum;
}

/*
 * 通用 Source Stepping 外层循环。
 * 电路本身只需在残差函数中将独立源写成 lambda * source_value。
 */
bool dc_source_stepping_solve(
    const DcProblem *problem,
    const DcSolverOptions *options,
    double *x,
    DcStepCallback step_callback,
    void *user_data,
    int *total_newton_iterations)
{
    const int n = problem->dimension;
    double lambda = 0.0;
    double lambda_step = options->initial_lambda_step;
    double previous_solution[DC_MAX_UNKNOWNS];
    double previous_lambda = 0.0;
    bool have_previous_solution = false;
    int previous_step_iterations = 1000000;
    int iterations = 0;
    int total_iterations = 0;

    if (n <= 0 || n > DC_MAX_UNKNOWNS) {
        return false;
    }
    if (options->initial_lambda_step <= 0.0 ||
        options->minimum_lambda_step <= 0.0 ||
        options->maximum_lambda_step < options->minimum_lambda_step ||
        options->initial_lambda_step > options->maximum_lambda_step ||
        options->maximum_lambda_step >= 1.0 ||
        options->maximum_junction_voltage_step < 0.0 ||
        options->fast_newton_iteration_threshold < 0 ||
        options->slow_newton_iteration_threshold <=
            options->fast_newton_iteration_threshold ||
        options->lambda_step_growth_factor < 1.0 ||
        options->lambda_step_shrink_factor <= 0.0 ||
        options->lambda_step_shrink_factor >= 1.0 ||
        (options->source_step_policy != DC_SOURCE_STEP_LEGACY &&
         options->source_step_policy != DC_SOURCE_STEP_SECANT_PREDICTOR)) {
        return false;
    }

    /* lambda=0：独立源关闭，先求容易的辅助问题。 */
    memset(x, 0, sizeof(double) * (size_t)n);
    DcNewtonReport report;
    if (!dc_newton_solve_with_report(problem, options, lambda, x, &report)) {
        return false;
    }
    iterations = report.iterations;
    total_iterations += iterations;
    if (step_callback != NULL) {
        step_callback(lambda, x, iterations, user_data);
    }

    /* 外层延拓：每次只增加一段 lambda。 */
    while (lambda < 1.0 - 1.0e-15) {
        const double target = fmin(lambda + lambda_step, 1.0);
        double trial[DC_MAX_UNKNOWNS];
        bool used_secant_predictor = false;

        /* 继承上一个成功工作点，作为 target 的 Newton 初值。 */
        memcpy(trial, x, sizeof(double) * (size_t)n);

        if (options->source_step_policy == DC_SOURCE_STEP_SECANT_PREDICTOR &&
            have_previous_solution && lambda > previous_lambda &&
            previous_step_iterations <= 12) {
            const double scale = (target - lambda) / (lambda - previous_lambda);
            double predicted[DC_MAX_UNKNOWNS];
            double inherited_residual[DC_MAX_UNKNOWNS];
            double predicted_residual[DC_MAX_UNKNOWNS];
            for (int i = 0; i < n; ++i) {
                predicted[i] = x[i] + scale * (x[i] - previous_solution[i]);
            }
            problem->build_residual(problem->context, trial, target,
                                    inherited_residual);
            problem->build_residual(problem->context, predicted, target,
                                    predicted_residual);
            if (isfinite(infinity_norm(predicted_residual, n)) &&
                infinity_norm(predicted_residual, n) <
                    0.25 * infinity_norm(inherited_residual, n)) {
                memcpy(trial, predicted, sizeof(double) * (size_t)n);
                used_secant_predictor = true;
            }
        }

        bool converged = dc_newton_solve_with_report(
            problem, options, target, trial, &report);
        int rejected_predictor_iterations = 0;
        if (!converged && used_secant_predictor) {
            rejected_predictor_iterations = report.iterations;
            memcpy(trial, x, sizeof(double) * (size_t)n);
            converged = dc_newton_solve_with_report(
                problem, options, target, trial, &report);
        }

        if (converged) {
            iterations = report.iterations;
            memcpy(previous_solution, x, sizeof(double) * (size_t)n);
            previous_lambda = lambda;
            have_previous_solution = true;
            previous_step_iterations = iterations;
            lambda = target;
            memcpy(x, trial, sizeof(double) * (size_t)n);
            total_iterations += rejected_predictor_iterations + iterations;

            if (step_callback != NULL) {
                step_callback(lambda, x, iterations, user_data);
            }

            /* 原规则式步长控制：收敛快则增大，收敛慢则减小。 */
            if (iterations <= options->fast_newton_iteration_threshold) {
                lambda_step = fmin(lambda_step *
                                       options->lambda_step_growth_factor,
                                   options->maximum_lambda_step);
            } else if (iterations >=
                       options->slow_newton_iteration_threshold) {
                lambda_step *= options->lambda_step_shrink_factor;
            }
        } else {
            /* 失败时不接受 trial，保留旧解并缩小步长。 */
            iterations = report.iterations;
            total_iterations += rejected_predictor_iterations + iterations;
            previous_step_iterations = 1000000;
            lambda_step *= options->lambda_step_shrink_factor;
            if (lambda_step < options->minimum_lambda_step) {
                return false;
            }
        }
    }

    if (total_newton_iterations != NULL) {
        *total_newton_iterations = total_iterations;
    }
    return true;
}

bool dc_sequential_source_stepping_solve(
    const DcProblem *problem,
    const DcSolverOptions *options,
    double *x,
    DcSequentialStepCallback step_callback,
    void *user_data,
    int *total_newton_iterations)
{
    const int n = problem->dimension;
    int total_iterations = 0;
    DcNewtonReport report;

    if (n <= 0 || n > DC_MAX_UNKNOWNS ||
        problem->independent_source_count <= 0 ||
        problem->set_source_scale == NULL ||
        options->initial_lambda_step <= 0.0 ||
        options->minimum_lambda_step <= 0.0 ||
        options->maximum_lambda_step < options->minimum_lambda_step ||
        options->initial_lambda_step > options->maximum_lambda_step ||
        options->maximum_lambda_step >= 1.0 ||
        options->maximum_junction_voltage_step < 0.0 ||
        options->fast_newton_iteration_threshold < 0 ||
        options->slow_newton_iteration_threshold <=
            options->fast_newton_iteration_threshold ||
        options->lambda_step_growth_factor < 1.0 ||
        options->lambda_step_shrink_factor <= 0.0 ||
        options->lambda_step_shrink_factor >= 1.0 ||
        (options->source_step_policy != DC_SOURCE_STEP_LEGACY &&
         options->source_step_policy != DC_SOURCE_STEP_SECANT_PREDICTOR)) {
        return false;
    }

    /* 所有独立源关闭后，先求解零源辅助电路。 */
    for (int i = 0; i < problem->independent_source_count; ++i) {
        problem->set_source_scale(problem->context, i, 0.0);
    }
    memset(x, 0, sizeof(double) * (size_t)n);
    if (!dc_newton_solve_with_report(problem, options, 1.0, x, &report)) {
        return false;
    }
    total_iterations += report.iterations;
    if (step_callback != NULL) {
        step_callback(-1, 0.0, x, report.iterations, user_data);
    }

    for (int source = 0;
         source < problem->independent_source_count;
         ++source) {
        double lambda = 0.0;
        double lambda_step = options->initial_lambda_step;
        double previous_solution[DC_MAX_UNKNOWNS];
        double previous_lambda = 0.0;
        bool have_previous_solution = false;
        int previous_step_iterations = 1000000;

        while (lambda < 1.0 - 1.0e-15) {
            const double target = fmin(lambda + lambda_step, 1.0);
            double trial[DC_MAX_UNKNOWNS];
            bool used_secant_predictor = false;
            int rejected_predictor_iterations = 0;

            /* 仅增加当前源；此前完成的源保持开启，后续源保持关闭。 */
            problem->set_source_scale(problem->context, source, target);
            memcpy(trial, x, sizeof(double) * (size_t)n);

            if (options->source_step_policy == DC_SOURCE_STEP_SECANT_PREDICTOR &&
                have_previous_solution && lambda > previous_lambda &&
                previous_step_iterations <= 12) {
                const double scale =
                    (target - lambda) / (lambda - previous_lambda);
                double predicted[DC_MAX_UNKNOWNS];
                double inherited_residual[DC_MAX_UNKNOWNS];
                double predicted_residual[DC_MAX_UNKNOWNS];

                for (int i = 0; i < n; ++i) {
                    predicted[i] = x[i] + scale * (x[i] - previous_solution[i]);
                }
                problem->build_residual(problem->context, trial, 1.0,
                                        inherited_residual);
                problem->build_residual(problem->context, predicted, 1.0,
                                        predicted_residual);
                if (isfinite(infinity_norm(predicted_residual, n)) &&
                    infinity_norm(predicted_residual, n) <
                        0.25 * infinity_norm(inherited_residual, n)) {
                    memcpy(trial, predicted, sizeof(double) * (size_t)n);
                    used_secant_predictor = true;
                }
            }

            bool converged = dc_newton_solve_with_report(
                problem, options, 1.0, trial, &report);
            if (!converged && used_secant_predictor) {
                rejected_predictor_iterations = report.iterations;
                memcpy(trial, x, sizeof(double) * (size_t)n);
                converged = dc_newton_solve_with_report(
                    problem, options, 1.0, trial, &report);
            }

            if (converged) {
                const int iterations = report.iterations;
                memcpy(previous_solution, x, sizeof(double) * (size_t)n);
                previous_lambda = lambda;
                have_previous_solution = true;
                previous_step_iterations = iterations;
                lambda = target;
                memcpy(x, trial, sizeof(double) * (size_t)n);
                total_iterations += rejected_predictor_iterations + iterations;

                if (step_callback != NULL) {
                    step_callback(source, lambda, x, iterations, user_data);
                }

                if (iterations <= options->fast_newton_iteration_threshold) {
                    lambda_step = fmin(lambda_step *
                                           options->lambda_step_growth_factor,
                                       options->maximum_lambda_step);
                } else if (iterations >=
                           options->slow_newton_iteration_threshold) {
                    lambda_step *= options->lambda_step_shrink_factor;
                }
            } else {
                total_iterations += rejected_predictor_iterations +
                                    report.iterations;
                previous_step_iterations = 1000000;
                problem->set_source_scale(problem->context, source, lambda);
                lambda_step *= options->lambda_step_shrink_factor;
                if (lambda_step < options->minimum_lambda_step) {
                    return false;
                }
            }
        }
    }

    if (total_newton_iterations != NULL) {
        *total_newton_iterations = total_iterations;
    }
    return true;
}
