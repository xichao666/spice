#include "gmin_stepping.h"

#include <math.h>
#include <string.h>

/*
 * 通用 GMIN Stepping 外层循环。
 * 求解器先让每个节点对地并联较大的附加电导，再逐步减小该电导，最后取为零。
 */
bool dc_gmin_stepping_solve(
    const DcProblem *problem,
    const DcSolverOptions *options,
    double *x,
    DcGminStepCallback step_callback,
    void *user_data,
    int *total_newton_iterations)
{
    const int n = problem->dimension;
    double gmin = options->initial_gmin;
    double reduction_factor = options->gmin_reduction_factor;
    int total_iterations = 0;
    DcNewtonReport report;

    if (n <= 0 || n > DC_MAX_UNKNOWNS ||
        problem->set_gmin == NULL ||
        options->initial_gmin <= 0.0 ||
        options->minimum_gmin <= 0.0 ||
        options->initial_gmin < options->minimum_gmin ||
        options->gmin_reduction_factor <= 0.0 ||
        options->gmin_reduction_factor >= 1.0) {
        return false;
    }

    /* 先以较大的附加电导求解；真实独立源保持完整开启。 */
    memset(x, 0, sizeof(double) * (size_t)n);
    problem->set_gmin(problem->context, gmin);
    if (!dc_newton_solve_with_report(problem, options, 1.0, x, &report)) {
        problem->set_gmin(problem->context, 0.0);
        if (total_newton_iterations != NULL) {
            *total_newton_iterations = report.iterations;
        }
        return false;
    }
    total_iterations += report.iterations;
    if (step_callback != NULL) {
        step_callback(gmin, x, &report, user_data);
    }

    /* 每个成功点均作为更小 GMIN 问题的 Newton 初值。 */
    while (gmin > options->minimum_gmin * (1.0 + 1.0e-12)) {
        const double target = fmax(gmin * reduction_factor,
                                   options->minimum_gmin);
        double trial[DC_MAX_UNKNOWNS];

        memcpy(trial, x, sizeof(double) * (size_t)n);
        problem->set_gmin(problem->context, target);
        if (dc_newton_solve_with_report(problem, options, 1.0,
                                        trial, &report)) {
            gmin = target;
            memcpy(x, trial, sizeof(double) * (size_t)n);
            total_iterations += report.iterations;
            if (step_callback != NULL) {
                step_callback(gmin, x, &report, user_data);
            }
        } else {
            /* 目标过激进时，保留旧解并采用更温和的电导下降比例。 */
            total_iterations += report.iterations;
            problem->set_gmin(problem->context, gmin);
            reduction_factor = sqrt(reduction_factor);
            if (reduction_factor >= 1.0 - 1.0e-12) {
                problem->set_gmin(problem->context, 0.0);
                if (total_newton_iterations != NULL) {
                    *total_newton_iterations = total_iterations;
                }
                return false;
            }
        }
    }

    /* 最后撤掉全部附加电导，得到原始电路的 DC 工作点。 */
    {
        double trial[DC_MAX_UNKNOWNS];
        memcpy(trial, x, sizeof(double) * (size_t)n);
        problem->set_gmin(problem->context, 0.0);
        if (!dc_newton_solve_with_report(problem, options, 1.0,
                                         trial, &report)) {
            total_iterations += report.iterations;
            if (total_newton_iterations != NULL) {
                *total_newton_iterations = total_iterations;
            }
            return false;
        }
        memcpy(x, trial, sizeof(double) * (size_t)n);
        total_iterations += report.iterations;
        if (step_callback != NULL) {
            step_callback(0.0, x, &report, user_data);
        }
    }

    if (total_newton_iterations != NULL) {
        *total_newton_iterations = total_iterations;
    }
    return true;
}
