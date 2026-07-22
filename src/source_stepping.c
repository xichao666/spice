#include "source_stepping.h"

#include <math.h>
#include <string.h>

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
    int iterations = 0;
    int total_iterations = 0;

    if (n <= 0 || n > DC_MAX_UNKNOWNS) {
        return false;
    }

    /* lambda=0：独立源关闭，先求容易的辅助问题。 */
    memset(x, 0, sizeof(double) * (size_t)n);
    if (!dc_newton_solve(problem, options, lambda, x, &iterations)) {
        return false;
    }
    total_iterations += iterations;
    if (step_callback != NULL) {
        step_callback(lambda, x, iterations, user_data);
    }

    /* 外层延拓：每次只增加一段 lambda。 */
    while (lambda < 1.0 - 1.0e-15) {
        const double target = fmin(lambda + lambda_step, 1.0);
        double trial[DC_MAX_UNKNOWNS];

        /* 继承上一个成功工作点，作为 target 的 Newton 初值。 */
        memcpy(trial, x, sizeof(double) * (size_t)n);

        if (dc_newton_solve(problem, options, target, trial, &iterations)) {
            lambda = target;
            memcpy(x, trial, sizeof(double) * (size_t)n);
            total_iterations += iterations;

            if (step_callback != NULL) {
                step_callback(lambda, x, iterations, user_data);
            }

            /* 收敛快则加大步长，收敛慢则减小步长。 */
            if (iterations <= 5) {
                lambda_step =
                    fmin(lambda_step * 1.5, options->maximum_lambda_step);
            } else if (iterations >= 15) {
                lambda_step *= 0.5;
            }
        } else {
            /* 失败时不接受 trial，保留旧解并缩小步长。 */
            lambda_step *= 0.5;
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
