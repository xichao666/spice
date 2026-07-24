#include "homotopy.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

/* 最小固定点同伦：F(x)=x-2 在 lambda=1 有唯一根 x=2。 */
static void build_residual(const void *context, const double *x,
                           double lambda, double *residual)
{
    (void)context;
    (void)lambda;
    residual[0] = x[0] - 2.0;
}

static void build_jacobian(
    const void *context, const double *x, double lambda,
    double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{
    (void)context;
    (void)x;
    (void)lambda;
    memset(jacobian, 0, sizeof(double) * DC_MAX_UNKNOWNS * DC_MAX_UNKNOWNS);
    jacobian[0][0] = 1.0;
}

typedef struct {
    int count;
    double solution;
    double residual_norm;
    double maximum_lambda;
} TestResult;

static void record_path(int step, double arc_length, double lambda,
                        double arc_step, double local_error, const double *x,
                        void *user_data)
{
    TestResult *result = user_data;
    (void)step;
    (void)arc_length;
    (void)arc_step;
    (void)local_error;
    (void)x;
    if (lambda > result->maximum_lambda) result->maximum_lambda = lambda;
}

/* Homotopy 模块只应在原始 F=0 的 Newton 复核后调用此回调。 */
static void receive_solution(int index, const double *x, double residual_norm,
                             void *user_data)
{
    TestResult *result = user_data;
    (void)index;
    ++result->count;
    result->solution = x[0];
    result->residual_norm = residual_norm;
}

int main(void)
{
    const DcProblem problem = {
        .dimension = 1,
        .context = NULL,
        .build_residual = build_residual,
        .build_jacobian = build_jacobian,
        .limit_newton_step = NULL
    };
    const DcSolverOptions newton_options = dc_solver_default_options();
    const double starting_point[] = {0.0};
    DcHomotopyOptions options = dc_homotopy_default_options(starting_point);
    TestResult result = {0, 0.0, INFINITY, -INFINITY};
    int solution_count = 0;

    options.initial_arc_step = 0.01;
    options.maximum_arc_step = 0.05;
    options.maximum_path_steps = 10000;
    if (!dc_fixed_point_homotopy_solve(
            &problem, &newton_options, &options, record_path, receive_solution,
            &result, &solution_count) ||
        solution_count != 1 || result.count != 1 ||
        fabs(result.solution - 2.0) > 1.0e-10 ||
        result.residual_norm > 1.0e-10) {
        fprintf(stderr, "fixed-point homotopy test failed (max lambda %.9g, solutions %d).\n",
                result.maximum_lambda, solution_count);
        return 1;
    }
    puts("fixed-point homotopy test passed.");
    return 0;
}
