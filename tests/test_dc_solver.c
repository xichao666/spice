#include "dc_solver.h"
#include "source_stepping.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

/* 构造两元线性方程 x0+x1=3、x0-x1=1 的残差。 */
static void build_residual(const void *context, const double *x,
                           double lambda, double *residual)
{
    (void)context;
    (void)lambda;
    residual[0] = x[0] + x[1] - 3.0;
    residual[1] = x[0] - x[1] - 1.0;
}

/* 构造该方程组的常量 Jacobian。 */
static void build_jacobian(const void *context, const double *x,
                           double lambda,
                           double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{
    (void)context;
    (void)x;
    (void)lambda;
    memset(jacobian, 0, sizeof(double) * DC_MAX_UNKNOWNS * DC_MAX_UNKNOWNS);
    jacobian[0][0] = 1.0;
    jacobian[0][1] = 1.0;
    jacobian[1][0] = 1.0;
    jacobian[1][1] = -1.0;
}

typedef struct {
    double scale[2];
} SequentialContext;

/* A one-variable circuit model driven by two independently scaled sources. */
static void build_sequential_residual(
    const void *context,
    const double *x,
    double lambda,
    double *residual)
{
    const SequentialContext *sources = context;
    (void)lambda;
    residual[0] = x[0] - sources->scale[0] - 2.0 * sources->scale[1];
}

static void build_sequential_jacobian(
    const void *context,
    const double *x,
    double lambda,
    double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{
    (void)context;
    (void)x;
    (void)lambda;
    memset(jacobian, 0, sizeof(double) * DC_MAX_UNKNOWNS * DC_MAX_UNKNOWNS);
    jacobian[0][0] = 1.0;
}

static void set_sequential_source_scale(
    const void *context,
    int source_index,
    double scale)
{
    SequentialContext *sources = (SequentialContext *)context;
    sources->scale[source_index] = scale;
}

/* 验证通用 Newton 求解器可正确求得已知的两元解。 */
int main(void)
{
    const DcProblem problem = {
        .dimension = 2,
        .context = NULL,
        .build_residual = build_residual,
        .build_jacobian = build_jacobian,
        .limit_newton_step = NULL
    };
    const DcSolverOptions options = dc_solver_default_options();
    double x[DC_MAX_UNKNOWNS] = { 0.0 };
    DcNewtonReport report;

    if (!dc_newton_solve_with_report(&problem, &options, 1.0, x, &report) ||
        fabs(x[0] - 2.0) > 1e-12 || fabs(x[1] - 1.0) > 1e-12) {
        fputs("dc_solver test failed.\n", stderr);
        return 1;
    }

    SequentialContext sequential_context = { { 1.0, 1.0 } };
    const DcProblem sequential_problem = {
        .dimension = 1,
        .context = &sequential_context,
        .build_residual = build_sequential_residual,
        .build_jacobian = build_sequential_jacobian,
        .limit_newton_step = NULL,
        .independent_source_count = 2,
        .set_source_scale = set_sequential_source_scale
    };
    int sequential_total = 0;

    if (!dc_sequential_source_stepping_solve(
            &sequential_problem, &options, x, NULL, NULL, &sequential_total) ||
        fabs(x[0] - 3.0) > 1e-12 ||
        fabs(sequential_context.scale[0] - 1.0) > 1e-12 ||
        fabs(sequential_context.scale[1] - 1.0) > 1e-12) {
        fputs("sequential source stepping test failed.\n", stderr);
        return 1;
    }

    printf("dc_solver test passed in %d Newton iteration(s), %d line-search reduction(s).\n",
           report.iterations, report.line_search_reductions);
    return 0;
}
