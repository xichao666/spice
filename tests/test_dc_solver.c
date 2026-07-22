#include "dc_solver.h"

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

/* 验证通用 Newton 求解器可正确求得已知的两元解。 */
int main(void)
{
    const DcProblem problem = { 2, NULL, build_residual, build_jacobian };
    const DcSolverOptions options = dc_solver_default_options();
    double x[DC_MAX_UNKNOWNS] = { 0.0 };
    int iterations = 0;

    if (!dc_newton_solve(&problem, &options, 1.0, x, &iterations) ||
        fabs(x[0] - 2.0) > 1e-12 || fabs(x[1] - 1.0) > 1e-12) {
        fputs("dc_solver test failed.\n", stderr);
        return 1;
    }

    printf("dc_solver test passed in %d Newton iteration(s).\n", iterations);
    return 0;
}
