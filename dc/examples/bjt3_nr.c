/* 复用电路三模型，仅求 lambda=1 时从零初值出发的直接 Newton 解。 */
#define BJT3_SOURCE_STEPPING_NO_MAIN
#include "bjt3_source_stepping.c"

int main(void)
{
    const Bjt3Circuit circuit = bjt3_default_circuit();
    DcSolverOptions options = dc_solver_default_options();
    DcProblem problem;
    double x[DC_MAX_UNKNOWNS] = {0.0};
    Bjt3Values values;
    int iterations = 0;

    options.maximum_newton_iterations = 120;
    bjt3_make_problem(&circuit, &problem);
    puts("Direct Newton-Raphson test for BJT3 (lambda = 1, initial guess = 0):");

    if (!dc_newton_solve(&problem, &options, 1.0, x, &iterations)) {
        fputs("Direct Newton-Raphson did not converge.\n", stderr);
        return 1;
    }

    bjt3_evaluate(&circuit, x, &values);
    printf("Direct Newton-Raphson converged in %d iterations.\n", iterations);
    puts("\nFinal DC operating point:");
    for (int i = 0; i < BJT3_NODE_COUNT; ++i) {
        printf("X%-2d = % .9f V\n", i + 1, x[i]);
    }
    return 0;
}
