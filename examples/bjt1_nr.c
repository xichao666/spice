/*
 * 直接 Newton-Raphson 实验复用 BJT1 电路模型，
 * 但通过宏屏蔽 bjt1_source_stepping.c 中的 Source Stepping 主程序。
 */
#define BJT1_SOURCE_STEPPING_NO_MAIN
#include "bjt1_source_stepping.c"

int main(void)
{
    const Bjt1Circuit circuit = bjt1_default_circuit();
    const DcSolverOptions options = dc_solver_default_options();
    DcProblem problem;
    double x[DC_MAX_UNKNOWNS] = {0.0};
    Bjt1Values values;
    int iterations = 0;

    bjt1_make_problem(&circuit, &problem);
    puts("Direct Newton-Raphson test (lambda = 1, initial guess = 0):");

    /* 固定 lambda=1，从零初值直接求解原电路。 */
    if (!dc_newton_solve(&problem, &options, 1.0, x, &iterations)) {
        fputs("Direct Newton-Raphson did not converge.\n", stderr);
        return 1;
    }

    bjt1_evaluate(&circuit, x, 1.0, &values);
    printf("Direct Newton-Raphson converged in %d iterations.\n", iterations);
    puts("\nFinal DC operating point:");
    printf("V1 = %.9f V\nV2 = %.9f V\nV3 = %.9f V\nV4 = %.9f V\n",
           x[BJT1_V1], x[BJT1_V2], x[BJT1_V3], x[BJT1_V4]);
    printf("IC1 = %.9e A, IE1 = %.9e A\n", values.ic1, values.ie1);
    printf("IC2 = %.9e A, IE2 = %.9e A\n", values.ic2, values.ie2);
    return 0;
}
