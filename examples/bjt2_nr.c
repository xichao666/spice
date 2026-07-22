/* 复用 BJT2 电路模型，只执行固定 lambda=1 的直接 Newton-Raphson。 */
#define BJT2_SOURCE_STEPPING_NO_MAIN
#include "bjt2_source_stepping.c"

int main(void)
{
    const Bjt2Circuit circuit = bjt2_default_circuit();
    const DcSolverOptions options = dc_solver_default_options();
    DcProblem problem;
    double x[DC_MAX_UNKNOWNS] = {0.0};
    Bjt2Values values;
    int iterations = 0;

    bjt2_make_problem(&circuit, &problem);
    puts("Direct Newton-Raphson test for BJT2 (lambda = 1, initial guess = 0):");

    if (!dc_newton_solve(&problem, &options, 1.0, x, &iterations)) {
        fputs("Direct Newton-Raphson did not converge.\n", stderr);
        return 1;
    }

    bjt2_evaluate(&circuit, x, &values);
    printf("Direct Newton-Raphson converged in %d iterations.\n", iterations);
    puts("\nFinal DC operating point:");
    printf("X1 = %.9f V\nX2 = %.9f V\nX3 = %.9f V\nX4 = %.9f V\nX5 = %.9f V\n",
           x[BJT2_X1], x[BJT2_X2], x[BJT2_X3], x[BJT2_X4], x[BJT2_X5]);
    printf("IC1 = %.9e A, IE1 = %.9e A\n", values.ic1, values.ie1);
    printf("IC2 = %.9e A, IE2 = %.9e A\n", values.ic2, values.ie2);
    return 0;
}
