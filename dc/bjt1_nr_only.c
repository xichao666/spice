/*
 * Direct Newton-Raphson experiment for the BJT circuit.
 *
 * Reuse the solver, transistor model, residual, and Jacobian from the
 * source-stepping program.  Renaming its main() avoids two entry points.
 */
#define main source_stepping_program_main
#include "bjt1_source_stepping.c"
#undef main

int main(void)
{
    const Circuit circuit = {
        .rc1 = 2000.0,
        .rc2 = 1000.0,
        .r3 = 10000.0,
        .re = 100.0,
        .vcc = 10.0,
        .vin = 1.50,
        .is1 = 1.0e-16,
        .is2 = 1.0e-16,
        .alpha_f = 100.0 / 101.0,
        .alpha_r = 1.0 / 2.0,
        .vt = 1.0 / 38.78
    };
    double x[N] = {0.0, 0.0, 0.0, 0.0};
    BjtValues b;
    int iterations = 0;

    puts("Direct Newton-Raphson test (lambda = 1, initial guess = 0):");

    /* No source stepping: solve the original circuit directly at lambda = 1. */
    if (!newton_solve(&circuit, 1.0, x, &iterations)) {
        fputs("Direct Newton-Raphson did not converge.\n", stderr);
        return 1;
    }

    evaluate_bjts(&circuit, x, 1.0, &b);
    printf("Direct Newton-Raphson converged in %d iterations.\n", iterations);
    puts("\nFinal DC operating point:");
    printf("V1 = %.9f V\nV2 = %.9f V\nV3 = %.9f V\nV4 = %.9f V\n",
           x[V1], x[V2], x[V3], x[V4]);
    printf("IC1 = %.9e A, IE1 = %.9e A\n", b.ic1, b.ie1);
    printf("IC2 = %.9e A, IE2 = %.9e A\n", b.ic2, b.ie2);
    return 0;
}
