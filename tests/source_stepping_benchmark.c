#define BJT1_SOURCE_STEPPING_NO_MAIN
#define BJT2_SOURCE_STEPPING_NO_MAIN
#define BJT3_SOURCE_STEPPING_NO_MAIN

#include "../examples/bjt1_source_stepping.c"
#include "../examples/bjt2_source_stepping.c"
#include "../examples/bjt3_source_stepping.c"

#include <stdio.h>

typedef struct {
    double step;
    const char *label;
} StepCase;

static int run_bjt1(double step)
{
    const Bjt1Circuit circuit = bjt1_default_circuit();
    DcProblem problem;
    DcSolverOptions options = dc_solver_default_options();
    double x[DC_MAX_UNKNOWNS];
    int total = 0;

    bjt1_make_problem(&circuit, &problem);
    options.initial_lambda_step = step;
    options.maximum_lambda_step = step;
    return dc_source_stepping_solve(&problem, &options, x, NULL, NULL,
                                    &total) ? total : -1;
}

static int run_bjt2(double step)
{
    const Bjt2Circuit circuit = bjt2_default_circuit();
    DcProblem problem;
    DcSolverOptions options = dc_solver_default_options();
    double x[DC_MAX_UNKNOWNS];
    int total = 0;

    bjt2_make_problem(&circuit, &problem);
    options.initial_lambda_step = step;
    options.maximum_lambda_step = step;
    return dc_source_stepping_solve(&problem, &options, x, NULL, NULL,
                                    &total) ? total : -1;
}

static int run_bjt3(double step)
{
    const Bjt3Circuit circuit = bjt3_default_circuit();
    DcProblem problem;
    DcSolverOptions options = dc_solver_default_options();
    double x[DC_MAX_UNKNOWNS];
    int total = 0;

    bjt3_make_problem(&circuit, &problem);
    options.initial_lambda_step = step;
    options.maximum_lambda_step = step;
    options.maximum_newton_iterations = 120;
    return dc_source_stepping_solve(&problem, &options, x, NULL, NULL,
                                    &total) ? total : -1;
}

static int run_bjt1_adaptive(double initial_step)
{
    const Bjt1Circuit circuit = bjt1_default_circuit();
    DcProblem problem;
    DcSolverOptions options = dc_solver_default_options();
    double x[DC_MAX_UNKNOWNS];
    int total = 0;

    bjt1_make_problem(&circuit, &problem);
    options.initial_lambda_step = initial_step;
    return dc_source_stepping_solve(&problem, &options, x, NULL, NULL,
                                    &total) ? total : -1;
}

int main(void)
{
    const StepCase cases[] = {
        {0.02, "0.02"}, {0.05, "0.05"}, {0.10, "0.10"},
        {0.125, "0.125"}, {0.20, "0.20"}, {0.25, "0.25"},
        {0.333333333333, "1/3"}, {0.50, "0.50"}
    };

    puts("constant source-step benchmark (-1 means the path failed)");
    puts(" step       BJT1   BJT2   BJT3");
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
        printf(" %-8s %5d  %5d  %5d\n", cases[i].label,
               run_bjt1(cases[i].step), run_bjt2(cases[i].step),
               run_bjt3(cases[i].step));
    }

    puts("\nadaptive BJT1 benchmark (legacy success rule, maximum step = 0.50)");
    puts(" initial step      total Newton iterations");
    printf(" 0.25              %d\n", run_bjt1_adaptive(0.25));
    printf(" 0.50              %d\n", run_bjt1_adaptive(0.50));
    return 0;
}
