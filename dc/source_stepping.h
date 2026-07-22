#ifndef SOURCE_STEPPING_H
#define SOURCE_STEPPING_H

#include "dc_solver.h"

/* 每个成功的 lambda 点都可通过该回调输出或保存结果。 */
typedef void (*DcStepCallback)(
    double lambda,
    const double *x,
    int newton_iterations,
    void *user_data);

/* lambda 从 0 步进到 1；每一步以之前的解作为 Newton 初值。 */
bool dc_source_stepping_solve(
    const DcProblem *problem,
    const DcSolverOptions *options,
    double *x,
    DcStepCallback step_callback,
    void *user_data,
    int *total_newton_iterations);

#endif
