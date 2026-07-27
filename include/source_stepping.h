#ifndef SOURCE_STEPPING_H
#define SOURCE_STEPPING_H

#include "dc_solver.h"

/* 每个成功的 lambda 点都可通过该回调输出或保存结果。 */
typedef void (*DcStepCallback)(
    double lambda,
    const double *x,
    int newton_iterations,
    void *user_data);

/* 顺序 Source Stepping 中，每个成功步进点调用的回调函数。 */
typedef void (*DcSequentialStepCallback)(
    int source_index,
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

/*
 * 顺序打开已注册的独立源。
 * 已完成的源保持完整值，后续源保持关闭状态。
 */
bool dc_sequential_source_stepping_solve(
    const DcProblem *problem,
    const DcSolverOptions *options,
    double *x,
    DcSequentialStepCallback step_callback,
    void *user_data,
    int *total_newton_iterations);

#endif
