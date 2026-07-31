#ifndef GMIN_STEPPING_H
#define GMIN_STEPPING_H

#include "dc_solver.h"

/* 每个成功 GMIN 步进点都可通过该回调输出或保存结果。 */
typedef void (*DcGminStepCallback)(
    double gmin,
    const double *x,
    const DcNewtonReport *newton_report,
    void *user_data);

/*
 * 从较大的 GMIN 附加电导逐步减小到零。
 * 每个成功工作点作为下一步 Newton 的初值，最终恢复原始电路。
 */
bool dc_gmin_stepping_solve(
    const DcProblem *problem,
    const DcSolverOptions *options,
    double *x,
    DcGminStepCallback step_callback,
    void *user_data,
    int *total_newton_iterations);

#endif
