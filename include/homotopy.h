#ifndef HOMOTOPY_H
#define HOMOTOPY_H

#include "dc_solver.h"

/*固定点同伦： H(x, lambda) = (1-lambda) G (x-a) + lambda F(x).
 * 令 G 为对角阵；diagonal_scale 为 NULL 时，所有对角元都取 uniform_scale。
 * lambda=0 时 x=a 是精确已知起点；lambda=1 时 H=F，因此路径与 lambda=1 平面的交点是原 DC 工作点。
 */
typedef struct {
    const double *starting_point;       /* 固定点 a，长度为问题维度。 */
    const double *diagonal_scale;       /* 可选的 G 对角元，长度为问题维度。 */
    double uniform_scale;                /* diagonal_scale 为 NULL 时的 Gii。 */
    double initial_arc_step;             /* 初始弧长步长 h。 */
    double minimum_arc_step;             /* 预测误差过大时的 h 下限。 */
    double maximum_arc_step;             /* 允许的 h 上限。 */
    double local_error_tolerance;        /* predictor-corrector 局部误差容差。 */
    double crossing_tolerance;           /* |lambda-1| 交点判定容差。 */
    int maximum_path_steps;              /* 防止无界路径跟踪。 */
} DcHomotopyOptions;

typedef void (*DcHomotopyPathCallback)(
    int step, double arc_length, double lambda, double arc_step,
    double local_error, const double *x, void *user_data);

typedef void (*DcHomotopySolutionCallback)(
    int solution_index, const double *x, double original_residual_norm,
    void *user_data);

/* 固定点同伦的路径积分与原始 Newton 复核统计。 */
typedef struct {
    int accepted_path_steps;
    int crossing_newton_calls;
    int crossing_newton_iterations;
} DcHomotopyReport;

DcHomotopyOptions dc_homotopy_default_options(const double *starting_point);

/* 使用伪弧长 predictor-corrector 跟踪固定点同伦路径。 */
bool dc_fixed_point_homotopy_solve(
    const DcProblem *problem, const DcSolverOptions *newton_options,
    const DcHomotopyOptions *homotopy_options,
    DcHomotopyPathCallback path_callback,
    DcHomotopySolutionCallback solution_callback, void *user_data,
    int *solution_count);

/* 与基础接口相同，但额外返回路径点和 lambda=1 复核统计。 */
bool dc_fixed_point_homotopy_solve_with_report(
    const DcProblem *problem, const DcSolverOptions *newton_options,
    const DcHomotopyOptions *homotopy_options,
    DcHomotopyPathCallback path_callback,
    DcHomotopySolutionCallback solution_callback, void *user_data,
    int *solution_count, DcHomotopyReport *report);

#endif
