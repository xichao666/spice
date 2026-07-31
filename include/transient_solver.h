#ifndef TRANSIENT_SOLVER_H
#define TRANSIENT_SOLVER_H

#include "dc_solver.h"

#include <stdbool.h>

/* 供瞬态控制器选择的积分公式阶数。 */
typedef enum {
    TRANSIENT_GEAR_1 = 1,
    TRANSIENT_GEAR_2 = 2
} TransientMethod;

/*
 * 传递给电路模型的本步离散化参数。
 * 电容、电感等动态元件据此建立当前时间点的 companion model：
 *
 * dx/dt = derivative_scale * x(n+1)
 *       + history_scale[0] * x(n)
 *       + history_scale[1] * x(n-1)
 */
/* 调用 transient_gear_solve() 时的时间推进与失败重试限制。 */
typedef struct {
    double time;
    double step;
    double previous_step;
    TransientMethod method;
    double derivative_scale;
    double history_scale[2];
} TransientStepInfo;

/* transient_gear_solve() 返回的执行统计，供调用方诊断步进行为。 */
typedef struct {
    double start_time;
    double stop_time;
    double initial_step;
    double maximum_step;
    double minimum_step;
    int maximum_gear_order;
    int maximum_retries_per_point;
    double newton_failure_shrink_factor;
} TransientSolverOptions;

typedef struct {
    int accepted_steps;
    int rejected_newton_steps;
    int gear1_steps;
    int gear2_steps;
    int total_newton_iterations;
    double final_time;
} TransientReport;

/*
 * 由电路模型实现的试算生命周期接口。
 * prepare_step 为当前时间点准备动态 stamp；commit_step 提交已接受解；
 * discard_step 丢弃失败试算的临时缓存。它们共同保证拒绝步不会污染历史。
 */
typedef bool (*TransientPrepareStepFunction)(
    void *context, const TransientStepInfo *step);
typedef void (*TransientCommitStepFunction)(
    void *context, const double *solution, double time);
typedef void (*TransientDiscardStepFunction)(void *context);
typedef bool (*TransientInitializeStateFunction)(
    void *context, const double *solution, double time);

/* 将现有 DcProblem 与动态元件的试算生命周期回调组合成瞬态问题。 */
typedef struct {
    DcProblem newton_problem;
    TransientInitializeStateFunction initialize_state;
    TransientPrepareStepFunction prepare_step;
    TransientCommitStepFunction commit_step;
    TransientDiscardStepFunction discard_step;
} TransientProblem;

/* 返回一组可复制后修改的 Gear 瞬态默认配置。 */
TransientSolverOptions transient_solver_default_options(void);

/*
 * 根据本步和前一步步长生成 Gear1/Gear2 导数系数。
 * 成功时填充 info，失败表示方法或时间步参数不合法。
 */
bool transient_gear_make_step_info(
    double time,
    double step,
    double previous_step,
    TransientMethod method,
    TransientStepInfo *info);

/*
 * 从 initial_and_final_solution 给出的初始状态推进到 stop_time。初始化回调
 * 必须将该初始状态同步到动态元件的历史缓存。
 * 成功后同一数组被覆盖为终态，report 返回接受/拒绝步统计。
 */
bool transient_gear_solve(
    const TransientProblem *problem,
    const DcSolverOptions *newton_options,
    const TransientSolverOptions *options,
    double *initial_and_final_solution,
    TransientReport *report);

#endif
