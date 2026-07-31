#include "transient_solver.h"

#include <math.h>
#include <string.h>

/* 瞬态外层的保守默认值：Gear2 可用，但不会主动增大步长。 */
TransientSolverOptions transient_solver_default_options(void)
{
    return (TransientSolverOptions) {
        .start_time = 0.0,
        .stop_time = 1.0,
        .initial_step = 1.0e-3,
        .maximum_step = 1.0e-3,
        .minimum_step = 1.0e-12,
        .maximum_gear_order = 2,
        .maximum_retries_per_point = 12,
        .newton_failure_shrink_factor = 0.5
    };
}

bool transient_gear_make_step_info(
    double time,
    double step,
    double previous_step,
    TransientMethod method,
    TransientStepInfo *info)
{
    /* 将积分公式归一化为动态元件可直接消费的 alpha/beta 系数。 */
    if (info == NULL || !isfinite(time) || !isfinite(step) || step <= 0.0) {
        return false;
    }

    *info = (TransientStepInfo) {
        .time = time,
        .step = step,
        .previous_step = previous_step,
        .method = method
    };

    if (method == TRANSIENT_GEAR_1) {
        /* Gear1 即后向欧拉，只依赖最后一个已接受点。 */
        info->derivative_scale = 1.0 / step;
        info->history_scale[0] = -1.0 / step;
        return true;
    }

    if (method == TRANSIENT_GEAR_2 &&
        isfinite(previous_step) && previous_step > 0.0) {
        /* 变步长 BDF2：系数由相邻步长比 r = h/h_prev 决定。 */
        const double ratio = step / previous_step;
        const double denominator = 1.0 + ratio;

        if (!isfinite(ratio) || denominator <= 0.0) return false;

        info->derivative_scale = (1.0 + 2.0 * ratio) / (denominator * step);
        info->history_scale[0] = -(1.0 + ratio) / step;
        info->history_scale[1] = ratio * ratio / (denominator * step);
        return isfinite(info->derivative_scale) &&
               isfinite(info->history_scale[0]) &&
               isfinite(info->history_scale[1]);
    }

    return false;
}

static bool valid_options(const TransientSolverOptions *options)
{
    /* 在分配历史和进入时间循环前拒绝无法保证推进语义的配置。 */
    return options != NULL &&
           isfinite(options->start_time) && isfinite(options->stop_time) &&
           options->stop_time >= options->start_time &&
           isfinite(options->initial_step) && options->initial_step > 0.0 &&
           isfinite(options->maximum_step) && options->maximum_step > 0.0 &&
           isfinite(options->minimum_step) && options->minimum_step > 0.0 &&
           options->maximum_step >= options->minimum_step &&
           (options->maximum_gear_order == 1 ||
            options->maximum_gear_order == 2) &&
           options->maximum_retries_per_point >= 0 &&
           options->newton_failure_shrink_factor > 0.0 &&
           options->newton_failure_shrink_factor < 1.0;
}

bool transient_gear_solve(
    const TransientProblem *problem,
    const DcSolverOptions *newton_options,
    const TransientSolverOptions *options,
    double *initial_and_final_solution,
    TransientReport *report)
{
    /* 该函数负责时间推进、历史管理与失败回滚；Newton 本身交给 DC 内核。 */
    if (report == NULL) return false;
    *report = (TransientReport) { 0 };

    if (problem == NULL || newton_options == NULL || !valid_options(options) ||
        initial_and_final_solution == NULL ||
        problem->initialize_state == NULL || problem->prepare_step == NULL ||
        problem->commit_step == NULL ||
        problem->newton_problem.context == NULL ||
        problem->newton_problem.dimension <= 0 ||
        problem->newton_problem.dimension > DC_MAX_UNKNOWNS) {
        return false;
    }

    const int n = problem->newton_problem.dimension;
    /* 三份解分别代表已接受点、次新历史点和当前可丢弃的候选点。 */
    double accepted[DC_MAX_UNKNOWNS] = { 0.0 };
    double previous[DC_MAX_UNKNOWNS] = { 0.0 };
    double trial[DC_MAX_UNKNOWNS] = { 0.0 };
    memcpy(accepted, initial_and_final_solution, (size_t)n * sizeof(double));

    if (!problem->initialize_state((void *)problem->newton_problem.context,
                                   accepted, options->start_time)) {
        return false;
    }

    double time = options->start_time;
    report->final_time = time;
    double proposed_step = fmin(options->initial_step, options->maximum_step);
    double previous_step = 0.0;
    int accepted_history_count = 1;
    bool force_gear1 = false;

    while (time < options->stop_time) {
        /* 本实现固定采用建议步长，仅在终点截断或 Newton 失败时改变它。 */
        const double remaining = options->stop_time - time;
        double step = fmin(proposed_step, remaining);
        int retries = 0;

        for (;;) {
            /* 历史不足或失败重试时降为 Gear1；其余接受点使用 Gear2。 */
            const bool may_use_gear2 =
                !force_gear1 && options->maximum_gear_order >= 2 &&
                accepted_history_count >= 2 && previous_step > 0.0;
            const TransientMethod method = may_use_gear2
                ? TRANSIENT_GEAR_2 : TRANSIENT_GEAR_1;
            TransientStepInfo step_info;

            if (!transient_gear_make_step_info(
                    time + step, step, previous_step, method, &step_info) ||
                !problem->prepare_step(
                    (void *)problem->newton_problem.context, &step_info)) {
                /* prepare 失败也不能留下试算缓存。 */
                if (problem->discard_step != NULL) {
                    problem->discard_step((void *)problem->newton_problem.context);
                }
                return false;
            }

            if (method == TRANSIENT_GEAR_2) {
                /* 用最近两个接受解线性外推，为本步 Newton 提供初值。 */
                const double ratio = step / previous_step;
                for (int i = 0; i < n; ++i) {
                    trial[i] = accepted[i] + ratio * (accepted[i] - previous[i]);
                }
            } else {
                memcpy(trial, accepted, (size_t)n * sizeof(double));
            }

            DcNewtonReport newton_report;
            const bool converged = dc_newton_solve_with_report(
                &problem->newton_problem, newton_options, 1.0, trial,
                &newton_report);
            report->total_newton_iterations += newton_report.iterations;

            if (converged) {
                /* 只有 Newton 收敛后才推进全局历史并提交器件历史。 */
                memcpy(previous, accepted, (size_t)n * sizeof(double));
                memcpy(accepted, trial, (size_t)n * sizeof(double));
                problem->commit_step(
                    (void *)problem->newton_problem.context, accepted,
                    time + step);
                time += step;
                previous_step = step;
                proposed_step = step;
                accepted_history_count = 2;
                force_gear1 = false;
                ++report->accepted_steps;
                if (method == TRANSIENT_GEAR_1) ++report->gear1_steps;
                else ++report->gear2_steps;
                break;
            }

            /* 失败候选不进入历史；缩步后从 accepted 状态重新试算。 */
            if (problem->discard_step != NULL) {
                problem->discard_step((void *)problem->newton_problem.context);
            }
            ++report->rejected_newton_steps;
            ++retries;
            step *= options->newton_failure_shrink_factor;
            force_gear1 = true;

            if (retries > options->maximum_retries_per_point ||
                step < options->minimum_step) {
                return false;
            }
        }
    }

    memcpy(initial_and_final_solution, accepted, (size_t)n * sizeof(double));
    report->final_time = time;
    return true;
}
