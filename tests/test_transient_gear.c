#include "transient_solver.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

/* 一状态 RC 阶跃模型：同时作为电容 Gear stamp 的集成测试夹具。 */
typedef struct {
    double resistance;
    double capacitance;
    double source_voltage;
    double accepted_voltage[2];
    TransientStepInfo step;
    int discard_count;
    double reject_steps_larger_than;
} RcCircuit;

static void rc_residual(const void *context, const double *x,
                        double lambda, double *residual)
{
    /* KCL = 电阻支路电流 + 按当前 Gear 系数离散的电容电流。 */
    const RcCircuit *circuit = context;
    (void)lambda;
    /* 测试专用开关：以 NaN 模拟过大时间步下的 Newton 失败。 */
    if (circuit->step.step > circuit->reject_steps_larger_than) {
        residual[0] = NAN;
        return;
    }
    residual[0] = (x[0] - circuit->source_voltage) / circuit->resistance +
        circuit->capacitance * (circuit->step.derivative_scale * x[0] +
            circuit->step.history_scale[0] * circuit->accepted_voltage[0] +
            circuit->step.history_scale[1] * circuit->accepted_voltage[1]);
}

static void rc_jacobian(const void *context, const double *x, double lambda,
                        double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{
    const RcCircuit *circuit = context;
    (void)x;
    (void)lambda;
    memset(jacobian, 0, sizeof(double) * DC_MAX_UNKNOWNS * DC_MAX_UNKNOWNS);
    jacobian[0][0] = 1.0 / circuit->resistance +
        circuit->capacitance * circuit->step.derivative_scale;
}

static bool rc_prepare(void *context, const TransientStepInfo *step)
{
    RcCircuit *circuit = context;
    circuit->step = *step;
    return true;
}

static bool rc_initialize(void *context, const double *x, double time)
{
    RcCircuit *circuit = context;
    (void)time;
    circuit->accepted_voltage[0] = x[0];
    circuit->accepted_voltage[1] = x[0];
    return true;
}

static void rc_commit(void *context, const double *x, double time)
{
    /* 仅在接受点滚动电容电压历史。 */
    RcCircuit *circuit = context;
    (void)time;
    circuit->accepted_voltage[1] = circuit->accepted_voltage[0];
    circuit->accepted_voltage[0] = x[0];
}

static void rc_discard(void *context)
{
    RcCircuit *circuit = context;
    ++circuit->discard_count;
}

/* 一状态 RL 阶跃模型：验证电感电流历史项的离散化。 */
typedef struct {
    double resistance;
    double inductance;
    double source_voltage;
    double accepted_current[2];
    TransientStepInfo step;
} RlCircuit;

static void rl_residual(const void *context, const double *x,
                        double lambda, double *residual)
{
    const RlCircuit *circuit = context;
    (void)lambda;
    residual[0] = circuit->inductance *
        (circuit->step.derivative_scale * x[0] +
         circuit->step.history_scale[0] * circuit->accepted_current[0] +
         circuit->step.history_scale[1] * circuit->accepted_current[1]) +
        circuit->resistance * x[0] - circuit->source_voltage;
}

static void rl_jacobian(const void *context, const double *x, double lambda,
                        double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{
    const RlCircuit *circuit = context;
    (void)x;
    (void)lambda;
    memset(jacobian, 0, sizeof(double) * DC_MAX_UNKNOWNS * DC_MAX_UNKNOWNS);
    jacobian[0][0] = circuit->inductance * circuit->step.derivative_scale +
        circuit->resistance;
}

static bool rl_prepare(void *context, const TransientStepInfo *step)
{
    ((RlCircuit *)context)->step = *step;
    return true;
}

static bool rl_initialize(void *context, const double *x, double time)
{
    RlCircuit *circuit = context;
    (void)time;
    circuit->accepted_current[0] = x[0];
    circuit->accepted_current[1] = x[0];
    return true;
}

static void rl_commit(void *context, const double *x, double time)
{
    RlCircuit *circuit = context;
    (void)time;
    circuit->accepted_current[1] = circuit->accepted_current[0];
    circuit->accepted_current[0] = x[0];
}

static TransientProblem rc_problem(RcCircuit *circuit)
{
    return (TransientProblem) {
        .newton_problem = {
            .dimension = 1,
            .context = circuit,
            .build_residual = rc_residual,
            .build_jacobian = rc_jacobian
        },
        .initialize_state = rc_initialize,
        .prepare_step = rc_prepare,
        .commit_step = rc_commit,
        .discard_step = rc_discard
    };
}

static TransientProblem rl_problem(RlCircuit *circuit)
{
    return (TransientProblem) {
        .newton_problem = {
            .dimension = 1,
            .context = circuit,
            .build_residual = rl_residual,
            .build_jacobian = rl_jacobian
        },
        .initialize_state = rl_initialize,
        .prepare_step = rl_prepare,
        .commit_step = rl_commit
    };
}

static bool test_coefficients(void)
{
    /* 同时锁定等步长和 r = 0.5 的变步长 BDF2 系数。 */
    TransientStepInfo info;
    if (!transient_gear_make_step_info(1.0, 0.2, 0.2,
                                       TRANSIENT_GEAR_2, &info) ||
        fabs(info.derivative_scale - 7.5) > 1e-14 ||
        fabs(info.history_scale[0] + 10.0) > 1e-14 ||
        fabs(info.history_scale[1] - 2.5) > 1e-14) return false;

    if (!transient_gear_make_step_info(1.0, 0.1, 0.2,
                                       TRANSIENT_GEAR_2, &info) ||
        fabs(info.derivative_scale - 13.333333333333334) > 1e-12 ||
        fabs(info.history_scale[0] + 15.0) > 1e-14 ||
        fabs(info.history_scale[1] - 1.6666666666666667) > 1e-12) return false;
    return true;
}

static bool test_rc_gear2_is_more_accurate(const DcSolverOptions *newton)
{
    /* 在相同步长下，二阶 Gear2 应优于一阶 Gear1。 */
    const double final_time = 1.0e-3;
    const double exact = 1.0 - exp(-10.0);
    double gear1_solution[DC_MAX_UNKNOWNS] = { 0.0 };
    double gear2_solution[DC_MAX_UNKNOWNS] = { 0.0 };
    RcCircuit gear1 = { .resistance = 1000.0, .capacitance = 1.0e-7,
                        .source_voltage = 1.0, .reject_steps_larger_than = INFINITY };
    RcCircuit gear2 = gear1;
    TransientProblem gear1_problem = rc_problem(&gear1);
    TransientProblem gear2_problem = rc_problem(&gear2);
    TransientSolverOptions options = transient_solver_default_options();
    options.stop_time = final_time;
    options.initial_step = 2.0e-5;
    options.maximum_step = 2.0e-5;
    options.minimum_step = 1.0e-9;
    options.maximum_gear_order = 1;
    TransientReport report;

    if (!transient_gear_solve(&gear1_problem, newton, &options,
                              gear1_solution, &report)) return false;
    options.maximum_gear_order = 2;
    if (!transient_gear_solve(&gear2_problem, newton, &options,
                              gear2_solution, &report)) return false;

    return fabs(gear2_solution[0] - exact) < fabs(gear1_solution[0] - exact) &&
           fabs(gear2_solution[0] - exact) < 2.0e-4;
}

static bool test_rl(const DcSolverOptions *newton)
{
    /* 与 RL 阶跃解析解比较，并确认起步后确实进入 Gear2。 */
    RlCircuit circuit = { .resistance = 10.0, .inductance = 1.0e-3,
                          .source_voltage = 2.0 };
    double solution[DC_MAX_UNKNOWNS] = { 0.0 };
    TransientSolverOptions options = transient_solver_default_options();
    options.stop_time = 1.0e-3;
    options.initial_step = 1.0e-5;
    options.maximum_step = 1.0e-5;
    options.minimum_step = 1.0e-9;
    TransientReport report;
    const double exact = 0.2 * (1.0 - exp(-10.0));
    TransientProblem problem = rl_problem(&circuit);

    return transient_gear_solve(&problem, newton, &options,
                                solution, &report) &&
           fabs(solution[0] - exact) < 2.0e-4 && report.gear2_steps > 0;
}

static bool test_nonzero_initial_condition(const DcSolverOptions *newton)
{
    /* 初始电容电压必须在第一步被写入 companion model 的历史项。 */
    RcCircuit circuit = { .resistance = 1000.0, .capacitance = 1.0e-7,
                          .source_voltage = 0.0,
                          .reject_steps_larger_than = INFINITY };
    double solution[DC_MAX_UNKNOWNS] = { 1.0 };
    TransientSolverOptions options = transient_solver_default_options();
    options.stop_time = 1.0e-5;
    options.initial_step = 1.0e-5;
    options.maximum_step = 1.0e-5;
    options.maximum_gear_order = 1;
    TransientReport report;
    TransientProblem problem = rc_problem(&circuit);

    return transient_gear_solve(&problem, newton, &options,
                                solution, &report) &&
           solution[0] > 0.8 && solution[0] < 1.0;
}

static bool test_newton_rollback(const DcSolverOptions *newton)
{
    /* 验证拒绝步调用 discard，且提交历史仍与最终接受解一致。 */
    RcCircuit circuit = { .resistance = 1000.0, .capacitance = 1.0e-7,
                          .source_voltage = 1.0,
                          .reject_steps_larger_than = 5.0e-5 };
    double solution[DC_MAX_UNKNOWNS] = { 0.0 };
    TransientSolverOptions options = transient_solver_default_options();
    options.stop_time = 2.0e-4;
    options.initial_step = 1.0e-4;
    options.maximum_step = 1.0e-4;
    options.minimum_step = 1.0e-8;
    TransientReport report;
    TransientProblem problem = rc_problem(&circuit);

    return transient_gear_solve(&problem, newton, &options,
                                solution, &report) &&
           report.rejected_newton_steps == 1 && circuit.discard_count == 1 &&
           fabs(circuit.accepted_voltage[0] - solution[0]) < 1e-14;
}

int main(void)
{
    DcSolverOptions newton = dc_solver_default_options();
    newton.maximum_component_step = 1.0;

    if (!test_coefficients() || !test_rc_gear2_is_more_accurate(&newton) ||
        !test_rl(&newton) || !test_nonzero_initial_condition(&newton) ||
        !test_newton_rollback(&newton)) {
        fputs("transient Gear test failed.\n", stderr);
        return 1;
    }
    puts("transient Gear tests passed.");
    return 0;
}
