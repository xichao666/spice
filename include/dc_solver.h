#ifndef DC_SOLVER_H
#define DC_SOLVER_H

#include <stdbool.h>

/* 稠密求解器支持的最大未知量数；更大电路可在后续改为稀疏矩阵。 */
#define DC_MAX_UNKNOWNS 32

/* 电路模型必须提供的残差函数：计算 F(x, lambda)。 */
typedef void (*DcResidualFunction)(
    const void *context,
    const double *x,
    double lambda,
    double *residual);

/* 电路模型必须提供的 Jacobian 函数：计算 dF/dx。 */
typedef void (*DcJacobianFunction)(
    const void *context,
    const double *x,
    double lambda,
    double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS]);

/* Optional circuit-specific limiter for a Newton correction vector. */
typedef void (*DcStepLimiterFunction)(
    const void *context,
    const double *x,
    double *delta,
    double maximum_junction_voltage_step);

/*
 * 通用电路问题描述。
 * context 指向具体电路参数；两个函数指针由具体电路模型实现。
 */
typedef struct {
    int dimension;
    const void *context;
    DcResidualFunction build_residual;
    DcJacobianFunction build_jacobian;
    DcStepLimiterFunction limit_newton_step;
} DcProblem;

typedef enum {
    DC_SOURCE_STEP_LEGACY,
    DC_SOURCE_STEP_SECANT_PREDICTOR
} DcSourceStepPolicy;

/* Newton 和 Source Stepping 的数值控制参数。 */
typedef struct {
    int maximum_newton_iterations;
    double residual_tolerance;
    double voltage_tolerance;
    double maximum_component_step;
    double maximum_junction_voltage_step;
    double initial_lambda_step;
    double minimum_lambda_step;
    double maximum_lambda_step;
    int fast_newton_iteration_threshold;
    int slow_newton_iteration_threshold;
    double lambda_step_growth_factor;
    double lambda_step_shrink_factor;
    DcSourceStepPolicy source_step_policy;
} DcSolverOptions;

/* Newton 求解过程的迭代和线搜索统计信息。 */
typedef struct {
    int iterations;
    int line_search_reductions;
    double final_residual_norm;
} DcNewtonReport;

/* 返回本示例使用的一组默认数值参数。 */
DcSolverOptions dc_solver_default_options(void);

/* 固定 lambda，在给定初值 x 下执行 Newton-Raphson。 */
bool dc_newton_solve(
    const DcProblem *problem,
    const DcSolverOptions *options,
    double lambda,
    double *x,
    int *iteration_count);

/* 固定 lambda 执行 Newton-Raphson，并返回完整的求解统计信息。 */
bool dc_newton_solve_with_report(
    const DcProblem *problem,
    const DcSolverOptions *options,
    double lambda,
    double *x,
    DcNewtonReport *report);

#endif
