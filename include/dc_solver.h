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

/*
 * 通用电路问题描述。
 * context 指向具体电路参数；两个函数指针由具体电路模型实现。
 */
typedef struct {
    int dimension;
    const void *context;
    DcResidualFunction build_residual;
    DcJacobianFunction build_jacobian;
} DcProblem;

/* Newton 和 Source Stepping 的数值控制参数。 */
typedef struct {
    int maximum_newton_iterations;
    double residual_tolerance;
    double voltage_tolerance;
    double maximum_component_step;
    double initial_lambda_step;
    double minimum_lambda_step;
    double maximum_lambda_step;
} DcSolverOptions;

/* 返回本示例使用的一组默认数值参数。 */
DcSolverOptions dc_solver_default_options(void);

/* 固定 lambda，在给定初值 x 下执行 Newton-Raphson。 */
bool dc_newton_solve(
    const DcProblem *problem,
    const DcSolverOptions *options,
    double lambda,
    double *x,
    int *iteration_count);

#endif
