# 通用 DC Source Stepping 框架

本目录将求解算法和具体电路模型分开。

```text
dc_solver.h / dc_solver.c
    通用 DC 求解内核：Newton-Raphson、阻尼、线搜索、稠密高斯消元。

source_stepping.h / source_stepping.c
    通用延拓外层：lambda 从 0 到 1 的自适应 Source Stepping。

netlist_dc_solver.c
    主程序：读取 R、V、Q（NPN）和 .model 网表，自动装配 MNA 方程。

examples/
    bjt1、bjt2、bjt3 的手写残差/Jacobian 示例，仅保留作学习和验证。
```

## 编译

在终端中运行：

```powershell
gcc -std=c11 -O2 -Wall -Wextra -Wpedantic `
    dc_solver.c source_stepping.c netlist_dc_solver.c `
    -lm -o netlist_dc_solver.exe
```

## 从网表直接读取电路

`netlist_dc_solver.c` 是一个简化网表 parser。当前支持 `R`、`V`、`Q`（NPN）及 `.model ... NPN IS/BF/BR`；不支持 MOS、受控源、子电路和其他复杂 SPICE 语法。

```powershell
gcc -std=c11 -O2 -Wall -Wextra -Wpedantic dc_solver.c source_stepping.c netlist_dc_solver.c -lm -o netlist_dc_solver.exe
.\netlist_dc_solver.exe ..\netlist\Schmitt3_netlist\Netlist.txt --temp 26
.\netlist_dc_solver.exe ..\netlist\Schmitt3_netlist\Netlist.txt --nr --temp 26
```

电路一、二、三都可以由同一可执行程序运行：

```powershell
.\netlist_dc_solver.exe ..\netlist\Schmitt1_netlist\Netlist.txt --temp 26
.\netlist_dc_solver.exe ..\netlist\Schmitt2_netlist\Netlist.txt --temp 26
.\netlist_dc_solver.exe ..\netlist\Schmitt3_netlist\Netlist.txt --temp 26
```

## 手写示例

`examples/` 内的专用程序不再是日常运行入口，它们展示了如何手动写出某个电路的残差和 Jacobian，并可用于验证 parser 的自动装配结果。编译电路三示例：

```powershell
gcc -std=c11 -O2 -Wall -Wextra -Wpedantic dc_solver.c source_stepping.c examples\bjt3_source_stepping.c -lm -o bjt3_source_stepping.exe
```

## 不使用 parser 的扩展方式

不需要修改 `dc_solver.c` 或 `dc_solver.h`。新建一个电路模型模块，例如：

```text
my_circuit_model.h
my_circuit_model.c
```

电路模型模块需要完成以下工作：

1. 定义电路参数结构体；
2. 指定未知量个数 `dimension`，当前上限为 `DC_MAX_UNKNOWNS`；
3. 实现残差函数 `F(x, lambda)`；
4. 实现 Jacobian 函数 `J(x, lambda)`；
5. 通过 `DcProblem` 把这两个函数和电路参数传给通用求解器。

核心接口为：

```c
typedef struct {
    int dimension;
    const void *context;
    DcResidualFunction build_residual;
    DcJacobianFunction build_jacobian;
} DcProblem;
```

其中 `context` 指向你的电路参数结构体。Source Stepping 会把同伦参数 `lambda` 传入残差和 Jacobian 函数。独立源应在电路模型中按 `lambda * source_value` 进行缩放。

调用方式：

```c
DcSolverOptions options = dc_solver_default_options();
DcProblem problem;
double x[DC_MAX_UNKNOWNS];

my_circuit_make_problem(&circuit, &problem);
dc_source_stepping_solve(&problem, &options, x, callback, NULL, NULL);
```

如果需要统计 Source Stepping 过程中所有成功步进点的 Newton 总迭代次数，可传入整数地址：

```c
int total_newton_iterations = 0;
dc_source_stepping_solve(
    &problem, &options, x, callback, NULL, &total_newton_iterations);
```

对于不使用 Source Stepping 的直接 Newton-Raphson：

```c
dc_newton_solve(&problem, &options, 1.0, x, &iterations);
```

当前框架仍采用稠密矩阵和固定最大维度，适合小规模电路。后续替换为稀疏矩阵时，主要修改 `dc_solver.c` 中的线性方程求解部分即可。
