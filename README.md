# SPICE DC 工作点求解学习项目

这是一个用 C 编写的电路 DC 工作点求解学习项目。当前支持：

- 直接 Newton-Raphson；
- Source Stepping；
- R、V、NPN/PNP BJT 的简化网表读取；
- 改进节点分析（MNA）；
- Ebers-Moll BJT 模型；
- 三个施密特触发器/Chua 电路示例。

## 目录结构

```text
spice/
├── include/       # 公共接口：DcProblem、Newton 与 Source Stepping
├── src/           # 求解器实现、网表读取与自动 MNA 方程构建
├── examples/      # BJT1、BJT2、BJT3 的手写方程验证程序
├── tests/         # 求解器基础测试
├── netlist/       # Schmitt1、2、3 的网表和 PSpice 工程文件
├── learning/      # 学习笔记
├── CMakeLists.txt # CMake 构建配置
└── README.md
```

## 模块职责

| 文件或目录 | 职责 |
| --- | --- |
| `include/dc_solver.h` | 定义通用 DC 问题、Newton 选项和求解接口。 |
| `include/source_stepping.h` | 定义 Source Stepping 接口。 |
| `src/dc_solver.c` | 稠密高斯消元、阻尼 Newton 和线搜索。 |
| `src/source_stepping.c` | 自适应 `lambda: 0 → 1` 源步进。 |
| `src/netlist_dc_solver.c` | 读取网表，自动构造 MNA 残差与 Jacobian。 |
| `examples/` | 用手写残差/Jacobian 验证通用求解器。 |
| `tests/test_dc_solver.c` | 验证 Newton 与线性方程求解的基本功能。 |

Newton 的扩展接口 `dc_newton_solve_with_report()` 会返回 `DcNewtonReport`：Newton 迭代次数、线搜索回退次数和最终残差范数。默认 Source Stepping 使用规则式步长控制：Newton 次数不超过 5 次时将步长乘 1.5，达到 15 次时将步长减半。可选的 `--secant-predictor` 会利用两个已收敛工作点预测下一初值；预测残差未显著降低时自动退回默认初值。

## 使用 CMake 构建

安装 CMake 后，在项目根目录执行：

```powershell
cmake -S . -B build -G "MinGW Makefiles"
cmake --build build
ctest --test-dir build --output-on-failure
```

网表程序的运行示例：

```powershell
.\build\netlist_dc_solver.exe .\netlist\Schmitt1_netlist\Netlist.txt --temp 26
.\build\netlist_dc_solver.exe .\netlist\Schmitt3_netlist\Netlist.txt --nr --temp 26
```

## 不使用 CMake 时的手动编译

当前本地 GCC 可用时，可以在项目根目录执行：

```powershell
gcc -std=c11 -O2 -Wall -Wextra -Wpedantic -Iinclude `
    src\dc_solver.c src\source_stepping.c src\netlist_dc_solver.c `
    -lm -o netlist_dc_solver.exe
```

运行：

```powershell
.\netlist_dc_solver.exe .\netlist\Schmitt2_netlist\Netlist.txt --temp 26
```

## Source Stepping 参数试验

通用网表入口支持在不修改源代码的情况下调整步长策略：

```powershell
.\netlist_dc_solver.exe .\netlist\Schmitt1_netlist\Netlist.txt --temp 26 `
    --initial-step 0.50 --max-step 0.50 `
    --slow-threshold 18 --shrink-factor 0.50
```

可调整的参数包括：

| 选项 | 含义 |
|---|---|
| `--initial-step` | 初始 λ 步长。 |
| `--max-step` | λ 步长上限。 |
| `--fast-threshold` | Newton 次数不超过该值时，步长增大。 |
| `--slow-threshold` | Newton 次数达到该值时，步长缩小。 |
| `--growth-factor` | 快速收敛时的步长放大倍率。 |
| `--shrink-factor` | 慢收敛或失败时的步长缩小倍率。 |
| `--junction-step` | BJT 结电压变化上限；`0` 表示关闭该后备保护。 |

默认参数保持偏稳健的 `initial=0.25`、`max=0.50`、快速阈值 `5`、慢速阈值 `15`。在 Schmitt1、Schmitt2、Schmitt3 三个 BJT 网表上的第一轮测试中，`initial=0.50`、`slow-threshold=18`、`shrink-factor=0.50` 分别得到 59、59、66 次总 Newton 迭代；因此它可作为这组三个电路的可选实验配置。它尚未替换通用默认值，因为还需要在更多类型的网表上验证稳定性。

对于网表中的 BJT，求解器还提供了结电压后备限制。普通的阻尼 Newton 与线搜索先照常执行；只有它们无法降低残差时，程序才会限制 B-E 和 B-C 结电压的修正量并重新尝试。默认结电压上限为 `0.10 V`。这种设计避免在正常收敛时额外增加迭代次数，同时为强指数非线性或困难初值提供保护。

## 阶段总结

项目当前的阶段性完成内容、求解器能力、性能实验和后续计划见：[summary/dc_summary.md](summary/dc_summary.md)。

## 当前限制与后续方向

当前实现使用固定最大维度和稠密矩阵，适合小规模学习电路。后续可将 `src/dc_solver.c` 的稠密线性方程求解替换为稀疏 LU，而网表读取、电路方程构建和 Source Stepping 的外层结构可以保持不变。
