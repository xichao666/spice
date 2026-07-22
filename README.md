# SPICE DC 工作点求解学习项目

这是一个用 C11 编写的电路 DC 工作点求解学习项目。当前支持：

- 直接 Newton-Raphson；
- Source Stepping；
- R、V、NPN BJT 的简化网表读取；
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

## 当前限制与后续方向

当前实现使用固定最大维度和稠密矩阵，适合小规模学习电路。后续可将 `src/dc_solver.c` 的稠密线性方程求解替换为稀疏 LU，而网表读取、电路方程构建和 Source Stepping 的外层结构可以保持不变。
