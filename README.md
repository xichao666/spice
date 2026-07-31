# SPICE DC 工作点求解学习项目

这是一个使用 C11 编写的 SPICE 风格 DC 工作点求解学习项目。项目围绕模拟电路的 DC 工作点分析，逐步实现并验证改进节点分析（MNA）、Newton-Raphson、Source Stepping、BJT Ebers-Moll 模型、网表解析与同伦延拓。

它是一个面向学习和实验的小型求解器，不是 ngspice、PSpice 等完整工业级 SPICE 的替代品。

## 当前阶段说明

已发布的版本标签为 `v1`、`v2` 和当前准备发布的 `v3`。下面的阶段名称用于说明代码能力的演进。

### v3：网表与收敛算法升级

- 同时 Source Stepping、顺序 Source Stepping 和 GMIN Stepping；
- 每次实验输出总 Newton 次数、每步迭代次数、线搜索回退、残差与运行时间；
- 支持 NPN / PNP BJT、独立电压源 `V`、独立电流源 `I`、DC 电容开路；
- 支持 `.subckt` / `X` 的层次展开，内部节点使用实例路径隔离；
- 支持数字节点、文字节点、BJT 实例参数和常见单位后缀；

### 基础教学阶段

- 针对 Schmitt1、Schmitt2、Schmitt3 手写残差与 Jacobian；
- 使用四变量 BJT 电路验证直接 Newton-Raphson；
- 实现最初的同时 Source Stepping；
- 对照 MATLAB / PSpice 工作点学习多解与初值敏感性。

### 通用 DC 求解阶段

- 抽象出通用 `DcProblem` 接口；
- 分离通用 Newton 求解器与具体电路模型；
- 增加简化网表 Parser 和自动 MNA 方程构建；
- 支持 NPN / PNP Ebers-Moll BJT、温度设置、步长控制和求解统计。

### 当前扩展阶段

- 同时 Source Stepping 的自适应步长、失败回退和受残差保护的 secant predictor；
- 顺序 Source Stepping，可按指定顺序逐个开启独立电压源；
- 固定点同伦与伪弧长 predictor-corrector 的独立实验模块；
- 保存同伦路径数据，并提供绘图脚本和学习笔记。

## 当前支持范围

### 网表器件与模型

| 项目 | 当前支持情况 |
| --- | --- |
| 电阻 `R` | 支持。 |
| 独立电压源 `V` | 支持；可读取普通 DC 数值以及 `DC`、`SIN`、`PULSE`、`PWL` 形式的 DC 初值。 |
| BJT `Q` | 支持 NPN / PNP Ebers-Moll 模型。 |
| `.model NPN/PNP` | 支持当前子集：`IS`、`BF`、`BR`。 |
| BJT 实例参数 | 支持 `area`。 |
| `.model` 参数续行 | 支持以 `+` 开头的模型参数续行。 |
| 独立电流源 `I` | 支持；与电压源一起参与 Source Stepping。 |
| 理想电容 `C` | 支持 DC 开路处理；不参与 DC MNA 方程。 |
| `.subckt` / `X` | 支持当前文件内的层次展开；内部节点按实例路径命名。 |
| 受控源、MOSFET、二极管、电感 | 当前网表 DC 求解器未实现。 |
| `.TRAN`、`.AC`、`.DC` | 未实现；当前项目只进行 DC 工作点分析。 |
| `.include` / `.lib`、参数化网表 | 未实现。 |

不要把“能够解析一部分 SPICE 风格网表”理解为“兼容完整 SPICE 网表语法或全部器件模型”。

## 求解流程

### 线性与非线性电路的共同框架

```text
网表 / 手写电路模型
        ↓
残差 F(x, lambda) 与 Jacobian J(x)
        ↓
Newton-Raphson 内层迭代
        ↓
直接 NR / 同时 Source Stepping / 顺序 Source Stepping
        ↓
DC 工作点、残差与 Newton 统计信息
```

对固定的 `lambda`，Newton-Raphson 求解：

```text
J(xk) * delta_x = -F(xk)
x(k+1) = x(k) + delta_x
```

通用 Newton 求解器包含：

- 稠密高斯消元与部分主元选取；
- 单个节点电压变化限制；
- BJT 结电压后备限制；
- 回溯线搜索；
- 残差和电压更新量双重收敛判断。

### 三种 DC 求解模式

| 模式 | 启用方式 | 主要用途 | 特点 |
| --- | --- | --- | --- |
| 直接 Newton-Raphson | `--nr` | 初值较好或用于比较 | 速度可能较快，但对初值敏感。 |
| 同时 Source Stepping | 默认 | 常规 DC 工作点求解 | 所有独立源共用一个 `lambda` 同时缩放。 |
| 顺序 Source Stepping | `--sequential-sources` | 同时步进失败、研究上电路径或多解 | 各独立源按指定顺序逐个开启，通常更稳健但不一定更快。 |
| GMIN Stepping | `--gmin-stepping` | 减弱浮空节点、病态矩阵或强非线性引起的收敛困难 | 为每个非地节点并联临时对地电导，再逐步撤除。 |

### 同时 Source Stepping

所有独立电压源同时按比例开启：

```text
Vk(lambda) = lambda * Vk
lambda: 0 -> 1
```

电路从“所有源关闭”逐步变为完整电路。每个成功工作点作为下一个 `lambda` 点的 Newton 初值。

### 顺序 Source Stepping

顺序模式为每个独立源维护独立缩放比例。例如四个源可按如下路径开启：

```text
(0, 0, 0, 0)
    -> (1, 0, 0, 0)
    -> (1, 1, 0, 0)
    -> (1, 1, 1, 0)
    -> (1, 1, 1, 1)
```

不同源顺序可能得到不同 DC 工作点；对于多稳态或强反馈电路，某些顺序也可能无法继续追踪到最终状态。这是延拓路径相关性，不代表原电路一定没有 DC 解。

### 固定点同伦与伪弧长实验

`src/homotopy.c` 与相关测试程序用于学习固定点同伦和伪弧长 predictor-corrector。它与默认网表入口的 Source Stepping 是相互独立的实验模块。

伪弧长方法将 `lambda` 和电路未知量一起看作扩展变量，并通过“同伦方程 + 切平面约束”进行校正。与仅按 `lambda` 单调增加的 Source Stepping 相比，它可以处理解曲线的折返，但实现和参数控制更复杂。

## 目录说明

```text
spice/
├── include/                 # 公开接口
│   ├── dc_solver.h          # DcProblem、Newton 选项和求解接口
│   ├── source_stepping.h    # 同时与顺序 Source Stepping 接口
│   ├── gmin_stepping.h      # GMIN Stepping 接口
│   └── homotopy.h           # 固定点同伦接口
├── src/                     # 求解器、Parser 与方程实现
│   ├── dc_solver.c
│   ├── source_stepping.c
│   ├── gmin_stepping.c
│   ├── netlist_dc_solver.c
│   └── homotopy.c
├── examples/                # BJT1、BJT2、BJT3 的手写方程示例
├── tests/                   # Newton、Source Stepping、同伦测试与实验
├── netlist/                 # Schmitt1、Schmitt2、Schmitt3 网表和 PSpice 文件
├── learning/                # 学习笔记、代码讲解、同伦路径数据和绘图工具
├── summary/                 # 阶段性总结
├── CMakeLists.txt
└── README.md
```

## 构建

### Windows：项目自带 w64devkit 工具链

在 PowerShell 中进入项目根目录后，先临时加入工具链：

```powershell
$env:Path = "$PWD\.toolchain\w64devkit\w64devkit\bin;" + $env:Path
```

确认环境：

```powershell
gcc --version
cmake --version
```

使用 CMake 构建：

```powershell
cmake -S . -B build -G "MinGW Makefiles"
cmake --build build
ctest --test-dir build --output-on-failure
```

旧的 `examples/` 手写 BJT 程序仍保留作学习资料，但尚未全部迁移到当前
`DcProblem` 回调接口，因此默认不参与构建。当前推荐从通用网表入口运行实验。

主要生成程序包括：

```text
build/netlist_dc_solver.exe
build/bjt1_source_stepping.exe
build/bjt1_nr.exe
build/bjt2_source_stepping.exe
build/bjt2_nr.exe
build/bjt3_source_stepping.exe
build/bjt3_nr.exe
```

### 手动编译网表程序

```powershell
$env:Path = "$PWD\.toolchain\w64devkit\w64devkit\bin;" + $env:Path

gcc -std=c11 -O2 -Wall -Wextra -Wpedantic -Iinclude `
    src\dc_solver.c src\source_stepping.c src\netlist_dc_solver.c `
    -lm -o netlist_dc_solver.exe
```

## 运行

### 默认：同时 Source Stepping

```powershell
.\build\netlist_dc_solver.exe .\netlist\Schmitt1_netlist\Netlist.txt --temp 26
```

不带 `--nr` 和 `--sequential-sources` 时，程序默认使用同时 Source Stepping。

### 直接 Newton-Raphson

```powershell
.\build\netlist_dc_solver.exe .\netlist\Schmitt1_netlist\Netlist.txt --nr --temp 26
```

### GMIN Stepping

```powershell
.\build\netlist_dc_solver.exe .\netlist\Schmitt1_netlist\Netlist.txt --temp 26 `
    --gmin-stepping
```

该模式保持原始独立源全开，从较大的 `gmin` 开始，向每个非地节点写入一条临时的对地电导。每个成功工作点作为下一档更小 `gmin` 的 Newton 初值，最后令 `gmin=0` 并再次求解，因此输出的是原始电路而不是被修改电路的工作点。

### 输出可比较的实验数据

为单次运行加上 `--report-dir`，程序会在已存在的目录写入：

- `dc_steps.csv`：每个成功步进点的参数、Newton 次数、线搜索回退次数、残差、该点耗时与全部未知量；
- `dc_summary.csv`：总 Newton 次数、成功点所用 Newton 次数、失败探测所用 Newton 次数、总耗时、最终残差和本次参数。

```powershell
New-Item -ItemType Directory -Force .\result\schmitt1-gmin
.\build\netlist_dc_solver.exe .\netlist\Schmitt1_netlist\Netlist.txt `
    --temp 26 --gmin-stepping --report-dir .\result\schmitt1-gmin
```

若要对 `netlist/` 下的全部 `Netlist.txt` 批量比较直接 Newton、同时/顺序 Source Stepping 和 GMIN Stepping，可执行：

```powershell
powershell -ExecutionPolicy Bypass -File .\tools\run_dc_experiments.ps1
```

它会生成 `result/dc_experiments/dc_algorithm_comparison.csv`。顺序 Source Stepping 的默认源顺序若在某个多稳态电路上失败，汇总表会保留对应的失败状态；可再单独指定 `--source-order` 做路径实验。

可使用标准库 Python 将该批量 CSV 生成两张 SVG 图：

```powershell
python .\tools\plot_dc_experiments.py .\result\dc_experiments
```

- `dc_algorithm_comparison.svg`：各电路、各算法的总 Newton 次数与运行时间；红色柱表示该路径未收敛。
- `dc_step_iterations.svg`：每个电路中，各延拓算法在每个成功步的 Newton 次数变化。

固定点同伦/伪弧长实验也可输出同类型的路径与汇总数据：

```powershell
New-Item -ItemType Directory -Force .\result\homotopy
.\build\reference_homotopy_experiments.exe .\result\homotopy
```

其中每个电路会生成 `*_path.csv`、`*_solutions.csv` 和 `*_summary.csv`。路径表包含弧长、`lambda`、弧长步长、局部误差和逐步耗时；汇总表包含接受路径点数、解数量、总耗时和 `lambda` 覆盖范围。

### 直接 Newton-Raphson

```powershell
.\build\netlist_dc_solver.exe .\netlist\Schmitt1_netlist\Netlist.txt --nr --temp 26
```

### 同时 Source Stepping + secant predictor

```powershell
.\build\netlist_dc_solver.exe .\netlist\Schmitt2_netlist\Netlist.txt --temp 26 `
    --secant-predictor
```

`--secant-predictor` 利用前两个成功步进点预测下一个初值。程序会先比较预测点与普通继承初值的残差，只有预测明显更好时才采用它。

### 顺序 Source Stepping

```powershell
.\build\netlist_dc_solver.exe .\netlist\Schmitt1_netlist\Netlist.txt --temp 26 `
    --sequential-sources --secant-predictor
```

### 指定电源开启顺序

`--source-order` 必须列出网表中全部独立电压源，每个源名称只出现一次：

```powershell
.\build\netlist_dc_solver.exe .\netlist\Schmitt3_netlist\Netlist.txt --temp 26 `
    --sequential-sources --source-order VCC1,VCC2,V2,V1 `
    --secant-predictor
```

该选项只改变电源开启顺序，不改变电路的连接关系。

## Source Stepping 参数

| 选项 | 含义 |
| --- | --- |
| `--initial-step value` | 初始 `lambda` 步长。 |
| `--max-step value` | `lambda` 步长上限。 |
| `--fast-threshold count` | Newton 次数不超过该值时，下一步增大步长。 |
| `--slow-threshold count` | Newton 次数达到该值时，下一步减小步长。 |
| `--growth-factor value` | 快速收敛时的步长放大倍率。 |
| `--shrink-factor value` | 慢收敛或失败时的步长缩小倍率。 |
| `--junction-step value` | BJT 结电压修正上限；`0` 表示关闭该后备保护。 |
| `--secant-predictor` | 启用受残差保护的割线预测初值。 |
| `--sequential-sources` | 启用顺序 Source Stepping。 |
| `--source-order names` | 指定顺序模式的独立电压源开启顺序。 |
| `--gmin-stepping` | 启用 GMIN Stepping。 |
| `--gmin-initial siemens` | 初始附加对地电导，默认 `1e-3 S`。 |
| `--gmin-min siemens` | 最后一档非零附加电导，默认 `1e-12 S`。 |
| `--gmin-factor value` | 每次成功后 `gmin` 的缩小倍率，默认 `0.1`。 |

默认策略使用初始步长 `0.25`、最大步长 `0.50`。若某步 Newton 迭代次数不超过 `5`，下一步将扩大；若达到 `15` 次或本步失败，下一步将缩小。失败点不会被接受，程序会从上一个成功工作点重新尝试。

## 测试与验证

构建后可运行：

```powershell
ctest --test-dir build --output-on-failure
```

测试内容包括：

- 通用 Newton 与稠密线性方程求解；
- 同时 Source Stepping 基础流程；
- 顺序 Source Stepping 的多源缩放过程；
- 固定点同伦与参考同伦实验。

已完成的电路实验包括：

- Schmitt1、Schmitt2、Schmitt3 的手写方程求解；
- 三个电路通过网表入口进行 DC 求解；
- 与 MATLAB / PSpice 工作点进行过对照；
- 验证直接 Newton、同时步进和顺序步进会因初值或路径不同而得到不同工作点；
- 验证 Schmitt3 的顺序步进对电源顺序敏感，`VCC1,VCC2,V2,V1` 是一个可行顺序。



## 学习资料

- [阶段总结](summary/dc_summary.md)
- [学习笔记目录](learning/note/README.md)
- [网表 Parser 代码讲解](learning/网表Parser代码讲解.md)
- [电路三 Chua 代码逐行解析](learning/电路三-Chua代码逐行解析.md)

## 已知限制与后续方向

| 当前限制 | 后续方向 |
| --- | --- |
| 使用固定最大维度和稠密矩阵 | 引入稀疏矩阵存储和稀疏 LU 分解。 |
| 网表器件与语法范围有限 | 扩展电流源、受控源、MOSFET、二极管与子电路。 |
| BJT 为简化 Ebers-Moll 子集 | 研究更完整的模型与参数温度处理。 |
| 顺序步进需要人工指定顺序 | 研究源顺序启发式策略。 |
| Source Stepping 难以越过折返 | 继续完善伪弧长 predictor-corrector。 |
| 多解依赖当前路径发现 | 研究多路径追踪与 branch switching。 |

更详细的阶段性进展见 [summary/dc_summary.md](summary/dc_summary.md)。
