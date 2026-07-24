# DC 求解器项目阶段总结

## 1. 项目目标

本项目的目标是实现一个用于非线性模拟电路直流工作点分析的小型 SPICE 风格求解器，并在实现过程中学习：

- 改进节点分析（MNA）；
- Newton-Raphson 非线性方程求解；
- Source Stepping 同伦/延拓方法；
- BJT 的 Ebers-Moll 模型；
- 网表解析、数值稳定性与求解器性能优化。

DC 工作点的数学目标是求解：

$$F(x)=0$$

其中 $x$ 包含节点电压与独立电压源的支路电流。

---

## 2. 已完成的通用求解框架

项目目前已按功能分为：

```text
include/      求解器公开接口
src/          通用求解器、Source Stepping、网表解析器
examples/     三个电路的教学型、写死拓扑示例
netlist/      Schmitt1、Schmitt2、Schmitt3 网表及 PSpice 文件
tests/        基础求解器与参数对比测试
learning/     学习笔记
```

核心模块如下。

| 文件 | 主要职责 |
|---|---|
| `include/dc_solver.h` | 定义 `DcProblem`、求解参数和 Newton 报告接口。 |
| `src/dc_solver.c` | 稠密高斯消元、Newton-Raphson、阻尼和线搜索。 |
| `src/source_stepping.c` | 外层 $\lambda$ 延拓、自适应步长、失败回退与割线预测。 |
| `src/netlist_dc_solver.c` | 读取网表，构造 MNA 残差与 Jacobian，并调用通用求解器。 |

`DcProblem` 将具体电路与通用算法分离。具体电路只需要提供：

1. 残差函数 $F(x,\lambda)$；
2. Jacobian 矩阵 $J=\partial F/\partial x$；
3. 可选的 Newton 修正限制函数。

这样，求解器不需要知道电路中具体有哪些元件。

---

## 3. 三个 BJT 电路的 DC 求解

已完成 Schmitt1、Schmitt2、Schmitt3 三个非线性 BJT 电路的求解。

对于电路 1，特别确认了以下正确拓扑：

```text
R3 连接 x1 与 x4
x3 是 T2 集电极节点
x4 是 T2 基极节点
x3 与 x4 不直接相连
x1 与 x3 不直接相连
```

三个电路均可以使用：

- 直接 Newton-Raphson；
- Source Stepping；
- 网表通用求解器。

求得的工作点已与 MATLAB / PSpice 参考结果进行过对照。由于器件模型、温度处理和数值容差不完全相同，存在很小数值差异是正常的。

---

## 4. 网表 Parser 的完成情况

通用网表程序 `netlist_dc_solver` 已支持：

- 电阻 `R`；
- 独立电压源 `V`；
- BJT 实例 `Q`；
- `.model ... NPN`；
- `.model ... PNP`；
- `.model` 的 `+` 续行；
- BJT 的 `area` 参数；
- `DC`、`SIN`、`PULSE`、`PWL` 电压源的 DC 初值。

NPN 与 PNP 共用同一套 Ebers-Moll 方程结构，只通过极性参数区分。对于 NPN，极性为 $+1$；对于 PNP，极性为 $-1$。

当前 Parser 的限制包括：

- 不支持 MOSFET；
- 不支持子电路 `.subckt` 展开；
- 不支持电流源和受控源；
- 不支持完整的 Gummel-Poon BJT 模型；
- 使用固定最大维度和稠密矩阵。

---

## 5. Newton-Raphson 与 Source Stepping

固定 $\lambda$ 时，Newton-Raphson 求解：

$$J(x^{(k)})\Delta x^{(k)}=-F(x^{(k)})$$

并更新：

$$x^{(k+1)}=x^{(k)}+\alpha\Delta x^{(k)}$$

其中 $\alpha$ 是线搜索得到的阻尼系数。

Source Stepping 将独立源缩放为：

$$V_{source}(\lambda)=\lambda V_{source}$$

并逐步增大：

$$\lambda:0\rightarrow1$$

每一个新的 $\lambda$ 点都把上一个已收敛工作点作为 Newton 初值。

默认规则式步长策略是：

```text
Newton 次数 <= 5：步长乘 1.5
Newton 次数 >= 15：步长乘 0.5
其他情况：步长不变
Newton 失败：保留旧解，步长乘 0.5 后重试
```

---

## 6. 已完成的数值稳定性改进

### 6.1 节点电压修正限制

程序限制单个未知量的最大修正量，避免 Newton 一步跳到远离物理解的位置。

默认设置为：

$$|\Delta x_i|\leq0.20V$$

### 6.2 BJT 结电压后备限制

BJT 电流包含指数项，因此真正敏感的是结电压：

$$V_{BE}=V_B-V_E$$

以及：

$$V_{BC}=V_B-V_C$$

项目已实现 B-E 和 B-C 结电压的后备限制。

该限制器不会在每一次 Newton 修正时强制执行，而是采用以下策略：

```text
先尝试普通阻尼 Newton 与线搜索
若不能降低残差
再限制 ΔVBE 与 ΔVBC
重新尝试该 Newton 修正方向
```

这样可避免普通收敛路径被过度限制。命令行参数为：

```powershell
--junction-step 0.10
```

传入 `0` 可关闭该后备保护。

---

## 7. Source Stepping 性能实验

### 7.1 步长策略实验

对三个 Schmitt 网表测试了不同初始步长。结果表明，过小初始步长虽然更保守，但会增加中间延拓点数量，从而增加总 Newton 次数。

| 初始步长 | Schmitt1 | Schmitt2 | Schmitt3 |
|---:|---:|---:|---:|
| 0.05 | 85 | 85 | 90 |
| 0.125 | 69 | 72 | 76 |
| 默认 0.25 | 65 | 68 | 73 |
| 0.50 | 60 | 62 | 67 |

对这三个电路，较大的首步在总 Newton 次数上更有优势；但该结论不能直接推广到所有电路，因此默认参数仍保持偏稳健的设置。

### 7.2 小步启动、提前加速实验

测试了：

```text
initial-step = 0.125
fast-threshold = 12
growth-factor = 1.5
```

结果为：

| 电路 | 总 Newton 次数 |
|---|---:|
| Schmitt1 | 67 |
| Schmitt2 | 67 |
| Schmitt3 | 71 |

该策略优于固定小步，但总体上仍不如较大首步。

---

## 8. 受保护割线预测器

默认 Source Stepping 使用前一个解作为下一点初值：

$$x^{(0)}_{k+1}=x_k$$

项目已增加可选的 `--secant-predictor`。它使用两个历史工作点构造预测：

$$x_{pred}=x_k+\frac{\lambda_{k+1}-\lambda_k}{\lambda_k-\lambda_{k-1}}(x_k-x_{k-1})$$

为保证稳定性，预测器具有保护机制：

1. 必须已有两个成功的延拓点；
2. 前一步 Newton 不能过于困难；
3. 预测点残差必须显著低于直接继承旧解的残差；
4. 预测初值求解失败时，回退到直接继承旧解；
5. 两种初值都失败时，才缩小 $\lambda$ 步长。

三个 Schmitt 电路的比较结果：

| 电路 | 默认方式 | `--secant-predictor` | 改善 |
|---|---:|---:|---:|
| Schmitt1 | 65 | 37 | 43.1% |
| Schmitt2 | 68 | 52 | 23.5% |
| Schmitt3 | 73 | 41 | 43.8% |

最终 DC 工作点保持一致，说明割线预测改善的是 Newton 初值，而不是改变了最终求得的电路解。

在 VS Code 终端中可运行：

```powershell
.\build-check\netlist_dc_solver.exe .\netlist\Schmitt1_netlist\Netlist.txt --temp 26 --secant-predictor
.\build-check\netlist_dc_solver.exe .\netlist\Schmitt2_netlist\Netlist.txt --temp 26 --secant-predictor
.\build-check\netlist_dc_solver.exe .\netlist\Schmitt3_netlist\Netlist.txt --temp 26 --secant-predictor
```

---

## 9. 当前结论

目前最有效的性能优化是受保护割线预测器，而不是单纯调整 $\lambda$ 步长。

其原因是：步长调整只能决定“走多远”；预测器直接改善每一个新延拓点的 Newton 初值，使 Newton 从更接近真实工作点的位置开始。

---

## 10. 后续计划

后续可以按照以下顺序继续扩展：

1. 在 `data` 目录中更多可解析 BJT 网表上验证 ；
2. 实现顺序 Source Stepping，例如先打开 $V_{CC}$，再打开输入源；
3. 实现 GMIN stepping；
4. 实现 pseudo-transient continuation；
5. 支持电流源、MOSFET、受控源和子电路
