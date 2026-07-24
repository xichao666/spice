# DC 求解器项目阶段总结

## 项目目标

本项目实现一个面向非线性模拟电路 DC 工作点分析的小型 SPICE 风格求解器，并学习 MNA、Newton-Raphson、Source Stepping、BJT 模型、网表解析和数值优化。

DC 工作点的目标是求解：

$$F(x)=0$$

其中 $x$ 包含节点电压和独立电压源支路电流。

## 已完成的框架

项目采用以下职责划分：

```text
include/      求解器公开接口
src/          通用 DC 求解器、Source Stepping、网表解析器
examples/     三个写死拓扑的教学电路
netlist/      Schmitt1、Schmitt2、Schmitt3 网表与 PSpice 文件
tests/        基础测试与参数对比程序
learning/     学习笔记和辅助工具
```

`DcProblem` 将具体电路与通用求解算法分离。电路模型提供残差 $F(x,\lambda)$、Jacobian $J=\partial F/\partial x$，并可提供可选的 Newton 修正限制函数。

## 已完成的电路与网表能力

- 已完成 Schmitt1、Schmitt2、Schmitt3 三个 BJT 电路的 DC 求解；
- 支持直接 Newton-Raphson 和 Source Stepping；
- 已确认电路 1 的关键关系：$R3$ 连接 $x1$ 与 $x4$，而非 $x1$ 与 $x3$；
- 已实现通用网表求解器 `netlist_dc_solver`；
- 支持电阻、电压源、BJT、`.model NPN/PNP`、模型续行、`area` 参数和多种电压源 DC 初值；
- NPN、PNP 共用极性统一的 Ebers-Moll 方程。

当前 Parser 尚不支持 MOSFET、受控源、电流源、`.subckt` 展开和完整 Gummel-Poon 模型。

## Newton 与 Source Stepping

固定 $\lambda$ 时，Newton 迭代求解：

$$J(x^{(k)})\Delta x^{(k)}=-F(x^{(k)})$$

并更新：

$$x^{(k+1)}=x^{(k)}+\alpha\Delta x^{(k)}$$

其中 $\alpha$ 由回溯线搜索决定。

Source Stepping 将独立源写为：

$$V_{source}(\lambda)=\lambda V_{source}$$

然后从 $\lambda=0$ 逐步推进到 $\lambda=1$。默认规则为：Newton 次数不超过 5 次时步长乘 1.5，达到 15 次时步长减半，失败时保留上一个解并缩小步长重试。

## 数值稳定性改进

- 实现稠密高斯消元和部分主元选取；
- 实现节点电压修正限制、阻尼 Newton 和回溯线搜索；
- 实现 BJT B-E、B-C 结电压后备限制；
- 结电压限制仅在普通 Newton 与线搜索均不能降低残差时触发，避免正常收敛路径额外变慢；
- 通过 `--junction-step` 参数调整结电压上限，传入 `0` 可关闭该保护。

## Source Stepping 性能实验

对三个 Schmitt 网表的步长实验表明，缩小初始步长会增加中间延拓点，通常不能降低总 Newton 次数。

| 初始步长 | Schmitt1 | Schmitt2 | Schmitt3 |
|---:|---:|---:|---:|
| 0.05 | 85 | 85 | 90 |
| 0.125 | 69 | 72 | 76 |
| 默认 0.25 | 65 | 68 | 73 |
| 0.50 | 60 | 62 | 67 |

大首步在当前三组电路上更快，但由于尚未覆盖更多网表，通用默认策略仍保持较稳健的参数。

## 受保护割线预测器

默认方式直接继承前一个工作点：

$$x_{k+1}^{(0)}=x_k$$

`--secant-predictor` 使用两个历史点预测下一工作点：

$$x_{pred}=x_k+\frac{\lambda_{k+1}-\lambda_k}{\lambda_k-\lambda_{k-1}}(x_k-x_{k-1})$$

预测器要求预测残差明显优于继承初值；若预测求解失败，会回退到继承初值。三个 Schmitt 电路的结果如下：

| 电路 | 默认方式 | `--secant-predictor` |
|---|---:|---:|
| Schmitt1 | 65 | 37 |
| Schmitt2 | 68 | 52 |
| Schmitt3 | 73 | 41 |

最终工作点一致，说明割线预测改善的是 Newton 初值和收敛效率，并未改变电路求得的解。

## 后续计划

1. 在更多可解析 BJT 网表上验证割线预测器；
2. 实现顺序 Source Stepping；
3. 加入 GMIN stepping 与 pseudo-transient；
4. 扩展 Parser；
5. 使用稀疏矩阵与稀疏 LU 替换当前稠密实现。
