# Source Stepping 与 Gmin Stepping

这两个方法都不是新的非线性求解器。它们构造一串容易到困难的辅助 DC 问题；每个点仍调用 `newton()` 校正，并把上一个成功解作为下一次初值。

## 1. Source Stepping

若原方程能分为：

$$F(x)=G(x)-b=0,$$

其中 $b$ 收集独立源的贡献，则定义：

$$\boxed{H_s(x,\lambda)=G(x)-\lambda b=0,\quad0\le\lambda\le1.}$$

端点：

```text
lambda = 0: 所有独立源关闭的辅助电路
lambda = 1: 原始电路 F(x)=0
```

例如二极管电路：

$$H_s(v,\lambda)=\frac{\lambda V_s-v}{R}-I_S(e^{v/V_T}-1).$$

固定 $\lambda$ 时：

$$\frac{\partial H_s}{\partial x}=J_F(x),$$

因此内层直接复用原 Newton。电压源的 MNA 约束应变成 $v_p-v_n-\lambda V=0$，独立电流源注入应变成 $\lambda I$。

### 1.1 可靠的外层循环

```text
lambda = 0; step = initial_step
x = newton(Hs(., 0), initial_guess)

while lambda < 1:
    target = min(lambda + step, 1)
    trial = copy(x)
    if newton(Hs(., target), trial) succeeds:
        commit x = trial, lambda = target
        adapt step from corrector iteration count
    else:
        discard trial
        step *= 0.5
        fail only if step < min_step
```

这里 `trial` 必须与已接受 `x` 分开：失败点不是可信的 continuation 起点。快速收敛可放大步长，接近最大 Newton 次数则缩小下一步。

### 1.2 不能解决什么

Source Stepping 不修复浮置节点和理想源矛盾，也不能强制路径的 $\lambda$ 单调。如果路径在 $(x,\lambda)$ 空间折返，持续增大 $\lambda$ 的算法会失效；这属于 arclength 的问题。

### 1.3 受保护割线预测：改进 Newton 初值而不改变步长规则

本项目默认仍使用规则式 Source Stepping：Newton 次数不超过 5 次时将步长乘 1.5，达到 15 次时将步长减半，失败时将步长减半。这个规则只决定下一步的 $\lambda$，不改变 Newton 的初值选择。

原始初值继承为：

$$x_{trial}=x_k$$

其中，$x_k$ 是当前已接受的 $\lambda_k$ 工作点。受保护割线预测额外保存前一个已接受点 $(\lambda_{k-1},x_{k-1})$，在目标 $\lambda_{k+1}$ 处构造：

$$x_{predict}=x_k+\frac{\lambda_{k+1}-\lambda_k}{\lambda_k-\lambda_{k-1}}(x_k-x_{k-1})$$

它把前两个工作点之间的变化趋势看成局部直线，并沿该方向预测下一个点。随后仍由 Newton 求解 $H_s(x,\lambda_{k+1})=0$，因此它是 predictor-corrector 的简化形式：割线负责预测，Newton 负责校正。

不能无条件相信预测点。强非线性电路中，解曲线可能快速弯曲、接近折返，甚至存在多个工作点。程序采用三层保护：

1. 只有前一个成功点的 Newton 次数不超过 12 次，才允许使用割线预测；次数较多说明附近路径不平稳。
2. 在目标 $\lambda$ 下分别计算继承初值和预测初值的残差。只有满足：

$$\left\|F(x_{predict},\lambda_{target})\right\|_\infty<0.25\left\|F(x_k,\lambda_{target})\right\|_\infty$$

才采用预测点。也就是预测必须把残差至少降至原初值的四分之一。
3. 即使预测点残差较小，Newton 也可能失败。这时不立即缩小步长，而是回退到 $x_k$，在相同的 $\lambda_{target}$ 重新执行一次 Newton；只有两次尝试都失败，才把步长减半。

伪代码为：

```text
trial = x_k

if 有前两个成功点 且 前一步足够平稳:
    predict = 割线外推(x_{k-1}, x_k, lambda_target)
    if 残差(predict) < 0.25 * 残差(x_k):
        trial = predict
        used_predictor = true

if Newton(trial, lambda_target) 失败 且 used_predictor:
    trial = x_k
    再次执行 Newton(trial, lambda_target)

if 仍失败:
    step = step / 2
else:
    接受 trial
    仍按原规则调整下一步 step
```

程序中默认使用原策略；命令行增加 `--secant-predictor` 时启用此方案：

```powershell
.\build-check\netlist_dc_solver.exe .\data\Benchmark\bjt\rca.sp --secant-predictor --temp 27
```

在当前可比较的 BJT 数据集上，统计所有实际执行过的 Newton 迭代（包括预测失败后的回退尝试）得到：

| 网表 | 原策略 | 受保护割线预测 | 变化 |
| --- | ---: | ---: | ---: |
| `astabl.sp` | 36 | 28 | -22.2% |
| `rca.sp` | 93 | 52 | -44.1% |
| `schmitecl.sp` | 37 | 28 | -24.3% |
| `bias.sp` | 335 | 331 | -1.2% |
| 合计 | 501 | 439 | -12.4% |

四个网表均收敛到相同工作点。`astabl.sp` 两种结果的最大节点电压差约为 $2\times10^{-9}\mathrm{V}$，属于数值舍入差异。这个方案仍是以 $\lambda$ 为自然参数的 continuation；若路径出现真正的折返点，仍需要伪弧长延拓处理。

## 2. Gmin Stepping

Gmin 在每个非地节点临时并联对地电导。原方程变为：

$$\boxed{H_g(x,g)=F(x)+gDx=0.}$$

$D$ 选择节点电压变量。例如 $x=[v_1,v_2,i_V]^T$ 时：

$$D=\mathrm{diag}(1,1,0).$$

因为 $i_V$ 是 MNA 的电压源支路电流，不应被错误解释为“对地漏电流”。

对 Jacobian：

$$\boxed{J_{H_g}=J_F+gD.}$$

从大的 $g_{start}$ 开始，逐步降至 $g_{final}$：

```text
g = g_start
solve H_g(x,g)=0
repeat:
    g_try = max(g / ratio, g_final)
    use last accepted x as Newton initial guess
    accept if successful; otherwise choose a less aggressive reduction
```

几何递减适合跨越多个数量级，如 $10^{-2},10^{-3},\ldots,10^{-12}\ \mathrm S$。

## 3. Gmin 有效的原因与风险

- 给浮置/近浮置节点 DC 泄放通路；
- 给 Jacobian 节点行增加正对角项，常能改善条件；
- 大 $g$ 会把节点温和地拉向地，形成更容易的起点。

但它改变了辅助电路的工作点。最终答案必须检查原残差 $F(x)$，并观察解对 $g_{final}$ 是否稳定。若 $g$ 越小解仍剧烈漂移，可能是原电路没有良定义 DC 解、存在多分支，或参数尚未足够小。

## 4. 两者的区别

| 项目 | Source | Gmin |
|---|---|---|
| 变化对象 | 独立激励 | 节点对地导通 |
| 起点 | 零源问题 | 强泄放问题 |
| 显式 Jacobian 修改 | 通常没有 | $+gD$ |
| 主要用途 | 降低激励造成的跳变 | 改善 DC 约束与病态 |

固定点同伦里的 `G(x-a)` 是**同伦起点缩放**，不要把它与 $gDx$ 混为同一个 “g”。前者通常作用于全部同伦方程，后者对应具体节点的物理辅助 stamp。

