# DC 工作点与弧长同伦学习笔记

这组笔记的目标不是只会调用某一个求解函数，而是能从网表和器件模型出发，写出非线性 MNA 方程，并理解不同 DC 收敛策略怎样共用同一个 Newton 内核。

```text
网表 / 元件模型
  -> MNA 未知量 x
  -> 原始残差 F(x)、Jacobian J_F(x)
  -> Newton corrector
  -> 路径控制：直接 Newton / Source / Gmin / Homotopy / 弧长续接
  -> 以原始 F(x) 验证、去重并记录 DC 解
```

## 文件顺序

| 文件 | 学习主题 | 完成后应获得的能力 |
|---|---|---|
| [newton.md](newton.md) | 非线性 DC 与 Newton-Raphson | 从 $F(x)=0$ 推到 $J\Delta x=-F$，实现可靠内层校正器 |
| [continuation.md](continuation.md) | Source/Gmin Stepping | 用连续辅助问题生成可靠初值，设计外层回退与步长控制 |
| [method.md](method.md) | 固定点同伦 | 从任意已知点 $a$ 构造通向原方程的同伦，并处理缩放 |
| [arclength.md](arclength.md) | 伪弧长预测—校正 | 追踪折返路径、记录每次穿越 $\lambda=1$ 的解 |
| [implementation.md](implementation.md) | 工程结构与实验 | 把算法落实成 Parser、MNA、LUP、测试和诊断日志 |

## 统一符号

| 符号 | 含义 |
|---|---|
| $x\in\mathbb R^n$ | MNA 未知量：节点电压和理想电压源电流 |
| $F(x)$ | 原始电路残差；原 DC 工作点满足 $F(x)=0$ |
| $J_F(x)$ | $F$ 对 $x$ 的 Jacobian |
| $\lambda$ | 同伦/源延拓参数；目标通常位于 $\lambda=1$ |
| $g$ | Gmin Stepping 中的辅助对地电导 |
| $a$ | 固定点同伦的已知起点 |
| $z=[\lambda;x]$ | 弧长续接的增广未知量 |
| $s,\Delta s$ | 路径弧长与一步弧长 |

本组始终区分三件事：

1. **原始问题**：$F(x)=0$；只有它定义真实 DC 工作点。
2. **辅助问题**：例如 $H(x,\lambda)=0$ 或 $H_g(x,g)=0$；它只用于找到原问题的解。
3. **数值策略**：Newton 限制、线搜索、步长调整、去重和日志；它们不应悄悄改变原始电路方程。

## 最终验收

- 能为电阻、电流源、电压源和一个非线性器件写出 $F,J$ 的 stamp。
- 能解释 Source、Gmin、固定点 Homotopy 的辅助方程和端点。
- 能写出伪弧长的增广残差与 Jacobian，并说明折返点为何需要它。
- 能在每个 $\lambda=1$ 交点以原始 $F$ 做 Newton 复核，并区分不同 DC 分支。

