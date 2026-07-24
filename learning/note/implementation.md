# DC 求解器实现、验证与改动记录

## 1. 推荐目录构成

```text
include/dcsolve/     数据结构、公共接口
src/                 Parser、Circuit、MNA、Nonlinear、LUP、Homotopy、main
examples/            .cir 网表示例
tests/               LUP、stamp、Newton、continuation 测试
learning/note/       本组学习笔记
```

这比“一个巨大 main.c”更利于验证：Parser 的错误、MNA 的编号错误、器件 Jacobian 错误与续接算法错误可以分开定位。

## 2. 建议的数据与接口

```c
typedef struct {
    int dimension;
    const void *context;
    void (*residual)(const void *context, const double *x, double *f);
    void (*jacobian)(const void *context, const double *x, double *j);
} DcProblem;
```

原始 `DcProblem` 只负责 $F,J_F$。Source、Gmin、Homotopy 应作为 wrapper 提供新的 residual/Jacobian，而不是污染每个器件函数。这样可验证：

```text
base problem        -> F, JF
source wrapper      -> Hs, Hs_x
gmin wrapper        -> Hg, Hg_x
fixed-point wrapper -> H, Hx, Hlambda
arclength wrapper   -> R, Jaug
```

## 3. MNA 关键约定

- 参考地不进入未知量；非地节点使用稳定编号。
- 电阻从节点 $p$ 到 $q$ 的导纳 stamp 为 $+g$ 对角、$-g$ 非对角。
- 电流源方向必须文档化；同一方向约定决定 KCL 残差符号。
- 理想电压源添加一个支路电流未知量和一条 KVL 行。
- 非线性器件每个 Newton 迭代都根据当前 $x$ 重新 stamp 残差与 Jacobian。

对同伦专用行重排，保留标准 MNA $F$ 不动，只在 `homotopy_residual/jacobian` wrapper 内产生重排版本。原始校验、普通 Newton 和输出应仍使用标准方程。

## 4. LUP 与数值检查

每个 Newton/Corrector 都依赖线性求解。LUP 至少需要部分主元选择：选择当前列绝对值最大的可用主元；若小于阈值，返回“奇异/近奇异”而不是除零。测试应覆盖：

1. 已知 $A x=b$ 的精确小矩阵；
2. 需要换行才能正确求解的矩阵；
3. 奇异矩阵的可诊断失败；
4. 求解后复算 $\|Ax-b\|$。

稠密 LUP 适合学习和几十个未知量；大型 SPICE 必须转向稀疏数据结构、填充控制和稀疏 LU，但接口层次无需改变。

## 5. 分层测试顺序

```text
Parser           -> 元件、节点、.MODEL、.START 的解析
Linear MNA       -> R/I/V 的 A*x=b
Nonlinear stamp  -> F/J 与数值差分一致
Newton           -> 温和二极管/BJT 的工作点
Continuation     -> Source/Gmin 的接受、回退、终点
Homotopy arc     -> 切线、corrector、lambda=1 交点、去重
```

每一步成功后才进入下一层。若多解例子失败，先把问题缩小为一个器件的 $F,J$；不要先怀疑 arclength。

## 6. 必要的运行日志

| 日志 | 用途 |
|---|---|
| Newton: iter、$||F||$、$||\Delta x||$、线搜索 $\alpha$ | 查 Jacobian 和局部收敛 |
| Source: $\lambda$、尝试步长、接受/拒绝 | 查外层路径 |
| Gmin: $g$、$||H_g||$、$||F||$ | 区分辅助收敛和原题收敛 |
| Arclength: $s,\lambda,\Delta s$、切线、corrector iter | 查折返与方向 |
| Solution: $x,||F||$、去重距离、发现步号 | 证明多解记录可靠 |

日志应当可以导出 CSV，用图观察而不是只看最终一行文本。

## 7. 当前阶段的完成定义

第一组真正完成不是“代码能打印一个电压”，而是：同一套非线性 MNA 能被直接 Newton、Source、Gmin 与固定点伪弧长调用；每种方法的失败可解释；每个报告的 DC 解都经原 $F(x)$ 复核；多解的路径、交点和去重规则可复现。

