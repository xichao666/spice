# 伪弧长延拓与多解 DC 记录

## 1. 从 lambda 参数改为弧长参数

令：

$$z(s)=\begin{bmatrix}\lambda(s)\\x(s)\end{bmatrix},\qquad H(x(s),\lambda(s))=0.$$

对 $s$ 求导：

$$H_\lambda\frac{d\lambda}{ds}+H_x\frac{dx}{ds}=0.$$

切向量 $t=[t_\lambda;t_x]=dz/ds$ 满足：

$$\begin{bmatrix}H_\lambda&H_x\end{bmatrix}t=0.$$

它只有方向而没有长度，因此还要以缩放度量归一化：

$$t^TWt=1.$$

这里的 $W$ 应和固定点同伦的变量尺度一致；不能把无量纲 $\lambda$、伏特和安培直接不加区分地做欧氏距离。

## 2. 切线如何计算与定向

在远离奇异点时，固定 $t_\lambda=1$，解：

$$H_xt_x=-H_\lambda,$$

再归一化即可。但 $H_x$ 接近奇异时这个办法失效，实际 arclength 实现应求增广 Jacobian 的零空间，例如用 QR/SVD 得到 $[H_\lambda\ H_x]$ 的零向量。

$t$ 与 $-t$ 都合法。设前一切线为 $t_{old}$，若：

$$t^TWt_{old}<0,$$

就翻转 $t\leftarrow-t$，保证路径持续向同一方向走。起点没有旧切线时，可选择 $t_\lambda>0$ 作为初始方向，但后续不应强制它一直为正。

## 3. Predictor

已接受点 $z_k$ 和单位切线 $t_k$ 给出：

$$\boxed{z_{pred}=z_k+\Delta s\,t_k.}$$

$\Delta s$ 是路径距离，不是 $\Delta\lambda$。预测点一般不在精确曲线上；其职责是给 corrector 一个靠近正确分支的初值。

## 4. Corrector：增广 Newton 方程

$H=0$ 只有 $n$ 条方程而 $z$ 有 $n+1$ 个未知量。加入通过预测点、法向为当前切线的超平面：

$$c(z)=t_k^TW(z-z_{pred})=0.$$

完整残差：

$$R(z)=\begin{bmatrix}H(x,\lambda)\\c(z)\end{bmatrix}.$$

完整 Jacobian：

$$\boxed{J_{aug}=
\begin{bmatrix}
H_\lambda&H_x\\
t_\lambda w_\lambda&t_x^TW_x
\end{bmatrix}.}$$

Newton 校正为：

$$J_{aug}\Delta z=-R(z),\qquad z\leftarrow z+\Delta z.$$

正确的实现同时检查 $\|H\|$、$|c|$ 和缩放后的 $\|\Delta z\|$。只检查 $H$ 会漏掉伪弧长约束未满足的情况。

## 5. 接受、拒绝与步长

```text
for each accepted point:
    predict with ds
    correct augmented system
    if corrector fails or error is too large:
        ds *= 0.5 and retry from last accepted point
    else:
        accept point, compute/orient next tangent
        enlarge ds if corrector was cheap; shrink if expensive
```

必须设置 `min_ds`、最大路径步数和最大 corrector 迭代数，防止无限重试。对高曲率区，正确的行为是少量小步，而不是一次 Newton 跳去远处分支。

## 6. lambda=1 的终局处理与去重

若相邻接受点的 $(\lambda_k-1)(\lambda_{k+1}-1)\le0$，路径可能穿越目标平面。处理流程：

1. 在该段上插值，得到接近 $\lambda=1$ 的 $x_{guess}$；线性插值可作起点，高阶 predictor 可用 Hermite 插值。
2. 固定 $\lambda=1$，以 $x_{guess}$ 对**原始** $F(x)=0$ 做 Newton。
3. 检查原始残差、变量容差和物理可行性。
4. 与既有解计算缩放距离；若大于去重阈值，记录为新解。

```text
solutions = []
if lambda plane is crossed:
    x = final_newton_on_original_F(interpolated_guess)
    if valid(x) and not duplicate(x, solutions):
        append x
```

注意：沿曲线擦碰 $\lambda=1$ 而不改变符号时，单纯符号检测会漏解；研究型实现可监视 $|\lambda-1|$ 的局部极小值，或对接近区间进一步求交。

## 7. 常见错误

| 错误 | 结果 |
|---|---|
| 只让 $\lambda$ 增加 | 折返点停止或漏分支 |
| 不做变量/方程缩放 | 步长和切线被某个量纲支配 |
| 忘记切线方向连续性 | 相邻步来回走 |
| 约束写成相对旧点而非预测点 | corrector 缺乏局部选根能力 |
| 首次 $\lambda=1$ 就停止 | 漏掉该同伦曲线上的其他交点 |
| 仅检查辅助 $H$ | 辅助问题收敛但原 DC 解未被复核 |

## 8. 最小可视化实验

每个接受点导出 `step,s,lambda,ds,corrector_iter,||H||` 与关键节点电压。画出 `lambda` 对某个节点电压的路径图；折返、多次穿越和步长收缩会一目了然。对 Schmitt、Chua 等多稳态例子，最终表格必须包含每个去重后的 $\lambda=1$ 解及其 $\|F\|_\infty$。

