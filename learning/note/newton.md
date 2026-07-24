# Newton 求非线性 DC 工作点

## 1. 问题从哪里来

DC 分析中，独立源是常数，电容在稳态开路、电感在稳态短路；但二极管、BJT、MOS 等器件的电流仍是端口电压的非线性函数。因此电路并非一般的 $Ax=b$，而是：

$$F(x)=0.$$

在 MNA 中，$x$ 不只是节点电压。若电路有 $n_v$ 个非地节点和 $n_s$ 个理想电压源：

$$x=[v_1,\ldots,v_{n_v},i_{V_1},\ldots,i_{V_{n_s}}]^T.$$

前 $n_v$ 行通常是 KCL 残差，后 $n_s$ 行是电压源 KVL 约束。对所有未知量同时求零，才是完整 DC 工作点。

## 2. 一维二极管例子

```text
 Vs --- R --- o v --- diode --- GND
```

取流入节点的电阻电流减去流出节点的二极管电流：

$$F(v)=\frac{V_s-v}{R}-I_S\left(e^{v/V_T}-1\right).$$

工作点满足 $F(v)=0$。对 $v$ 求导：

$$F'(v)=-\frac1R-\frac{I_S}{V_T}e^{v/V_T}.$$

残差的整体符号可以反过来定义，但 `F` 与 `J=dF/dx` 必须使用同一约定；否则 Newton 会朝错误方向走。

## 3. Newton 推导

当前猜测为 $x_k$。在它附近一阶展开：

$$F(x_k+\Delta x)\approx F(x_k)+J_F(x_k)\Delta x.$$

令局部线性模型的残差为零：

$$\boxed{J_F(x_k)\Delta x=-F(x_k)}$$

$$\boxed{x_{k+1}=x_k+\Delta x.}$$

一维时就是 $\Delta v=-F(v_k)/F'(v_k)$。多变量时，每次迭代都要通过 LUP 或稀疏 LU 解一个线性系统。Newton 的“快”来自根附近的二次收敛；它并不表示从任意初值都可靠。

## 4. 从 Jacobian 看非线性器件 stamp

对二极管 $i_D(v)=I_S(e^{v/V_T}-1)$，在当前 $v_k$：

$$i_D(v)\approx i_D(v_k)+g_D(v_k)(v-v_k),$$

$$g_D(v_k)=\frac{I_S}{V_T}e^{v_k/V_T}.$$

整理：

$$i_D(v)\approx g_Dv+I_{eq},\qquad I_{eq}=i_D(v_k)-g_Dv_k.$$

所以每一次 Newton 都把二极管暂时替代为“微分电导 $g_D$ + 等效电流源 $I_{eq}$”，组成一张线性电路求解。这个 stamp 视角和 $J\Delta x=-F$ 完全等价。

## 5. 一个通用 corrector 的职责

```text
newton(problem, parameter, x):
  repeat at most max_iter times:
      F = residual(problem, parameter, x)
      J = jacobian(problem, parameter, x)
      fail if F/J contains NaN, Inf or J cannot be factored
      solve J * dx = -F
      dx = device_limit(x, dx)
      alpha = line_search(x, dx)       # optional but strongly recommended
      x_new = x + alpha * dx
      accept only after convergence checks
```

`parameter` 可以是 Source 的 $\lambda$、Gmin 的 $g$，或同伦的 $\lambda$。这使所有外层方法共用同一个 corrector。

### 5.1 收敛不能只看一个数字

残差小表示电路方程近似成立；修正小表示迭代点不再明显移动。建议同时检查：

$$\|W_FF(x)\|_\infty\le \tau_F,$$

$$\|W_x\Delta x\|_\infty\le \tau_x.$$

$W_F,W_x$ 是按方程、变量量纲设置的缩放。只检查步长会把“线搜索已把步长压得很小”误判为成功；只检查未缩放残差则会混淆 A、V 和不同量级的支路电流。

## 6. 失败模式与预防

| 现象 | 典型原因 | 处理 |
|---|---|---|
| `exp` 溢出 | PN 结修正过大 | 指数保护、pn-junction voltage limiting |
| 两点之间振荡 | 局部线性模型太差 | 阻尼或残差下降线搜索 |
| LU 奇异 | 浮置节点、矛盾理想源 | 拓扑检查、主元选择、Gmin 辅助 |
| 收敛到错误分支 | 多解、初值不同 | 记录初值/路径；使用延拓追踪分支 |
| 残差不降 | Jacobian 符号/偏导错误 | 数值差分验证 Jacobian |

数值差分是最重要的调试工具之一：

$$\frac{\partial F_i}{\partial x_j}\approx
\frac{F_i(x+h e_j)-F_i(x-h e_j)}{2h}.$$

它不替代解析 Jacobian，但能快速发现电流方向、节点编号或偏导漏项。

