# 瞬态求解器流程：从初值到接受的时间点

## 1. 求解器的职责

完整瞬态求解器反复执行：选择试算时间点、生成动态伴随模型、进行 Newton、估计误差、接受或回退。它复用 DC 求解器的非线性残差、Jacobian 和线性求解器；新增的是时间、历史和步长控制。

## 2. 初始化

程序需要准备：

```text
电路拓扑与器件模型
积分方法、容差、最大 Newton 次数
tstop、TMAX、输出时间序列、已知断点
初始状态与动态元件历史
```

默认初始状态来自 DC 工作点；若用户明确指定初值并要求跳过工作点，则直接采用该状态。电容初始电压和电感初始电流必须与状态初始化一致。

## 3. 每一步的预测—校正—接受

```text
accepted state at tn
      ↓  选择 htrial，不能跨越断点或 tstop
预测 xn+1
      ↓
用 htrial 建立电容/电感/互感伴随模型
      ↓
在 tn+1 进行 Newton 迭代
      ↓
Newton 失败？──是→ 回退、缩小 h、重试
      ↓ 否
估计归一化 LTE
      ↓
LTE 超差？────是→ 回退、缩小 h、重试
      ↓ 否
接受：提交状态、更新历史、输出、选择下一 h
```

## 4. 当前时间点的方程装配

每次 Newton 迭代都应在当前猜测 `x(k)` 下：

```text
1. stamp 电阻、受控源及 t_n+1 时刻的独立源；
2. stamp 当前 h 下储能元件的 Geq 和历史项；
3. 线性化二极管/BJT/MOS 等非线性器件；
4. 构造 F(x(k)) 与 J(x(k))；
5. 解 J Δx = -F 并更新。
```

`h` 改变时，储能 stamp 会改变；不能只沿用旧矩阵数值。

## 5. Newton 收敛判断

应至少检查更新量，并最好检查残差。变量容差通常采用相对加绝对形式：

```text
|Δvi| <= VNTOL + RELTOL max(|vi,new|,|vi,old|)
|Δii| <= ABSTOL + RELTOL max(|ii,new|,|ii,old|)
```

仅看更新量可能造成假收敛，例如矩阵病态时步长很小而 KCL 残差仍大。

## 6. 状态与回滚设计

建议显式区分：

```text
accepted_state：xn 与所有已提交器件历史
trial_state：   tn+1 的候选解
newton_state：  trial 内部迭代值
history：       多步法的已接受时间/状态序列
```

重试必须从 `accepted_state` 重新开始；不要复制失败候选中的电容电压、电感电流或 Gear 历史。

## 7. 断点、输出与内部步长

内部积分点由误差控制决定，输出间隔不应被误认为固定积分步长。无论提出的 `h` 多大，都要截断到下一个已知 PULSE/PWL/开关断点。输出点可通过内插获得，但动态断点必须真实落点。

## 8. 一个实现化伪代码

```text
state = initialize_from_DC_or_user_IC()
time = t0
h = initial_step

while time < tstop:
    htrial = limit_by_stop_breakpoint_and_TMAX(time, h)
    trial = predict(state, history, htrial)

    if not newton_transient(time+htrial, htrial, state, trial):
        h = shrink(htrial)
        continue

    r = normalized_LTE(state, trial, history, htrial)
    if r > 1:
        h = step_from_error(htrial, r)
        continue

    commit(trial, state, history)
    time += htrial
    write_requested_output(time, state)
    h = proposed_next_step(htrial, r)
```

## 9. 常见实现错误

* 将输出步长误作积分步长，导致错误的时间控制。
* Newton 收敛后立即提交，忽略 LTE。
* 步长变化却不更新伴随模型。
* 重试时保留失败点历史。
* 让一步跨越电源边沿或开关时刻。
* 用同一个数组同时存 accepted 与 Newton 临时状态。

正确的瞬态程序，本质上是可靠的“试算—验证—提交”事务循环。
