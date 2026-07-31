# 第五章：瞬态分析学习笔记

本目录把第五章按“一个瞬态求解器如何工作”的链路拆成独立主题，而非按课本小节编号命名。阅读顺序如下：

1. [瞬态分析基础](transient_foundations.md)：明确要求解的数学问题，以及 DC 与瞬态的关系。
2. [数值积分方法](numerical_integration.md)：向前欧拉、向后欧拉、梯形法和多步法。
3. [储能元件伴随模型](companion_models.md)：把电容、电感和互感变成当前时间点的代数模型。
4. [误差估计与时间步控制](error_and_step_control.md)：决定一个试算时间点是否可接受。
5. [稳定性与刚性电路](stability_and_stiffness.md)：解释隐式积分为何适合 SPICE。
6. [Gear 与后向差分公式](gear_bdf.md)：理解 BDF/Gear 的历史点与刚性处理。
7. [瞬态求解器流程](transient_solver.md)：把上述模块组合成可实现的程序。
8. [振荡器瞬态分析实例](oscillator_examples.md)：将算法用于考比兹振荡器和多谐振荡器。

整章主线为：

```text
F(x, x_dot, t) = 0
        ↓ 数值积分
F_n+1(x_n+1) = 0
        ↓ Newton-Raphson
得到试算解
        ↓ LTE 与收敛检查
接受并保存历史，或回退并缩小步长
```

文中的 `n` 表示已接受的时间点编号；`k` 表示固定时间点内部的 Newton 迭代编号。二者不可混淆。
