# SFU 固定点 Homotopy 复现实验

本实验严格参考根页面的固定点同伦：

```text
H(x, λ) = (1 − λ)G(x − a) + λF(x)
```

其中 G=10⁻³I，a 使用三个参考 MATLAB `Main.m` 给出的起始向量。路径由弧长 s 参数化；每步计算 `DH=[Hλ Hx]` 的切线，使用 Euler predictor 与 Adams-Moulton 梯形 corrector 控制局部误差。穿过 λ=1 时，以插值点为初值回到原始 F(x)=0 做 Newton 复核，并按变量最大差去重。

## 文件与边界

- `src/homotopy.c`：通用固定点同伦路径算法。
- `tests/reference_homotopy_experiments.c`：根页面三套参考 MNA 方程：8 维 Schmitt 1、7 维 Schmitt 2、18 维 Chua。
- 不修改 `examples/`。参考 MNA 与项目早期简化 BJT 示例的维度不同，故独立存放在测试中。

为避免手抄大型符号 Jacobian 出错，参考实验的 Jacobian 用中心差分构造；残差方程、变量顺序、器件参数和起始向量逐项来自根页面直接下载的 MATLAB 包。

## 已复现结果

| 电路 | MNA 维度 | 参考工作点数 | 本实验复现数 |
|---|---:|---:|---:|
| Schmitt Trigger 1 | 8 | 3 | 3 |
| Schmitt Trigger 2 | 7 | 3 | 3 |
| 四管 Chua | 18 | 9 | 9 |

每个解的完整 MNA 向量和 `||F||∞` 由实验程序打印。CTest 将 3/3/9 作为回归断言：任一数量变化即测试失败。

运行：

```powershell
cmake --build build-homotopy --target reference_homotopy_experiments
.\build-homotopy\reference_homotopy_experiments.exe
```

## 路径图输出

实验程序可选地把每一个被接受的弧长步导出为 CSV；绘图脚本不会删减路径点。每个子图是一项 MNA 未知量随 λ 的变化；红色虚线为 λ=1，红点为经原始方程 Newton 复核的工作点。

```powershell
New-Item -ItemType Directory -Force learning\note\homotopy-path-data, learning\note\assets | Out-Null
.\build-homotopy\reference_homotopy_experiments.exe .\learning\note\homotopy-path-data
D:\anaconda\python.exe .\learning\tools\plot_homotopy_paths.py .\learning\note\homotopy-path-data .\learning\note\assets
```

将得到 `schmitt1_homotopy_paths.png`、`schmitt2_homotopy_paths.png` 和 `chua_homotopy_paths.png`，并同时生成可编辑文字的 SVG、PDF 以及 600 dpi TIFF 版本。

参考：<https://www.sfu.ca/~ljilja/cnl/projects/Homotopy/>。
