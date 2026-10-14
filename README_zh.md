# 高维 AFT 模型的因子调整平滑秩估计

本仓库提供右删失高维加速失效时间（AFT）模型的 R 研究代码，包括训练内 PCA、部分 L1 惩罚的高斯平滑 Gehan 目标、完整坐标 KKT 诊断和冻结预测映射。

**作者：** Weineng Zhu、Yang Ju · **代码版本：** 0.17.0 · [English](README.md)

## 快速开始

使用 [renv.lock](renv.lock) 记录的 **R 4.5.2**。在新检出的仓库根目录依次运行：

```sh
Rscript --vanilla scripts/restore_dependencies.R
Rscript --vanilla scripts/check_environment.R
Rscript --vanilla scripts/run_checks.R
Rscript --vanilla examples/fixed_array_demo.R
```

依赖恢复用于首次准备环境；已有锁定版本的环境可从环境检查开始。演示使用固定数组和有更新上限的小例，无需下载数据。完整操作见[复现说明](docs/REPRODUCIBILITY.md)。

## 代码内容

| 内容 | 入口 |
|---|---|
| 训练内中心化、PCA 与冻结变换 | `aft_preprocess()` |
| 平滑目标、完整梯度与 Hessian 向量积 | `aft_loss()`、`aft_hvp()` |
| 原近端梯度求解器 | `aft_solve()` |
| 非重启回溯 FISTA | `aft_solve_fista_experimental()` |
| 固定活跃符号面的 Newton 精修 | `aft_polish_active_face()` |
| 独立终点复核与条件性精修调度 | `aft_solve_cv_v012()` |
| 新样本冻结预测 | `aft_predict_map()` |
| 折内预处理与验证算术 | `p3a_prepare_split()`、`p3a_validation_blocked()` |

[算法与接口](docs/ALGORITHMS.md) 说明目标、参数排列、停止条件和返回字段；[审阅指南](docs/REVIEW_GUIDE.md) 对应源码和测试。

在 R 中可从任意工作目录加载：

```r
source("/path/to/checkout/scripts/load_project.R")
config <- aft_public_config()
```

共同适配器的下游计算应使用 `result$standard`。精修成功后，其 `coefficients`、`gradient`、`beta`、`b_pred` 和 `fitted_scores` 均从最终新系数重建。原始 FISTA 和精修记录分别保存在 `result$fista`、`result$polish`。

## 文件与研究材料

`R/`、`optim/`、`diagnostics/` 和 `tuning/` 保存算法及验证组件；`config/` 保存可移植控制参数；`examples/` 和 `tests/` 提供固定小例与检查；`provenance/` 记录归档来源及 SHA256；`docs/` 提供审阅和复现文档。

本次发布面向代码审阅和使用者自备训练数据的计算。论文、研究队列、大型拟合对象和研究轨迹保存在作者本地归档中。仓库检查复现固定计算案例；完整研究结果的重建还需要归档中的输入和执行记录。

## 引用与贡献

软件引用信息见 [CITATION.cff](CITATION.cff)，问题反馈与修改要求见 [CONTRIBUTING.md](CONTRIBUTING.md)。

作者尚未为本代码指定再使用许可证，保留全部权利；如需再使用授权，请通过仓库联系作者。依赖包遵循各自许可证。
