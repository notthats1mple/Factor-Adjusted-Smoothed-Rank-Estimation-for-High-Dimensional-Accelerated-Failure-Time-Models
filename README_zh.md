# 高维加速失效时间模型的因子调整平滑秩估计

本仓库提供 v0.16 论文、完整证明附录、图件和 LaTeX 编译材料，供同行阅读与审阅。

论文研究因子结构协变量下的右删失 AFT 稀疏系数估计，将估计因子坐标与原模型的稀疏斜率相联系，在明确的模型和设计条件下建立有限样本界，并说明数值精度如何接入系数分析。计算诊断与训练内交叉验证使用六个固定合成队列。研究算法源码、训练输入和大型拟合对象保存在研究档案中；本仓库收录论文及其编译材料。

**作者：** Weineng ZHU、Yang JU。**状态：** v0.16 未发表作者审阅稿，采用 *Lifetime Data Analysis* 投稿格式。作者声明及通信作者身份仍待确认，所附 PDF 保留这一稿件状态。

## 阅读入口

| 文件 | 内容 |
|---|---|
| [主文](manuscript/main.pdf) | 模型、估计量、理论、计算与数值研究 |
| [Online Resource 1](manuscript/ESM_1.pdf) | 完整证明、数值证据说明与实现细节 |
| [连续阅读版](manuscript/complete.pdf) | 主文与附录合并阅读 |

以上 PDF 从 v0.16 作者包原样复制。原件身份与仓库调整见 [来源说明](docs/PROVENANCE.md)。

[审阅指南](docs/REVIEW_GUIDE.md) 对应主要论证和适合反馈的位置；[复现说明](docs/REPRODUCIBILITY.md) 说明编译、引文检查和已有核验记录。

## 编译

在仓库根目录运行：

```sh
sh manuscript/build.sh
```

需要 TeX Live 或 MacTeX 中的 pdfLaTeX，以及 Python 3。交叉引用导出仅使用 Python 标准库，正常编译读取随包的参考文献文件。Springer 模板文件及独立字号覆盖的来源见 [来源说明](docs/PROVENANCE.md)。

脚本编译 `main.tex`、`supplement.tex`、`complete.tex`，并将 `supplement.pdf` 复制为 `ESM_1.pdf`；执行内容限于文档构建。PDF 核查和图件格式处理所需的额外工具列在 [复现说明](docs/REPRODUCIBILITY.md)。

## 引用与反馈

引用本版本时，请使用 [CITATION.cff](CITATION.cff) 中的未发表论文条目。反馈请标明文档、节或公式标签，并给出具体问题或建议修改。涉及数值结果时，请注明场景、方法和调参值，以便对应原始证据。


## 使用许可

作者尚未为论文和项目脚本指定再使用许可证。Springer 模板保留上游使用条件，详见[第三方材料说明](docs/THIRD_PARTY.md)。
