# Course materials map

| Source | Role | Sections | Current destination |
|---|---|---|---|
| `课程课件B/chap-0.tex` | 课程介绍 Beamer 主文件 | 课程信息、目标、内容与考核 | `00-course-introduction.Rmd` / `en/00-course-introduction.Rmd` |
| `课程课件B/chap-1.tex` | 第 1 章 Beamer 主文件 | 六个一级主题 | `01-probability-theory.Rmd` / `en/01-probability-theory.Rmd` |
| `课程课件B/tex/chap-1/1-1_ST.tex` | 集合论片段 | 样本空间、事件、集合运算、划分 | 中英文第 1 章 |
| `课程课件B/tex/chap-1/1-2_PT.tex` | 概率论片段 | $\sigma$-代数、Kolmogorov 公理 | 中英文第 1 章 |
| `课程课件B/tex/chap-1/1-3_CP.tex` | 条件概率片段 | 条件概率、Bayes 公式、独立性 | 中英文第 1 章 |
| `课程课件B/tex/chap-1/1-4_RV.tex` | 随机变量片段 | 随机变量及三次抛硬币例题 | 中英文第 1 章 |
| `课程课件B/tex/chap-1/1-5_DF.tex` | 分布函数片段 | cdf 定义与性质 | 中英文第 1 章 |
| `课程课件B/tex/chap-1/1-6_DM.tex` | pmf/pdf 片段 | 质量函数、密度函数、Logistic 分布 | 中英文第 1 章 |
| `课程课件B/pic/chap-1.R` | 配套绘图代码 | Logistic cdf | `01-probability-theory.Rmd` 的可执行 chunk |
| `课程课件B/pic/fig1_6_1.png` | 原课件配图 | Logistic 区间概率 | `images/ch01/logistic-density-interval.png` |
| `statistical_inference_code/week_1/random variables and univariate distributions.rmd` | 配套 R 示例 | 三次抛硬币 | `01-probability-theory.Rmd` 的可执行 chunk |
| `课程课件B/chap-2.tex` 与 `tex/chap-2/*.tex` | 第 2 章 Beamer 主文件与四个片段 | 随机变量变换、期望、矩与 MGF、积分号下求导 | 中英文第 2 章 |
| `课程课件B/pic/fig_2-1-1.png` | 原课件配图 | 正弦平方变换 | `images/ch02/fig_2-1-1.png` |
| `Statistical_Inference_Casella_Berger.pdf` | 教材概念主线 | 全书 | 仅归档于 `legacy/private-sources/`，不进入公开输出 |

源目录共有 51 个 `.tex` 文件，其中 11 个是章节主文件，38 个是非空章节片段，另有两个位于 `tex/chap-1/` 的零字节占位文件。完整的逐章覆盖状态见 `source_coverage.md`。

第 0–2 章已完成中英文迁移；第 3–10 章已建立中英文目标 Rmd 文件，尚未迁移正文、代码或图片。
