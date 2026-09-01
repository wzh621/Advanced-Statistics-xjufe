# LaTeX source coverage

本清单是后续迁移的最低内容边界。状态“已覆盖”表示相应非空源片段的教学内容已经进入中英文讲义；“待迁移”表示当前只有章节骨架。章节主文件主要负责 Beamer 元数据、宏包、主题和 `\input` 关系，也必须在迁移时检查其中未下沉到片段文件的正文。

| 章节 | 源文件 | 中文目标 | 英文目标 | 状态 |
|---|---|---|---|---|
| 0 | `chap-0.tex` | `00-course-introduction.Rmd` | `en/00-course-introduction.Rmd` | 已覆盖 |
| 1 | `chap-1.tex`; `1-1_ST.tex`; `1-2_PT.tex`; `1-3_CP.tex`; `1-4_RV.tex`; `1-5_DF.tex`; `1-6_DM.tex` | `01-probability-theory.Rmd` | `en/01-probability-theory.Rmd` | 已覆盖；源缺失与修正见迁移记录 |
| 1 空文件 | `tex/chap-1/2-3_MM.tex`; `tex/chap-1/2-4_DU.tex` | 不产生正文 | 不产生正文 | 已核实为零字节 |
| 2 | `chap-2.tex`; `2-1_DF.tex`; `2-2_EV.tex`; `2-3_MM.tex`; `2-4_DU.tex` | `02-transformations-expectations.Rmd` | `en/02-transformations-expectations.Rmd` | 已覆盖；公式修正见迁移记录 |
| 3 | `chap-3.tex`; `3-1_DD.tex`; `3-2_CD.tex`; `3-3_EF.tex`; `3-4_DU.tex`; `pic/fig_distr.png` | `03-distribution-families.Rmd` | `en/03-distribution-families.Rmd` | 已覆盖；分布公式及推导修正见迁移记录 |
| 4 | `chap-4.tex`; `4-1_JM.tex`; `4-2_CD.tex`; `4-3_BT.tex`; `4-4_HM.tex`; `4-5_MD.tex`; `4-6_IN.tex`; `MD.png`; `CD.png`; `501003.pdf`; `501004.pdf`; `p180.jpg`; `p190.jpg` | `04-multiple-random-variables.Rmd` | `en/04-multiple-random-variables.Rmd` | 已覆盖；公式及源图处理见迁移记录 |
| 5 | `chap-5.tex`; `5-1_BC.tex`; `5-2_SR.tex`; `5-3_OS.tex`; `5-4_Co.tex`; `chapter-5.md` | `05-random-samples.Rmd` | `en/05-random-samples.Rmd` | 已覆盖；截断补稿及公式修正见迁移记录 |
| 6 | `chap-6.tex`; `6-1_In.tex`; `6-2_SP.tex`; `6-3_LP.tex`; `6-4_EP.tex` | `06-data-reduction.Rmd` | `en/06-data-reduction.Rmd` | 待迁移 |
| 7 | `chap-7.tex`; `7-1_In.tex`; `7-2_MF.tex`; `7-3_ME.tex` | `07-point-estimation.Rmd` | `en/07-point-estimation.Rmd` | 待迁移 |
| 8 | `chap-8.tex`; `8-1_In.tex`; `8-2_FT.tex`; `8-3_ET.tex` | `08-hypothesis-testing.Rmd` | `en/08-hypothesis-testing.Rmd` | 待迁移 |
| 9 | `chap-9.tex`; `9-1_In.tex`; `9-2_FT.tex`; `9-3_ET.tex` | `09-interval-estimation.Rmd` | `en/09-interval-estimation.Rmd` | 待迁移 |
| 10 | `chap-10.tex`; `10-1_PE.tex`; `10-2_Ro.tex`; `10-3_HT.tex`; `10-4_IE.tex` | `10-asymptotic-evaluations.Rmd` | `en/10-asymptotic-evaluations.Rmd` | 待迁移 |

## 每章完成条件

- 主文件和全部片段中的标题、定义、定理、引理、例题、证明、公式、表格、图形与教学说明均有明确去向。
- 中英文版内容范围一致，公式和 R 方法一致。
- 原始 R 代码已转为可执行且标签唯一的 knitr chunk，或明确记录无法运行的原因。
- 原稿错误、缺失条件及兼容性修正已写入 `migration_notes.md`。
- 两个语言版本都能从干净 R 会话构建，并完成图形文字和页面视觉检查。
