# 高等统计学在线课程

本项目正在把原有 LaTeX/Beamer 课件和配套统计代码迁移为中英文两套 Bookdown 在线讲义。

## 在线阅读

- [课程双语入口](https://wzh621.github.io/Advanced-Statistics-xjufe/)
- [中文版](https://wzh621.github.io/Advanced-Statistics-xjufe/zh/)
- [English edition](https://wzh621.github.io/Advanced-Statistics-xjufe/en/)

GitHub 仓库：<https://github.com/wzh621/Advanced-Statistics-xjufe>

## 构建

在项目根目录执行：

```r
source("build.R")
```

构建脚本依次生成中文版和英文版：

- `docs/zh/`：中文版；
- `docs/en/`：English edition；
- `docs/index.html`：双语入口。
- 封面包含“王芝皓副教授”作者卡片，并链接到中英文作者简介页。

`site/index.html` 是可直接在本地打开的封面源文件，其按钮指向 `docs/zh/index.html` 和 `docs/en/index.html`。构建时 `build.R` 会自动把链接改写为适合 GitHub Pages 的 `zh/index.html` 和 `en/index.html`。

图形通过 `showtext::showtext_auto()` 渲染文字。构建环境应安装 `showtext`、`sysfonts` 以及至少一种中日韩字体；推荐 `Noto Sans CJK SC`。

## GitHub Pages

网站构建产物保存在 `docs/`。GitHub Pages 使用 `main` 分支的 `/docs` 目录发布；每次更新课程后，应先运行 `Rscript build.R` 并通过站点质量检查，再提交生成结果。

## 当前进度

- Bookdown 项目骨架已建立；
- 第 0–10 章已完成中英文迁移，并持续补全证明、例题与教学实验；
- 第 7–9 章新增 Rao--Blackwell 方差比较、功效与 p 值实验、置信区间覆盖比较，以及贯穿三章的产品缺陷率案例；
- 每个主要小节使用接近课堂投影的独立教学画布；
- 桌面端教学画布与双语封面统一采用 16:9；超长内容在画布内部滚动，移动端恢复自然高度；
- 桌面端方格内部滚动时，小节标题固定在方格顶部，帮助读者持续确认当前主题；
- 后续各章必须逐项覆盖源 `.tex` 中的定义、定理、例题、证明、公式、图形和教学说明；
- 原始材料在原位置与 `legacy/` 归档副本中同时保留。
