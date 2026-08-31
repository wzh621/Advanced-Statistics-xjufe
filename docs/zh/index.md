---
title: "高等统计学"
subtitle: "Advanced Statistics"
author: "王芝皓 副教授"
institute: "新疆财经大学统计与数据科学学院"
date: "2026"
lang: zh-CN
site: bookdown::bookdown_site
documentclass: book
bibliography:
  - references.bib
  - packages.bib
description: "高等统计学课程在线讲义，依据原有课程课件与配套代码整理。"
---



# 前言 {-}

**新疆财经大学统计与数据科学学院**

本讲义是“高等统计学（Advanced Statistics）”课程的在线版本。课程以概率论为基础，依次讨论随机变量及其分布、随机样本、数据约简、点估计、假设检验、区间估计与渐近理论。

当前版本处于逐章迁移阶段：第 1 章已经根据原 Beamer 课件整理为可阅读、可执行的 Bookdown 章节；其余章节仅建立文件骨架，尚未迁移正文。每个主要小节采用接近课堂投影的独立“教学画布”，但仍保留在线讲义所需的完整推导与上下文。

## 使用说明 {-}

- 数学公式使用 MathJax 渲染。
- 图形使用 `showtext_auto()` 和自动探测的中日韩字体渲染中文。
- R 示例均从项目根目录以相对路径运行。
- 原始 LaTeX、代码、图片和教材资料保存在 `legacy/` 的归档副本中；原位置的材料也继续保留。
- 迁移中发现的原稿疑点和必要修正记录在 `migration_notes.md`。
- 右上角的语言入口可切换至英文版。
