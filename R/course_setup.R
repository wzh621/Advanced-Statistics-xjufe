# Shared rendering setup for both language editions.

course_font_candidates <- c(
  "Noto Sans CJK SC",
  "Noto Sans SC",
  "Source Han Sans SC",
  "WenQuanYi Micro Hei",
  "Microsoft YaHei",
  "PingFang SC",
  "Arial Unicode MS"
)

course_find_cjk_font <- function() {
  available <- sysfonts::font_families()
  matched <- course_font_candidates[course_font_candidates %in% available]

  if (length(matched)) {
    return(matched[[1]])
  }

  font_table <- sysfonts::font_files()
  pattern <- paste(
    c("NotoSansCJK", "NotoSansSC", "SourceHanSans", "WenQuanYi", "YaHei", "PingFang", "Arial Unicode"),
    collapse = "|"
  )
  matched_file <- font_table[grepl(pattern, font_table$file, ignore.case = TRUE), , drop = FALSE]

  if (nrow(matched_file)) {
    font_path <- file.path(matched_file$path[[1]], matched_file$file[[1]])
    sysfonts::font_add("course-cjk", regular = font_path)
    return("course-cjk")
  }

  warning(
    "No preferred CJK font was found. Falling back to the system sans font; ",
    "install Noto Sans CJK SC for portable Chinese figures."
  )
  "sans"
}

course_setup <- function() {
  required <- c("knitr", "showtext", "sysfonts")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing)) {
    stop(
      "Missing packages required for course rendering: ",
      paste(missing, collapse = ", "),
      "."
    )
  }

  course_font <- course_find_cjk_font()

  # autoshow: route text through showtext for every graphics device opened by knitr.
  showtext::showtext_auto(enable = TRUE)
  showtext::showtext_opts(dpi = 144)

  options(course.plot.family = course_font)
  knitr::opts_chunk$set(
    echo = TRUE,
    message = FALSE,
    warning = FALSE,
    fig.align = "center",
    fig.retina = 2,
    dpi = 144,
    out.width = "90%"
  )

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::theme_set(ggplot2::theme_minimal(base_family = course_font))
  }

  invisible(course_font)
}

course_plot_family <- function() {
  getOption("course.plot.family", "sans")
}

