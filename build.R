required <- c("bookdown", "rmarkdown", "knitr", "showtext", "sysfonts")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  stop(
    "Missing packages: ",
    paste(missing, collapse = ", "),
    ". Install them before building the course."
  )
}

render_edition <- function(directory, input = "index.Rmd") {
  previous <- setwd(directory)
  on.exit(setwd(previous), add = TRUE)

  bookdown::render_book(
    input = input,
    output_format = "bookdown::gitbook",
    envir = new.env(parent = globalenv())
  )
}

render_edition(".")
render_edition("en")

dir.create("docs", recursive = TRUE, showWarnings = FALSE)
landing_page <- readLines("site/index.html", encoding = "UTF-8", warn = FALSE)
landing_page <- gsub(
  "../docs/zh/index.html",
  "zh/index.html",
  landing_page,
  fixed = TRUE
)
landing_page <- gsub(
  "../docs/en/index.html",
  "en/index.html",
  landing_page,
  fixed = TRUE
)
landing_page <- gsub(
  "../docs/zh/author-bio.html",
  "zh/author-bio.html",
  landing_page,
  fixed = TRUE
)
landing_page <- gsub(
  "../docs/en/author-bio.html",
  "en/author-bio.html",
  landing_page,
  fixed = TRUE
)
writeLines(landing_page, "docs/index.html", useBytes = TRUE)
file.copy("site/landing.css", "docs/landing.css", overwrite = TRUE)
