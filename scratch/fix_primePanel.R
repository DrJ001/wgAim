f <- "C:/jtaylor/Research/Computing/R Packages/Claude/wgaim/R/primePanel.R"
x <- readLines(f, encoding = "UTF-8")
x[106] <- "#' @return An object of class \\code{c(\"wgPanel\", \"wgCross\")} -- a list with:"
writeLines(x, f)
cat("Done. Line 106 now reads:\n", x[106], "\n")
