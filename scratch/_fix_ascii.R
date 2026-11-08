base  <- "C:/jtaylor/Research/Computing/R Packages/Claude/wgaim/"
files <- list.files(paste0(base, "R"), pattern = "\\.R$", full.names = FALSE)
files <- paste0("R/", files)

subs <- list(
    c("\u2014", "--"),       # em dash
    c("\u00b1", "+/-"),      # plus-minus sign
    c("\u00d7", "x"),        # multiplication sign
    c("\u00b2", "^2"),       # superscript 2 (h^2)
    c("\u00b7", "."),        # middle dot
    c("\u2019", "'"),        # right single quotation mark
    c("\u2018", "'"),        # left single quotation mark
    c("\u201c", '"'),        # left double quotation mark
    c("\u201d", '"'),        # right double quotation mark
    c("\u2192", "->"),       # rightwards arrow
    c("\u2190", "<-"),       # leftwards arrow
    c("\u2248", "~="),       # almost equal
    c("\u00e9", "e"),        # e with acute accent
    c("\u03b2", "beta"),     # greek small beta
    c("\u03a3", "Sigma"),    # greek capital sigma
    c("\u2264", "<="),       # less-than or equal
    c("\u2265", ">="),       # greater-than or equal
    c("\u00b0", "deg"),      # degree sign
    c("\u00e0", "a"),        # a with grave
    c("\u2026", "..."),      # horizontal ellipsis
    c("\u2212", "-")         # minus sign (U+2212, distinct from ASCII hyphen)
)

for (f in files) {
    path <- paste0(base, f)
    raw  <- readBin(path, "raw", n = file.info(path)$size)
    txt  <- rawToChar(raw)
    for (s in subs) txt <- gsub(s[1], s[2], txt, fixed = TRUE)
    writeBin(charToRaw(txt), path)

    remaining <- which(grepl("[^\x01-\x7F]", readLines(path, warn = FALSE)))
    if (length(remaining)) {
        cat(sprintf("STILL DIRTY: %s  lines: %s\n", f,
                    paste(remaining, collapse = ", ")))
    } else {
        cat(sprintf("CLEAN: %s\n", f))
    }
}
