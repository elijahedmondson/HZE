factor_ <- function(x, levels=NULL, labels=levels, na.last=NA) {

        if (is.factor(x)) return(x)
        if (is.null(levels)) levels <- sort(unique.default(x), na.last=na.last)
        suppressWarnings(f <- fmatch(x, levels, nomatch=if (isTRUE(na.last)) length(levels) else NA_integer_))
        levels(f) <- as.character(labels)
        class(f) <- "factor"
        f

}
