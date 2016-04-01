#' @title Convert a Large list of data.frames into a GRangesList
#' @author Elijah F Edmondson, \email{elijah.edmondson@@gmail.com}
#' @export

GRconvert = function(file) {
        load(file = file)
        result = qtl
        library(GenomicRanges)
        library(BSgenome.Mmusculus.UCSC.mm10)
        library(DOQTL)

        chrs = c(1:19, "X")
        qtl = GRangesList(GRanges("list", length(result)))

        for(i in 1:length(chrs)) {
                print(i)
                qtl[[i]] <- GRanges(seqnames = Rle(result[[i]]$CHR),
                                    ranges = IRanges(start = result[[i]]$POS, width = 1),
                                    p.value = result[[i]]$pv)
        } # for(i)
        save(qtl, file.prefix, file = paste0(file.prefix, "_GR_QTL.Rdata"))
        return(qtl)

}
