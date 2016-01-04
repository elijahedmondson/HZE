#' @title Plot function (w/ binning to average markers and max LOD)
#' @export

GRconvert = function(result) {
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

}
