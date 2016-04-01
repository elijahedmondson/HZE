#' @title MGI features of HZE QTL Support Intervals.
#' @author Elijah F Edmondson, \email{elijah.edmondson@@gmail.com}
#' @details Jan. 2, 2016
#' @param file: the location of your file
#' @param chr: the chromosome to evaluate
#' @param LODcutoff: the LOD score threshold
#' @param type: "all", "gene", "pseudogenic_transcript", "pseudogenic_exon",
#' "pseudogene", "match", "match-part", "transcript", "exon",
#' "mRNA", "five_prime_UTR", "start_codon", "CDS", "stop_codon",
#' three_prime_UTR", "pseudogenic_mRNA", "pseudogenic_start_codon",
#' "pseudogenic_CDS", "pseudogenic_stop_codon", "pseudogenic_five_prime_UTR",
#' "pseudogenic_three_prime_UTR", "sequence_feature"
#' @export

maxLOD = function(file, chr, type = "gene", LODcutoff = 6) {
        library(DOQTL)
        print("Loading file...")
        load(file)
        top <- max(-log10(result[[chr]]$pv))
        print(paste0("Maximum LOD score on Chr ", chr, " is ", top, ","))
        print(paste("located at",
                    result[[chr]]$POS[which(-log10(result[[chr]]$pv) >= top)], "bp."))
        max.LOD.position <- result[[chr]]$POS[which(-log10(result[[chr]]$pv) > LODcutoff)]
        print(paste("QTL interval for LOD >", LODcutoff, ": CHROMOSOME", chr, ",",
                    min(max.LOD.position), "-", max(max.LOD.position)))
        start = min(max.LOD.position)
        end = max(max.LOD.position)
        print("Loading information...")
        mgi = get.mgi.features(chr = chr, start = start, end = end, type = type, source = "MGI")
        print(mgi$Name)
        return(mgi)
}
