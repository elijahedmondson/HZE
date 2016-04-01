#' @title Genome wide association mapping using CoxPH1
#' @author Elijah F Edmondson, \email{elijah.edmondson@@gmail.com}
#' @details Dec. 29, 2015
#' @param pheno: data.frame containing phenotypes in columns and samples in
#'                   rows. Rownames must contain sample IDs.
#' @param pheno.col: the phenotype column to map. surv = Surv(pheno$cat.days, pheno$cataract)
#' @param probs: 3D numeric array containing the haplotype probabilities
#'                   for each sample. Samples in rows, founders in columns,
#'                   markers in slices. Samnple IDs must be in rownames. Marker
#'                   names must be in dimnames(probs)[[3]].
#' @param K: List of kinship matrices, one per chromosome in markers.
#' @param addcovar: data.frame of additive covariates to use in the mapping.
#'                      Sample IDs must be in rownames.
#' @param markers: data.frame containing at least 3 columns with marker names,
#'                     chr, Mb postion.
#' @param snp.file: character string containing the full path to the Sanger
#'                      SNP file containing genomic positions and SNPs. This file
#'                      is created using condense.sanger.snps().
#' @export


GRSD.coxph1 = function(pheno, pheno.col, probs, K, addcovar, markers, snp.file,
                      outdir = "~/Desktop/files/", tx = c("Gamma", "HZE", "Unirradiated")){
        begin <- Sys.time()
        begin
        # COVARIATES #

        samples = intersect(rownames(pheno), rownames(probs))
        samples = intersect(samples, rownames(addcovar))
        samples = intersect(samples, rownames(K[[1]]))
        stopifnot(length(samples) > 0)
        print(paste("A total of", length(samples), tx, "samples are complete."))

        pheno = pheno[samples,,drop = FALSE]
        addcovar = addcovar[samples,,drop = FALSE]
        probs = probs[samples,,,drop = FALSE]

        file.prefix = paste(tx, pheno.col, sep = "_")
        print(file.prefix)
        plot.title = paste(tx, pheno.col, sep = " ")
        print(plot.title)

        # COX PH MODEL #
        surv = Surv(pheno$days, pheno[,pheno.col])
        fit = survfit(surv ~ addcovar)
        plot(fit, col = 1:2, las = 1, main = plot.title)
        legend("bottomleft", col = 1:2, lty = 1, legend = c("female", "male"))
        mod = coxph(surv ~ addcovar)
        text(x = 25, y = 0.15, labels = paste("p =", format(anova(mod)[2,4],
                                                            digits = 2)), adj = 0)
        #print(paste(round(100*(sum(trait) / length(samples)), digits = 1),
        #            "% display the", pheno.col, "phenotype in the", tx, "group."))

        # REGRESSION MODEL #

        for(i in 1:length(K)) {
                K[[i]] = K[[i]][samples, samples]
        } # for(i)

        chrs = c(1:19, "X")
        data = vector("list", length(chrs))
        names(data) = chrs
        for(i in 1:length(chrs)) {

                rng = which(markers[,2] == chrs[i])
                data[[i]] = list(probs = probs[,,rng], K = K[[i]],
                                 markers = markers[rng,])

        } # for(i)

        rm(probs, K, markers)

        setwd = outdir

        # MAPPING ANALYSES #

        result = vector("list", length(data))
        names(result) = names(data)
        print(paste("Mapping with", length(samples), tx, "samples..."))

        for(i in 1:19) {
                print(paste("CHROMOSOME", i))
                result[[i]] = GRSDcoxph(data[[i]], pheno, surv, addcovar, tx)
        } #for(i)

        print("X CHROMOSOME")
        result[["X"]] = GRSDcoxph.xchr(data[["X"]], pheno, surv, addcovar, tx)

        print(paste(round(difftime(Sys.time(), begin, units = 'hours'), digits = 2),
                    "hours elapsed during mapping."))

        # Convert to GRangesList for storage
        chrs = c(1:19, "X")
        qtl = GRangesList(GRanges("list", length(result)))

        for(i in 1:length(chrs)) {
                print(i)
                qtl[[i]] <- GRanges(seqnames = Rle(result[[i]]$CHR),
                                    ranges = IRanges(start = result[[i]]$POS, width = 1),
                                    p.value = result[[i]]$pv)
        } # for(i)

        save(qtl, file.prefix, file = paste0(file.prefix, "_QTL.Rdata"))

}
