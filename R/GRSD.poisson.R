#' @title Genome wide association mapping with glm(family = poisson).
#' @author Elijah F Edmondson, \email{elijah.edmondson@@gmail.com}
#' @details Dec. 29, 2015
#' @param pheno: data.frame containing phenotypes in columns and samples in
#'                   rows. Rownames must contain sample IDs.
#' @param pheno.col: the phenotype column to map.
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


GRSD.poisson = function(pheno, pheno.col, probs, K, addcovar, markers, snp.file,
                      outdir = "~/Desktop/", tx = c("Gamma", "HZE", "Unirradiated")){
        begin <- Sys.time()
        # COVARIATES #

        samples = intersect(rownames(pheno), rownames(probs))
        samples = intersect(samples, rownames(addcovar))
        samples = intersect(samples, rownames(K[[1]]))
        stopifnot(length(samples) > 0)
        print(paste("A total of", length(samples), tx, "samples are complete."))

        pheno = pheno[samples,,drop = FALSE]
        addcovar = addcovar[samples,,drop = FALSE]
        probs = probs[samples,,,drop = FALSE]



        # DEFINE TRAIT #

        file.prefix = paste(tx, pheno.col, "poisson", sep = "_")

        plot.title = paste(tx, pheno.col, sep = " ")
        print(plot.title)

        trait = pheno[,pheno.col]
        print(table(trait))
        print(paste(round(100*(sum(trait) / length(samples)), digits = 1),
                    "% display the", pheno.col, "phenotype in the", tx, "group."))

        # LOGISTIC REGRESSION MODEL #

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

        setwd(outdir)



        # MAPPING ANALYSES #

        result = vector("list", length(data))
        names(result) = names(data)
        print(paste("Mapping with", length(samples), tx, "samples."))

        for(i in 1:19) {
                print(i)
                result[[i]] = GRSDpoisson(data[[i]], pheno, pheno.col, tx)
        } #for(i)

        print("X")
        result[["X"]] = GRSDpoisson.xchr(data[["X"]])

        print(paste(round(difftime(Sys.time(), begin, units = 'hours'), digits = 2),
                    "hours elapsed during mapping."))

}
