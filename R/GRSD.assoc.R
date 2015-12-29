#' @title Genome wide association mapping for binary phenotypes.
#' @author Elijah F. Edmondson <elijah.edmondson@gmail.com>
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


GRSD.assoc = function(pheno, pheno.col, probs, K, addcovar, markers, snp.file,
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

        file.prefix = paste(tx, pheno.col, sep = "_")

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
        print(paste("Mapping with", length(samples), tx, "samples..."))

        for(i in 1:19) {
                print(paste("CHROMOSOME", i))
                result[[i]] = GRSDbinom(data[[i]], pheno, pheno.col, addcovar, tx)
        } #for(i)

        print("X CHROMOSOME")
        result[["X"]] = GRSDbinom.xchr(data[["X"]], pheno, pheno.col, addcovar, tx)

        print(paste(round(difftime(Sys.time(), begin, units = 'hours'), digits = 2),
                    "hours elapsed during mapping."))


        # PLOTTING
        plotter <- Sys.time()

        setwd(outdir)
        files = dir(pattern = file.prefix)
        files = files[files != paste0(file.prefix, ".Rdata")]
        png.files = grep("png$", files)
        if(length(png.files) > 0) {
                files = files[-png.files]
        }
        num = gsub(paste0("^", file.prefix, "_chr|\\.Rdata$"), "", files)
        files = files[order(as.numeric(num))]

        data = vector("list", length(files))
        names(data) = num[order(as.numeric(num))]
        print("Plotting...")
        for(i in 1:length(files)) {

                load(files[i])
                data[[i]] = pv
                data[[i]][,6] = -log10(data[[i]][,6])

        } # for(i)

        num.snps = sapply(data, nrow)
        chrs = c(1:19, "X")

        xlim = c(0, sum(num.snps))
        ylim = c(0, max(sapply(data, function(z) { max(z[,6]) })))

        # PLOT ALL CHROMOSOMES #

        chrlen = get.chr.lengths()[1:20]
        chrsum = cumsum(chrlen)
        chrmid = c(1, chrsum[-length(chrsum)]) + chrlen * 0.5
        names(chrmid) = names(chrlen)

        png(paste0(file.prefix, "_QTL.png"), width = 2600, height = 1200, res = 200)
        plot(-1, -1, col = 0, xlim = c(0, max(chrsum)), ylim = ylim, xlab = "",
             ylab = "-log10(p-value)", las = 1, main = plot.title, xaxt = "n")
        for(i in 1:length(data)) {
                print(paste("Plotting chromosome", i))
                pos = data[[i]][,3] * 1e-6 + c(0, chrsum)[i]
                points(pos, data[[i]][,6], col = c("black", "grey50")[i %% 2 + 1],
                       pch = 20)
        } # for(i)
        mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.5)
        dev.off()


        save(result, file.prefix, file = paste0(file.prefix, "_QTL.Rdata"))

        print(paste(round(difftime(Sys.time(), plotter, units = 'hours'), digits = 2),
                    "hours elapsed during plotting."))

}
