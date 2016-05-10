#' @title QTL Confidence Interval with Resample Model Averaging for Cox Regression
#'
#' For experimental cross data on n individuals, one makes n draws, with replacement,
#' from the observed individuals to form a new data set in which some individuals are
#' omitted and some appear multiple times. An estimate of QTL location is calculated
#' with these new data, and the process is repeated many times. A ∼95% confidence
#' interval for the location of the QTL is obtained as the interval containing 95%
#' of the estimated locations from the bootstrap replicates. Visscher et al. (1996)
#'
#'@param Window, determines the evaluation window around each peak (Mb)
#'
#' Arguments:
#'
#'
#' @export


HS.cox.bootstrap = function(perms, chr, pheno, pheno.col, days.col, probs, K, addcovar,
                              markers, snp.file, outdir = "~/Desktop/", peakMB, window = 8000000,
                              tx = "", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")
{
        begin <- Sys.time()
        file.prefix = paste0("Bootstrap.", tx, ".", pheno.col, "(Chr", chr, ")")
        plot.title = paste("Bootstrap ", tx, " ", pheno.col, " ", "(Chr ", chr, "):", sep = "")
        print(paste(plot.title, Sys.time()))

        load(file = paste0(sanger.dir, "probs", chr, ".Rdata"))

        samples = intersect(rownames(pheno), rownames(probs))
        samples = intersect(samples, rownames(addcovar))
        samples = intersect(samples, rownames(K[[1]]))
        stopifnot(length(samples) > 0)

        pheno = pheno[samples,,drop = FALSE]
        addcovar = addcovar[samples,,drop = FALSE]
        probs = probs[samples,,,drop = FALSE]

        #samples2 = markers$SNP_ID[which(markers$Chr == chr & markers$Mb_NCBI38 > (peakMB - 4000000) & markers$Mb_NCBI38 < (peakMB + 4000000))]
        samples2 = markers$SNP_ID[which(markers$Chr == chr)]
        probs = probs[,,samples2, drop = FALSE]

        K = K[[chr]][samples,samples,drop = FALSE]

        setwd(outdir)

        ##

        permutations = matrix(1, nrow = perms, ncol = 5, dimnames = list(1:perms, c("LOD", "min", "max", "average", "#Markers")))

        rm(samples, samples2)
        for(p in 1:perms) {
                LODtime = Sys.time()
                print(p)

                repeat {
                        phenoperm = data.frame(pheno[sample(nrow(pheno), replace = TRUE), ],
                                               check.names = FALSE, check.rows = FALSE)
                        sanger.dir = sanger.dir
                        ### FROM DG ###

                        samples = sub("\\.[0-9]$", "", rownames(phenoperm))
                        probsperm = probs[samples,,]
                        rownames(probsperm) = make.unique(rownames(probsperm))


                        Kperm = K
                        Kperm = K[samples,samples,drop = FALSE]
                        rownames(Kperm) = make.unique(rownames(Kperm))


                        ### Move the model into the loop REGRESSION MODEL ###

                        Kperm = Kperm[samples, samples]
                        chrs = chr
                        data = vector("list", length(chrs))
                        names(data) = chrs

                        rng = which(markers[,2] == chrs)
                        data = list(probsperm = probsperm, Kperm = Kperm,
                                    markers = markers[rng,])

                        result = vector("list", length(data))
                        names(result) = names(data)
                        rm(probsperm, Kperm)



                        ### WORK HORSE ###
                        ### ONLY RETURN POS WITHIN 8 MB WINDOW AROUND PEAK ###


                        result = GRSD.coxph4perms(data, chr = chr, pheno = phenoperm, pheno.col, days.col, addcovar, tx, sanger.dir)

                        top = max(-log10(result$pv))
                        MAX.LOD = result$POS[which(-log10(result$pv) == top)]
                        MegaBase = (min(MAX.LOD) + max(MAX.LOD))/2000000

                        print(paste((round(MegaBase, digits = 2)), "Mb: LOD", round(top, digits = 2)))
                        if (MAX.LOD > (peakMB - (window/2)) & MAX.LOD < (peakMB + (window/2))) break

                }



                print(paste0("Accepted locus: ", (round(MegaBase, digits = 2)), " Mb"))



                # Save the locus.
                permutations[p,] = c(top, min(MAX.LOD), max(MAX.LOD), ((min(MAX.LOD) + max(MAX.LOD))/2), length(MAX.LOD))
                print(paste(round(difftime(Sys.time(), LODtime, units = 'mins'), digits = 2),
                            "minutes..."))

        }

        print(paste(round(difftime(Sys.time(), begin, units = 'hours'), digits = 2),
                    "hours elapsed during analysis"))

        quant = quantile(permutations[,4], c(0.025,0.975))
        print(paste("95% Confidence Interval for QTL:", min(quant), "-", max(quant)))
        print(paste("Interval =", (round((((max(quant) - min(quant))/1000000)), digits = 2)), "Mb"))

return(permutations)


} #HS.cox.bootstrap
