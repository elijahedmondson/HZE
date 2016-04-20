#' For experimental cross data on n individuals, one makes n draws, with replacement,
#' from the observed individuals to form a new data set in which some individuals are
#' omitted and some appear multiple times. An estimate of QTL location is calculated
#' with these new data, and the process is repeated many times. A ∼95% confidence
#' interval for the location of the QTL is obtained as the interval containing 95%
#' of the estimated locations from the bootstrap replicates. Visscher et al. (1996)
#'
#' Arguments:
#'
#'
#' @export


HS.assoc.bootstrap = function(perms, chr, pheno, pheno.col, probs, K, addcovar,
                 markers, snp.file, outdir = "~/Desktop/",
                 tx = "", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")
{
                begin <- Sys.time()
                file.prefix = paste(tx, pheno.col, sep = "_")
                plot.title = paste(tx, pheno.col, sep = " ")
                print(paste(plot.title, "Bootstrap Analysis:", Sys.time()))

                samples = intersect(rownames(pheno), rownames(probs))
                samples = intersect(samples, rownames(addcovar))
                samples = intersect(samples, rownames(K[[1]]))
                stopifnot(length(samples) > 0)

                pheno = pheno[samples,,drop = FALSE]
                addcovar = addcovar[samples,,drop = FALSE]
                probs = probs[samples,,,drop = FALSE]



                # LOGISTIC REGRESSION MODEL #
                for(i in 1:length(K)) {
                        K[[i]] = K[[i]][samples, samples]
                } # for(i)

                chrs = c(chr)
                data = vector("list", length(chr))
                names(data) = chrs

                rng = which(markers[,2] == chrs[chr])
                data[[chr]] = list(probs = probs[,,rng], K = K[[chr]],
                                         markers = markers[rng,])


                data = pheno[data(nrow(data), replace = TRUE), ]

                rm(probs, K, markers)

                setwd(outdir)

                ##
                result = vector("list", length(data))
                names(result) = names(data)

                result = GRSDbinom.permsfast(data[[chr]], pheno = pheno, pheno.col = "Thyroid", addcovar, tx, sanger.dir)


                permutations = matrix(1, nrow = perms, ncol = 2, dimnames = list(1:perms, c("min", "max")))
                sanger.dir = sanger.dir
                for(p in 1:perms) {
                        LODtime = Sys.time()
                        print(p)

                        phenoperm = data.frame(pheno[sample(nrow(pheno), replace = TRUE), ], row.names = pheno$rownames, check.names = FALSE, check.rows = FALSE)

                        # FROM DG #
                        samples = sub("\\.[0-9]$", "", rownames(phenoperm))
                        probsperm = probs[samples,,]
                        rownames(probsperm) = make.unique(rownames(probsperm))


                        result = GRSDbinom.permsfast(data[[chr]], pheno = phenoperm, pheno.col = "pheno.col", addcovar, tx, sanger.dir)

                        top <- max(-log10(result$pv))
                        MAX.LOD = result$POS[which(-log10(result$pv) >= top)]

                        print(paste0("Maximum LOD score on Chr ", chr, " is ", top, ","))
                        print(paste("located between", max(MAX.LOD), "and ", min(MAX.LOD), "bp."))

                        #if(chr == "X") {
                                result = GRSDbinom.xchr.permsfast(data[["X"]], pheno = phenonew, pheno.col, addcovar, tx, sanger.dir)
                                min.x.pv = min(result$pv)
                        }

                        # Save the locus.
                        permutations[p,] = c(min(MAX.LOD), max(MAX.LOD))
                        print(paste(round(difftime(Sys.time(), LODtime, units = 'mins'), digits = 2),
                                    "minutes..."))

                }
                print(paste(round(difftime(Sys.time(), begin, units = 'hours'), digits = 2),
                            "hours elapsed during analysis"))

                return(permutations)
                QUANTILE = quantile(permutations,c(0.025,0.975))
                print("95% Confidence Interval for QTL:")
                print(paste(QUANTILE))
                return(QUANTILE)

}
