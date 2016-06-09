#' @title Effect size and odds ratio for SNP profiles
#' @author Elijah F Edmondson, \email{elijah.edmondson@@gmail.com}
#' @export

get.effect.size.coxPH = function(pheno = pheno, pheno.col, days.col = days, probs = probs, sdp.file = "~/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz",
                                 markers, threshold = 5.05, dir = "/Users/elijah/Desktop/R/QTL/WD/2.\ Binomial\ Mapping/")
{
        load("/Users/elijah/Desktop/R/QTL/WD/hs.colors.Rdata")
        #Enter the directory of QTL files
        files <- (Sys.glob(paste0(dir,"*.Rdata")))

        # Create a matrix of SDPs.
        sdp.mat = matrix(as.numeric(intToBits(1:2^8)), nrow = 32)
        sdp.mat = sdp.mat[8:1,]
        dimnames(sdp.mat) = list(LETTERS[1:8], 1:2^8)

        #helper function from DG
        get.genotype = function(chr, pos, snp, markers, probs) {

                # Convert the SNP to numbers.
                snp = unlist(snp)
                names(snp) = make.names(sub("_", ".", names(snp)))
                strains = make.names(hs.colors[,2])

                # Get the slices from the haplotype probs matrix.
                markers = markers[markers[,1] %in% dimnames(probs)[[3]],]
                probs = probs[,,dimnames(probs)[[3]] %in% markers[,1]]
                markers = markers[markers[,2] == chr,]
                probs = probs[,,markers[,1]]
                markers = markers[max(which(markers[,3] < pos)):min(which(markers[,3] > pos)),]

                # Get the probs for these markers.
                probs = probs[,,markers[,1], drop = FALSE]
                probs = apply(probs, 1:2, mean)

                # Multiply the two matrices and return the result.
                return(probs %*% snp)

        } # get.genotype()

        for(j in 1:length(files)){
                #Create the matrix in which all data will be stored
                EFFECT = matrix(0, nrow = 20, ncol = 10,
                                dimnames = list(1:20, c("PHENOTYPE", "CHR", "SNP", "LOD", "Logrank P",
                                                        "BB Hazard", "BB Hazard Ratio", "SE BB HR", "Pr(>|z|)", "Rsquare")))
                load(file = files[j])
                print(files[j])
                for(i in 1:19) {
                        tryCatch({

                                # Determine most significant SNP and LOD score on the chromosome of interest.
                                SNP = qtl[[i]]@ranges@start[which.min(qtl[[i]]@elementMetadata@listData$p.value)]
                                LOD = -log10(qtl[[i]]@elementMetadata@listData$p.value[which(qtl[[i]]@ranges@start == SNP)])


                                # Run the loop for all significant loci
                                if(LOD > threshold) {
                                        print(i)
                                        # Read in the unique SDPs.
                                        tf = TabixFile(file = sdp.file)
                                        sdps = scanTabix(file = sdp.file, param = GRanges(seqnames = i, ranges = SNP))[[1]]
                                        sdps = strsplit(sdps, split = "\t")
                                        sdps = matrix(unlist(sdps), ncol = 3, byrow = T)
                                        chr  = sdps[1,1]
                                        pos  = as.numeric(sdps[,2])
                                        sdps = as.numeric(sdps[,3])

                                        geno = get.genotype(chr = chr,
                                                            pos = pos,
                                                            snp = sdp.mat[,sdps],
                                                            markers = markers,
                                                            probs = probs)
                                        geno = round(geno, digits = 1)
                                        geno = ifelse(geno < 0.25, "AA",
                                                ifelse(geno >=.25 & geno <= 0.75, "AB",
                                                ifelse(geno > .75, "BB",
                                                       NA)))


                                        # Fit the model.
                                        samples = intersect(rownames(pheno), rownames(probs))
                                        samples = intersect(samples, rownames(addcovar))
                                        samples = intersect(samples, rownames(geno))
                                        stopifnot(length(samples) > 0)
                                        pheno = pheno[samples,,drop = FALSE]
                                        addcovar = addcovar[samples,,drop = FALSE]
                                        geno = geno[samples,,drop = FALSE]
                                        addcovar = addcovar[samples,,drop = FALSE]
                                        #probs = probs[samples,,,drop = FALSE]

                                        ###### CoxPH model ######

                                        # Measures of explained variation, such as the coefficient of determination (R2) in linear models,
                                        # are helpful in assessing the explanatory power of a model. In survival analysis, these measures
                                        # help quantify the ability of prognostic factors to predict a patient's time until death. As in
                                        # linear models, covariates in Cox regression may be statistically significant but still have very
                                        # little predictive power. In the censored data setting, the definition of such a measure is not
                                        # straightforward; several measures of explained variation have been proposed. The most popular of
                                        # these is the generalized R-squared, calculated as 1-exp((χLR2)/n), where (χLR2) is the chi-square
                                        # statistic for the likelihood ratio test for the overall model, and n is the total number of
                                        # patients. Although the generalized R-squared is commonly recommended for the Cox model, its
                                        # sensitivity to the proportion of censored values is not often mentioned. In fact, the expected
                                        # value of R-squared decreases substantially as a function of the percent censored, with early
                                        # censoring having a greater impact than later censoring. Simulations show that complete data
                                        # R-squared values from the Cox model are very close to those from a similar linear model. However,
                                        # average R-squared values can decrease by 20% or more (e.g., R-squared from 0.5 to 0.4) with heavy
                                        # censoring (e.g., 50% censoring) compared to complete data. Simulation results will be presented,
                                        # and alternatives to the generalized R-squared will be discussed. SCHEMPER 1996

                                        surv = Surv(pheno[,days.col], pheno[,pheno.col])
                                        fit = survfit(surv ~ geno)
                                        mod = coxph(surv ~ geno)

                                        png(paste0(pheno.col, "_chr", chr,".png"), width = 2000,
                                            height = 1600, res = 250)
                                        plot(fit, col = 1:3, las = 1, main = paste0(pheno.col, ": chr ", chr, " bp ", SNP))
                                        legend("bottomleft", col = 1:3, lty = 1, legend = c("AA", "AB", "BB"))
                                        text(x = 5, y = 0.25, labels = paste("BB Hazard =", format(mod$coefficients[2], digits = 4),
                                                                             "(HR =", format(exp(mod$coefficients[2]), digits = 4), ")",
                                                                             "P.value =", format(anova(mod)[2,4], digits = 2)), adj = 0)
                                        dev.off()




                                        EFFECT[i,] = c(pheno.col, chr, SNP, LOD, format(anova(mod)[2,4]),
                                                       mod$coefficients[2], summary(mod)$coefficients["genoBB","exp(coef)"],
                                                       summary(mod)$coefficients["genoBB","se(coef)"],
                                                       summary(mod)$coefficients["genoBB","Pr(>|z|)"], summary(mod)$rsq[1])
                                        print(EFFECT[i,])
                                }
                }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})


                }
                write.csv(EFFECT, file = paste0(files[j], "QTL", ".csv"))
                print(EFFECT)
                rm(qtl, EFFECT)
        }


}#get.effect.size.coxPH()
