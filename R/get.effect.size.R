#' @title Effect size and odds ratio for SNP profiles
#' @author Elijah F Edmondson, \email{elijah.edmondson@@gmail.com}
#' @export

get.effect.size = function(pheno = pheno, pheno.col, probs = probs, sdp.file = "~/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz",
                           markers, threshold = 5.05, dir = "/Users/elijah/Desktop/R/QTL/WD/2.\ Binomial\ Mapping/")
{
        library(Rsamtools)

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
                EFFECT = matrix(0, nrow = 20, ncol = 15,
                                dimnames = list(1:20, c("PHENOTYPE", "CHR", "SNP", "LOD", "ODDS",
                                                        "2.5% ODDS", "97.5% ODDS", "ANOVA Pr(>Chi)",
                                                        "AIC", "CoxSnell", "Nagelkerke", "McFadden",
                                                        "Tjur","sqPearson", "sqD")))
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
                                        
                                        mod0 = glm(pheno[,pheno.col] ~ addcovar, family = binomial(logit))
                                        #mod0
                                        mod1 = glm(pheno[,pheno.col] ~ addcovar + geno[,1], family = binomial(logit))
                                        #mod1
                                        #summary(mod1)
                                        
                                        
                                        ANOVA = anova(mod0,mod1,test = "Chisq")
                                        
                                        #print("Odds of developing tumor as a function of genotype:")
                                        odds = exp(coef(mod1))
                                        #odds
                                        oddsCI = exp(confint.default(mod1))
                                        #oddsCI
                                        
                                        #There are several ways of calculating (pseudo) R-squared values for logistic regression models,
                                        #with no consensus about which is best. The RsqGLM function, now included in the modEvA package,
                                        #calculates those of McFadden (1974), Cox & Snell (1989), Nagelkerke (1991), Tjur (2009), and
                                        #the squared Pearson correlation between observed and predicted values.
                                        R2 = RsqGLM(model = mod1)
                                        
                                        #Linear models come with an R-squared value that measures the proportion of variation that the
                                        #model accounts for. The R-squared is provided with summary(model) in R. For generalized linear
                                        #models (GLMs), the equivalent is the amount of deviance accounted for D-squared (Guisan &
                                        #Zimmermann 2000), but this value is not normally provided with the model summary. The Dsquared
                                        #function, now included in the modEvA package (Barbosa et al. 2014), calculates it. There is also
                                        #an option to calculate the adjusted D-squared, which takes into account the number of observations
                                        #and the number of model parameters, thus allowing direct comparison among different models
                                        #(Weisberg 1980, Guisan & Zimmermann 2000).
                                        D2 = Dsquared(model = mod1)
                                        
                                        EFFECT[i,] = c(pheno.col, chr, SNP, LOD, odds[3], oddsCI[3], oddsCI[6], ANOVA$`Pr(>Chi)`[2],
                                                       mod1$aic, R2$CoxSnell, R2$Nagelkerke, R2$McFadden, R2$Tjur, R2$sqPearson, D2)
                                        print(EFFECT[i,])
                                        rm(LOD, oddsCI, ANOVA, mod1, mod0, R2, D2, SNP)
                                }
                        }, error=function(e){print(ERROR)})

                        
                }
                write.csv(EFFECT, file = paste0(files[j], "QTL", ".csv"))
                print(EFFECT)
                rm(qtl, EFFECT)
        }


}
