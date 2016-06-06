#' @title Effect size and odds ratio for SNP profiles
#' @author Elijah F Edmondson, \email{elijah.edmondson@@gmail.com}
#' @export

get.effect.size = function(pheno = pheno, pheno.col, chr, probs = probs, sdp.file = "~/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz",
                           markers, dir = "/Users/elijah/Desktop/R/QTL/WD/2.\ Binomial\ Mapping/")
{
        library(Rsamtools)

        #Create the matrix in which all data will be stored
        EFFECT = matrix(0, nrow = 1000, ncol = 8, dimnames = list(1:1000, c("PHENOTYPE", "CHR", "SNP", "LOD", "Effect Size", "ODDS", "2.5% ODDS", "97.5% ODDS")))
        #Enter the directory of QTL files
        files <- (Sys.glob(paste0(dir,"*.Rdata")))
        # Create a matrix of SDPs.
        sdp.mat = matrix(as.numeric(intToBits(1:2^8)), nrow = 32)
        sdp.mat = sdp.mat[8:1,]
        dimnames(sdp.mat) = list(LETTERS[1:8], 1:2^8)

        for(j in 1:length(files)){

                load(file = files[j])
                print(files[j])

                qtl = qtl[[chr]]
                qtl = as.data.frame(qtl)
                qtl = qtl[match(markers$Mb_NCBI38, qtl$start, nomatch=0),]
                stopifnot(length(qtl) > 0)

                # Get the SNP at the minimum p-value.
                SNP = qtl$start[which.min(qtl$p.value)]
                LOD = -log10(qtl$p.value[which(qtl$start == SNP)])

                # Read in the unique SDPs.
                tf = TabixFile(file = sdp.file)
                sdps = scanTabix(file = sdp.file, param = GRanges(seqnames = chr, ranges = SNP))[[1]]
                sdps = strsplit(sdps, split = "\t")
                sdps = matrix(unlist(sdps), ncol = 3, byrow = T)
                chr  = sdps[1,1]
                pos  = as.numeric(sdps[,2])
                sdps = as.numeric(sdps[,3])

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
                geno = geno[samples,,drop = FALSE]
                addcovar = addcovar[samples,,drop = FALSE]
                probs = probs[samples,,,drop = FALSE]


                mod0 = glm(pheno[,pheno.col] ~ addcovar, family = binomial(logit))
                mod0

                mod1 = glm(pheno[,pheno.col] ~ addcovar + geno[,1], family = binomial(logit))
                mod1
                anova(mod0,mod1,test = "Chisq")

                summary(mod1)
                print("Odds of developing tumoras a function of genotype:")
                odds = exp(coef(mod1))
                oddsCI = exp(confint.default(mod1))


                # Get RSS / SST for the genotypes.
                effect.size = anova(mod)[2,2] / sum(anova(mod)[,2])

                #perc.var = 100 * (1.0 - (ss / ss.null))

                EFFECT[i,] = c(pheno.col, chr, SNP, LOD, effect.size, odds[3], oddsCI[3], oddsCI[6])
                print(EFFECT)
                write.csv(LODmat, file = paste0(files[j], "QTL", ".csv"))
                rm(qtl)
        }


}
