#' @title Permutation analysis for continuous traits in HZE ion exposed genetically heterogeneous mice.
#' @author Elijah F Edmondson, \email{elijah.edmondson@@gmail.com}
#' @export

Scanone.assoc.perms = function(perms, pheno = pheno, pheno.col, probs, K = K, tx = "",
                               addcovar = addcovar, markers = MM_snps, sdp.file = sdp.file, ncl = 4) {
        begin <- Sys.time()
        print(paste(tx, pheno.col, "Permutation Analysis:", Sys.time(), sep = " "))
        females = which(pheno$sex == "0")
        males = which(pheno$sex == "1")
        sdp.file = sdp.file
        ncl = ncl
        permutations = matrix(1, nrow = perms, ncol = 2, dimnames = list(1:perms, c("A", "X")))

        for(p in 1:perms) {
                print(p)
                LODtime = Sys.time()
                new.order = rep(0, length(pheno[,pheno.col]))
                new.order[females] = sample(females)
                new.order[males] = sample(males)

                X.frame = rep("X", length(new.order))

                row.names = data.frame(X.frame, new.order)
                row.names = apply(row.names, 1, paste, collapse="")

                phenoperm = data.frame(row.names = row.names, sex = as.numeric(pheno$sex == "1"),
                                       "2" = as.numeric(pheno[,pheno.col]))

                min.a.pv = 1

                qtl = scanone.assoc(pheno = phenoperm, pheno.col = 2, probs, K = K, addcovar,
                                    markers = MM_snps, sdp.file = sdp.file, ncl = 4)

                min.a.pv = min(min.a.pv, min(qtl$`1`@elementMetadata$p.value),
                               min(qtl$`2`@elementMetadata$p.value),
                               min(qtl$`3`@elementMetadata$p.value),
                               min(qtl$`4`@elementMetadata$p.value),
                               min(qtl$`5`@elementMetadata$p.value),
                               min(qtl$`6`@elementMetadata$p.value),
                               min(qtl$`7`@elementMetadata$p.value),
                               min(qtl$`8`@elementMetadata$p.value),
                               min(qtl$`9`@elementMetadata$p.value),
                               min(qtl$`10`@elementMetadata$p.value),
                               min(qtl$`11`@elementMetadata$p.value),
                               min(qtl$`13`@elementMetadata$p.value),
                               min(qtl$`14`@elementMetadata$p.value),
                               min(qtl$`15`@elementMetadata$p.value),
                               min(qtl$`16`@elementMetadata$p.value),
                               min(qtl$`17`@elementMetadata$p.value),
                               min(qtl$`18`@elementMetadata$p.value),
                               min(qtl$`19`@elementMetadata$p.value))
                min.x.pv = min(qtl$`X`@elementMetadata$p.value)

                permutations[p,] = c(-log10(min.a.pv), -log10(min(min.x.pv)))

                print(paste("Max random autosomal LOD was", -log10(min.a.pv), "X chr", -log10(min.x.pv)))
                print(paste(round(difftime(Sys.time(), LODtime, units = 'mins'), digits = 2),
                            "minutes..."))


        }
        print(paste(round(difftime(Sys.time(), begin, units = 'hours'), digits = 2),
                    "hours elapsed during analysis"))

        return(permutations)


} #Scanone.assoc.perms
