#' @title Permutation analysis for continuous traits in HZE ion exposed genetically heterogeneous mice.
#' @author Elijah F. Edmondson <elijah.edmondson@gmail.com>
#' @export

Scanone.assoc.perms = function(perms, chr = 1:19, Xchr = FALSE,
                           pheno, pheno.col, probs, K, addcovar,
                           markers, snp.file, outdir = "~/Desktop/files",
                           tx = "", ncl = 4) {
        begin <- Sys.time()


        females = which(pheno$sex == "0")
        males = which(pheno$sex == "1")

        perms = matrix(1, nrow = perms, ncol = 2, dimnames = list(1:perms, c("A", "X")))

        for(p in 1:perms) {
                print(p)
                LODtime = Sys.time()
                new.order = rep(0, length(trait))
                new.order[females] = sample(females)
                new.order[males] = sample(males)

                log.perm = trait[new.order]
                trait = log.perm

                phenonew = data.frame(cbind("sex" = pheno$sex, trait))

                min.a.pv = 1

                qtl = scanone.assoc(pheno = phenonew, pheno.col = 11, probs = model.probs, K = K, addcovar, markers = MM_snps, sdp.file = sdp.file, ncl = 4)

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

                perms[p,] = c(-log10(min.a.pv), -log10(min(qtl$`X`@elementMetadata$p.value)))

                print(paste("Max random autosomal LOD was", -log10(min.a.pv), "X chr", -log10(min.x.pv)))
                print(paste(round(difftime(Sys.time(), LODtime, units = 'mins'), digits = 2),
                            "minutes..."))


        }
        print(paste(round(difftime(Sys.time(), begin, units = 'hours'), digits = 2),
                    "hours elapsed during analysis"))

        save(permutations, file.prefix, file = paste0(file.prefix, "_perms.Rdata"))
        return(permutations)


} #Scanone.assoc.perms
