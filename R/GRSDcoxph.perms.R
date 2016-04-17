#' @title Permutation analysis for CoxPH model in HZE.
#' @author Elijah F Edmondson, \email{elijah.edmondson@@gmail.com}
#' @export

GRSDcoxph.perms = function(perms, chr = 1:19, Xchr = FALSE,
                           pheno, pheno.col, days.col, probs, K, addcovar,
                           markers, snp.file, outdir = "/home/ubuntu",
                           tx = "", sanger.dir = "home/ubuntu/HS.sanger.files/") {
        begin <- Sys.time()
        begin

        samples = intersect(rownames(pheno), rownames(probs))
        samples = intersect(samples, rownames(addcovar))
        samples = intersect(samples, rownames(K[[1]]))
        stopifnot(length(samples) > 0)

        pheno = pheno[samples,,drop = FALSE]
        addcovar = addcovar[samples,,drop = FALSE]
        probs = probs[samples,,,drop = FALSE]

        # DEFINE surv #

        file.prefix = paste(tx, pheno.col, sep = "_")

        plot.title = paste(tx, pheno.col, sep = " ")
        print(paste(plot.title, "Permutation Analysis:", Sys.time()))

        # COX PH MODEL #
        surv = Surv(pheno[,days.col], pheno[,pheno.col])
        fit = survfit(surv ~ addcovar)


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

        ##
        result = vector("list", length(data))
        names(result) = names(data)

        females = which(pheno$sex == "0")
        males = which(pheno$sex == "1")

        permutations = matrix(1, nrow = perms, ncol = 2, dimnames = list(1:perms, c("A", "X")))
        sanger.dir = sanger.dir
        for(p in 1:perms) {
                LODtime = Sys.time()
                print(p)
                new.order = rep(0, length(surv))
                new.order[females] = sample(females)
                new.order[males] = sample(males)
                log.perm = pheno[,pheno.col][new.order]


                pheno["new.col"] <- log.perm

                #surv = Surv(pheno[,days.col], pheno[,new.col])

                min.a.pv = 1

                for(i in 1:length(chr)) {
                        result = GRSD.coxph4perms(data[[i]], pheno = pheno, pheno.col = "new.col", days.col = days.col, addcovar, tx, sanger.dir)
                        min.a.pv = min(min.a.pv, min(result$pv))
                } #for(i)

                min.a.pv = min(min.a.pv, min(result$pv))
                min.x.pv = 1

                if(Xchr) {
                        result = GRSDcoxph.xchr4perms(data[["X"]], pheno = pheno, pheno.col = "new.col", days.col, addcovar, tx, sanger.dir)
                        min.x.pv = min(result$pv)
                }

                # Save the minimum p-values.
                permutations[p,] = c(-log10(min.a.pv), -log10(min.x.pv))
                print(paste("Max randomized LOD was", -log10(min.a.pv), -log10(min.x.pv)))
                print(paste(round(difftime(Sys.time(), LODtime, units = 'mins'), digits = 2),
                            "minutes..."))

        }
        print(paste(round(difftime(Sys.time(), begin, units = 'hours'), digits = 2),
                    "hours elapsed during permutation analysis"))

        save(permutations, file.prefix, file = paste0(file.prefix, "_perms.Rdata"))
        return(permutations)


} #GRSDcoxph.perms
