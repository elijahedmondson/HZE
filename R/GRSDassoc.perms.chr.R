
GRSDassoc.perms.chr = function(result, perms, chr, Xchr = TRUE,
                           pheno, pheno.col, addcovar, tx) {
        for(p in 1:perms) {
                print(p)
                new.order = rep(0, length(trait))
                new.order[females] = sample(females)
                new.order[males] = sample(males)

                log.perm = trait[new.order]
                trait = log.perm

                min.a.pv = 1


                print(chr)
                result = GRSDbinom(data[[chr]])
                min.a.pv = min(min.a.pv, min(result$pv))

                if(Xchr) {
                        print("X")
                        result = GRSDbinom.xchr(data[["X"]])
                        min.x.pv = min(result$pv)
                }
        # Save the minimum p-values.
        perms[p,] = c(-log10(min.a.pv), -log10(min.x.pv))

} # for(p)

