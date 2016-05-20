#' Association mapping with survival outcomes using the CoxPH model.
#' @author Elijah F Edmondson, \email{elijah.edmondson@@gmail.com}
#'
#' Performs association mapping in multiparent mouse populations.
#' @export

GRSDcoxph.xchr = function(obj, pheno, pheno.col, surv, addcovar, tx, sanger.dir) {

        chr = obj$markers[1,2]

        strains = sub("/", "_", hs.colors[,2])

        file.prefix = paste("CoxPH", tx, pheno.col, sep = "_")
        plot.title = paste("CoxPH", tx, pheno.col, sep = " ")

        load(file = paste0(sanger.dir, "X.Rdata"))

        null.mod = coxph(surv ~ addcovar)
        null.ll = logLik(null.mod)
        pv = rep(0, nrow(sanger))

        # Get the unique SDPs between each pair of markers and
        # calculate the COXPH LOD.

        # CoxPH function.
        coxph.fxn = function(snp.rng, local.probs) {

                # Get the SDPs.
                sdp.nums = sanger[snp.rng,] %*% 2^(7:0)
                sdps2keep = which(!duplicated(sdp.nums))
                cur.sdps = sanger[snp.rng,,drop = FALSE][sdps2keep,,drop = FALSE]
                unique.sdp.nums = sdp.nums[sdps2keep]
                m = match(sdp.nums, unique.sdp.nums)

                # Multiply the SDPs by the haplotype probabilities.
                cur.alleles = tcrossprod(cur.sdps, local.probs)
                cur.ll = rep(null.ll, nrow(cur.sdps))

                # Check for low allele frequencies and remove SDPs with too
                # few samples carrying one allele.
                sdps.to.use = which(rowSums(cur.alleles) > 2.0)

                sex.col = which(colnames(addcovar) == "sex")
                if(length(sex.col) != 1) {
                        stop("One of the columns of addcovar MUST be named 'sex'.")
                } # if(length(sex.col) != 1)

                # Run the Cox PH model at each unique SDP.
                for(j in sdps.to.use) {

                        ###DMG
                        ### Put the logistic regression or linear model here.

                        # For the X chromosome we map with sex as an interactive
                        # covariate with genotype.
                        mod = coxph(surv ~ addcovar + cur.alleles[j,] +
                                            addcovar[,sex.col] * cur.alleles[j,])
                        cur.ll[j] = logLik(mod)

                } # for(j)

                # This is the LRS.
                cur.ll = cur.ll - null.ll

                # Return the results.
                cur.ll[m]

        } # coxph.fxn()

        # SNPs before the first marker.
        snp.rng = which(sanger.hdr$POS <= obj$markers[1,3])
        if(length(snp.rng) > 0) {

                pv[snp.rng] = coxph.fxn(snp.rng, obj$probs[,,1])

        } # if(length(snp.rng) > 0)

        # SNPs between Markers.
        for(i in 1:(nrow(obj$markers)-1)) {

                snp.rng = which(sanger.hdr$POS > obj$markers[i,3] &
                                        sanger.hdr$POS <= obj$markers[i+1,3])

                if(length(snp.rng) > 0) {

                        # Take the mean of the haplotype probs at the surrounding markers.
                        pv[snp.rng] = coxph.fxn(snp.rng, (obj$probs[,,i] +
                                                                  obj$probs[,,i+1]) * 0.5)

                } # if(length(snp.rng) > 0)

        } # for(i)

        # SNPs after the last marker.
        snp.rng = which(sanger.hdr$POS > obj$markers[nrow(obj$markers),3])
        if(length(snp.rng) > 0) {

                pv[snp.rng] = coxph.fxn(snp.rng, obj$probs[,,nrow(obj$markers)])

        } # if(length(snp.rng) > 0)

        # Convert LRS to p-values using the chi-squared distribution.
        # Note that we have more degrees of freedom in the model.
        pv = pchisq(2 * pv, df = 1, lower.tail = FALSE)
        pv = data.frame(sanger.hdr, pv, stringsAsFactors = FALSE)

        save(pv, file = paste0(file.prefix, "_chr", chr, ".Rdata"))

        png(paste0(file.prefix, "_chr", chr,".png"), width = 2000,
            height = 1600, res = 200)
        plot(as.numeric(pv[,3]) * 1e-6, -log10(pv[,6]), pch = 20)
        mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr", chr))
        dev.off()

        # Return the positions and p-values.
        return(pv)

} # GRSDcoxph.xchr()
