#' @title Association mapping on autosomal chromosomes with binary variable outcomes optimized for speed by
#' pre-formating sanger snps to sdp files.
#'
#' @author Elijah F Edmondson, \email{elijah.edmondson@@gmail.com}
#' Performs association mapping in multiparent mouse populations.
#' @export

GRSDbinom.permsfast = function(obj, pheno, pheno.col, addcovar, tx, sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/") {
        chr = obj$markers[1,2]

        setwd(outdir)

        file.prefix = paste(tx, pheno.col, sep = "_")

        plot.title = paste(tx, pheno.col, sep = " ")

        strains = sub("/", "_", hs.colors[,2])

        load(file = paste0(sanger.dir, chr, ".Rdata"))

        null.mod = glm(pheno[,pheno.col] ~ addcovar, family = binomial(logit))
        #null.mod = glm(trait ~ addcovar, family = poisson(link = "log"))
        null.ll = logLik(null.mod)
        pv = rep(0, nrow(sanger))

        glm.fxn = function(snp.rng, local.probs) {

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
                sdps.to.use = which(rowSums(cur.alleles) > 1.0)

                # Run the model at each unique SDP.
                for(j in sdps.to.use) {


                        full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = binomial(logit))
                        #full.mod = glm(trait ~ addcovar + cur.alleles[j,], family = poisson(link = "log"))
                        cur.ll[j] = logLik(full.mod)

                } # for(j)

                # This is the LRS.
                cur.ll = cur.ll - null.ll

                # Return the results.
                cur.ll[m]

        } # glm.fxn()

        # SNPs before the first marker.
        snp.rng = which(sanger.hdr$POS <= obj$markers[1,3])
        if(length(snp.rng) > 0) {

                pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,1])

        } # if(length(snp.rng) > 0)

        # SNPs between Markers.
        for(i in 1:(nrow(obj$markers)-1)) {

                snp.rng = which(sanger.hdr$POS > obj$markers[i,3] &
                                        sanger.hdr$POS <= obj$markers[i+1,3])

                if(length(snp.rng) > 0) {

                        # Take the mean of the haplotype probs at the surrounding markers.
                        pv[snp.rng] = glm.fxn(snp.rng, (obj$probs[,,i] +
                                                                obj$probs[,,i+1]) * 0.5)

                } # if(length(snp.rng) > 0)

        } # for(i)

        # SNPs after the last marker.
        snp.rng = which(sanger.hdr$POS > obj$markers[nrow(obj$markers),3])
        if(length(snp.rng) > 0) {

                pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,nrow(obj$markers)])

        } # if(length(snp.rng) > 0)

        # Convert LRS to p-values using the chi-squared distribution.
        pv = pchisq(2 * pv, df = 1, lower.tail = FALSE)
        pv = data.frame(sanger.hdr, pv, stringsAsFactors = FALSE)



        # Return the positions and p-values.
        return(pv)


} # GRSDbinom.permsfast()
