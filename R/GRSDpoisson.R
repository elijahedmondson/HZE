#' @title Association mapping on autosomal chromosomes with binary variable outcomes.
#'
#' Performs association mapping in multiparent mouse populations.
#' @export

GRSDpoisson = function(obj, pheno, pheno.col, tx) {

        chr = obj$markers[1,2]

        setwd(outdir)

        file.prefix = paste(tx, pheno.col, "poisson", sep = "_")

        plot.title = paste(tx, pheno.col, sep = " ")

        strains = sub("/", "_", hs.colors[,2])

        hdr = scanVcfHeader(snp.file)
        gr = GRanges(seqnames = chr, range = IRanges(start = 0,
                                                     end = 200e6))
        param = ScanVcfParam(geno = c("GT", "FI"), fixed = "ALT",
                             samples = strains[strains != "C57BL_6J"], which = gr)
        sanger = readVcf(file = snp.file, genome = "mm10", param = param)

        # Keep high quality SNPs (quality == 1)
        sanger = sanger[rowSums(geno(sanger)$FI, na.rm = TRUE) == 7]

        # Keep polymorphic SNPs.
        keep = which(rowSums(geno(sanger)$GT == "0/0", na.rm = TRUE) < 7)
        sanger = sanger[keep]
        rm(keep)

        # We have to do some work to extract the alternate allele.
        alt = CharacterList(fixed(sanger)$ALT)
        alt = unstrsplit(alt, sep = ",")


        sanger.hdr = data.frame(ID = names(rowRanges(sanger)), CHR = as.character(seqnames(sanger)),
                                POS = start(sanger), REF = as.character(fixed(sanger)$REF),
                                ALT = alt, stringsAsFactors = FALSE)
        rm(alt)


        sanger = cbind(geno(sanger)$GT[,1:4,drop = FALSE],
                       "C57BL_6J" = "0/0",
                       geno(sanger)$GT[,5:7,drop = FALSE])

        sanger = (sanger != "0/0") * 1

        # Make the MAF between 1/8 and 4/8.
        flip = which(rowSums(sanger) > 4)
        sanger[flip,] = 1 - sanger[flip,,drop = FALSE]
        rm(flip)

        #null.mod = glm(pheno[,pheno.col] ~ addcovar, family = binomial(logit))
        null.mod = glm(pheno[,pheno.col] ~ addcovar, family = poisson(link = "log"))
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


                        #full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = binomial(logit))
                        full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = poisson(link = "log"))
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

        save(pv, file = paste0(file.prefix, "_chr", chr, ".Rdata"))

        png(paste0(file.prefix, "_chr", chr,".png"), width = 2000,
            height = 1600, res = 200)
        plot(as.numeric(pv[,3]) * 1e-6, -log10(pv[,6]), pch = 20)
        mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr", chr))
        dev.off()

        # Return the positions and p-values.
        return(pv)

} # GRSDpoisson()
