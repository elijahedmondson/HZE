#' @title Genome wide association plots with separate colors
#'
#' param bin.width = size of loop, in which maximum LOD value for a set of bin of SNPs is recorded
#' param Red == HZE, Blue == Gamma, Purple == All Irradiated, Green == Unirradiated, Black == All mice
#'
#' author Elijah F Edmondson, elijah.edmondson@@gmail.com
#' @export

plot.hs.color.qtl = function(qtl, bin.width = 1000, color = "black", ...) {

        new.qtl = NULL
        for(chr in 1:length(qtl)) {
                library(BSgenome.Mmusculus.UCSC.mm10)
                print(chr)

                # Create 100 SNP bins.
                brks = cut(x = 1:length(qtl[[chr]]), breaks = length(qtl[[chr]]) / bin.width)
                # Split up the SNP positions and get the mean.
                pos = split(start(qtl[[chr]]), brks)
                pos = sapply(pos, mean)
                # Split up the p-values and get the max.
                pv = split(mcols(qtl[[chr]])$p.value, brks)
                pv = sapply(pv, min)

                # Make a single new GRanges object to return.
                gr = GRanges(seqnames = seqnames(qtl[[chr]])[1],
                             ranges = IRanges(start = pos, width = 1), p.value = pv)

                if(chr == 1) {
                        new.qtl = gr
                } else {
                        new.qtl = c(new.qtl, gr)
                } # else

        } # for(chr)

        # Get the chromosome lengths.
        chrlen = seqlengths(BSgenome.Mmusculus.UCSC.mm10)
        names(chrlen) = sub("^chr", "", names(chrlen))
        chrlen = chrlen[seqlevels(new.qtl)] * 1e-6

        # Add the chr lengths to the chromosomes for plotting.
        # Switch positions to genome Mb.
        gmb = start(new.qtl) * 1e-6
        for(chr in 2:length(chrlen)) {

                wh = which(seqnames(new.qtl) == names(chrlen)[chr])
                gmb[wh] = gmb[wh] + sum(chrlen[1:(chr - 1)])

        } # for(chr)

        # Get chromosome mid-points for plotting the Chr name.
        chrmid = (chrlen / 2) + cumsum(c(1, chrlen[-length(chrlen)]))

        # Make the plot.
        if(color == "purple") {
                col = rep(rgb(.4, .2, .6), length(new.qtl))
                even.chr = which(seqnames(new.qtl) %in% (1:10 * 2))
                col[even.chr] = rgb(.549, .325, .776)
                plot(gmb, -log10(new.qtl$p.value), pch = 20, xaxt = "n",
                     col = col, las = 1, xlab = "", ylab = "-log10(p-value)", ...)
                mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)

        }
        if(color == "red") {
                col = rep(rgb(.859, .169, .239), length(new.qtl))
                even.chr = which(seqnames(new.qtl) %in% (1:10 * 2))
                col[even.chr] = rgb(.902, .424, .475)
                plot(gmb, -log10(new.qtl$p.value), pch = 20, xaxt = "n",
                     col = col, las = 1, xlab = "", ylab = "-log10(p-value)", ...)
                mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)
        }
        if(color == "blue") {
                col = rep(rgb(.122, .471, .706), length(new.qtl))
                even.chr = which(seqnames(new.qtl) %in% (1:10 * 2))
                col[even.chr] = rgb(.255, .624, .871)
                plot(gmb, -log10(new.qtl$p.value), pch = 20, xaxt = "n",
                     col = col, las = 1, xlab = "", ylab = "-log10(p-value)", ...)
                mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)
        }
        if(color == "green") {
                col = rep(rgb(.2, .627, .173), length(new.qtl))
                even.chr = which(seqnames(new.qtl) %in% (1:10 * 2))
                col[even.chr] = rgb(.325, .808 ,.294)
                plot(gmb, -log10(new.qtl$p.value), pch = 20, xaxt = "n",
                     col = col, las = 1, xlab = "", ylab = "-log10(p-value)", ...)
                mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)

        }
        if(color == "black") {
                col = rep(rgb(0,0,0), length(new.qtl))
                even.chr = which(seqnames(new.qtl) %in% (1:10 * 2))
                col[even.chr] = rgb(0.7,0.7,0.7)
                plot(gmb, -log10(new.qtl$p.value), pch = 20, xaxt = "n",
                     col = col, las = 1, xlab = "", ylab = "-log10(p-value)", ...)
                mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)

        }
        return(new.qtl)

} # plot.hs.color.qtl
