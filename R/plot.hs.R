#' @title Plot function (w/ binning to average markers and max LOD)
#' @export

plot.hs = function(qtl, ...) {

        # Get the chromosome lengths.
        chrlen = seqlengths(BSgenome.Mmusculus.UCSC.mm10)
        names(chrlen) = sub("^chr", "", names(chrlen))
        chrlen = chrlen[seqlevels(qtl)] * 1e-6

        # Add the chr lengths to the chromosomes for plotting.
        # Switch positions to genome Mb.
        gmb = start(qtl) * 1e-6
        for(chr in 2:length(chrlen)) {

                wh = which(seqnames(qtl) == names(chrlen)[chr])
                gmb[wh] = gmb[wh] + sum(chrlen[1:(chr - 1)])

        } # for(chr)

        # Get chromosome mid-points for plotting the Chr name.
        chrmid = (chrlen / 2) + cumsum(c(1, chrlen[-length(chrlen)]))

        col = rep(rgb(0,0,0), length(qtl))
        even.chr = which(seqnames(qtl) %in% (1:10 * 2))
        col[even.chr] = rgb(0.7,0.7,0.7)
        plot(gmb, -log10(qtl$p.value), pch = 20, xaxt = "n",
             col = col, las = 1, xlab = "", ylab = "-log10(p-value)", ...)
        mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)


} # plot.hs
