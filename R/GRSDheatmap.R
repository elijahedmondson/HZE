library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(DOQTL)
setwd("/Users/elijah/Desktop/R/QTL/WD/Heatmap/")
load(file = "~/Desktop/R/QTL/WD/__")

pv = split(mcols(result[[1]])$pv, brks)
pv = sapply(pv, min)
pos = split(start(result[[chr]]), brks)
pos = sapply(pos, mean)
gr = GRanges(seqnames = seqnames(result[[chr]])[1],
             ranges = IRanges(start = pos, width = 1), p.value = pv)









# Plot function (w/ binning to average markers and max LOD)
plot.hs.qtl = function(qtl, bin.width = 1000, ...) {

        new.qtl = NULL
        for(chr in 1:length(qtl)) {

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
        col = rep(rgb(0,0,0), length(new.qtl))
        even.chr = which(seqnames(new.qtl) %in% (1:10 * 2))
        col[even.chr] = rgb(0.7,0.7,0.7)
        plot(gmb, -log10(new.qtl$p.value), pch = 20, xaxt = "n",
             col = col, las = 1, xlab = "", ylab = "-log10(p-value)", ...)
        mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)

        return(new.qtl)

} # plot.hs.qtl

setwd("/Users/elijahedmondson/Desktop/R/QTL/WD/Heatmap/")

load(file ="/Users/elijah/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Ectoderm.Rdata")
qtl <- AM.qtl
rm(AM.qtl)
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.Ectoderm.heatmap.Rdata")
rm(qtl, qtl.smaller)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Endoderm.Rdata")
qtl <- AM.qtl
rm(AM.qtl)
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.Endoderm.heatmap.Rdata")
rm(qtl, qtl.smaller)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Mesoderm.Rdata")
qtl <- AM.qtl
rm(AM.qtl)
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.Mesoderm.heatmap.Rdata")
rm(qtl, qtl.smaller)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/Gamma/AMQTL.Ectoderm.Rdata")
qtl <- AM.qtl
rm(AM.qtl)
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.Ectoderm.heatmap.Rdata")
rm(qtl, qtl.smaller)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/Gamma/AMQTL.Endoderm.Rdata")
qtl <- AM.qtl
rm(AM.qtl)
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.Endoderm.heatmap.Rdata")
rm(qtl, qtl.smaller)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/Gamma/AMQTL.Mesoderm.Rdata")
qtl <- AM.qtl
rm(AM.qtl)
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.Mesoderm.heatmap.Rdata")
rm(qtl, qtl.smaller)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/Background/AMQTL.Ectoderm.Rdata")
qtl <- AM.qtl
rm(AM.qtl)
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Background.Ectoderm.heatmap.Rdata")
rm(qtl, qtl.smaller)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/Background/AMQTL.Endoderm.Rdata")
qtl <- AM.qtl
rm(AM.qtl)
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Background.Endoderm.heatmap.Rdata")
rm(qtl, qtl.smaller)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/Background/AMQTL.Mesoderm.Rdata")
qtl <- AM.qtl
rm(AM.qtl)
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Background.Mesoderm.heatmap.Rdata")
rm(qtl, qtl.smaller)




#Load all files for combination
load(file="/Users/elijahedmondson/Desktop/R/QTL/WD/Heatmap/Background.Ectoderm.heatmap.Rdata")
Background.Ectoderm <- qtl.smaller
rm(qtl.smaller)

load(file="/Users/elijahedmondson/Desktop/R/QTL/WD/Heatmap/Background.Endoderm.heatmap.Rdata")
Background.Endoderm <- qtl.smaller
rm(qtl.smaller)

load(file="/Users/elijahedmondson/Desktop/R/QTL/WD/Heatmap/Background.Mesoderm.heatmap.Rdata")
Background.Mesoderm <- qtl.smaller
rm(qtl.smaller)

load(file="/Users/elijahedmondson/Desktop/R/QTL/WD/Heatmap/HZE.Ectoderm.heatmap.Rdata")
HZE.Ectoderm <- qtl.smaller
rm(qtl.smaller)

load(file="/Users/elijahedmondson/Desktop/R/QTL/WD/Heatmap/HZE.Endoderm.heatmap.Rdata")
HZE.Endoderm <- qtl.smaller
rm(qtl.smaller)

load(file="/Users/elijahedmondson/Desktop/R/QTL/WD/Heatmap/HZE.Mesoderm.heatmap.Rdata")
HZE.Mesoderm <- qtl.smaller
rm(qtl.smaller)

load(file="/Users/elijahedmondson/Desktop/R/QTL/WD/Heatmap/Gamma.Ectoderm.heatmap.Rdata")
Gamma.Ectoderm <- qtl.smaller
rm(qtl.smaller)

load(file="/Users/elijahedmondson/Desktop/R/QTL/WD/Heatmap/Gamma.Endoderm.heatmap.Rdata")
Gamma.Endoderm <- qtl.smaller
rm(qtl.smaller)

load(file="/Users/elijah/Desktop/R/QTL/WD/Heatmap/Gamma.Mesoderm.heatmap.Rdata")
Gamma.Mesoderm <- qtl.smaller
rm(qtl.smaller)

combine <- cbind(seqnames=as.character(result[[1:2]]$ID), -log10(result[[1]]$pv) -log10(result[[2]]$pv))

#Combining the columns
combined <- cbind(seqnames=as.character(seqnames))

combined <- cbind(seqnames=as.character(seqnames(HZE.Mesoderm)),
                  -log10(Background.Ectoderm$p.value),
                  -log10(Background.Endoderm$p.value),
                  -log10(Background.Mesoderm$p.value),
                  -log10(HZE.Ectoderm$p.value),
                  -log10(HZE.Endoderm$p.value),
                  -log10(HZE.Mesoderm$p.value),
                  -log10(Gamma.Ectoderm$p.value),
                  -log10(Gamma.Endoderm$p.value),
                  -log10(Gamma.Mesoderm$p.value))
head(combined)
heatmap(combined, Rowv = NA)

combined <- cbind(-log10(Background.Ectoderm$p.value),
                  -log10(Background.Endoderm$p.value),
                  -log10(Background.Mesoderm$p.value),
                  -log10(HZE.Ectoderm$p.value),
                  -log10(HZE.Endoderm$p.value),
                  -log10(HZE.Mesoderm$p.value),
                  -log10(Gamma.Ectoderm$p.value),
                  -log10(Gamma.Endoderm$p.value),
                  -log10(Gamma.Mesoderm$p.value))

##with heatmap.2()
mypalette <- colorRampPalette(c("green", "yellow", "red"))(n = 299)

heatmap.2(t(combined), Colv=NA, col=mypalette,
          labCol=NA, sepwidth = 5, trace = "row", tracecol = "black",
          RowSideColors = c(
                  rep("gray", 3),
                  rep("blue", 3),
                  rep("black", 3)))


par(lend = 1)
legend(.75, 1.04, legend = c("Unirradiated", "HZE", "Gamma"),
       col = c("gray", "blue", "black"), lty= 1, lwd = 10)


##Plotting 3 QTL maps for comparison##
layout(matrix(3:1, 3, 1))
DOQTL:::plot.scanone.assoc(HZE, chr = 14, bin.size = 100, main = "HZE Ion")
DOQTL:::plot.scanone.assoc(Gamma, chr = 14, bin.size = 100, main = "Gamma ray")
DOQTL:::plot.scanone.assoc(Background, chr = 14, bin.size = 100, main = "Unirradiated")
