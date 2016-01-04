GRconvert = function(result) {
        library(GenomicRanges)
        library(BSgenome.Mmusculus.UCSC.mm10)
        library(DOQTL)

        chrs = c(1:19, "X")
        qtl = GRangesList(GRanges("list", length(result)))

        for(i in 1:length(chrs)) {
                print(i)
                qtl[[i]] <- GRanges(seqnames = Rle(result[[i]]$CHR),
                                    ranges = IRanges(start = result[[i]]$POS, width = 1),
                                    p.value = result[[i]]$pv)
        } # for(i)

}




load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/Background/AMQTL.Mesoderm.Rdata")
qtl <- AM.qtl
rm(AM.qtl)
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Background.Mesoderm.heatmap.Rdata")
rm(qtl, qtl.smaller)

load(file="/Users/elijah/Desktop/R/QTL/WD/Heatmap/Gamma.Mesoderm.heatmap.Rdata")
Gamma.Mesoderm <- qtl.smaller
rm(qtl.smaller)

#Combining the columns

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
