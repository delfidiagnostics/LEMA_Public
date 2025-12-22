### Fragment weighted GC correction (As in Mathios 2021)
gcCorrectTarget <- function(fragments, ref, bychr=TRUE){
    fragments[, gc := round(gc, 2)]
    if(bychr) {
        DT.gc <- fragments[,.(n=.N), by=.(gc, chr)]
        DT.gc <- DT.gc[gc >= .20 & gc <= .80]
        DT.gc <- DT.gc[order(gc, chr)]
    } else {
        DT.gc <- fragments[,.(n=.N), by=gc]
        DT.gc <- DT.gc[gc >= .20 & gc <= .80]
        DT.gc <- DT.gc[order(gc)]
    }
# setkey(mediandt, gc, seqnames)

    if(bychr) {
        setkey(DT.gc, gc, chr)
        setkey(ref, gc, chr)
    } else {
        setkey(DT.gc, gc)
        setkey(ref, gc)
    }
#     DT.gc <- DT.gc[ref][order(chr, gc)]
    DT.gc <- DT.gc[ref]
    DT.gc[,w:=target/n]
    if(bychr) {
        fragments[DT.gc, on= .(chr, gc), weight := i.w]
    }
    else fragments[DT.gc, on= .(gc), weight := i.w]
    fragments <- fragments[!is.na(weight)]
    fragments[,weight := weight * .N/sum(weight)]
    fragments[]
}

### Aggregate fragments into bins
binFrags <- function(fragments, bins,
                     chromosomes=paste0("chr",c(1:22, "X", "Y"))) {
    bins <- bins[chr %in% paste0("chr", 1:22)]
    fragments <- fragments[chr %in% paste0("chr", 1:22)]
    setkey(bins, chr, start, end)
    setkey(fragments, chr, start, end)
    fragbins <- foverlaps(fragments[chr %in% chromosomes],
                          bins, type="within", nomatch=NULL)
    bins2 <- fragbins[,.(arm=unique(arm), gc=gc[1], map=map[1],
                         cov = sum(weight),
                         short = sum(weight[w >= 100 & w <= 150]),
                         long = sum(weight[w > 150 & w <= 250]),
                         amplitude = (sum(weight[w==134])+sum(weight[w==145]))/2,
                         mediansize = as.integer(median(w)),
                         frag.gc = mean(fraggc)),
            by=.(chr, start, end)]
    setkey(bins2, chr, start, end)
    bins2 <- bins2[bins]
    bins2 <- bins2[is.na(i.gc), which(grepl("cov", colnames(bins2))):=0]
    bins2[,`:=`(gc=i.gc, map=i.map, arm=i.arm)]
    bins2[,which(grepl("^i.", colnames(bins2))):=NULL]
    bins2[, bin:=1:.N]
    setcolorder(bins2, c("chr", "start", "end", "bin"))
    bins2[]
}
