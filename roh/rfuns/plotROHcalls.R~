# Rscript plotROHcalls.R inrohsfile outpng autosomebed famfile
source("/home/genis/impala/analyses/goatMapV3/roh/roh_call/rohrfuns/rohfuns.R")


args <- commandArgs(trailingOnly=T)

f <- args[1] # plink .hom file with rohs
outpng<- args[2] # name of output png file
g <- args[3] # bed file with size of autosomes
fam <- args[4] # paht to fam file, assumed that individual ids are first column
mergedist <- 5e5

autosome <- read.table(g)
autosome <- autosome[autosome$V3 > 5e6,]

autosome_len <- sum(as.numeric(autosome$V3))

roh <- read.table(f, h=T)
roh <- roh[roh$CHR %in% autosome$V1,]

samples <- as.character(unique(roh$FID))
chrs <- sort(unique(roh$CHR))


getROHproportion <- function(roh, s, autosome_len){

    x <- roh[roh$IID==s,]
    d <- c(
        sum(x$KB[x$KB >= 0.5e3 & x$KB < 1e3] * 1e3),
        sum(x$KB[x$KB >= 1e3 &x$KB < 2e3] * 1e3),
        sum(x$KB[x$KB >= 2e3 & x$KB < 5e3] * 1e3),
        sum(x$KB[x$KB >= 5e3 & x$KB < 10e3] * 1e3),
        sum(x$KB[x$KB >= 10e3] * 1e3)
    )

    return(d/autosome_len)
}


l <- apply(as.matrix(expand.grid(chrs, samples)),1, function(x) combineROH(roh[roh$CHR==x[1] & roh$FID == x[2],], maxdist=mergedist))
roh <- do.call("rbind", l)


ids <- read.table(fam, h=F, stringsAsFactors=F)$V1

m <- sapply(ids, getROHproportion, roh=roh, autosome_len=autosome_len)


#outpng <- "/newHome/genis/mydata/impala/analyses/goatMapV3/roh/roh_call/plots/impalarohcallsgoatmap.png"

bitmap(outpng, h=5,w=6,res=300)
par(oma=c(2,0,0,6.5))
barplot(m, names.arg=ids, las=2, col=RColorBrewer::brewer.pal("Set2", n=5), ylab="Genome fraction", main="ROH distribution", cex.lab=1.2, cex.main=1.4)

legend(x=ncol(m) * 1.2, y=max(colSums(m)) * 0.75, legend=c("0.5-1", "1-2","2-5","5-10",">10"), fill = RColorBrewer::brewer.pal("Set2", n=5), title="ROH length (Mbp)", bty="n", xpd=NA, cex=1.2)

#text(y=-0.08, x=c(1.5,6), labels=c("Central", "North"), xpd=NA, cex=1.5)
dev.off()
