# Rscript plotROHgenomeROHcalls.R inplink inrohs outpng sampleid faifile

source("rohfuns.R")
source("http://popgen.dk/albrecht/open/online.R")
library(scales)

args <- commandArgs(trailingOnly=T)

inplink <- args[1]
inrohs <- args[2]
outpng <- args[3]
s <- args[4]
faifile <- args[5]
mergedist <- 5e5


cat("running with following parameters:\n")
cat("inplink:", inplink, "\n")
cat("inrohs:", inrohs, "\n")
cat("outpng:", outpng, "\n")
cat("s:", s, "\n")


# load data
rohs <- read.table(inrohs, stringsAsFactors=F, h=T)

cat("loading plink file...\n")
dat <- plink(plinkFile = inplink)

fai <- read.table(faifile, stringsAsFactors=F)

chrs <- fai$V1[fai$V1 %in% unique(dat$bim$V1) & fai$V2>5e6]
chrlen <- fai$V2[fai$V1 %in% unique(dat$bim$V1) & fai$V2>5e6]
#chrs <- fai$V1[fai$V1 %in% unique(dat$bim$V1)]
#chrlen <- fai$V2[fai$V1 %in% unique(dat$bim$V1)]


rohs <- do.call("rbind", lapply(chrs, function(x) combineROH(rohs[rohs$CHR == x,], maxdist=mergedist)))


xlim <- c(0, max(chrlen))
ylim <- c(0, length(chrs) * 2)


nsnps <- 100

cols <- c( "#0095EF", "#3C50B1", "#6A38B3", "#A224AD", "#F31D64", "#FE433C")

cat("finished reading plink, will make plot!\n")
bitmap(outpng, w=12, h=6, res=300)
# make empty plot to draw chromosomes
plot(1, type="n",
     yaxt="n",
     ylab="",
     xlim=xlim, ylim=ylim,
     xlab="Position (Mbps)",
     main=paste("Runs of Homozygosity sample", s))

axis(2, at = length(chrs)*2- 0:(length(chrs)-1) *2 - 1.5, labels=chrs, las=2)
#for(i in 1:length(chrs)){
for(i in 1:length(chrs)){
    
    chrom <- chrs[i]

    k <- dat$bim$V1 == chrom
    g <- dat$geno[,k]
    pos <- dat$bim[dat$bim$V1==chrom, 4]
    
    r <- rohs[rohs$CHR == chrom,]
    
    wins <-  ceiling(seq_along(1:length(g))/nsnps)

    nhets <- tapply(g, wins, function(x) sum(x==1, na.rm=T))
    plotcols <- cols[ifelse(nhets > 5, 5,nhets)+1]

    start <-  tapply(pos, wins, function(x) x[1])
    end <-  tapply(pos, wins, function(x) x[length(x)])

    i <- i - 1

    rect(xleft=start, xright=end,
         ytop=length(chrs)*2-i*2 - 1,
         ybottom=length(chrs)*2-(i + 1)*2,
         col=plotcols, border=plotcols)


    if(nrow(r) > 0){
        segments(x0 = r$POS1, x1 = r$POS2,
                 y0 = length(chrs)*2-i*2 - 0.5,
                 y1 = length(chrs)*2-i*2 - 0.5,
                 lwd=1)

        segments(x0=c(r$POS1, r$POS2), x1=c(r$POS1, r$POS2),
                 y0 = length(chrs)*2-i*2 - 0.75,
                 y1 = length(chrs)*2-i*2 - 0.25,
                 lwd=1)
    }
}
dev.off()

cat("wrote plot to", outpng, "\n")
