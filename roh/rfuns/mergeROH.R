#Rscript mergeROH.R inrohsfile mergedist outrohsfile



whereDir <- function(){
    # function to get directory where scripts are, so apohFuns.R can be sourced when run from any folderfrom outside. Assumes apoh.R and apohFuns.R are in the same folder

    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file"
    match <- grep(needle, cmdArgs)
    tf <- unlist(strsplit(cmdArgs[match], "="))[2]
    d <- dirname(tf)
    return(d)
}

d <- whereDir()
source(paste(d, "rohfuns.R", sep="/"))

args <- commandArgs(trailingOnly=T)

f <- args[1] # plink .hom file with rohs
mergedist <- as.integer(args[2])
outf <- args[3]

roh <- read.table(f, h=T)
#roh <- roh[roh$CHR,]

samples <- as.character(unique(roh$FID))
chrs <- sort(unique(roh$CHR))

l <- do.call("rbind", apply(as.matrix(expand.grid(chrs, samples)),1, function(x) combineROH(roh[roh$CHR==x[1] & roh$FID == x[2],], maxdist=mergedist)))


write.table(l, outf, col.names=T, row.names=F, quote=F, sep="\t")
