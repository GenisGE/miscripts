# funciton go merge rohs within a certain distance of each other (default 100 kb)
# d is a dataframe from plink .hom format, should be given data for a single chromosome and sample to make sure it only
# merges rohs that should be merged, funciton does not check
combineROH <- function(d, maxdist = 1e5){


    newd <- data.frame(
        FID=character(0),
        IID=character(0),
        CHR=character(0),
        POS1=integer(0),
        POS2=integer(0),
        KB=numeric(0)
    )

    if(nrow(d) == 0) return(newd)
    if(nrow(d) == 1) {
 newd <- rbind(newd, data.frame(FID=d$FID[1],
                                          IID=d$IID[1],
                                          CHR=d$CHR[1],
                                          POS1=d$POS1[1],
                                          POS2=d$POS2[1],
                                          KB=(d$POS2[1] - d$POS1[1]) / 1e3
                                          )
                          )
return(newd)
}

    dists <- d[2:nrow(d),"POS1"] - d[1:(nrow(d)-1), "POS2"]

    i1 <- 1

    while(i1 < nrow(d)){

        i2 <- i1

        while(dists[i2] < maxdist){

            i2 <- i2 + 1

            if(i2 > length(dists)){

                newd <- rbind(newd, data.frame(FID=d$FID[i1],
                                              IID=d$IID[i1],
                                              CHR=d$CHR[i1],
                                              POS1=d$POS1[i1],
                                              POS2=d$POS2[i2],
                                              KB=(d$POS2[i2] - d$POS1[i1]) / 1e3
                                              )
                              )
                break
                }
        }

        if(i2 > length(dists)) break

        newd <- rbind(newd, data.frame(FID=d$FID[i1],
                                IID=d$IID[i1],
                                CHR=d$CHR[i1],
                                POS1=d$POS1[i1],
                                POS2=d$POS2[i2],
                                KB=(d$POS2[i2] - d$POS1[i1]) / 1e3
                                )
                      )

        i1 <- i1 + 1

        if(i1 == nrow(d)){

            newd <- rbind(newd, data.frame(FID=d$FID[i1],
                                          IID=d$IID[i1],
                                          CHR=d$CHR[i1],
                                          POS1=d$POS1[i1],
                                          POS2=d$POS2[i1],
                                          KB=(d$POS2[i1] - d$POS1[i1]) / 1e3
                                          )
                          )

            }

    }

    return(newd)
}





## function to get proprotion of genome (genome where roh could have been called is given as autosome_len)
## for a certain sample, binned by roh lenght. hardcoded roh lenght bins, might want to change
getROHproportion <- function(roh, s, autosome_len){

    x <- roh[roh$IID==s,]
    
    d <- c(
        sum(x$KB[x$KB >= 0.5e3 & x$KB < 1e3] * 1e3),
        sum(x$KB[x$KB >= 1e3 &x$KB < 2e3] * 1e3),
        sum(x$KB[x$KB >= 2e3 & x$KB < 5e3] * 1e3),
        # sum(x$KB[x$KB >= 5e3] * 1e3)
        sum(x$KB[x$KB >= 5e3 & x$KB < 10e3] * 1e3),
        sum(x$KB[x$KB >= 10e3] * 1e3)
    )

    return(d/autosome_len)
}
