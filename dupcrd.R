#!/usr/bin/env Rscript
# dupcrd.R duplicate CSV read ... uses the output of csvrdm.c
# args <- commandArgs(trailingOnly = TRUE)
# numargs <- length(args)
# enumargs <- 1 # expected numargs
# if(numargs != enumargs) {
#     print("This script processes the output of csvrdm ... needs one argument, a filename.")
#     stop("Stopping right here")
# }

# csv <- read.csv("sigu.csv")
# csv <- read.csv("sigcsvrd.csv")
# this next one is where I corrected for av3->vsz < av4->vsz and av3->vsz>1, which was a problem.
# csv <- read.csv("sig22.csv")
# unfortunately a blight with the above was csvrdm retain a typeo GpgGrp instead of CpgGrp so
# yes, it was this that caused fewer genes in dupcrd.R
csv <- read.csv("sig222.csv")

doon <- paste0(csv$Genename, "__", csv$CpgGrp)
doon2 <- paste0(csv$Genename, "__", csv$CpgGrp, "__", csv$Relation_to_Island)
d3 <- paste0(csv$Genegrp, "__", csv$Relation_to_Island)
# the duplicated annotation provide no extra information.
w <- which(!duplicated(doon))
w2 <- which(!duplicated(doon2))
w2m <- setdiff(w2,w)
w2mi <- which(!w2 %in% w)
csvw <- csv[w,]
csvw <- csvw[order(csvw$Genename),]
csvw2 <- csv[w2,]
csvw2 <- csvw2[order(csvw2$Genename),]
# first one ignores Relation_to_Island
write.csv(csvw, "sigcdd_.csv", quote=F, row.names=F)
# second csv does not.
write.csv(csvw2, "sigcdd3_.csv", quote=F, row.names=F)
