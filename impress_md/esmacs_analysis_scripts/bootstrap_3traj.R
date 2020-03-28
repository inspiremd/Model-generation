library(boot)
args <- commandArgs(trailingOnly = TRUE)
file_com <- paste('g-com-tot-',args[1],'-rep-all.dat',sep="")
file_rec <- paste('g-rec-tot-',args[1],'-rep-all.dat',sep="")
file_lig <- paste('g-lig-tot-',args[1],'-rep-all.dat',sep="")
gcom <-read.table(file_com)
grec <-read.table(file_rec)
glig <-read.table(file_lig)
num_reps <- as.integer(args[2])
sim.comb <- NULL
for (n in 1:5000) {
    sim.comb <- rbind(sim.comb,mean(sample(gcom[,1],num_reps,replace=TRUE,prob=NULL),na.rm=TRUE) - mean(sample(grec[,1],num_reps,replace=TRUE,prob=NULL),na.rm=TRUE) - mean(sample(glig[,1],num_reps,replace=TRUE,prob=NULL),na.rm=TRUE))
}
file_out <- paste('resample-avg-',args[1],'.dat',sep="")
write.table(sim.comb, file=file_out, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
