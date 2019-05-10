#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
infile	<- args[1]

if (length(args)<1) {stop("[USAGE] Rscript plot_TPM.r infile")}

a<-read.table(infile, header=T)
pdfile<-paste(infile, ".pdf", sep = "")

pdf(pdfile)
plot(a, type="l")
dev.off()
