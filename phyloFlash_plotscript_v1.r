#!/usr/bin/Rscript
## plotscript.r by Brandon

## Prerequisites: APE package in R
## Usage: Rscript plotscript.r --args <treefile> <histofile>
## Output: <treefile>.pdf, <treefile>.png, <histofile>.pdf, <histofile>.png

## Get arguments from command line
args <- commandArgs(trailingOnly=TRUE)
treefile <- args[2]
histofile <- args[3]

## Plot tree from Newick .tree file (guide tree produced by MAFFT) and output PDF plot
library(ape)							# ape package for phylogenetics
thetree <- read.tree(file=treefile)

pdf(file=paste(treefile,"pdf",sep="."),width=11,height=8)
plot(thetree,type="phylogram",no.margin=TRUE,font=1,cex=0.5)
dev.off()

png(file=paste(treefile,"png",sep="."),width=1100,height=800)
plot(thetree,type="phylogram",no.margin=TRUE,font=1,cex=0.5)
dev.off()


## Plot insert size histogram
histo <- read.table(file=histofile,header=F,sep="\t",skip=1)
names(histo) <- c("InsertSize", "Count")

xx <- c(histo$InsertSize, rev(histo$InsertSize))		# Define polygon for solid plot, using trick from 
yy <- c(rep(0,nrow(histo)),rev(histo$Count))			# http://earlh.com/blog/2009/07/28/filled-line-plots-graphs-in-r-part-10-in-a-series/

pdf(file=paste(histofile,"pdf",sep="."),width=11,height=8)
plot(histo, xlab="Insert size", ylab="Counts", col="grey", type="l")
polygon(xx, yy, border=NA, col="grey")
dev.off()

png(file=paste(histofile,"png",sep="."),width=1100,height=800)
plot(histo, xlab="Insert size", ylab="Counts", col="grey", type="l")
polygon(xx, yy, border=NA, col="grey")
dev.off()
