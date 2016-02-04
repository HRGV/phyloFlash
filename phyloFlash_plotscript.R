#!/usr/bin/env Rscript
## plotscript.r by Brandon Seah (kbseah@mpi-bremen.de)

## Prerequisites: APE package in R
## Usage: Rscript plotscript.r --args <treefile> <histofile>
## Output: <treefile>.pdf, <treefile>.png, <histofile>.pdf, <histofile>.png

## Get arguments from command line
args <- commandArgs(trailingOnly=TRUE)
treefile <- args[2]
histofile <- args[3]

if (treefile != "NULL") {   # If no Newick tree file was generated (when -skip_emirge and -skip_spades options activated), skip this step
    ## Plot tree from Newick .tree file (guide tree produced by MAFFT) and output PDF plot
    library(ape)							# ape package for phylogenetics
    thetree <- read.tree(file=treefile)
    
    pdf(file=paste(treefile,"pdf",sep="."),width=11,height=8)
    plot(thetree,type="phylogram",no.margin=TRUE,font=1,cex=0.5)
    dev.off()
    
    png(file=paste(treefile,"png",sep="."),width=1100,height=800)
    plot(thetree,type="phylogram",no.margin=TRUE,font=1,cex=0.5)
    dev.off()
}

## Plot insert size histogram
histo <- read.table(file=histofile,header=F,sep="\t",skip=1)
names(histo) <- c("InsertSize", "Count")
# "Untabulate" the tabulated insert size counts
histvals <- as.vector( # Convert to vector
                      unlist( # "Flatten" list to atomic elements
                             apply(histo, # Table of counts of insert sizes
                                   1, # Margin - by row
                                   function(x) rep(x[1],x[2]) # Repeat each value by number of counts
                                   )
                             )
                      ) 
# Export PDF version
pdf(file=paste(histofile,"pdf",sep="."),width=11,height=8)
hist(histvals,col="grey",border="grey",xlab="Insert size (bp)",main="Insert size histogram")
dev.off()
# Export PNG version
png(file=paste(histofile,"png",sep="."),width=1100,height=800)
hist(histvals,col="grey",border="grey",xlab="Insert size (bp)",main="Insert size histogram")
dev.off()
