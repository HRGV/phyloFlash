#!/usr/bin/env Rscript
## plotscript.r by Brandon Seah (kbseah@mpi-bremen.de)

## Prerequisites: APE package in R
## Usage: Rscript plotscript.r --args <treefile> <histofile> <idhistofile> <decimalcomma>
## Output: <treefile>.pdf, <treefile>.png, <histofile>.pdf, <histofile>.png

## Get arguments from command line
args         <- commandArgs(trailingOnly=TRUE)
treefile     <- args[2] # Tree from MAFFT
histofile    <- args[3] # Insert size histogram
idhistofile  <- args[4] # Mapping %id histogram
decimalcomma <- args[5] # Use decimal comma?

# Skip drawing tree if no Newick tree file was generated
# i.e. -skip_emirge and -skip_spades options in phyloFlash.pl
if (treefile != "NULL") {
    ## Plot tree from Newick .tree file (guide tree produced by MAFFT) and output PDF plot
    library(ape) # ape package for phylogenetics
    thetree <- read.tree(file=treefile)
    
    pdf(file=paste(treefile,"pdf",sep="."),width=11,height=8)
    plot(thetree,type="phylogram",no.margin=TRUE,font=1,cex=0.5)
    dev.off()
    
    png(file=paste(treefile,"png",sep="."),width=1100,height=800)
    plot(thetree,type="phylogram",no.margin=TRUE,font=1,cex=0.5)
    dev.off()
}

# If -decimalcomma option used in phyloFlash.pl
if (decimalcomma == 1) {
    dec <- ","
} else {
    dec <- "."
}

## Plot insert size histogram if not running in SE mode
if (histofile != "SEmode") {
    histo <- try(
        read.table(file=histofile,
                        header=F,
                        sep="\t", # Tab-separated table
                        dec=dec, 
                        comment.char="#")
    )
    if (class(histo) == "try-error") {
        quit()
    }
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
    pdf(file=paste(histofile,"pdf",sep="."),
        width=11,
        height=8)
    hist(histvals,
         col="grey",
         border="grey",
         xlab="Insert size (bp)",
         main="Insert size histogram"
         )
    dev.off()
    # Export PNG version
    png(file=paste(histofile,"png",sep="."),
        width=660,
        height=480
        )
    hist(histvals,
         col="grey",
         border="grey",
         xlab="Insert size (bp)",
         main="Insert size histogram"
         )
    dev.off()
}

## Plot percent-identity histogram
idhisto <- try(
    read.table (file=idhistofile,
                       header=F,
                       sep="\t",
                       dec=dec, # Detect decimal separator for this locale
                       comment.char="#")
)
if (class(idhisto) != "try-error") {
idhistvals <- as.vector(
                      unlist(
                             apply(idhisto,
                                   1,
                                   function(x) rep(x[1],x[2])
                                   )
                             )
                     )

# Export PDF
pdf(file=paste(idhistofile,"pdf",sep="."),
    width=11,
    height=8)
hist(idhistvals,
     breaks=100,
     col="grey",
     border="grey",
     xlab="Percentage identity",
     main="Histogram of read identity to reference"
     )
dev.off()
# Export PNG
png(file=paste(idhistofile,"png",sep="."),
    width=660,
    height=480)
hist(idhistvals,
     breaks=100,
     col="grey",
     border="grey",
     xlab="Percentage identity",
     main="Histogram of read identity to reference"
     )
dev.off()
}
