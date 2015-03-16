#!/usr/bin/env Rscript
#
#  phyloFlash_heatmap.R
#
#  Copyright (C) 2014-2015 Elmar Pruesse <elmar@pruesse.net>
#
#  This script generates heatmaps from phyloFlash output files.
#  Run as ./phyloFlash_heatmap.R to see options. 
#
#  LICENCE
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.
#  If not, see <http://www.gnu.org/licenses/>.

library(optparse);

load_libraries <- function() {
    library(methods);
    library(grid);
    library(ggplot2);
    library(reshape2);
    library(ggdendro);
    library(gtable);
}



defaultColorScheme <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9",
                        "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61",
                        "#F46D43", "#D73027", "#A50026");

####  helper functions

### logging

pf_logLevel <- 2;
pf_setLogLevel <- function(x) { assign("pf_logLevel", x, .GlobalEnv); }
msg <- function(...,lvl=2) {
    # levels:
    # 0 error   (always)
    # 1 warning 
    # 2 notice  
    # 3 info    
    # 4 debug

    if (lvl=="error") {
        stop(...);
    } else if (pf_logLevel >= lvl) {
        cat(...,'\n');
    }
}
err   <- function(...) msg(lvl=0,...);
warn  <- function(...) msg(lvl=1,...);
info  <- function(...) msg(lvl=3,...);
debug <- function(...) msg(lvl=4,...);

### ggplot helpers

# extract just the legend from a ggplot object
g_get_legend <- function(a.gplot) {
    # from http://stackoverflow.com/questions/11883844/
    # "inserting-a-table-under-the-legend-in-a-ggplot2-histogram"
    tmp <- ggplot_gtable(ggplot_build(a.gplot));
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box");
    legend <- tmp$grobs[[leg]];
    return(legend);
}

g_get_map <- function(a.gplot) {
    gb <- ggplotGrob(a.gplot);
    return(gtable_filter(gb,pattern="axis-l|panel"));
}

# hide the legend in a ggplot object
g_hide_legend <- function(a.gplot) {
    return (a.gplot + theme(legend.position = "none"));
}

g_hide_x_axis <- function(a.gplot) {
    return (a.gplot + theme(axis.text.x=element_blank())
#            + theme(plot.margin =unit(c(0,0,0,0),"null"))
            );
    return 
}
            
# prepare a pure dendrogram plot from a dendro_data object
# @param bool vertical   If true, plot has leaves as rows.
# @param bool labels     If true, plot includes labels
g_make_dendro_plot <- function(dendro, vertical=TRUE, labels=TRUE) {
    ddata <- dendro_data(dendro, type="rectangle");
    
    # Unexpanded, the dendrogram will span maximum width. That is,
    # the outer leaf nodes will end at the upper and lower corner
    # of the heapmap. To get them to end up at the center of the
    # heatmap rows/columns, we expand the axis. Expanding by .5
    # yields extra space to the left and right of the dendrogram
    # equivalent to the width of the dendrogram. Given n dendrogram
    # leaves, we have (n-1) spaces between the leaves. We need one
    # extra space, so we need an expansion factor of  .5/(n-1):
    expandFactor <- 0.5/(length(ddata$labels$label)-1);

    # plot the dendrogram without any labels or ticks or spaces
    # for ticks.
    p <- ggplot() +
         geom_segment(data = segment(ddata),
                      aes(x=x, y=y, xend=xend, yend=yend)) +
         theme_dendro() +
         labs(x = NULL, y = NULL) +
         scale_y_continuous(expand=c(0,0)) +
         theme(axis.ticks.length = unit(0,"null"),
               axis.ticks.margin = unit(0,"null")
               );

    # flip if vertical and add 1 mm space on the outer
    # edge to make sure the outmost connection is visible
    if (vertical) {
        p <- p + coord_flip() +
             theme(plot.margin   = unit(c(0,1,0,0),"mm"))
    } else {
        p <- p +
             theme(plot.margin   = unit(c(1,0,0,0),"mm"))
    }

    # add the appropriate scale configuration
    # if labels are to be shown, pull them from the
    # dendro_data object and decide whether the text
    # should be vertical or horizontal
    if (labels) {
        p <- p + scale_x_continuous(
            expand = c(expandFactor,0),
            breaks = 1:length(ddata$labels$label),
            labels = ddata$labels$label);
        if (vertical) {
            p <- p + theme(axis.text.y = element_text(angle = 0, hjust = 1));
        } else {
            p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1));
        }
    } else {
        p <- p + scale_x_continuous(expand=c(expandFactor,0));
    }
    
    return (p);             
}


# load phyloFlash data
read.phyloFlash <- function(files) {
    ### Extract library names from command line
    # remove .phyloFlash...csv 
    libs <- lapply(files[grepl("phyloFlash.*csv$",files)],
                   function(x) sub("\\.phyloFlash.*\\.csv","", x));
    libs <- unique(libs);
    
    if (length(libs) == 0)  {
        err("No phyloFlash output CSVs found on command line.")
    }

    msg("Selected phyloFlash libraries:");
    msg(unlist(libs), fill=77);
    
    msg("Loading CSV files...");
    for (lib in libs) {
        fileName <- paste(lib, ".phyloFlash.NTUabundance.csv", sep="");
        
        info("Reading: ",fileName);
        fileData <- read.csv(fileName);
       
        # rename "reads" column to lib name
        colnames(fileData)[colnames(fileData)=="reads"] = lib;
       
        # merge into single dataframe
        if (!exists("NTUcounts")) {
            NTUcounts <- fileData;
        } else {
            NTUcounts <- merge(x=NTUcounts, y=fileData, by="NTU", all=TRUE);
        }
    }
    
    ntu_names      <- NTUcounts$NTU;
    sample_names   <- colnames(NTUcounts[,-1]);
    
    NTUcounts  <- as.matrix(NTUcounts[,-1]);
    rownames(NTUcounts) <- ntu_names;
    
    # turn NA into 0
    NTUcounts[is.na(NTUcounts)] <- 0;

    return (NTUcounts);
}

shorten_taxnames <- function(data) {
    # convert NTU name column into rownames
    shortnames <- lapply(rownames(data),
                         function(x) gsub(".*;(.*;.*)","\\1", x));
    rownames(data) <- shortnames;

    return(data);
}

### remove taxa observed rarely
merge_low_counts <- function(data, thres=50, other="Other") {
    msg("A total of",nrow(data),"taxa were observed.",
        "Merging taxa with <", thres, "observations into \"",
        other, "\".");

    ndata <- data[rowSums(data) >= thres,];
    
    odata <- colSums(data[rowSums(data) < thres,]);
    odata <- matrix(odata, ncol=length(odata),
                    dimnames=list(c(other), names(odata)));

    data <- rbind(ndata, odata);

    msg("Removed taxa with <",thres,"observations.",
        nrow(data),"taxa left.");

    return(data);
}

scale_to_percent <- function(mat) {
    return(scale(mat, center=FALSE, scale=colSums(mat)) * 100);
}

cluster <- function(mat,method="ward") {
    return(as.dendrogram(hclust(dist(mat), method)));
}

g_make_heatmap <- function(mat, colorScheme=defaultColorScheme) {
    ## factorize dims
    matNames <- attr(mat, "dimnames");
    df <- as.data.frame(mat);
    colnames(df) <- matNames[[2]];
    df$y.variable <- matNames[[1]];
    df$y.variable <- with(df, factor(y.variable,levels=y.variable,ordered=TRUE));

    x=mean(mat);
    
    mat <- melt(df, id.vars="y.variable");
    
    heatMapPlot <- ggplot(mat, aes(variable,y.variable)) +
       geom_tile(aes(fill=value)) +
       #scale_fill_gradientn(colours = colorScheme) +
       #scale_fill_gradientn(colours = rainbow(3)) +
       #scale_fill_gradient2(low="blue", mid="green", high="red", midpoint = 30)+
       scale_fill_gradient(low="white", high="steelblue")+
       labs(x = NULL, y = NULL) +
       scale_x_discrete(expand=c(0,0)) +
       scale_y_discrete(expand=c(0,0)) +
       theme(axis.text.x = element_text(angle=90));
#                  theme(plot.margin = unit(c(0,0,1,1), "cm"));
    
    
    return(heatMapPlot);
}


plot.phyloFlash <- function(pdata) {
    row_dendro <- cluster(pdata);
    col_dendro <- cluster(t(pdata));

    ## reorder to match clustering
    pdata <- pdata [order.dendrogram(row_dendro),
                    order.dendrogram(col_dendro)];

    heatmap    <- g_make_heatmap(pdata);
    ntuTree    <- g_make_dendro_plot(row_dendro, TRUE, FALSE);
    sampleTree <- g_make_dendro_plot(col_dendro, FALSE, FALSE);

    g <- g_make_grid(heatmap, ntuTree, sampleTree);
    g <- g_add_to_grid(g, heatmap, ntuTree);
    g <- g_add_to_grid(g, heatmap, ntuTree);
    return(g);
}

g_make_grid <- function(heatMap, ntuTree, sampleTree) {
    ## start with the raw heatmap (no legend)
    g <- ggplotGrob(g_hide_legend(heatMap));

    ## add a column to the right
    g <- gtable_add_cols(g, unit(5, "cm"));
    ## put in the dendrogram for the NTUs
    g <- gtable_add_grob(g, ggplotGrob(ntuTree),
                         t=3, b=3, l=6, r=6);
    
    ## add a row above 
    g <- gtable_add_rows(g, unit(5, "cm"),0);
    ## put the dendrogram with the sample clustering in
    g <- gtable_add_grob(g, ggplotGrob(sampleTree),
                         t=1, b=1, l=4, r=4);
    ## put the legend in the upper right corner
    g <- gtable_add_grob(g, g_get_legend(heatMap),
                         t=1, b=1, l=6, r=6);

    return (g);
}

g_add_to_grid <- function(g, heatMap, ntuTree) {
    ## add some space
    g <- gtable_add_rows(g, unit(.2,"lines"),nrow(g)-3);

    ## add another heatmap
    g <- gtable_add_rows(g, unit(1, "null"),nrow(g)-3);
    ro <- nrow(g)-3;
    g <- gtable_add_grob(g, g_get_map(heatMap),
                         t=ro, b=ro, l=3, r=4);
    g <- gtable_add_grob(g, ggplotGrob(ntuTree),
                         t=ro, b=ro, l=6, r=6);
  
    return (g);
}    

pF_main <- function() {
    options <- list(
        make_option(
            c("-v", "--verbose"),
            action="store_true",
            default=FALSE,
            help="Be more talkative"
            ),
        make_option(
            c("-q", "--quiet"),
            action="store_true",
            default=FALSE,
            help="Be less talkative"
            ),
        make_option(
            "--min-ntu-count",
            default=50,
            type="integer",
            help="Sum NTUs with less counts in pseudo NTU \"Other\". Default %default."
            ),
        make_option(
            "--no-split-euks",
            action="store_false",
            help="Do not show Eukaryotes in separate plot section"
            ),
        make_option(
            "--out",
            default="out.png",
            help="Name of output file. Must end in .png or .pdf. Default %default."
            ),
        make_option(
            "--out-size",
            default="1024x768",
            help="Size of output graphic. Default %default"
            )
        );
    
    parser <- OptionParser(
        option_list=options,
        usage="usage: %prog [options] [files]",
        description="
Generates a heatmap plot from multiple phyloFlash result sets.

Files:
        A list of files and/or directories that will be searched
        for phyloFlash results."
        );

    conf <- parse_args(parser, positional_arguments = TRUE);

    # fixme: unparsed "-*" type options

    if (length(conf$args)==0) {
        print_help(parser);
        quit(status=2);
    }

    #str(conf$options);
    
    # set loglevel
    if (conf$options$quiet) {
        pf_setLogLevel(1);
    } else if (conf$options$verbose) {
        pf_setLogLevel(3);
    }

    load_libraries();

    pdata <- read.phyloFlash(conf$args);
    pdata <- shorten_taxnames(pdata);
    pdata <- merge_low_counts(pdata, thres=conf$options$"min-ntu-count");
    pdata <- scale_to_percent(pdata);

    g <- plot.phyloFlash(pdata);
    
    png(file="out.png", width=1280,height=1024);
    #svg(file="out.svg");
    #png(file="out.png",width=1280,height=1280);
    #pdf(file="out.pdf");
    grid.newpage();
    grid.draw(g);
    dev.off();

    png(file="out2.png", width=1280,height=1024);
    grid.newpage();
    gtable_show_layout(g);
    dev.off();

    


    
    msg("Brief summary of counts:");
    if (pf_logLevel >= 2) { ### "cat(summary(...))" does
        print(summary(pdata));
    }

    print(colSums(pdata));
    
}
    
# if we are run as a script from the cmdline
if (!interactive()) {
    pF_main();
} else {
    load_libraries();
}
