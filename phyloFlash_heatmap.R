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


### prerequisite packages

required.packages = c("optparse", "methods", "grid", "ggplot2", "reshape2",
                       "ggdendro", "gtable");
check_libraries <- function() {
    missing.packages <- required.packages[!(required.packages %in% installed.packages())];
    if (length(missing.packages)) {
        msg("Additional packages required: ", missing.packages);
        if (options("repos")[[1]] == "@CRAN@") {
            options(repos = "https://cran.rstudio.com/")
        }
        install.packages(missing.packages);
    }
}

load_libraries <- function() {
    library(methods);
    library(grid);
    library(ggplot2);
    library(reshape2);
    library(ggdendro);
    library(gtable);
}

### logging

## more useful errors on cmdline (traceback)
pf_debug <- FALSE;
pf_setDebug <- function(x) {
    debugCode = quote({
        dump.frames();
        cat(paste("  ", 1L:length(last.dump), ": ",
                  names(last.dump), sep = ""),"",
            sep = "\n", file=stderr())
    });
    if (x) {
        options(warn=2, keep.source=TRUE, error = debugCode);
    } else {
        options(warn=0, keep.source=TRUE, error=NULL);
    }
    assign("pf_debug", x, .GlobalEnv);
}


pf_logLevel <- 2;
pf_setLogLevel <- function(x) { assign("pf_logLevel", x, .GlobalEnv); }
msg <- function(...,lvl=2) {
    # levels:
    # 0 error   (always)
    # 1 warning 
    # 2 notice  
    # 3 info    
    # 4 debug

    if (lvl=="0") {
        if (pf_debug) {
            stop(...);
        } else {
            cat("ERROR: ", ...,'\n');
            quit();
        }
    } else if (pf_logLevel >= lvl) {
        cat(...,'\n');
    }
}
err   <- function(...) msg(lvl=0,...);
warn  <- function(...) msg(lvl=1,...);
info  <- function(...) msg(lvl=3,...);
debug <- function(...) msg(lvl=4,...);

## workaround for not working rbind(gtable...)
## adapted from http://stackoverflow.com/questions/24234791
rbind_max <- function(...,size=grid::unit.pmax){
    bind2 <- function (x, y) {
        stopifnot(ncol(x) == ncol(y))
        if (nrow(x) == 0) return(y)
        if (nrow(y) == 0) return(x)
        y$layout$t <- y$layout$t + nrow(x)
        y$layout$b <- y$layout$b + nrow(x)
        x$layout <- rbind(x$layout, y$layout)
        x$heights <- gtable:::insert.unit(x$heights, y$heights)
        x$rownames <- c(x$rownames, y$rownames)
        if (is.function(size)) {
            x$widths <- do.call(size, list(x$widths, y$widths))
        }
        x$grobs <- append(x$grobs, y$grobs)
        x
    }
    Reduce(bind2, list(...))
}

cbind_max <- function(...,size=grid::unit.pmax){
    bind2 <- function (x, y) {
        stopifnot(nrow(x) == nrow(y))
        if (ncol(x) == 0) return(y)
        if (ncol(y) == 0) return(x)
        y$layout$l <- y$layout$l + ncol(x)
        y$layout$r <- y$layout$r + ncol(x)
        x$layout <- rbind(x$layout, y$layout)
        x$widths <- gtable:::insert.unit(x$widths, y$widths)
        x$colnames <- c(x$colnames, y$colnames)
        if (is.function(size)) {
            x$heights <- do.call(size, list(x$heights, y$heights))
        }
        x$grobs <- append(x$grobs, y$grobs)
        x
    }
    Reduce(bind2, list(...))
}

arg_select <- function(select, ...) {
    data <- list(...);
    if (is.numeric(select)) {
        data[select];
    } else {
        from <- d<-sapply(substitute(list(...)), deparse)[-1]
        data <- data[match(select,from)];
    }
}

cbind_select <- function(select, ..., size=grid::unit.pmax) {
    do.call(cbind_max, c(arg_select(select, ...), list(size=size)));
}

rbind_select <- function(select, ..., size=grid::unit.pmax) {
    do.call(rbind_max, c(arg_select(select, ...), list(size=size)));
}

# extract a grob from a ggplot/gtable
g_get <- function(obj, pat) {
    if (is.ggplot(obj)) obj <- ggplotGrob(obj);
    if (!is.grob(obj)) err("not a grob?!");
    return (gtable_filter(obj,pattern=pat));
}

# prepare a pure dendrogram plot from a dendro_data object
# @param bool axis       1:4=below,left,above,right
# @param bool labels     If true, plot includes labels
g_make_dendro_plot <- function(dendro, axis, labels=TRUE) {
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

    # The axis is placed as follows: 1=below, 2=left,
    # 3=above and 4=right.
    trans <- ifelse(axis < 3, "reverse", "identity");

    # plot the dendrogram without any labels or ticks or spaces
    # for ticks.

    p <- ggplot() +
         geom_segment(data = segment(ddata),
                      aes(x=x, y=y, xend=xend, yend=yend)) +
         theme_dendro() +
         labs(x = NULL, y = NULL) +
         scale_y_continuous(expand=c(0,0), trans=trans) +
         theme(axis.ticks.length = unit(0,"null"),
               axis.ticks.margin = unit(0,"null")
               );

    # flip if vertical and add 1 mm space on the outer
    # edge to make sure the outmost connection is visible
    
    if (axis %% 2 == 0) {
        p <- p + coord_flip()
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
        if (axis %% 2 == 0) {
            p <- p + theme(axis.text.y = element_text(angle = 0, hjust = 1));
        } else {
            p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1));
        }
    } else {
        p <- p + scale_x_continuous(expand=c(expandFactor,0));
    }
    
    return (invisible(p));
}

# makes a ggplot heatmap from a matrix
g_make_heatmap <- function(mat, highcol, angle=90, hjust=0,vjust=0.6) {
    ## factorize dims
    matNames <- attr(mat, "dimnames");
    df <- as.data.frame(mat);
    colnames(df) <- matNames[[2]];
    df$y.variable <- matNames[[1]];
    df$y.variable <- with(df, factor(y.variable,levels=y.variable,ordered=TRUE));

    mat <- melt(df, id.vars="y.variable");

    breaks <- function(x) axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE);
    
    heatMapPlot <- ggplot(mat, aes(variable,y.variable)) +
       geom_tile(aes(fill=value)) +
       scale_fill_gradient(low="white", high=highcol, trans="log",
                           breaks=breaks,
                           limits=c(1,max(mat$value)),
                           na.value="white") +
       labs(x = NULL, y = NULL) +
       scale_x_discrete(expand=c(0,0)) +
       scale_y_discrete(expand=c(0,0)) +
       theme(axis.text.x = element_text(angle=angle, hjust=hjust,vjust=vjust),
             axis.ticks.length = unit(0,"null"),
             legend.title=element_blank());

    return(invisible(heatMapPlot));
}

# makes a grob containing a row of strings from a string vector
gtable_text_row <- function(strvec) {
    grobs <- lapply(strvec, function(str) {
        textGrob(str,gp=gpar(fontsize=8))
    })
    gt <- gtable(heights=unit(1,"lines"), widths=unit(0,"null"));

    if (length(strvec) == 0) return(gt);

    g <- gtable_row("textrow", grobs);
    gt <- gtable_add_grob(gt, g, t=1, l=1);

    invisible(gt)
}

# loads phyloFlash output files into R
read.phyloFlash <- function(files=".",sampleNameFromMeta=TRUE) {
    if (length(files) == 1) {
        files <- list.files(pattern=files);
    }
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
    for (lib in libs) tryCatch({
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

        # read meta-data
        fileName <- paste(lib, ".phyloFlash.report.csv", sep="");
        info("Reading: ", fileName);
        fileData <- read.csv(fileName, col.names=c("key",lib));
        if (!exists("MetaData")) {
            MetaData <- fileData;
        } else {
            MetaData <- merge(x=MetaData, y=fileData, by="key");
        }
    }, warning = function(e) {
        err("Error in read.phyloFlash:\n",e$message);
    });

    pfData <- list(); # result is a list

    ## turn key column into row names, transpose
    pfData$meta <- MetaData[,-1];
    rownames(pfData$meta) <- MetaData[,1];
    pfData$meta <- data.frame(t(pfData$meta));

    ## turn data.frame into matrix, get dimnames right
    ntu_names      <- NTUcounts$NTU;
    sample_names   <- colnames(NTUcounts[,-1]);
    NTUcounts  <- as.matrix(NTUcounts[,-1]);
    rownames(NTUcounts) <- ntu_names;
    if (sampleNameFromMeta) {
        colnames(NTUcounts) <- pfData$meta$library.name
    } else {
        colnames(NTUcounts) <- sample_names;
    }
    NTUcounts[is.na(NTUcounts)] <- 0;     # turn NA into 0
    pfData$data <- list(NTUcounts);

    return (pfData);
}

# shortens taxnames to last to group names
shorten_taxnames <- function(data) {
    if (is.list(data)) return (lapply(data, shorten_taxnames));
    # convert NTU name column into rownames
    shortnames <- lapply(rownames(data),
                         function(x) gsub(".*;(.*;.*)","\\1", x));
    rownames(data) <- shortnames;

    return(data);
}

# splits matrix using regex (returns list)
split_by_name <- function(data, re_list) {
    names  <- rownames(data[[1]]);
    groups <- rep(0, length(names));
    i <- 1;
    for (pat in re_list) {
        groups[grepl(pat, names) & groups==0] = i;
        i <- i+1;
    }
    sd <- split(data.frame(data), groups);
    return (lapply(sd, as.matrix));
}

### remove taxa observed rarely
merge_low_counts <- function(data, thres=50, other="Other") {
    if (is.list(data)) {
        return(lapply(data, merge_low_counts, thres, other));
    }

    msg("A total of",nrow(data),"taxa were observed.",
        "Merging taxa with <", thres, "observations into \"",
        other, "\".");

    ndata <- data[rowSums(data) >= thres,];

    if (length(ndata) == length(data)) {
        msg("No taxa to merge");
        return(data);
    }
    
    odata <- colSums(data[rowSums(data) < thres,]);
    odata <- matrix(odata, ncol=length(odata),
                    dimnames=list(c(other), names(odata)));

    data <- rbind(ndata, odata);

    msg("Removed taxa with <",thres,"observations.",
        nrow(data),"taxa left.");

    return (data);
}

# scales matrix columns to percent
scale_to_percent <- function(mat) {
    if (is.list(mat))
        return (lapply(mat, scale_to_percent));
    return (scale(mat, center=FALSE, scale=colSums(mat)) * 100);
}

# cluster, create dendrograms and reorder data
cluster <- function(pf, method="ward") {
    mkdendro <- function(mat) {
        return(as.dendrogram(hclust(dist(mat), method)));
    }

    ## re-join data if list
    joined = do.call(rbind, pf$data);
    ## create horizontal clusters
    pf$col_dendro <- mkdendro(t(joined));
    ## re-order meta-data
    pf$meta <- pf$meta[order.dendrogram(pf$col_dendro),];

    ## create vertical clusters
    pf$row_dendro <- lapply(pf$data, mkdendro);
    ## re-order data matrices
    rorder <- function(mat, dendr) {
        return (mat[order.dendrogram(dendr),
                    order.dendrogram(pf$col_dendro)]);
    }

    pf$data <- mapply(rorder, pf$data, pf$row_dendro, SIMPLIFY=FALSE);

    return(pf);
}

# creates a plot from a phyloFlash "object"
plot.phyloFlash <- function(pf,
                            row.order=c("tree","map","chao","labels"),
                            col.order=c("labels","map","tree"),
                            map.colors=c("steelblue", "indianred")
) {
    ## turn arguments into vectors (workaround)
    row.order = strsplit(paste(collapse=",",row.order),",")[[1]];
    col.order = strsplit(paste(collapse=",",col.order),",")[[1]];
    map.colors = strsplit(paste(collapse=",",map.colors),",")[[1]];
    
    ## get number of maps and number of rows per map
    nmaps <- length(pf$data);
    nrows <- sapply(pf$data,nrow);

    ## some empty tables
    zero <- gtable(widths=unit(0,"null"),heights=unit(0,"null"));
    zero1 <- gtable(widths=unit(1,"null"),heights=unit(0,"null"));
    
    ## render the heapmaps
    gg_heatmaps <- mapply(g_make_heatmap, pf$data, map.colors, SIMPLIFY=FALSE);

    ## extract maps column
    gr_heatmaps <- lapply(gg_heatmaps, g_get, "panel");
    gr_heatmaps <- do.call(rbind_max, gr_heatmaps);

    ## extract labels column
    gr_labels   <- lapply(gg_heatmaps, g_get, "axis-l");
    gr_labels   <- do.call(rbind_max, gr_labels);

    ## extract legends
    gr_legends <- lapply(gg_heatmaps, function(x) {
        g_get(g_get(x, "guide-box")$grobs[[1]], "guides") })
    gr_legends <-  do.call(cbind_max, gr_legends)
    
    gr_legend <- gtable_add_grob(zero, gr_legends , t=1, l=1)
    gr_legend$heights = max(gr_legends$heights)
    gr_legend$widths = sum(gr_legends$widths)

    # extract sample labels
    gr_sample_labels  <- g_get(gg_heatmaps[[1]], "axis-b");

    ## render trees over taxa
    axis     <- ifelse(match("tree", col.order) < match("map", col.order), 2, 4);
    gr_trees <- lapply(pf$row_dendro, g_make_dendro_plot, axis=axis);
    gr_trees <- lapply(gr_trees, g_get, "panel");
    gr_trees <- do.call(rbind_max, gr_trees);

    ## scale heights by number of rows
    gr_heatmaps$heights <- gr_heatmaps$heights * (nrows/sum(nrows));
    gr_labels$heights   <- gr_labels$heights   * (nrows/sum(nrows));
    gr_trees$heights    <- gr_trees$heights    * (nrows/sum(nrows));

    # render tree over samples
    axis     <- ifelse(match("tree",row.order) < match("map", row.order), 3, 1);
    gr_sampleTree     <- g_get(g_make_dendro_plot(pf$col_dendro, axis=axis), "panel");
    gr_sampleTree$heights <- unit(0.1,"null")

    # make chao line
    chao <- pf$meta$NTU.Chao1.richness.estimate;
    chao <- round(as.numeric(as.character(chao)))
    gr_chao_grob <- textGrob("Chao1",x=unit(.99,"npc"),just="right",gp=gpar(fontsize=8))
    gr_chao_lab <- gtable_add_grob(zero1,gr_chao_grob,t=1,l=1,r=1,b=1);
   
    ## handle ordering / component selection
    ## columns:
    corder <- match(col.order, c("labels","map","tree"));
    tree   <- cbind_select(corder, zero,        gr_sampleTree,         gr_legend);
    chao   <- cbind_select(corder, gr_chao_lab, gtable_text_row(chao), zero);
    labels <- cbind_select(corder, zero,        gr_sample_labels,      zero);
    map    <- cbind_select(corder, gr_labels,   gr_heatmaps,           gr_trees, size=1);
    map    <- gtable_add_row_space(map, unit(.2,"lines"));

    ## rows
    g <- rbind_select(row.order, tree, map, chao, labels);

    ## add some spacing
    g <- gtable_add_row_space(g, unit(.1,"lines"));
    g <- gtable_add_col_space(g, unit(.1,"lines"));
    g <- gtable_add_padding(g, unit(.3,"lines"));
    
    return (invisible(g));
}

pF_main <- function() {
    require(optparse);

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
            c("-d", "--debug"),
            action="store_true",
            default=FALSE,
            help="Show debug messages"
            ),
        make_option(
            c("-n", "--min-ntu-count"),
            default=50,
            type="integer",
            help="Sum NTUs with less counts in pseudo NTU \"Other\". Default %default."
            ),
        make_option(
            "--no-split",
            action="store_true",
            default=FALSE,
            help="Do not split heatmap"
            ),
        make_option(
            c("-t", "--split-regex"),
            default="Eukaryota",
            type="character",
            help="Split heatmap using this regex on taxa. Default '%default'",
            ),
        make_option(
            c("-l", "--long-taxnames"),
            action="store_true",
            default=FALSE,
            help="Do not shorten taxa names to last two groups",
            ),
        make_option(
            c("-a", "--absolute"),
            action="store_true",
            default=FALSE,
            help="Do not scale columns to percentages"
            ),
        make_option(
            c("-m","--hclust-method"),
            default="ward",
            help="Use this method for hclust clustering. Can be:
                ward, single, complete, average, mcquitty, median or centroid.
                Default is %default."
            ),
        make_option(
            c("-r","--rows"),
            default="tree,map,chao,labels",
            help="Component rows, in order, to render (separated by commas).
                Valid terms are: tree, map, chao and labels.
                Default is %default."
            ),
        make_option(
            c("-c", "--cols"),
            default="labels,map,tree",
            help="Component columns, in order, to render (separated by commas).
                Valid terms are: labels, map and tree.
                Default is %default."
            ),
        make_option(
            c("--colors"),
            default="steelblue,indianred",
            help="Colors for heatmaps. Default is %default."
        ),
        make_option(
            c("-o","--out"),
            default="out.png",
            help="Name of output file. Must end in .png or .pdf. Default is %default."
            ),
        make_option(
            c("--aa"),
            default="gray",
            help="Type of anti-aliasing to use for PNG output. Can be one of default,
                none, gray, or subpixel. Default is %default."
            ),
        make_option(
            c("-s", "--out-size"),
            default="autoXauto",
            help="Size of output graphic in pixels (e.g. 100x100). Assumes 72 DPI for
                PDF. Using \"auto\" for a dimension will attempt to guess at suitable
                size. Default %default"
            ),
        make_option(
            c("--library-name-from-file"),
            action="store_true",
            default=FALSE,
            help="Use thee filename to derive library name instead of parsing ...report.csv"
            )
        );
    
    parser <- OptionParser(
        option_list=options,
        usage="usage: %prog [options] [files]",
        description="
Generates a heatmap plot from multiple phyloFlash result sets. For more control,
source this file from R.

Files:
        A list of files and/or directories that will be searched
        for phyloFlash results."
        );

    conf <- parse_args(parser, positional_arguments = TRUE);

    if (length(conf$args)==0) {
        print_help(parser);
        quit(status=2);
    }
    
    ## set loglevel
    if (conf$options$quiet) {
        pf_setLogLevel(1);
    } else if (conf$options$verbose) {
        pf_setLogLevel(3);
    }

    ## set debug mode
    pf_setDebug(conf$options$debug);

    info("Loading libraries");
    load_libraries();

    pf      <- read.phyloFlash(conf$args,
                               sampleNameFromMeta = !conf$options$"library-name-from-file");

    ## split by domain
    if (!conf$options$"no-split") {
        msg("Splitting data according to regex ", conf$options$"split-regex");
        pat <- strsplit(conf$options$"split-regex",",")[[1]];
        pf$data <- split_by_name(pf$data, pat);
    }
    if (!conf$options$"long-taxnames") {
        pf$data <- shorten_taxnames(pf$data);
    }

    pf$data <- merge_low_counts(pf$data,
                                thres=conf$options$"min-ntu-count");
    
    if (!conf$options$absolute) {
        msg("Rescaling counts to percentages");
        pf$data <- scale_to_percent(pf$data);
    }

    msg("Clustering...");
    pf      <- cluster(pf, method=conf$options$"hclust-method");

    ntaxa   <- sum(sapply(pf$data, nrow))
    nsample <- ncol(pf$data[[1]])
    msg("Found", ntaxa, "taxa and", nsample, "samples.")

    msg("Creating plot...");
    g       <- plot.phyloFlash(pf,
                               row.order=conf$options$rows,
                               col.order=conf$options$cols,
                               map.colors=conf$options$colors);

    msg(paste(sep="","Printing plot to \"", conf$options$out, "\"..."));

    outdim = strsplit(conf$options$"out-size","[Xx]")[[1]];
    if (outdim[1] == "auto") {
        labelwidth <- max(nchar(unlist(sapply(pf$data, rownames)))) * 5;
        width <- 80 + labelwidth + nsample * 25;
    } else {
        width=as.integer(outdim[1])
    }
    if (outdim[2] == "auto") {
        labelwidth <- max(nchar(unlist(sapply(pf$data, colnames)))) * 5;
        height <- 120 + labelwidth + ntaxa * 10;
    } else {
        height=as.integer(outdim[2])
    }

    switch(strsplit(conf$options$out, "[.]")[[1]][-1],
           png = png(file = conf$options$out,
               width=width, height=height,
               antialias=conf$options$aa),
           svg = svg(file = conf$options$out,
               width=width, height=height),
           pdf = pdf(file = conf$options$out,
               width=width/72, height=height/72)
           );

    grid.newpage();
    grid.draw(g);
    dev.off();

    invisible(1);
}
    
# if we are run as a script from the cmdline
if (!interactive()) {
    check_libraries();

    pF_main();
} else {
    load_libraries();
    msg("Loaded phyloFlash R functions. Example usage:
 pf      <- read.phyloFlash()
 pf$data <- split_by_name(pf$data, \"Euk\")
 pf$data <- shorten_taxnames(pf$data)
 pf$data <- scale_to_percent(pf$data)
 pf      <- cluster(pf, method=\"ward\")
 g       <- plot.phyloFlash(pf)
");
}
