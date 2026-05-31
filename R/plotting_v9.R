#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import GenomicRanges
#' @import dplyr
#' @import tidyr
#' @import S4Vectors
#' @import pheatmap
#' @import IRanges
#' @import ggplot2
#' @import RColorBrewer
#' @import BiocFileCache
#' @import matrixStats
#' @import gridExtra
#' @import patchwork
#' @import rtracklayer
#' @importFrom Chicago defaultSettings
#' @importFrom Chicago setExperiment
#' @importFrom Chicago readAndMerge
#' @importFrom Chicago chicagoPipeline
#' @importFrom Chicago exportResults
#' @importFrom utils read.delim str write.table
#' @importFrom graphics abline barplot boxplot legend lines par rect text polygon
#' @importFrom grDevices dev.off palette.colors pdf png rgb
#' @importFrom stats density
#' @importFrom stringr str_remove
#' @importFrom matrixStats colMedians
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices adjustcolor boxplot.stats
## usethis namespace: end
NULL

## global variables:
utils::globalVariables(c(
    "grl",
    "resfolder",
    "m",
    "d",
    "i",
    "a",
    "zoom",
    "oeids",
    "rmapgr",
    "baited_genes",
    "L",
    "type",
    "troubleshooting",
    "cd",
    "bg",
    "scalingL",
    "id",
    "bait",
    "Lm",
    "make_plot"
))


## Central plotting fnction:


plotInteractions <- function(
      L,
      id,
      k,
      zoom = 200000,
      rmapgr,
      ylim = NULL,
      show.legend = TRUE,
      d = NULL,
      preselected = FALSE,
      name = NULL,
      intervals = NULL,
      show.legend.intervals = TRUE,
      xlim = NULL,
      col = NULL,
      colintervals = NULL,
      colints = NULL,
      troubleshooting = FALSE,
      lwd = 2,
      lty = NULL,
      show.legend.outside = FALSE,
      cex.intervals = 0.7,
      cex.legend = 1
) {
    #' @title Plot Lineplots of Interactions
    #' @description Plots lineplots of normalized Capture-C read counts around one view point (bait) from one or many samples.
    #' Optional annotation with genomic intervals of choice and significant interactions.
    #' @param L List with the ChicagoData tables
    #' @param id View point (baitID), corresponds to a fragment id in rmapgr.
    #' @param k Smoothing factor across k restriction fragments (running average).
    #' @param zoom Default: 200000; Area around the bait in +/-bp to plot.
    #' @param rmapgr GenomicRanges object containing restriction fragments from the reference genome. 
    #' IDs must be saved in a column called 'id'.
    #' @param ylim optional; Plotting limit of the y axis.
    #' @param show.legend TRUE(default)/FALSE; should the sample legend be shown inside the plot area?
    #' @param d Default: FALSE; interactions table to highlight significant interactions as transparent rectangles (shading).
    #' By default, if d is given, highlights any interaction. If different categories of interactions should be highlighted, 
    #' these are specified in a column called 'clusters_refined'. Example usage would be 'lost', 'gained' and 'stable' interactions. 
    #' Up to 8 categories can be defined which then receive a different color each.
    #' @param preselected , Default: FALSE; are interactions in d preselected for the correct baitID?
    #' Useful in cases when ids in the interactions table do not match the ids in L.
    #' @param name optional; Name of the view point, Default: name = baitID.
    #' @param intervals optional; Intervals to be plotted, provided as a named GenomicRangesList, 
    #' useful for annotation of the plot with ChIP-Seq data.
    #' @param show.legend.intervals TRUE(default)/FALSE; should the legend for intervals be shown? 
    #' @param xlim optional; Plotting limit of the x axis, as xlim=c(30000,200000), overrides zoom.
    #' @param col Custom colors of the lines. Default: 'Okabe-Ito' with a max of 8 different samples. 
    #' @param colintervals Custom colors for interval annotation. 
    #' @param colints Custom colors for significant interaction shading.
    #' @param troubleshooting if TRUE will return messages showing how far the process is, Default: FALSE
    #' @param lwd line width, Default 2pt
    #' @param lty line type, Default: 1
    #' @param show.legend.outside If TRUE, shows the legend outside of the plot. 
    #' Note that the total plot area then becomes wider. Default: FALSE
    #' @param cex.intervals Defines the size of the intervals legend. Default: 0.7
    #' @param cex.legend Defines the size of the samples legend. Default: 1
    #' @examples
    #' # example code
    #'
    #' #define directories
    #' extdata <- system.file("extdata", package="PostChicago")
    #' chicagoOutputDir <- file.path(extdata,'results')
    #' designDir <- file.path(extdata,'designDir')
    #'
    #' #create baited_genes and rmapgr:
    #' baited_genes=baitmap2baited_genes(designDir,save=FALSE)
    #' rmapgr=rmap2rmapgr(designDir,save=FALSE)
    #'
    #' #load list with ChicagoData objects:
    #' cdlist <- loadCdList(resultsDir = chicagoOutputDir, baited_genes = baited_genes)
    #' cdlist=cdlist[-grep('rep',names(cdlist))]
    #'
    #' #standard plot:
    #' col=c('indianred1','seagreen')
    #' name='Hoxb3'
    #' id=baited_genes[baited_genes$genename==name]$re_id
    #' k=15
    #' ylim=c(0,200)
    #' plotInteractions(cdlist,id,k,ylim=ylim,show.legend = TRUE,name=name, rmapgr=rmapgr,col=col)
    #'
    #' #vary line type and width:
    #' plotInteractions(cdlist,id,k,ylim=ylim,show.legend = TRUE,name=name, rmapgr=rmapgr,col=col,
    #'     lwd=1,lty=c(1,2))
    #'
    #' #change zoom around the bait:
    #' plotInteractions(cdlist,id,k,ylim=ylim,show.legend = TRUE,name=name, rmapgr=rmapgr,col=col,
    #'     zoom=500000)
    #'
    #' #change zoom the number of restriction fragments over which the running mean is calculated:
    #' plotInteractions(cdlist,id,k=61,ylim=ylim,show.legend = TRUE,name=name, rmapgr=rmapgr,col=col,
    #'     zoom=500000)
    #'
    #' @returns Line plots of Capture-C reads surrounding one view point.
    #' @export

    # check for optional variables:

    # name to use in the plot title:

    if (!is.null(name)) {
        if (length(grep(",", name)) > 0 || length(grep("\n", name)) > 0) {
            ## only add comma before ID if name doesn't contain comma or a line break!
            main <- paste0(name, "id=", id)
        } else {
            main <- paste0(name, ",id=", id)
        }
    } else {
        main <- paste0("id=", id)
    }

    # optional zoom use in the plot title:

    if (!is.na(zoom)) {
        main <- paste0(main, ",TSS+/-", zoom / 1e+06, "Mb,k=", k)
    } else {
        main <- paste0(main, ",k=", k)
    }

    if (!is.null(d)) {
        if ("clusters_refined" %in% names(d) == FALSE) {
            d$clusters_refined <- "interaction"
        }
        clust <- unique(unlist(strsplit(d$clusters_refined, split = ",")))
    }

    if (is.null(col)) {
        colors <- palette.colors(n = length(L), recycle = TRUE)
    } else {
        colors <- col
    }

    if (is.null(lty)) {
        lty <- rep(1, length(L))
    }

    if (troubleshooting) {
        message(paste("getting bait..."))
    }
    bait <- getbait(id, rmapgr = rmapgr)

    if (troubleshooting) {
        message(paste("bait retrieved..."))
    }

    if (is.null(xlim)) {
        xlim <- c(
            GenomicRanges::start(bait) - zoom,
            GenomicRanges::end(bait) + zoom
        )
    }

    if (troubleshooting) {
        message(paste("Bait: \n", print(bait)))
    }
    
    ##Note: Added AF 31.05. 
    if (show.legend.outside) {
      par(mfrow=c(1,2))
    }

    LL <- list() ## make a list with bed files for all samples
    legend <- vector()
    all <- vector()

    for (i in 1:length(L)) {
        # for every interaction

        if (troubleshooting) {
            message(paste("sample", i))
        }
        samplename <- names(L)[i]
        legend <- c(legend, samplename)
        if (troubleshooting) {
            message(paste("getting normbed..."))
        }

        bed <- getnormbed(
            id,
            data = L[[i]],
            samplename = names(L)[i],
            L = L,
            same.chromosome = TRUE,
            rmapgr = rmapgr
        )

        if (troubleshooting) {
            message(paste("normbed retrieved, coverage..."))
            message(str(bed)) # Check structure of bed
        }

        bed <- bed2cov(bed, bait, rmapgr = rmapgr)

        if (troubleshooting) {
            message(paste("cov retrieved, getrunmean..."))
            message(str(bed)) # Check structure of bed
            message(paste0(
                "Bed is between ",
                min(GenomicRanges::start(bed)),
                " - ",
                max(GenomicRanges::start(bed))
            ))
            message(paste0("Bait is ", GenomicRanges::start(bait)))
        }

        bed <- getrunmeanbed(bed, bait, xlim, k)
        if (troubleshooting) {
            message(paste("runmean retrieved"))
        }

        LL <- c(LL, list(bed))
        all <- c(all, as.vector(bed$x))
    }
    if (troubleshooting) {
        message(paste("list made"))
    }

    ## make plots:

    # check for graph specifications
    if (is.null(ylim)) {
        ylim <- max(LL[[1]]$x, na.rm = TRUE)

        if (length(LL) > 1) {
            for (j in 2:length(LL)) {
                ylim <- c(ylim, max(LL[[j]]$x, na.rm = TRUE))
            }
        }

        ylim <- as.numeric(c(0, max(ylim, na.rm = TRUE) * 1.2))
    }

    if (!is.null(intervals)) {
        len <- length(intervals)
        ylim <- c(-ylim[2] / 20 * len, ylim[2])
    }

    if (troubleshooting) {
        message(paste("xlab doing, Bait: \n", print(bait)))
    }

    ylab <- "normalised read counts (PRPP)"
    xlab <- paste(
        as.character(GenomicRanges::seqnames(bait)),
        ":",
        xlim[1],
        "-",
        xlim[2],
        sep = ""
    )

    if (troubleshooting) {
        message(paste("xlab done, Bait: \n", print(bait)))
    }

    bed <- LL[[1]]

    plot(
        (GenomicRanges::end(bed) + GenomicRanges::start(bed)) / 2,
        bed$x,
        type = "l",
        lwd = lwd,
        xlim = xlim,
        col = colors[1],
        ylab = ylab,
        xlab = xlab,
        main = main,
        ylim = ylim,
        lty = lty[1]
    )
    for (j in 2:length(LL)) {
        bed <- LL[[j]]
        lines(
            (GenomicRanges::end(bed) + GenomicRanges::start(bed)) / 2,
            bed$x,
            type = "l",
            xlim = xlim,
            col = colors[j],
            lwd = lwd,
            lty = lty[j]
        )
    }

    ## plot bait:
    b <- getbaitplus(id, rmapgr = rmapgr)
    # rect(xleft = GenomicRanges::start(b), xright = GenomicRanges::end(b), ybottom = 0, ytop = ylim[2], col = "green",
    # border = NA)
    baitposition <- mean(c(GenomicRanges::start(b), GenomicRanges::end(b)))
    x <- c(
        baitposition - 0.02 * (xlim[2] - xlim[1]),
        baitposition + 0.02 * (xlim[2] - xlim[1]),
        baitposition
    )
    y <- c(0, 0, 0.075 * max(ylim))
    polygon(x, y, col = "darkgrey", border = NA)

    ## plot signficant interactions

    if (!is.null(d)) {
        if (preselected == FALSE) {
            ib <- d[d$baitID == id, ]
        } else {
            ib <- d
        }
        ibb <- GenomicRanges::GRanges(
            ib$seqnames_otherEnd,
            IRanges::IRanges(ib$start_otherEnd, ib$end_otherEnd)
        )
        values(ibb) <- ib
        ib <- ibb

        # cols = list(rgb(0, 0, 0, alpha = 0.1), rgb(1, 0, 0, alpha = 0.1), rgb(0, 1, 0, alpha = 0.1),
        #            rgb(0, 0, 1, alpha = 0.1),  rgb(1, 1, 0, alpha = 0.1), rgb(1, 0, 1, alpha = 0.1), rgb(0, 1, 1, alpha = 0.1), rgb(0, 0.5, 0.5, alpha = 0.1))
        # 01Feb: SK this change also removes transparancy
        # 22Feb: SK addding back transparency
        ## 24.3.: AF: transparency increased again
        if (is.null(colints)) {
            cols <- adjustcolor(palette.colors(n = 8, palette = "Okabe-Ito"), 0.2)
        } else {
            cols <- adjustcolor(colints, 0.2)
        }
        clust <- clust[order(clust)]
        # message('plotting interaction types:', paste(clust))

        for (j in 1:length(clust)) {
            cl <- clust[j]
            col <- cols[j]

            lost <- ib[grep(cl, ib$clusters_refined)]

            if (length(lost) > 0) {
                rect(
                    xleft = GenomicRanges::start(lost),
                    xright = GenomicRanges::end(lost),
                    ybottom = -0,
                    ytop = (max(all) + 1),
                    col = col,
                    border = NA
                )
            }
        }

        legend(
            "topleft",
            legend = clust,
            bty = "n",
            cex = 0.8,
            fill = unlist(cols)[1:length(clust)],
            ncol = 1
        )
    }

    ## Plot intervals

    if (!is.null(intervals)) {
        abline(h = 0)

        if (is.null(colintervals)) {
            # does this work???  colintervals = palette.colors(length(intervals), recycle = TRUE)
            colintervals <- c(
                "grey",
                "pink",
                "purple",
                "green",
                "red",
                "blue",
                "orange",
                "black",
                "darkred",
                "darkgreen",
                "cornflowerblue"
            )
        }

        grl <- intervals
        col <- colintervals[1:length(grl)]

        for (i in 1:length(grl)) {
            gr <- grl[[i]]
            gr <- gr[
                as.character(GenomicRanges::seqnames(gr)) ==
                    as.character(GenomicRanges::seqnames(bait))
            ]
            if (length(gr) > 0) {
                st <- 0.05 * ylim[2]
                step <- st * i
                if (length(grep("TAD", names(grl)[i])) > 0) {
                    rect(
                        xleft = start(gr),
                        xright = end(gr),
                        ybottom = (-st - step),
                        ytop = (-st - step) + step / i,
                        col = col[i]
                    )
                } else {
                    rect(
                        xleft = start(gr),
                        xright = end(gr),
                        ybottom = (-st - step),
                        ytop = (-st - step) + step / i,
                        col = col[i],
                        border = NA
                    )
                }
            }
        }
        if (show.legend.intervals) {
            legend(
                "bottomright",
                fill = col,
                legend = names(grl),
                bty = "n",
                cex = cex.intervals
            )
        }
    }

    if (show.legend == TRUE) {
        legend(
            "topright",
            legend = legend,
            col = colors,
            lwd = lwd,
            bty = "n",
            ncol = 1,
            cex = cex.legend,
            lty = lty
        )
    }

    ## plot the legend outside of the plot, for both: intervals and samples
    if (show.legend.outside) {
        plot(0, 0, col = "white", axes = FALSE, xlab = "", ylab = "")
        legend(
            "topleft",
            ncol = 1,
            col = colors,
            lwd = lwd,
            legend = names(L),
            bty = "n",
            cex = cex.legend,
            lty = lty
        )
        if (!is.null(intervals)) {
            legend(
                "bottomleft",
                fill = colintervals,
                legend = names(grl),
                bty = "n",
                cex = cex.intervals
            )
        }
    }
}


## Plot boxplot quantifications from interval/interactionPeak quantifications:

boxplotsCapC <- function(annotatedCapCTable, reps = FALSE, col = NULL, show_outliers = TRUE,
    show_notch = TRUE, name = NULL) {
    #' @title Boxplots of summarized interactions data.
    #' @description
    #' Boxplots of normalized interaction data quantified from either interaction peaks or oneGeneOnePeak tables.
    #' Returns 4 boxplots, containing downsampled read and Chicago score counts, total and size-normalized values.
    #' Function assumes that replicate data are saved in a column with the name structure '_rep1'.
    #' @param annotatedCapCTable Table containing Capture-C quantification either over intervals or interaction peaks.
    #' Assumes the output from quantifyWithinPeaks().
    #' @param reps TRUE/FALSE(default); should replicate or summary data be plotted? Assumes replicate data are stored
    #' in columns names '_rep1', '_rep2' etc.
    #' @param col Custom plot color.
    #' @param show_outliers TRUE(default)/FALSE; Should outliers be shown?
    #' @param show_notch TRUE(default)/FALSE
    #' @param name optional, plot names
    #' @examples
    #' # example code
    #'
    #' ##Create an rmapgr object
    #' rmapgr <- GenomicRanges::GRanges(rep(1,100), IRanges::IRanges( 1:100,2:101  ))
    #' rmapgr$id=1:length(rmapgr)
    #'
    #' ##create interaction data table with two interactions
    #' data <- data.frame(intID=c('1;2', '1;3,1;5'),
    #'                    baitID = as.character(c(1,1)),
    #'                    otherEndID=as.character(c(2, '3,5')))
    #'
    #' ##create a ChicagoData object with two samples
    #' cdlist <- list(
    #'   sample1 = data.frame(baitID=rep(1,11),
    #'                        otherEndID=c(2:12),
    #'                        N=c(sample(1:100,10,replace=TRUE),1000),
    #'                        score=c(sample(0:10,10,replace=TRUE),100)),
    #'
    #'   sample2 = data.frame(baitID=rep(1,11),
    #'                        otherEndID=c(2:12),
    #'                        N=sample(90:110,11,replace=TRUE),
    #'                        score=sample(0:10,11,replace=TRUE))
    #' )
    #'
    #'
    #' ##create replicate ChicagoData object from cdlist:
    #' cdlist.reps <- c(cdlist,cdlist)
    #' f<-function(x){
    #'   x$N=x$N*sample((7:13)/10,11,replace=TRUE)
    #'   x$score=x$score*sample((7:13)/10,11,replace=TRUE)
    #'   return(x)
    #' }
    #' cdlist.reps=lapply(cdlist.reps,FUN=f)
    #' names(cdlist.reps)=paste0(rep(names(cdlist),2),'_rep',rep(1:2,each=2))
    #' cdlist.reps=cdlist.reps[order(names(cdlist.reps))]
    #'
    #' ##run quantifyWithinPeaks
    #' annotatedCapCTable <- quantifyWithinPeaks(data, cdlist, cdlist.reps, rmapgr=rmapgr)
    #'
    #' ##plots
    #' col=c('grey','red')
    #' boxplotsCapC(annotatedCapCTable, col=col)
    #' boxplotsCapC(annotatedCapCTable, col=col,show_notch=FALSE)
    #' boxplotsCapC(annotatedCapCTable, col=rep(col,each=2),reps=TRUE)
    #'
    #' @returns Boxplots of scores and reads of each experiment of a specific feature
    #' @export

    ## default colors
    if (is.null(col)) {
        col <- palette.colors(n = 8, palette = "Okabe-Ito")
    }

    ## helper ------------------------------------------------------------
    ## keep within the function, since otherwise too many other arguments would need to be defined!
    make_plot <- function(df, cols, sample_names, title, ylab) {
        plot_df <- df[, cols, drop = FALSE] %>%
            mutate(row_id = seq_len(n())) %>%
            pivot_longer(
                cols = -row_id,
                names_to = "sample",
                values_to = "value"
            )

        plot_df$sample <- factor(
            sample_names[match(plot_df$sample, cols)],
            levels = unique(sample_names)
        )

        ## dynamic upper y-limit
        if (show_outliers) {
            ymax <- max(plot_df$value, na.rm = TRUE) * 1.10
        } else {
            whisker_max <- tapply(
                plot_df$value,
                plot_df$sample,
                function(x) boxplot.stats(x)$stats[5] # upper whisker
            )

            ymax <- max(whisker_max, na.rm = TRUE) * 1.10
        }

        ggplot(plot_df, aes(x = sample, y = value, fill = sample)) +
            geom_boxplot(
                notch = show_notch,
                outlier.shape = if (show_outliers) 19 else NA,
                outlier.alpha = if (show_outliers) 0.5 else 0
            ) +
            scale_fill_manual(
                values = rep(col, length.out = length(unique(sample_names)))
            ) +
            scale_y_continuous(
                limits = c(0, ymax),
                expand = c(0, 0)
            ) +
            labs(
                x = NULL,
                y = ylab, # <- restored y-axis labels
                title = title
            ) +
            theme_bw(base_size = 12) +
            theme(
                legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(
                    angle = 45,
                    hjust = 1,
                    vjust = 1
                ),
                plot.title = element_text(
                    face = "bold",
                    hjust = 0.5
                )
            )
    }
    ## extract data ------------------------------------------------------
    reads <- annotatedCapCTable[, grep("N_downsampled$", names(annotatedCapCTable)), drop = FALSE]
    mean_reads <- annotatedCapCTable[, grep("N_downsampled_mean", names(annotatedCapCTable)), drop = FALSE]
    scores <- annotatedCapCTable[, grep("_score$", names(annotatedCapCTable)), drop = FALSE]
    mean_scores <- annotatedCapCTable[, grep("_score_mean", names(annotatedCapCTable)), drop = FALSE]

    rep_cols <- grep("rep", names(annotatedCapCTable))

    ## choose replicate/non-replicate columns ----------------------------
    if (!reps) {
        if (length(rep_cols) > 0) {
            reads <- reads[, -grep("rep", names(reads)), drop = FALSE]
            mean_reads <- mean_reads[, -grep("rep", names(mean_reads)), drop = FALSE]
            scores <- scores[, -grep("rep", names(scores)), drop = FALSE]
            mean_scores <- mean_scores[, -grep("rep", names(mean_scores)), drop = FALSE]
        }

        sample_names <- sapply(strsplit(names(reads), "_N_"), FUN = function(x) x[1])
    } else {
        if (length(rep_cols) == 0) {
            message(
                "Replicate data must be included in columns named <mysample>_rep<n>"
            )
            return(NULL)
        }

        reads <- reads[, grep("rep", names(reads)), drop = FALSE]
        mean_reads <- mean_reads[, grep("rep", names(mean_reads)), drop = FALSE]
        scores <- scores[, grep("rep", names(scores)), drop = FALSE]
        mean_scores <- mean_scores[, grep("rep", names(mean_scores)), drop = FALSE]

        # sample_names <- sapply(strsplit(names(reads), "_N_"), `[`, 1)
        # rep_names <- sapply(strsplit(names(reads), "_N_"), `[`, 2)
        sample_names <- sapply(strsplit(names(reads), "_N_"), FUN = function(x) x[1])
    }

    ## plots -------------------------------------------------------------
    p1 <- make_plot(
        reads,
        names(reads),
        sample_names,
        paste0(name, "reads"),
        "CapC downsampled_reads"
    )

    p2 <- make_plot(
        mean_reads,
        names(mean_reads),
        sample_names,
        paste0(name, "mean reads"),
        "CapC downsampled_reads size-norm"
    )

    p3 <- make_plot(
        scores,
        names(scores),
        sample_names,
        paste0(name, "scores"),
        "Chicago scores"
    )

    p4 <- make_plot(
        mean_scores,
        names(mean_scores),
        sample_names,
        paste0(name, "mean scores"),
        "Chicago scores size-norm"
    )

    ## 2x2 arrangement
    final_plot <- (p1 | p2) / (p3 | p4)

    return(final_plot)
}


## function to plot aggregate peaks from Capture-C analysis matrices:

plotAggregatePeaks <- function(
      Lm,
      mainprefix = NULL,
      col = NULL,
      ylim = NULL,
      lty = NULL,
      xlim = NULL,
      k = 3,
      ylab = "normalised reads",
      plotwhat = "all"
) {
    #' @title Plot Aggregate Interaction Profiles
    #' @description Takes a list of matrices and plots aggregate interaction signal, similar for metaprofiles in ChIP-Seq.
    #' @param Lm List of matrices containing interaction reads or scores.
    #' @param mainprefix Optional; text added to the plot title.
    #' @param col Custom colors for aggregate plots. Default: 4 colors.
    #' @param lty Optional; Line type for the aggregate plots
    #' @param ylim Optional
    #' @param xlim Optional; should be added in restriction fragments surrounding the interval/peak center. Example: c(-50,50).
    #' Default: All restriction fragments represented in the matrix.
    #' @param k Numeric; over how many fragments should the running mean be calculated? Default: 3
    #' @param ylab Optional; Default: 'normalised reads'
    #' @param plotwhat Defines what types of plots will be plotted. Options are 'all','mean','median'
    #' 'mean' and 'median' plot running averages over k fragments.
    #' Default: 'all': Four plots representing mean and median read counts both smoothed and unsmoothed over k fragments.
    #' @examples
    #' # example code
    #'
    #' ##create list of matrices:
    #'
    #' m1<-matrix(c(c(1:10,11,9:0)*sample((8:1.2)/10,21,replace=TRUE),
    #'             c(1:10,11,9:0)*sample((8:1.2)/10,21,replace=TRUE),
    #'             c(1:10,11,9:0)*sample((8:1.2)/10,21,replace=TRUE)),nrow=3)
    #' m2<-matrix(c(c(21:30,39,29:20)*sample((8:1.2)/10,21,replace=TRUE),
    #'             c(21:30,39,29:20)*sample((8:1.2)/10,21,replace=TRUE),
    #'             c(21:30,39,29:20)*sample((8:1.2)/10,21,replace=TRUE)),nrow=3)
    #' rownames(m1)=c('1;2','2;3','3;4')
    #' rownames(m2)=rownames(m1)
    #'
    #' Lm=list(m1=m1,m2=m2,m1_con=m1*0.2,m2_con=m2*0.1)
    #'
    #' ##plot median profiles::
    #' plotAggregatePeaks(Lm,plotwhat='median',k=3,ylab='reads',
    #'                   col=rep(c('red','blue'),2), lty=c(1,1,3,3))
    #'
    #' ##plot mean profiles::
    #' plotAggregatePeaks(Lm,plotwhat='mean',k=3,ylab='reads',
    #'                   col=rep(c('red','blue'),2), lty=c(1,1,3,3))
    #'
    #' @return Plot containing aggregate profiles of interactions.
    #' @export
    ## arguments: optional: xlim: Region (in restriction fragments) surrounding the peak summit to plot.
    ## Default: All restriction fragments represented in the matrix. plotwhat: What plots to plot? Options are
    ## 'all','mean','median' Default: four plots representing mean and median read counts both smoothed and
    ## unsmoothed.  Mean and median plot running averages over k fragments.  which: which of the plots to
    ## plot?

    ## removed: AF 25.08.2025 if (is.null(ylim)){ ylim=c(0,1.5) }
    if (is.null(col)) {
        col <- rep(c("blue", "red", "green", "grey"), 2)
    }
    if (is.null(lty)) {
        lty <- rep(c(1, 2), each = 4)
    }

    ## columns for the legend (up to 6 samples per column)
    ncol <- 1
    if (length(Lm) > 8) {
        ncol <- ceiling(length(Lm) / 8)
    }

    if (plotwhat == "all") {
        for (type in c("mean", "median")) {
            if (type == "mean") {
                main <- paste(mainprefix, "mean enrichment", sep = "")
            } else {
                main <- paste(mainprefix, "median enrichment", sep = "")
            }

            ## aggregate plots:
            Lx <- list()
            for (i in 1:length(Lm)) {
                if (type == "mean") {
                    Lx <- c(Lx, list(colMeans(Lm[[i]])))
                } else {
                    Lx <- c(Lx, list(matrixStats::colMedians(Lm[[i]])))
                }
            }
            names(Lx) <- names(Lm)

            Lxx <- Lx
            for (i in 1:length(Lxx)) {
                Lxx[[i]] <- runmean(Rle(Lxx[[i]]), k = k, endrule = "constant")
            }

            ## plot:

            dis <- (length(Lx[[1]]) - 1) / 2
            x <- -dis:dis

            if (!is.null(xlim)) {
                for (i in 1:length(Lx)) {
                    a <- Lx[[i]]
                    df <- data.frame(a, x)
                    Lx[[i]] <- df[df$x %in% xlim[1]:xlim[2], ]$a

                    a <- Lxx[[i]]
                    df <- data.frame(a, x)
                    Lxx[[i]] <- df[df$x %in% xlim[1]:xlim[2], ]$a
                }
                x <- xlim[1]:xlim[2]
            }

            ## Added AF 26.08.2025: Data-driven ylim:
            ylim <- c(0, max(unlist(Lx)) * 1.1)
            if (max(ylim) < 1) {
                ylim <- c(0, 1)
            }

            plot(
                x,
                Lx[[1]],
                type = "l",
                ylim = ylim,
                lwd = 2,
                ylab = ylab,
                xlab = "distance from summit [RE frags]",
                col = col[1],
                main = main
            )
            for (i in 1:length(Lx)) {
                lines(x, Lx[[i]], lwd = 2, col = col[i], lty = lty[i])
            }
            legend(
                "topright",
                cex = 0.8,
                bty = "n",
                lwd = 2,
                col = col,
                legend = names(Lm),
                lty = lty,
                ncol = ncol
            )
            legend(
                "topleft",
                bty = "n",
                cex = 0.8,
                legend = c(
                    paste("n_ints=", nrow(Lm[[1]]), sep = ""),
                    paste0(
                        "n_prom=",
                        length(unique(strsplit_string(row.names(Lm[[1]]), s2 = ";")[
                            1:length(strsplit_string(row.names(Lm[[1]]), s2 = ";")) %% 2 == 1
                        ]))
                    )
                )
            )

            plot(
                x,
                Lxx[[1]],
                type = "l",
                ylim = ylim,
                lwd = 2,
                ylab = ylab,
                xlab = "distance from summit [RE frags]",
                col = col[1],
                main = paste0(main, " k=", k)
            )
            for (i in 1:length(Lxx)) {
                lines(x, Lxx[[i]], lwd = 2, col = col[i], lty = lty[i])
            }
            legend(
                "topright",
                cex = 0.7,
                bty = "n",
                lwd = 2,
                col = col,
                legend = names(Lm),
                lty = lty,
                ncol = ncol
            )
            legend(
                "topleft",
                bty = "n",
                cex = 0.8,
                legend = c(
                    paste("n_ints=", nrow(Lm[[1]]), sep = ""),
                    paste0(
                        "n_prom=",
                        length(unique(strsplit_string(row.names(Lm[[1]]), s2 = ";")[
                            1:length(strsplit_string(row.names(Lm[[1]]), s2 = ";")) %% 2 == 1
                        ]))
                    )
                )
            )
        }
    } else {
        type <- plotwhat

        if (type == "mean") {
            main <- paste(mainprefix, "mean enrichment", sep = "")
        } else {
            main <- paste(mainprefix, "median enrichment", sep = "")
        }

        ## aggregate plots:
        Lx <- list()
        for (i in 1:length(Lm)) {
            if (type == "mean") {
                Lx <- c(Lx, list(colMeans(Lm[[i]])))
            } else {
                Lx <- c(Lx, list(matrixStats::colMedians(Lm[[i]])))
            }
        }
        names(Lx) <- names(Lm)

        Lxx <- Lx
        for (i in 1:length(Lxx)) {
            Lxx[[i]] <- runmean(Rle(Lxx[[i]]), k = k, endrule = "constant")
        }

        ## plot:

        dis <- (length(Lx[[1]]) - 1) / 2
        x <- -dis:dis

        if (!is.null(xlim)) {
            for (i in 1:length(Lx)) {
                a <- Lx[[i]]
                df <- data.frame(a, x)
                Lx[[i]] <- df[df$x %in% xlim[1]:xlim[2], ]$a

                a <- Lxx[[i]]
                df <- data.frame(a, x)
                Lxx[[i]] <- df[df$x %in% xlim[1]:xlim[2], ]$a
            }
            x <- xlim[1]:xlim[2]
        }

        ## Added AF 26.08.2025: Data-driven ylim:
        ylim <- c(0, max(unlist(Lx)) * 1.1)

        plot(
            x,
            Lxx[[1]],
            type = "l",
            ylim = ylim,
            lwd = 2,
            ylab = ylab,
            xlab = "distance from summit [RE frags]",
            col = col[1],
            main = paste0(main, " k=", k)
        )
        for (i in 1:length(Lxx)) {
            lines(x, Lxx[[i]], lwd = 2, col = col[i], lty = lty[i])
        }
        legend(
            "topright",
            cex = 0.7,
            bty = "n",
            lwd = 2,
            col = col,
            legend = names(Lm),
            lty = lty,
            ncol = ncol
        )
        legend(
            "topleft",
            bty = "n",
            cex = 0.8,
            legend = c(
                paste("n_ints=", nrow(Lm[[1]]), sep = ""),
                paste0(
                    "n_prom=",
                    length(unique(strsplit_string(row.names(Lm[[1]]), s2 = ";")[
                        1:length(strsplit_string(row.names(Lm[[1]]), s2 = ";")) %% 2 == 1
                    ]))
                )
            )
        )
    }
}


## AF, 24.03.2026: Changed name to makeQCplots, changed output in a single grid using gridExtra::grid.arrange

makeQCplots <- function(ints, folder = NULL, fontsize = 10) {
    #' @title Quality Control Plots
    #' @description
    #' Will create and save QC heatmaps of downsampled and raw data from the interactions matrix by reads
    #' and the scores of each replicate. 
    #' @param ints Interactions Table containing all significant interactions. Typically created using makeIntsTable().
    #' @param folder Default: Current folder; Path to the directory in which the plots are saved. 
    #' @param fontsize Default: 10.
    #' @examples
    #' # example code:
    #' extdata <- system.file("extdata", package="PostChicago")
    #' outputDir <- file.path(extdata,'postchicago')
    #'
    #' #read in interactions:
    #' ints=read.delim(paste0(outputDir,'/ints_all.txt'),stringsAsFactors=FALSE)
    #'
    #' #assess interactions:
    #' makeQCplots(ints, outputDir)
    #'
    #' @returns A pdf file containing QC heatmaps: pheatmaps_QC.pdf.
    #' @export

    if (is.null(folder)) {
        folder <- getwd()
    }

    m <- ints[, grep("N[.]", names(ints))]
    mraw <- log2(m[, -grep("downsampled", names(m))] + 1)
    mnorm <- log2(m[, grep("downsampled", names(m))] + 1)

    m2 <- ints[, grep("score", names(ints))]
    msc <- asinh(m2[, grep("rep", names(m2))])

    ## reads:
    plt1 <- pheatmap::pheatmap(
        mraw,
        show_rownames = FALSE,
        main = "log2(raw_reads+1)",
        border_color = NA,
        silent = TRUE,
        fontsize = fontsize,
        angle_col = 45
    )[[4]]

    plt2 <- pheatmap::pheatmap(
        mnorm,
        show_rownames = FALSE,
        main = "log2(downsampled_reads+1)",
        border_color = NA,
        silent = TRUE,
        fontsize = fontsize,
        angle_col = 45
    )[[4]]

    plt3 <- pheatmap::pheatmap(
        cor(mnorm, use = "pairwise.complete.obs"),
        main = "frags log2(downsampled_reads+1), \n cor pearson",
        border_color = NA,
        silent = TRUE,
        fontsize = fontsize,
        fontsize_row = 4,
        angle_col = 45
    )[[4]]

    plt4 <- pheatmap::pheatmap(
        cor(mnorm, method = "spearman", use = "pairwise.complete.obs"),
        main = "frags log2(downsampled_reads+1), \n cor spearman",
        border_color = NA,
        silent = TRUE,
        fontsize = fontsize,
        fontsize_row = 4,
        angle_col = 45
    )[[4]]

    ## scores:
    plt5 <- pheatmap::pheatmap(
        msc,
        show_rownames = FALSE,
        main = "asinh(scores)",
        border_color = NA,
        silent = TRUE,
        fontsize = fontsize,
        angle_col = 45
    )[[4]]

    plt6 <- pheatmap::pheatmap(
        cor(msc, use = "pairwise.complete.obs"),
        main = "frags asinh(scores), \n cor pearson",
        border_color = NA,
        silent = TRUE,
        fontsize = fontsize,
        fontsize_row = 4,
        angle_col = 45
    )[[4]]

    plt7 <- pheatmap::pheatmap(
        cor(msc, method = "spearman", use = "pairwise.complete.obs"),
        main = "frags asinh(scores), \n cor spearman",
        border_color = NA,
        silent = TRUE,
        fontsize = fontsize,
        fontsize_row = 4,
        angle_col = 45
    )[[4]]

    ## plot all plots, starting with global cor plots:
    g <- gridExtra::grid.arrange(
        plt3,
        plt4,
        plt1,
        plt2,
        plt6,
        plt7,
        plt5,
        nrow = 2
    )

    ## save:
    ggplot2::ggsave(
        paste0(folder, "/pheatmaps_QC.pdf"),
        g,
        width = 20,
        height = 10
    )
}


###########################################################
## Function to plot statistics from the ChicagoData lists
###########################################################

## helper functions:
nGenesSig <- function(x) {
    return(length(table(x[x$score >= 5, ]$baitID)))
}
nreads <- function(x) {
    return(sum(x$N))
}
fracReadsInSig <- function(x) {
    return(round(sum(x[x$score >= 5, ]$N) / sum(x$N) * 100, 2))
}

## AF: 8.2.26: Changed according to reviewer's comments, to adjust the formatting if only
## the bar plots are plotted
## Also changed the pdf file format to be longer (7/10) and cex=0.85 (from 0.75)

plotSigIntsStats <- function(
      L,
      col = NULL,
      lty = NULL,
      filename = NULL,
      plotDistribution = TRUE,
      plotExamples = FALSE
) {
    #' @title Plot Interaction Statistics
    #' @description Takes a list of ChicagoData tables and plots a per sample statistical summary.
    #' Includes the total number of reads mapped to baits, total genes with significant interactions,
    #' total unique and significant interactions, % of reads within interactions, as well as
    #' corresponding distributions per capture bait.
    #' @param L list of ChicagoData tables for which to extract statistics.
    #' @param col Custom colors for the plots. Default: Okabe-Ito, 8 colors.
    #' @param lty optional; vector defining line types.
    #' @param filename optional; if provided, the plots will be saved under filename;
    #' must have the extension '.pdf'
    #' @param plotDistribution if TRUE (default), distributions of reads and interactions per bait
    #' are plotted. 
    #' @param plotExamples Default: FALSE; if TRUE, example distributions will be plotted based on the data
    #' from Mahara et al. 
    #' @examples
    #' # example code:
    #' extdata <- system.file("extdata", package="PostChicago")
    #' dataPath <- file.path(extdata, 'dataPath')
    #' chicagoOutputDir <- file.path(extdata,'results')
    #' designDir <- file.path(extdata,'designDir')
    #' intDir <- file.path(extdata,'intervals')
    #' outputDir <- file.path(extdata,'postchicago')
    #' baited_genes=baitmap2baited_genes(designDir,save=FALSE)
    #' cdlist <- loadCdList(resultsDir = chicagoOutputDir, baited_genes = baited_genes)
    #' L <- cdlist[-grep('rep',names(cdlist))]
    #' plotSigIntsStats(L,plotDistribution=FALSE,plotExamples=TRUE)  ##this plots example distributions
    #' @returns Plots with statistics of interactions for each sample and bait.
    #' @export

    if (is.null(col)) {
        col <- palette.colors(palette = "Okabe-Ito")
    }

    if (is.null(lty)) {
        lty <- rep(1, length(L))
    }

    if (!is.null(filename)) {
        pdf(filename, width = 7, height = 10)
        par(mfrow = c(3, 3))
    } else if (plotDistribution || plotExamples) {
        par(mfrow = c(3, 3))
    } else {
        par(mfrow = c(2, 3))
    }

    ## a) total reads

    bp <- barplot(
        sapply(L, FUN = nreads),
        las = 3,
        main = "total reads mapped to baits",
        col = col,
        ylim = c(0, 1.4 * max(sapply(L, FUN = nreads)))
    )
    text(
        bp,
        sapply(L, FUN = nreads) * 1.2,
        labels = sapply(L, FUN = nreads),
        cex = 0.85
    )

    ## b) total genes with sig ints

    bp <- barplot(
        sapply(L, FUN = nGenesSig),
        las = 3,
        main = "total genes with significant ints",
        ylim = c(0, 1.4 * max(sapply(L, FUN = nGenesSig))),
        col = col
    )
    text(
        bp,
        sapply(L, FUN = nGenesSig) * 1.2,
        labels = sapply(L, FUN = nGenesSig),
        cex = 0.85
    )

    ## c) total unique ints

    bp <- barplot(
        sapply(L, FUN = function(x) nrow(x)),
        las = 3,
        main = "total unique ints",
        ylim = c(0, 1.4 * max(sapply(L, FUN = function(x) nrow(x)))),
        col = col
    )
    text(
        bp,
        sapply(L, FUN = function(x) nrow(x)) * 1.3,
        labels = sapply(L, FUN = function(x) nrow(x)),
        cex = 0.85
    )

    ## d) total sig ints

    bp <- barplot(
        sapply(L, FUN = function(x) nrow(x[x$score >= 5, ])),
        las = 3,
        main = "total significant ints",
        ylim = c(
            0,
            1.4 * max(sapply(L, FUN = function(x) nrow(x[x$score >= 5, ])))
        ),
        col = col
    )
    text(
        bp,
        sapply(L, FUN = function(x) nrow(x[x$score >= 5, ])) * 1.3,
        labels = sapply(L, FUN = function(x) nrow(x[x$score >= 5, ])),
        cex = 0.85
    )

    ## e) % reads within sig ints

    bp <- barplot(
        sapply(L, FUN = fracReadsInSig),
        las = 3,
        main = "% reads within sig ints",
        ylim = c(0, 1.4 * max(sapply(L, FUN = fracReadsInSig))),
        col = col
    )
    text(
        bp,
        sapply(L, FUN = fracReadsInSig) * 1.3,
        labels = sapply(L, FUN = fracReadsInSig),
        cex = 0.85
    )

    ## if distribution should not be plotted, plot examples

    if (plotDistribution) {
        ## f) distribution total reads per bait

        nreadsPerBait <- function(x) {
            return(sapply(unique(x$baitID), FUN = function(y) {
                sum(x[x$baitID == y, ]$N)
            }))
        }
        dat <- lapply(L, FUN = nreadsPerBait)

        plot(
            density(dat[[1]]),
            lwd = 2,
            xlab = "total reads per bait",
            main = "total reads per bait distribution",
            col = col[1],
            lty = lty[1],
            xlim = c(0, 1.25 * max(unlist(dat)))
        )
        for (i in 2:length(L)) {
            lines(density(dat[[i]]), col = col[i], lwd = 2, lty = lty[i])
        }
        legend(
            "topright",
            cex = 0.85,
            bty = "n",
            col = col[1:length(L)],
            legend = paste0(
                names(L),
                ",median=",
                round(sapply(dat, FUN = "median"), 0)
            ),
            lwd = 2,
            lty = lty
        )

        ## g) distribution total ints per bait

        dat <- lapply(L, FUN = function(x) table(x$baitID))
        plot(
            density(dat[[1]]),
            lwd = 2,
            xlab = "unique ints per bait",
            main = "unique ints per bait distribution",
            col = col[1],
            lty = lty[1],
            xlim = c(0, 1.25 * max(unlist(dat)))
        )
        for (i in 2:length(L)) {
            lines(density(dat[[i]]), col = col[i], lwd = 2, lty = lty[i])
        }
        legend(
            "topright",
            cex = 0.85,
            bty = "n",
            col = col[1:length(L)],
            legend = paste0(
                names(L),
                ",median=",
                round(sapply(dat, FUN = "median"), 0)
            ),
            lwd = 2,
            lty = lty
        )

        ## h) distribution sig ints per bait

        dat <- lapply(L, FUN = function(x) table(x[x$score >= 5, ]$baitID))
        plot(
            density(dat[[1]]),
            lwd = 2,
            ylim = c(0, 0.05),
            xlab = "significant ints per bait",
            main = "sig ints per bait distribution",
            col = col[1],
            lty = lty[1]
        )
        for (i in 2:length(L)) {
            lines(density(dat[[i]]), col = col[i], lwd = 2, lty = lty[i])
        }
        legend(
            "topright",
            cex = 0.85,
            bty = "n",
            col = col[1:length(L)],
            legend = paste0(
                names(L),
                ",median=",
                round(sapply(dat, FUN = "median"), 0)
            ),
            lwd = 2,
            lty = lty
        )

        ## i) distribution reads per sig int per bait (log2)

        nSigReadsPerIntPerBait <- function(x) {
            return(sapply(unique(x$baitID), FUN = function(y) {
                median(x[x$baitID == y & x$score >= 5, ]$N, na.rm = TRUE)
            }))
        }
        dat <- lapply(L, FUN = nSigReadsPerIntPerBait)
        ## AF, 22.4.26: Changed to make sure that NAs (meaning the bait
        ## has no significant ints!) are removed completely from the vector, not replaced with 0
        ## dat <- lapply(dat, FUN = function(x) replace(x, is.na(x), 0))
        dat <- lapply(dat, FUN = function(x) x[!is.na(x)])


        plot(
            density(dat[[1]]),
            lwd = 2,
            xlab = "Median reads in sig ints per bait",
            main = "Median reads in sig ints per bait distribution",
            col = col[1],
            lty = lty[1]
        )
        for (i in 2:length(L)) {
            lines(density(dat[[i]]), col = col[i], lwd = 2, lty = lty[i])
        }
        legend(
            "topright",
            cex = 0.85,
            bty = "n",
            col = col[1:length(L)],
            legend = paste0(
                names(L),
                ",median=",
                round(sapply(dat, FUN = "median"), 0)
            ),
            lwd = 2,
            lty = lty
        )
    }

    if (plotExamples) {
        ## define directories from the PostChicago package:
        extdata <- system.file("extdata", package = "PostChicago")
        chicagoOutputDir <- file.path(extdata, "results")

        lty <- c(1, 1)
        col <- col[1:2]

        ## f) distribution total reads per bait

        load(paste0(chicagoOutputDir, "/fnreadsPerBait"))

        plot(
            density(dat[[1]]),
            lwd = 2,
            xlab = "total reads per bait",
            main = "Example: \n total reads per bait distribution",
            col = col[1],
            lty = lty[1],
            xlim = c(0, 1.25 * max(unlist(dat)))
        )
        for (i in 2:length(dat)) {
            lines(density(dat[[i]]), col = col[i], lwd = 2, lty = lty[i])
        }
        legend(
            "topright",
            cex = 0.75,
            bty = "n",
            col = col[1:length(dat)],
            legend = paste0(
                names(dat),
                ",median=",
                round(sapply(dat, FUN = "median"), 0)
            ),
            lwd = 2,
            lty = lty
        )

        ## g) distribution total ints per bait

        load(paste0(chicagoOutputDir, "/gDistrTotIntsPerBait"))

        plot(
            density(dat[[1]]),
            lwd = 2,
            xlab = "unique ints per bait",
            main = "Example: \n unique ints per bait distribution",
            col = col[1],
            lty = lty[1],
            xlim = c(0, 1.25 * max(unlist(dat)))
        )
        for (i in 2:length(dat)) {
            lines(density(dat[[i]]), col = col[i], lwd = 2, lty = lty[i])
        }
        legend(
            "topright",
            cex = 0.75,
            bty = "n",
            col = col[1:length(dat)],
            legend = paste0(
                names(dat),
                ",median=",
                round(sapply(dat, FUN = "median"), 0)
            ),
            lwd = 2,
            lty = lty
        )

        ## h) distribution sig ints per bait

        load(paste0(chicagoOutputDir, "/hSigIntsPerBait"))

        plot(
            density(dat[[1]]),
            lwd = 2,
            ylim = c(0, 0.05),
            xlab = "significant ints per bait",
            main = "Example: \n sig ints per bait distribution",
            col = col[1],
            lty = lty[1]
        )
        for (i in 2:length(dat)) {
            lines(density(dat[[i]]), col = col[i], lwd = 2, lty = lty[i])
        }
        legend(
            "topright",
            cex = 0.75,
            bty = "n",
            col = col[1:length(dat)],
            legend = paste0(
                names(dat),
                ",median=",
                round(sapply(dat, FUN = "median"), 0)
            ),
            lwd = 2,
            lty = lty
        )

        ## i) distribution reads per sig int per bait (log2)

        load(paste0(chicagoOutputDir, "/iDistReadsPerSigIntPerBait"))

        plot(
            density(dat[[1]]),
            lwd = 2,
            xlab = "Median reads in sig ints per bait",
            main = "Example: \n Median reads in sig ints per bait distribution",
            col = col[1],
            lty = lty[1]
        )
        for (i in 2:length(dat)) {
            lines(density(dat[[i]]), col = col[i], lwd = 2, lty = lty[i])
        }
        legend(
            "topright",
            cex = 0.75,
            bty = "n",
            col = col[1:length(dat)],
            legend = paste0(
                names(dat),
                ",median=",
                round(sapply(dat, FUN = "median"), 0)
            ),
            lwd = 2,
            lty = lty
        )
    }
    if (!is.null(filename)) {
        dev.off()
    }
}
