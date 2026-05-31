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


colMeans <- function(m) {
    ncols <- ncol(m)
    v <- vector()
    for (i in 1:ncols) {
        a <- m[, i]
        num <- mean(a)
        v <- c(v, num)
    }
    return(v)
}


strsplit_string <- function(x, s1 = NULL, s2 = NULL) {
    #' @title String Cutter
    #' @description Cuts regular expressions at the beginning and/or the end of a string and only returns the middle string
    #' @param x String or vector of strings to be split
    #' @param s1 a regular expression before the desired string
    #' @param s2 a regular expression after the desired string
    #' @return x substring which is between s1 and s2
    #' @examples
    #' # basic usage of strsplit_string:
    #' strsplit_string('myfunnystring',s1='my',s2='string')
    #' @export

    a <- x
    if (!is.null(s1)) {
        a <- unlist(strsplit(a, split = s1))
        a <- a[1:length(a) %% 2 == 0]
    }
    if (!is.null(s2)) {
        a <- unlist(strsplit(a, split = s2))
    }
    return(a)
}


getnormbed <- function(
      id,
      data,
      samplename,
      L,
      same.chromosome = FALSE,
      p = 0,
      rmapgr,
      sk = NULL
) {
    #' title Get Normalised Bed From ChicagoData Objects
    #' description Normalises the read counts between samples in one Capture-C experiment by bait-mapped library size. Returns a
    #' GenomicRanges object containing normalized read counts for one bait.
    #' param id baitID
    #' param data ChicagoData, 'cd(at)x', for a specific sample
    #' param samplename name of the sample, has to be the same as the name of the samples in L
    #' param L list of ChicagoData tables, 'cd(at)x', for which to calculate scaling factors sk
    #' param same.chromosome should only reads from the same chromosome be taken?; Default: FALSE
    #' param p pseudocount to add before normalization; Default:0
    #' param rmapgr GenomicRanges object containing restriction fragments. IDs must be saved in a column called 'id'.
    #' param sk optional, scaling factors can be provided separately.
    #' return normalised read counts in Genomic Ranges format
    #' @noRd

    if (is.null(sk)) {
        sk <- getsk(L)
    }
    names(sk) <- names(L)
    s <- sk[samplename]

    a <- data[data$baitID == id, ]
    a <- a[order(a$otherEndID), ]
    bed <- rmapgr[rmapgr$id %in% a$otherEndID]
    S4Vectors::values(bed) <- cbind(
        as.data.frame(S4Vectors::values(bed)),
        as.data.frame(a)[, grep("N", names(a))]
    )
    summary(bed$id == a$otherEndID)
    bait <- IRanges::resize(rmapgr[rmapgr$id %in% a$baitID], 1, fix = "center")

    if (same.chromosome == TRUE) {
        bed <- bed[GenomicRanges::seqnames(bed) == GenomicRanges::seqnames(bait)]
    }
    bed$intID <- paste(bait$id, bed$id, sep = ";")

    ## normalize N using the scaling factors:
    bed$N <- (bed$N + p) * s

    return(bed)
}


bed2cov <- function(bed, bait, rmapgr) {
    #' title Readcounts to Genome Coverage
    #' description Converts a GenomicRanges object derived by getnormbed()
    #' (containing normalised read counts in a metadata column)
    #' into a continuous genome coverage GenomicRanges object that includes also 0 values
    #' param bed getnormbed()-derived object to be converted
    #' param bait baitID
    #' param rmapgr GenomicRanges object containing restriction fragments. IDs must be saved in a column called 'id'.
    #' return GenomicRanges object containing continuous genome coverage information for one bait
    #' @noRd

    bed <- bed[order(GenomicRanges::start(bed), decreasing = FALSE)]

    ## cis interactions: for runmean include also all frags that have 0 interactions:
    b <- rmapgr[GenomicRanges::seqnames(rmapgr) == GenomicRanges::seqnames(bait)]
    min <- min(
        bed[GenomicRanges::seqnames(bed) == GenomicRanges::seqnames(bait)]$id
    )
    max <- max(
        bed[GenomicRanges::seqnames(bed) == GenomicRanges::seqnames(bait)]$id
    )
    b <- b[b$id %in% c(min:max)]
    b$N <- 0
    b[b$id %in% bed$id]$N <- bed$N

    ## trans interactions separately

    # trans=bed[seqnames(bed) %in% seqnames(bait)==FALSE]
    trans <- bed[
        !as.character(GenomicRanges::seqnames(bed)) %in%
            as.character(GenomicRanges::seqnames(bait))
    ]
    S4Vectors::values(trans) <- S4Vectors::values(trans)[names(S4Vectors::values(
        b
    ))]

    ## combine cis and trans interactions:
    bed <- c(b, trans)

    return(bed)
}


getrunmeanbed <- function(bed, bait, xlim, k) {
    #' title Get Running Average Capture-C Data
    #' description Calculates runmean of reads in a GenomicRanges object containing continuous coverage
    #' information in the metadata column 'N', smoothing over k fragments. Used for plotting Capture-C profiles
    #' as line plots.
    #' Common usage:
    #' getrunmeanbed(bed2cov(getnrombed(id,...)))
    #' param bed GenomicRAnges object with normalised counts to be averaged
    #' param bait baitID
    #' param xlim area on the chromosome in (start,end) to be averaged
    #' param k amount of fragments to be averaged in one run
    #' return GenomicRanges object with running average
    #' @noRd

    x <- (S4Vectors::runmean(S4Vectors::Rle(bed$N), k = k, endrule = "constant"))

    ylab <- "% reads per promoter"
    xlab <- paste(
        as.character(GenomicRanges::seqnames(bait)),
        ":",
        xlim[1],
        "-",
        xlim[2],
        sep = ""
    )

    bed$x <- x
    bed <- bed[
        start(bed) > (xlim[1] - (0.1 * abs(xlim[1] - xlim[2]))) &
            end(bed) <
                (xlim[2]) +
                    (0.1 *
                        abs(
                            xlim[1] -
                                xlim[2]
                        ))
    ]
    return(bed)
}


getbait <- function(id, rmapgr) {
    #' title Get Bait
    #' description Returns a GenomicRanges object with the position of the bait
    #' param id baitID
    #' param rmapgr  GenomicRanges object containing restriction fragments. IDs must be saved in a column called 'id'.
    #' return GenomicRanges object with bait position
    #' @noRd

    bait <- IRanges::resize(rmapgr[rmapgr$id %in% id], 1, fix = "center")
    return(bait)
}


baitmap2baited_genes <- function(designDir = "designDir", save = FALSE) {
    #' @title Convert baitmap to baited_genes
    #' @description Searches for the Chicago baitmap file in designDir and converts to a GRanges object.
    #' Only executed if baited_genes.rds doesn't exist yet, otherwise loads the existing file
    #' @param designDir  Path to the directory with Chicago design files. Default: 'designDir'
    #' @param save If TRUE, saves baited_genes as baited_genes.rds in designDir. Default: FALSE
    #' @examples
    #' # Example running baitmap2baited_genes
    #' extdata <- system.file("extdata", package="PostChicago")
    #' designDir <- file.path(extdata,'designDir')
    #' baitmap2baited_genes(designDir,save=FALSE)
    #' @return Returns a GenomicRanges object with the positions of the baits, annotated with fragment IDs (re_id) and gene names
    #' @export

    f <- paste0(designDir, "/baited_genes.rds")

    ## only executed if baited_genes.rds doesn't exist yet, otherwise loads the existing file:
    if (f %in% list.files(designDir)) {
        message("file exists, loading...")
        bg <- readRDS(f)
    } else {
        # 22.Feb: SK replaced but not sure if it works
        baitmap <- list.files(designDir, full.names = TRUE, pattern = "baitmap$")
        # baitmap = paste(designDir, list.files(designDir)[grep("baitmap$", list.files(designDir))], sep = "/")
        baitmap <- read.delim(baitmap, header = FALSE, stringsAsFactors = FALSE)
        ## AF: 05.02.2026: updated from bed2gr:
        # bg = bed2gr(baitmap)
        bg <- makeGRangesFromDataFrame(
            baitmap,
            seqnames.field = "V1",
            start.field = "V2",
            end.field = "V3",
            ignore.strand = TRUE
        )
        bg$re_id <- baitmap$V4
        bg$genename <- baitmap$V5
        if (save) {
            saveRDS(bg, f)
        }
    }
    return(bg)
}


rmap2rmapgr <- function(designDir = "designDir", save = FALSE) {
    #' @title Convert .rmap to the GenomicRanges object 'rmapgr'
    #' @description Searches for .rmap file in designDir and converts to a GenomicRanges object rmapgr.
    #' Only executed if the GenomicRanges object doesn't exist yet, otherwise loads the existing file.
    #' @param designDir Default: 'designDir'; Path to the directory with Chicago design files. 
    #' @param save Default: FALSE; If TRUE, saves the converted object as .rds in the designDir,
    #' using the same name as the .rmap file.
    #' @examples
    #' # Example running rmap2rmapgr
    #' extdata <- system.file("extdata", package="PostChicago")
    #' designDir <- file.path(extdata,'designDir')
    #' rmap2rmapgr(designDir,save=FALSE)
    #' @return Returns a GenomicRanges object derived from .rmap.
    #' @export

    # 22.Feb: SK Replaced grep
    # f <- paste0(designDir, "/", list.files(designDir)[grep("rmap$", list.files(designDir))], ".rds")
    f <- list.files(designDir, full.names = TRUE, pattern = "rmap$")
    ## only executed if rmapgr doesn't exist yet, otherwise loads the existing file:
    ## if (paste0(list.files(designDir)[grep("rmap$", list.files(designDir))], ".rds") %in% list.files(designDir)){
    if (file.exists(paste0(f, ".rds"))) {
        message("file exists, loading...")
        rmapgr <- readRDS(paste0(f, ".rds"))
    } else {
        ## AF 19.03.2026 replaced:
        # rmap = paste(designDir, list.files(designDir)[grep("rmap$", list.files(designDir))], sep = "/")
        rmap <- read.delim(f, header = FALSE, stringsAsFactors = FALSE)
        ## AF: 05.02.2026: updated to have an S4 function:
        # rmapgr = bed2gr(rmap)
        rmapgr <- makeGRangesFromDataFrame(
            rmap,
            seqnames.field = "V1",
            start.field = "V2",
            end.field = "V3",
            ignore.strand = TRUE
        )
        rmapgr$id <- rmap$V4
        if (save) {
            saveRDS(rmapgr, paste0(f, ".rds"))
        }
    }
    return(rmapgr)
}


getbaitplus <- function(id, rmapgr) {
    #' title Get bait +/- 1
    #' description returns a gr object with the positons of the bait +/- 1 fragment
    #' param id baitID
    #' param rmapgr GenomicRanges object containing restriction fragments. IDs must be saved in a column called 'id'
    #' return GRanges Object of bait +/- 1 fragment
    #' @noRd

    min <- as.numeric(id) - 1
    max <- as.numeric(id) + 1
    bait <- GenomicRanges::reduce(rmapgr[rmapgr$id %in% min:max])
    return(bait)
}


getsk <- function(L) {
    #' title Get Scaling Factor
    #' description returns a vector with scaling factors
    #' for all cd(at)x objects in list of cd(at)x objects L
    #' scaling factors are then used to multiply N reads in cd(at)x
    #' objects with them to obtain %reads per sample per promoter (PRSP) * 1000
    #' param L List of all Chicago output file in (cd(at)x)
    #' returns scaling factor
    #' @noRd

    cov <- vector()
    nprom <- vector()
    for (i in 1:length(L)) {
        a <- L[[i]]
        b <- sum(a$N)
        cov <- c(cov, b)
        nprom <- c(nprom, length(unique(a$baitID)))
    }

    names(cov) <- names(L)
    sk <- 1 / cov * nprom * 100 * 1000
    names(sk) <- names(L)

    return(sk)
}


getnormscores <- function(
      id,
      data,
      samplename,
      L,
      same.chromosome = FALSE,
      bait,
      p = 0,
      rmapgr = rmapgr
) {
    ## returns a gr object with chicago scores for a specific bait arguments: id = baitID data=cd@x for a
    ## specific sample samplename = name of the sample, has to be the same as the name of the scaling factors
    ## sk L=list of cd@x objects for which to calculate sk same.chromosome: should only reads from the same
    ## chromosome be taken? Default:FALSE p:pseudocount to add before normalization, Default:0

    a <- data[data$baitID == id, ]
    a <- a[order(a$otherEndID), ]
    bed <- rmapgr[rmapgr$id %in% a$otherEndID]
    S4Vectors::values(bed) <- cbind(
        S4Vectors::values(bed),
        data.frame(N = as.data.frame(a)[, grep("score", names(a))])
    )
    summary(bed$id == a$otherEndID)
    bait <- resize(rmapgr[rmapgr$id %in% a$baitID], 1, fix = "center")
    if (same.chromosome == TRUE) {
        bed <- bed[GenomicRanges::seqnames(bed) == GenomicRanges::seqnames(bait)]
    }
    bed$intID <- paste(bait$id, bed$id, sep = ";")

    return(bed)
}


findOverlapsDf <- function(df1, df2, col1, col2) {
    #' title find Dataframe Overlaps
    #' description Function for mapping entire DFs or GR to each other based on IDs in cols:
    #' returns overlap matrix for df1 and df2 indices in column 1 and column 2, resp
    #' param df1 df or gr to be reordered
    #' param df2 second df or gr to be reordered
    #' param col1 character vector of length 1 or 3, based on which NUMERIC column(s) should the DFs be reordered?
    #' If one column is given, an artificial gr is created from that column,
    #' if three columns are given, automatic assumption that these are seqnames, start, end and these are used to create a gr
    #' param col2 character vector of length 1 or 3, based on which NUMERIC column(s) should the GRs be reordered?
    #' If one column is given, an artificial gr is created from that column,
    #' if three columns are given, automatic assumption that these are seqnames, start, end and these are used to create a gr
    #' returns Overlaps of Df1 and Df2
    #' @noRd

    if (length(col1) != 3) {
        ## get gr1 for overlaps:
        if (
            "data.frame" %in% class(df1) & data.table::is.data.table(df1) == FALSE
        ) {
            gr1 <- GenomicRanges::GRanges(
                rep(1, dim(df1)[1]),
                IRanges::IRanges(df1[, col1], df1[, col1])
            )
        } else if (data.table::is.data.table(df1)) {
            gr1 <- GenomicRanges::GRanges(
                rep(1, dim(df1)[1]),
                IRanges::IRanges(df1[[col1]], df1[[col1]])
            )
        } else {
            gr1 <- GenomicRanges::GRanges(
                rep(1, length(df1)[1]),
                IRanges::IRanges(
                    S4Vectors::values(df1)[, col1],
                    S4Vectors::values(df1)[, col1]
                )
            )
        }

        ## get gr2 for overlaps
        if (
            "data.frame" %in% class(df2) & data.table::is.data.table(df2) == FALSE
        ) {
            gr2 <- GenomicRanges::GRanges(
                rep(1, dim(df2)[1]),
                IRanges::IRanges(df2[, col2], df2[, col2])
            )
        } else if (data.table::is.data.table(df2)) {} else {
            gr2 <- GenomicRanges::GRanges(
                rep(1, length(df2)[1]),
                IRanges::IRanges(
                    S4Vectors::values(df2)[, col2],
                    S4Vectors::values(df2)[, col2]
                )
            )
        }
    }

    if (length(col1) == 3) {
        ## get gr1 for overlaps:
        if ("data.frame" %in% class(df1)) {
            gr1 <- GenomicRanges::GRanges(
                df1[, col1[1]],
                IRanges::IRanges(df1[, col1[2]], df1[, col1[3]])
            )
        } else {
            gr1 <- GenomicRanges::GRanges(
                S4Vectors::values(df1)[, col1[1]],
                IRanges::IRanges(
                    S4Vectors::values(df1)[
                        ,
                        col1[2]
                    ],
                    S4Vectors::values(df1)[, col1[3]]
                )
            )
        }

        ## get gr2 for overlaps
        if ("data.frame" %in% class(df2)) {
            gr2 <- GenomicRanges::GRanges(
                df2[, col2[1]],
                IRanges::IRanges(df2[, col2[2]], df2[, col2[3]])
            )
        } else {
            gr2 <- GenomicRanges::GRanges(
                S4Vectors::values(df2)[, col2[1]],
                IRanges::IRanges(
                    S4Vectors::values(df2)[
                        ,
                        col2[2]
                    ],
                    S4Vectors::values(df2)[, col2[3]]
                )
            )
        }
    }

    ov <- as.matrix(GenomicRanges::findOverlaps(gr1, gr2))

    return(ov)
}


annotateIntsRepscores <- function(L, ints) {
    #' @title Annotate Interactions with Replicate Scores
    #' @description annotates ints table with replicate scores; only works if scores were determined for each of the replicates
    #' @param L List of replicate cd objects
    #' @param ints interactions table from getInts
    #' @return Table of interactions with annotations
    #' @noRd

    ids <- ints$intID
    ints <- ints[order(ints$intID), ]

    for (i in 1:length(L)) {
        L[[i]] <- L[[i]][L[[i]]$baitID %in% ints$baitID, ]
        L[[i]]$intID <- paste(L[[i]]$baitID, L[[i]]$otherEndID, sep = ";")
        L[[i]] <- L[[i]][order(L[[i]]$intID), ]
        L[[i]] <- L[[i]][L[[i]]$intID %in% ids, ]
        ints[, paste(names(L)[i], "score", sep = "_")] <- 0
        ints[
            ints$intID %in% L[[i]]$intID,
            paste(names(L)[i], "score", sep = "_")
        ] <- L[[i]]$score
    }

    return(ints)
}


# Helper for makeIntsTable()

getInts <- function(
      L,
      baited_genes,
      repscores = FALSE,
      LL = NULL,
      scorecut = 5,
      readcut = 1,
      mode = "score",
      intervals = NULL
) {
    #' @title Create interactions data table
    #' @description
    #' Normalises, sorts, downsamples and saves ChicagoData of interactions
    #' @param L List with ChicagoData (cd(at)x)
    #' @param baited_genes Genomic Ranges object or data.frame, with columns 're_id' (bait restrictions fragments),
    #' 'genenames' (names of baited genes)
    #' @param repscores logical input, should replicate scores be added, if TRUE, LL (list of replicate ChicagoData) is required,
    #' optional: Default: FALSE
    #' @param LL list of ChicagoData objects (cd(at)x) of replicates of summarized interactions, optional: Default: FALSE
    #' @param scorecut score cutoff, optional: Default: 5
    #' @param readcut reads cutoff, optional: Default: 1
    #' @param mode used for cutoff when getting interactions, by 'score', by 'read' or 'both', optional: Default: 'score'
    #' @param intervals A genomic ranges object with intervals with which the otherEnds should overlap.
    #' Useful if mode='reads', to reduce the number of interactions in ints, optional: Default: NULL
    #' @returns 'int' variable, table of interactions data, saved .png files of quality control: library size both by reps and total
    #' @noRd

    ids <- vector()

    for (i in 1:length(L)) {
        # L[[i]]=L[[i]][L[[i]]$baitID %in% baited_genes$re_id,]
        L[[i]]$intID <- paste(L[[i]]$baitID, L[[i]]$otherEndID, sep = ";")

        if (mode == "score") {
            ids <- unique(c(
                ids,
                L[[i]][
                    L[[i]]$intID %in% ids == FALSE & L[[i]]$score >= scorecut,
                ]$intID
            ))
        } else if (mode == "read") {
            if (is.null(intervals)) {
                ids <- unique(c(
                    ids,
                    L[[i]][L[[i]]$intID %in% ids == FALSE & L[[i]]$N >= readcut, ]$intID
                ))
            } else {
                ids <- unique(c(
                    ids,
                    L[[i]][
                        L[[i]]$intID %in% ids == FALSE &
                            L[[i]]$N >= readcut &
                            L[[i]]$otherEndID %in%
                                oeids,
                    ]$intID
                ))
            }
        } else {
            ids <- unique(c(
                ids,
                L[[i]][
                    L[[i]]$intID %in% ids == FALSE &
                        L[[i]]$N >= readcut &
                        L[[i]]$score >= scorecut,
                ]$intID
            ))
        }
    }

    ids <- unique(ids)

    if (length(ids) > 0) {
        ints <- data.frame(intID = ids[order(ids)])

        for (i in 1:length(L)) {
            a <- L[[i]][L[[i]]$intID %in% ids, ]
            a <- a[order(a$intID), ]
            ints[, paste(names(L)[i], "score", sep = "_")] <- 0
            ints[, paste(names(L)[i], "N", sep = "_")] <- 0

            ## get pooled read numbers and scores
            ints[
                ints$intID %in% a$intID,
                paste(names(L)[i], "score", sep = "_")
            ] <- a$score
            ints[ints$intID %in% a$intID, paste(names(L)[i], "N", sep = "_")] <- a$N

            ## get replicate read numbers
            x <- names(a)
            x <- x[grep("N", x)]
            x <- x[grep("[.]", x)]
            if (length(x) > 0) {
                n <- as.numeric(unlist(strsplit(
                    x,
                    split = "[.]"
                ))[length(unlist(strsplit(x, split = "[.]")))])
                for (j in 1:n) {
                    ints[, paste0(names(L)[i], "_N.", j)] <- 0
                    ints[
                        ints$intID %in% a$intID,
                        paste0(names(L)[i], "_N.", j)
                    ] <- as.data.frame(a)[, paste0(
                        "N.",
                        j
                    )]
                }
            }
        }

        ints$intID <- as.character(ints$intID)
        ints$baitID <- as.numeric(unlist(strsplit(ints$intID, split = ";"))[
            1:length(unlist(strsplit(ints$intID, split = ";"))) %% 2 == 1
        ])
        ints$otherEndID <- as.numeric(unlist(strsplit(ints$intID, split = ";"))[
            1:length(unlist(strsplit(ints$intID, split = ";"))) %% 2 == 0
        ])
        ## AF: 5.2.26 updated the redundant gr2bed function:
        # oe = gr2bed(rmapgr[findOverlapsDf(ints, rmapgr, "otherEndID", "id")[, 2]])
        # b = gr2bed(rmapgr[findOverlapsDf(ints, rmapgr, "baitID", "id")[, 2]])
        oe <- as.data.frame(rmapgr[findOverlapsDf(
            ints,
            rmapgr,
            "otherEndID",
            "id"
        )[, 2]])[, 1:3]
        b <- as.data.frame(rmapgr[findOverlapsDf(ints, rmapgr, "baitID", "id")[
            ,
            2
        ]])[, 1:3]

        names(oe) <- paste(names(oe), "otherEnd", sep = "_")
        names(b) <- paste(names(b), "bait", sep = "_")
        names(oe)[1] <- "seqnames_otherEnd"
        names(b)[1] <- "seqnames_bait"
        ints <- cbind(ints, oe, b)
        row.names(ints) <- ints$intID

        ints <- ints[order(ints$otherEndID), ]
        ints <- ints[order(ints$baitID), ]

        ## normalize read counts of replicates by the read numbers in the respective libraries (written to
        ## vector v):
        m <- ints[, grep("[.]", names(ints))]

        v <- vector()
        for (i in 1:length(L)) {
            x <- names(L[[i]])
            x <- x[grep("N", x)]
            x <- x[grep("[.]", x)]
            if (length(x) > 0) {
                n <- as.numeric(unlist(strsplit(
                    x,
                    split = "[.]"
                ))[length(unlist(strsplit(x, split = "[.]")))])
                for (j in 1:n) {
                    v <- c(v, sum(as.data.frame(L[[i]])[, paste0("N.", j)]))
                    names(v)[length(v)] <- paste0(names(L)[i], "_N.", j)
                }
            }
        }
        libsize <- v

        ## downsample to smallest:
        if (
            length(summary(names(m) == names(libsize))[
                names(summary(names(m) == names(libsize))) == FALSE
            ]) ==
                0
        ) {
            mat <- m
            for (i in 1:ncol(m)) {
                m[, i] <- m[, i] / libsize[i] * min(libsize)
            }

            names(m) <- paste(names(m), "downsampled", sep = "_")
            ints <- cbind(ints, m)
        }

        ## normalize total read counts by the read numbers in the respective libraries (written to vecor v):

        m <- ints[, grep("_N$", names(ints))]

        v <- vector()
        for (i in 1:length(L)) {
            v <- c(v, sum(as.data.frame(L[[i]])[, "N"]))
            names(v)[length(v)] <- paste0(names(L)[i], "_N")
        }
        libsize <- v

        ## TO DO: Allow for separate normalization of grouped samples (samplegroups=NULL), such as Med13fl and
        ## Pcgf2fl samples norm to each othe ronly downsample to smallest:
        names(m) == names(libsize)
        mat <- m
        for (i in 1:ncol(m)) {
            m[, i] <- m[, i] / libsize[i] * min(libsize)
        }

        names(m) <- paste(names(m), "downsampled", sep = "_")
        ints <- cbind(ints, m)

        ## annotate ints with genenames:
        ints$genename <- "NA"
        i <- 0
        for (b in unique(ints$baitID)) {
            ## check that baited_genes is valid:
            nonames <- c(
                "seqnames",
                "ranges",
                "strand",
                "seqlevels",
                "seqlengths",
                "isCircular",
                "start",
                "end",
                "width",
                "element"
            )
            if (TRUE %in% (names(S4Vectors::values(baited_genes)) %in% nonames)) {
                names(S4Vectors::values(baited_genes)) -> v
                x <- as.data.frame(S4Vectors::values(baited_genes))
                x <- x[, v %in% nonames == FALSE]
                S4Vectors::values(baited_genes) <- x
            }
            n <- unique(as.character(as.data.frame(baited_genes)[
                as.data.frame(baited_genes)$re_id == b,
                grep(
                    "genename",
                    names(as.data.frame(baited_genes))
                )
            ]))

            if (length(n) > 1) {
                nl <- n[1]
                for (j in 2:length(n)) {
                    nl <- paste0(nl, ", ", n[j])
                }
                n <- nl
            }
            ints[ints$baitID == b, ]$genename <- n
        }

        if (repscores == TRUE) {
            ints <- annotateIntsRepscores(LL, ints)
        }

        return(ints)
    }
}


## Helper function for makeOneGeneOnePeak()

concatenate <- function(ids, separator = ",") {
    #' @title Concatenates strings
    #' @description Concatenates strings with a separator
    #' @noRd
    ids2 <- ids[1]
    if (length(ids) > 1) {
        for (i in 2:length(ids)) {
            ids2 <- paste(ids2, ids[i], sep = separator)
        }
    }
    return(ids2)
}


#############################################
#############################################
## Quantification in peaks:
#############################################
#############################################

## New: AF 29.4.2026

##################################################
## 2) Direct Quantification under interaction peaks
##################################################

## Helper functions:

## only obtain list of lists with all data tables one list for each interaction
## with these tables,w e will then be able to easily extract all the necessary data.
## May be faster!
obtainLists <- function(data, intlist, downsampled = FALSE) {
    lists <- lapply(1:nrow(data), FUN = function(x) {
        lapply(intlist, FUN = function(y) {
            (y[y$baitID == data[x, ]$baitID & y$otherEndID %in%
                min(as.numeric(unlist(strsplit(as.character(data[x, ]$otherEndID), split = ",")))):
                max(as.numeric(unlist(strsplit(as.character(data[x, ]$otherEndID), split = ",")))), ])
        })
    })

    names(lists) <- data$intID

    return(lists)
}

## downsample intlist
downsampleInts <- function(L) {
    nreads <- sapply(L, FUN = function(x) sum(x$N))
    downsamplingFacs <- min(nreads) / nreads * 0.9999
    Ldownsampled <- L
    for (i in 1:length(Ldownsampled)) {
        Ldownsampled[[i]]$N <- Ldownsampled[[i]]$N * downsamplingFacs[i]
    }
    return(Ldownsampled)
}


annotateIntsWithCounts <- function(data, intlist, downsampled = FALSE) {
    if (downsampled) {
        intlist <- downsampleInts(intlist)
    }

    lists <- obtainLists(data, intlist)

    ## total counts:
    all.counts <- do.call("rbind", lapply(lists, FUN = function(x) sapply(x, FUN = function(y) sum(y$N))))
    row.names(all.counts) <- names(lists)
    colnames(all.counts) <- paste0(colnames(all.counts), "_N")

    ## Only extracts scores if downsampled=FALSE
    if (downsampled) {
        colnames(all.counts) <- paste0(colnames(all.counts), "_downsampled")
        data <- cbind(data, all.counts)
    } else {
        ## scores:
        all.scores <- do.call("rbind", lapply(lists, FUN = function(x) sapply(x, FUN = function(y) sum(y$score))))
        row.names(all.scores) <- names(lists)
        colnames(all.scores) <- paste0(colnames(all.scores), "_score")

        ## combined:
        data <- cbind(data, all.counts, all.scores)
    }

    return(data)
}


saveBed <- function(bed, bedfile) {
    #' @title Save Bed File
    #' @description Saves a table without column names.
    #' @param bed Table in a bed format or a regular data.frame that should be saved without headers.
    #' @param bedfile Name of the file to which to save the table.
    #' @examples
    #' # basic usage of saveBed:
    #' bed <- data.frame(chr='chr1',start=1000, end=2000)
    #' saveBed(bed, 'myfile.bed')
    #' @return Saves a table without column names.
    #' @export

    write.table(
        bed,
        bedfile,
        sep = "\t",
        eol = "\n",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )
}


####################################################################################################
## Helpers for getMatrix reads
####################################################################################################

findOthers <- function(intID, zoom) {
    oe <- as.numeric(unlist(strsplit(as.character(intID), split = ";"))[2])
    minoe <- oe - zoom
    maxoe <- oe + zoom
    return(minoe:maxoe)
}

extractValues <- function(downsampledReads, bait, others) {
    df <- data.frame(others = others, Ndownsampled = 0)
    reads <- downsampledReads[
        downsampledReads$baitID == bait & downsampledReads$otherEndID %in% others,
        c("otherEndID", "N")
    ][
        order(
            downsampledReads[
                downsampledReads$baitID == bait &
                    downsampledReads$otherEndID %in% others,
                c("otherEndID", "N")
            ]$otherEndID
        ),
    ]
    df[df$others %in% reads$otherEndID, ]$Ndownsampled <- reads$N
    return(df$Ndownsampled)
}

annotateSampleMatrix <- function(intID, downsampledReads, zoom) {
    bait <- as.numeric(unlist(strsplit(as.character(intID), split = ";"))[1])
    others <- findOthers(intID, zoom)
    return(extractValues(downsampledReads, bait, others))
}

annotateControlMatrix <- function(intID, downsampledReads, zoom) {
    bait <- as.numeric(unlist(strsplit(as.character(intID), split = ";"))[1])
    otherEnd <- as.numeric(unlist(strsplit(as.character(intID), split = ";"))[2])
    ## remake intID to be distance-matched on the other side of the bait:
    intID <- paste(bait, (bait + bait - otherEnd), sep = ";")
    others <- findOthers(intID, zoom) ## no need to check chromosome, values will be 0 if different chromosome!
    return(extractValues(downsampledReads, bait, others))
}

####################################################################################################
## Helpers for getMatrix scores
####################################################################################################

extractValuesScores <- function(downsampledReads, bait, others) {
    df <- data.frame(others = others, Ndownsampled = 0)
    reads <- downsampledReads[
        downsampledReads$baitID == bait & downsampledReads$otherEndID %in% others,
        c("otherEndID", "score")
    ][
        order(
            downsampledReads[
                downsampledReads$baitID == bait &
                    downsampledReads$otherEndID %in% others,
                c("otherEndID", "score")
            ]$otherEndID
        ),
    ]
    df[df$others %in% reads$otherEndID, ]$Ndownsampled <- reads$score
    return(df$Ndownsampled)
}

annotateSampleMatrixScores <- function(intID, downsampledReads, zoom) {
    bait <- as.numeric(unlist(strsplit(as.character(intID), split = ";"))[1])
    others <- findOthers(intID, zoom)
    return(extractValuesScores(downsampledReads, bait, others))
}

annotateControlMatrixScores <- function(intID, downsampledReads, zoom) {
    bait <- as.numeric(unlist(strsplit(as.character(intID), split = ";"))[1])
    otherEnd <- as.numeric(unlist(strsplit(as.character(intID), split = ";"))[2])
    ## remake intID to be distance-matched on the other side of the bait:
    intID <- paste(bait, (bait + bait - otherEnd), sep = ";")
    others <- findOthers(intID, zoom) ## no need to check chromosome, values will be 0 if different chromosome!
    return(extractValuesScores(downsampledReads, bait, others))
}
