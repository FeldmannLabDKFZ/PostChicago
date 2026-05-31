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


## annotate interactions with peaks from your grl:

annotateInts <- function(ints, grl, minov = NULL, maxgp = NULL) {
    #' @title Annotate Interactions
    #' @description Annotates interaction tables containing CaptureC data with features (typically ChIP-Seq peaks)
    #' from a GRangesList object (grl). Annotation is based on feature overlap.
    #' @param ints a table with interactions derived from ChicagoData objects by makeIntsTable().
    #' Assumes six named columns:
    #' 'seqnames_bait','start_bait','end_bait' (coordinates of the baited restriction fragment) and
    #' 'seqnames_otherEnd','start_otherEnd','end_otherEnd' (coordinates of the interacting restriction fragment).
    #' @param grl GRangesList containing features with which ints should be annotated.
    #' Names of the grl correspond to the names that will be given to table columns.
    #' @param minov Default:NULL, minimum overlap (minoverlap in findOverlaps) between otherEnd or bait fragment
    #' and the feature.
    #' @param maxgp Default:NULL, maximum gap (maxgap in findOverlaps) between otherEnd or bait fragment and the feature.
    #' @examples
    #' # basic usage of annotateInts
    #' ints=data.frame(seqnames_bait='chr1', start_bait=1, end_bait=1000,
    #'                 seqnames_otherEnd='chr1', start_otherEnd=10000, end_otherEnd=20000)
    #' gr1=GenomicRanges::GRanges('chr1',IRanges::IRanges(10500,11000))
    #' gr2=GenomicRanges::GRanges('chr1',IRanges::IRanges(900,1100))
    #' grl=GenomicRanges::GRangesList(gr1,gr2)
    #' names(grl)=c('mytestpeaks1','mytestpeaks2')
    #' grl
    #' annotateInts(ints,grl)
    #' @return interactions file with annotations from GRanges file
    #' @export

    ## AF: 05.02.2026: Updated to an S4 function:
    # b = bed2gr(ints[, c("seqnames_bait", "start_bait", "end_bait")])
    b <- makeGRangesFromDataFrame(
        ints,
        seqnames.field = "seqnames_bait",
        start.field = "start_bait",
        end.field = "end_bait",
        ignore.strand = TRUE
    )

    for (i in 1:length(grl)) {
        iv <- grl[[i]]
        n <- names(grl)[i]
        # initialise variable values before
        S4Vectors::values(b)[, n] <- FALSE
        if (!is.null(minov)) {
            S4Vectors::values(b[as.matrix(GenomicRanges::findOverlaps(
                b,
                iv,
                minoverlap = minov
            ))[, 1]])[, n] <- TRUE
        } else if (!is.null(maxgp)) {
            S4Vectors::values(b[as.matrix(GenomicRanges::findOverlaps(
                b,
                iv,
                maxgap = maxgp
            ))[, 1]])[, n] <- TRUE
        } else {
            S4Vectors::values(b[as.matrix(GenomicRanges::findOverlaps(b, iv))[, 1]])[
                ,
                n
            ] <- TRUE
        }
    }

    names(S4Vectors::values(b)) <- paste0(names(S4Vectors::values(b)), "_bait")
    ## ensure there is no duplication of the same column names. If ints already contains info on interval
    ## overlaps, they will now be replaced:
    ints <- ints[, names(ints) %in% names(S4Vectors::values(b)) == FALSE]
    ints <- cbind(ints, as.data.frame(S4Vectors::values(b)))

    ## AF: 05.02.2026: Updated to an S4 function:
    # b = bed2gr(ints[, c("seqnames_otherEnd", "start_otherEnd", "end_otherEnd")])
    b <- makeGRangesFromDataFrame(
        ints,
        seqnames.field = "seqnames_otherEnd",
        start.field = "start_otherEnd",
        end.field = "end_otherEnd",
        ignore.strand = TRUE
    )

    for (i in 1:length(grl)) {
        iv <- grl[[i]]
        n <- names(grl)[i]
        S4Vectors::values(b)[, n] <- FALSE
        if (!is.null(minov)) {
            S4Vectors::values(b[as.matrix(GenomicRanges::findOverlaps(
                b,
                iv,
                minoverlap = minov
            ))[, 1]])[, n] <- TRUE
        } else if (!is.null(maxgp)) {
            S4Vectors::values(b[as.matrix(GenomicRanges::findOverlaps(
                b,
                iv,
                maxgap = maxgp
            ))[, 1]])[, n] <- TRUE
        } else {
            S4Vectors::values(b[as.matrix(GenomicRanges::findOverlaps(b, iv))[, 1]])[
                ,
                n
            ] <- TRUE
        }
    }

    names(S4Vectors::values(b)) <- paste0(
        names(S4Vectors::values(b)),
        "_otherEnd"
    )
    ## ensure there is no duplication of the same column names. If ints already contains info on interval
    ## overlaps, they will now be replaced:
    ints <- ints[, names(ints) %in% names(S4Vectors::values(b)) == FALSE]
    ints <- cbind(ints, as.data.frame(S4Vectors::values(b)))

    return(ints)
}


## AF 20.3.2026: Removed heatmaps to add to makeQC plots
makeIntsTable <- function(
      L,
      baited_genes,
      repscores = FALSE,
      LL = NULL,
      ngroups = 1,
      scorecut = 5,
      readcut = 1,
      mode = "score",
      intervals = NULL,
      outfolder = NULL,
      overwrite = FALSE,
      saveFile = TRUE
) {
    #' @title Create Interactions Table
    #' @description
    #' Integrates all significant interactions from multiple datasets, as supplied by the list of ChicagoData tables in L.
    #' The table is populated with scores, raw and downsampled reads from pooled and (optionally) replicate data.
    #' @param L list of ChicagoData objects containing summarized interactions.
    #' @param baited_genes GenomicRanges object with at least two metadata columns:
    #' 1) 're_id' (captured fragment ID (bait), corresponding to fragment numbers indicated in rmapgr),
    #' 2) genenames' (names of the corresponding baits, typically gene names).
    #' Can be created automatically from the Chicago baitmap table (see R Bioconductor package 'Chicago')
    #' using the function baitmap2baited_genes().
    #' @param repscores TRUE/FALSE (default): Should replicate scores be added? if TRUE, LL must be supplied,
    #' Default: FALSE
    #' @param LL optional: List of ChicagoData tables containing individual replicates of summarized interactions.
    #' Names are expected to follow the naming convention '_rep1', '_rep2',..., '_repn'
    #' @param ngroups Default: 1; number of groups in which unique view points (baits) should be processed.
    #' Efficient computing is achieved for a maximum of ca. 100-200 entries per group.
    #' @param scorecut Default: 5 (as recommended by Chicago); Chicago score cutoff for significant interactions.
    #' Sometimes useful to lower to 3.
    #' @param readcut Default: 1; Reads cutoff for significant interactions.
    #' @param mode Options: 'score', 'read' or 'both'; Default: 'score'; Which parameter should be used for the
    #' extraction of significant interactions?
    #' @param intervals optional: A GenomicRanges object with intervals with which the otherEnds should overlap.
    #' Useful if mode='reads', to reduce the number of interactions in ints.
    #' @param outfolder Default: current directory; Path to the directory in which the output is saved.
    #' @param overwrite Default: FALSE. If TRUE runs the function even if an ints_all.txt is present in the outputdir,
    #' otherwise loads ints_all.txt.
    #' @param saveFile Default: TRUE; saves ints_all.txt in outfolder.
    #' @examples
    #' # example code
    #'
    #' ##Create an rmapgr object
    #' rmapgr <- GenomicRanges::GRanges(rep(1,5), IRanges::IRanges( 1:5,2:6))
    #' rmapgr$id=1:5
    #'
    #' ##create baited_genes with one gene
    #' baited_genes <- GenomicRanges::GRanges(1,IRanges::IRanges(1,2))
    #' baited_genes$re_id=1
    #' baited_genes$genename='testgene'
    #'
    #' ##create a list of ChicagoData objects with two samples
    #' cdlist <- list(
    #' sample1 = data.frame(baitID=rep(1,4),
    #'                       otherEndID=c(2,3:5),
    #'                       N=c(2,5,10,8),
    #'                       N.1=c(2,5,10,8),
    #'                       N.2=c(3,4,9,7),
    #'                       score=c(5, 1, 2, 4)),
    #'
    #' sample2 = data.frame(baitID=rep(1,4),
    #'                       otherEndID=c(2,3:5),
    #'                       N=c(20,3,1,6),
    #'                       N.1=c(20,3,1,6),
    #'                       N.2=c(17,4,2,3),
    #'                       score=c(3, 5, 3, 5))
    #'       )
    #'
    #' ##run makeIntsTable
    #' ints = makeIntsTable(cdlist,baited_genes,saveFile=FALSE)
    #' head(ints)
    #'
    #' @returns Table of all significant interactions.
    #' @export

    ## checking table structures
    if (length(grep("genename", names(S4Vectors::values(baited_genes)))) < 1) {
        stop("baited genes lacks genenames column")
    }
    if (length(grep("re_id", names(S4Vectors::values(baited_genes)))) < 1) {
        stop("baited genes lacks re_id column")
    }
    if (length("id" %in% names(S4Vectors::values(baited_genes))) < 1) {
        stop("rmapgr lacks id column")
    }

    # 8Feb26: SK: changed output management & syntax
    if (is.null(outfolder) && saveFile) {
        outfolder <- getwd()
        message("Output saved here: ", outfolder)
    }

    ## Checking if ints_all.txt is already in the outfolder
    file <- "ints_all.txt"

    ## AF: 09.02.26: Added naming convention for LL after talking to Ezo
    ## Checking if LL has the correct naming convention containing 'rep':
    if (!is.null(LL)) {
        if (length(grep("rep", names(LL))) == 0) {
            stop(
                "Names of the replicate interactions list LL should be in the format .._rep1,.._rep2, etc."
            )
        }
    }

    if (!(file.exists(paste0(outfolder, "/", file))) || overwrite == TRUE) {
        # why is there a nested function??
        f <- function(
      y,
      scorecut = scorecut,
      readcut = readcut,
      mode = "score",
      intervals = NULL
        ) {
            filtered_L <- lapply(L, function(df) {
                # Check if baitID is in any of the arrays in y
                match_found <- sapply(y, function(arr) any(df$baitID %in% arr))

                # If any match is found in y, filter the rows
                if (any(match_found)) {
                    return(df[df$baitID %in% unlist(y[match_found]), ])
                } else {
                    return(NULL) # If no match found, return NULL (or empty data frame)
                }
            })

            # Clean up the list: remove NULL entries (if any)
            filtered_L <- filtered_L[!sapply(filtered_L, is.null)]
            return(getInts(
                filtered_L,
                baited_genes[baited_genes$re_id %in% unlist(y)],
                repscores,
                LL,
                scorecut = scorecut,
                readcut = readcut,
                mode = mode,
                intervals = intervals
            ))
        }

        ## split baitIDs into groups of ~100-200:
        v <- vector()
        for (i in 1:length(L)) {
            v <- c(v, unique(L[[i]]$baitID))
        }
        v <- unique(v)

        # check so that the samples measured are a multiple of the ngroups
        a <- split(baited_genes$re_id[baited_genes$re_id %in% v], 1:ngroups)

        ints <- do.call(
            "rbind",
            BiocGenerics::lapply(
                a,
                FUN = f,
                scorecut = scorecut,
                readcut = readcut,
                mode = mode,
                intervals = intervals
            )
        )
        row.names(ints) <- ints$intID

        ## save the file:
        if (saveFile) {
            write.table(
                ints,
                paste0(outfolder, "/ints_all.txt"),
                sep = "\t",
                eol = "\n",
                quote = FALSE,
                row.names = FALSE
            )
        }

        # dev.off()
    } else {
        message("ints_all.txt found in outfolder, loading...")
        ints <- read.delim(
            paste(outfolder, file, sep = "/"),
            stringsAsFactors = FALSE
        )
    }

    return(ints)
}


## Newest function,
## written on 4.5.2026 (AF)

makeOneGeneOnePeak <- function(
      grl,
      rmapgr,
      ints,
      maxgap = -1,
      folder = "oneGeneOnePeak"
) {
    #' @title Make one-Gene-one-Peak Interaction Tables.
    #' @description This function takes the interaction table ints (typically containing all significant interactions)
    #' and intersects all interactions from this table with intervals provided in the GRangesList grl. Each entry in grl
    #' receives a unique ID and is assigned only one interaction with each bait with which at least one significant
    #' interaction could be detected. OneGeneOnePeak tables for each interval set within grl are saved in folder.
    #' Within this table, start_otherEnd' and 'end_otherEnd' correspond to the start of the first and end of the last
    #' restriction fragment that overlap a particular interval. 'otherEndID' corresponds to all overlapped restriction
    #' fragments from rmapgr separated by commas ('1,2,3'). 'summit_re_id' corresponds to the restriction ID of the
    #' central restriction fragment within rmapgr that overlaps an interval.
    #' If two intervals overlap exactly the same set of restriction fragments from rmapgr, this will be visible
    #' in the output table column 'interval_id_<my_interval>' as two comma-separated ids.
    #' @param grl GRangesList containing features for which oneGeneOnePeak tables should be made.
    #' Names of the grl correspond to the names that will be given to table columns.
    #' @param rmapgr GenomicRanges object containing restriction fragments from the reference genome.
    #' IDs must be saved in a column called 'id'.
    #' @param ints Interactions table containing all significant interactions and typically derived by makeIntsTable().
    #' @param maxgap Default: -1, full overlap; Maximum gap between the otherEnd and an interval to count as overlap
    #' (see findOverlaps()).
    #' @param folder Output folder. Default: 'oneGeneOnePeak'
    #' @examples
    #' # example code
    #'
    #' ##create ints table with 4 interactions
    #' ints <- data.frame(intID=c('1;2', '1;3', '1;4','2;8'),
    #'                    baitID = as.character(c(1,1,1,2)),
    #'                    otherEndID=as.character(c(2,3,4,8)),
    #'                    seqnames_bait=as.character((paste0('chr',rep(1,4)))),
    #'                    start_bait=as.character(c(1,1,1,2)),
    #'                    end_bait=as.character(c(2,2,2,3)),
    #'                    seqnames_otherEnd=as.character((paste0('chr',rep(1,4)))),
    #'                    start_otherEnd=as.character(c(3,5,7,15)),
    #'                    end_otherEnd=as.character(c(4,6,8,16))
    #'                    )
    #'
    #' ##make rmapgr
    #' rmapgr <- GenomicRanges::GRanges(rep('chr1',8),
    #'               IRanges::IRanges( c(1,3,5,7,9,11,13,15),
    #'               c(2,4,6,8,10,12,14,16)  ) )
    #'               rmapgr$id=1:length(rmapgr)
    #'
    #' ##make grl with one interval type and two intervals
    #' grl <- GenomicRanges::GRangesList(testintervals=GenomicRanges::GRanges(rep('chr1',3),
    #'                    IRanges::IRanges(c(2,5,6), c(4,7,8))),
    #'                    testintervals2=GenomicRanges::GRanges(rep('chr1',2),
    #'                    IRanges::IRanges(c(5,11), c(10,13))),
    #'                    testintervals3 = GenomicRanges::GRanges('chr1',IRanges::IRanges(1,2)))
    #'
    #' ##run makeOneGeneOnePeak where intervals have to fully overlap
    #' ogopList <- makeOneGeneOnePeak(grl = grl,rmapgr = rmapgr, ints = ints, maxgap=-1)
    #'
    #' @returns A list with oneGeneOnePeak interactions.
    #' Saves a total of four files, for each of the overlapping intervals in grl:
    #' 1) An annotated OneGeneOnePeak interaction table as .txt.
    #' 2) An interaction table containing all information about interacting intervals,
    #' added to each row of the original interaction table ints,
    #' (not oneGeneOnePeak yet!): 'ints_<myintervaltype>_interacting.txt'
    #' 3) A reference table containing all interacting intervals, their ids (interval_id) and covered re_ids.
    #' In this table, intervals are summarized to one interval if they overlap exactly the same restriction fragments:
    #' '<myintervaltype>_interacting_rmapFrag_annotated.txt'
    #' 4) A reference table containing all intervals with interval_id and summit_re_id:
    #' '<myintervaltype>_withID.txt'
    #' oneGeneOnePeak tables can next be annotated with quantitative data using quantifyWithinPeaks().
    #' @export

    if (!file.exists(folder)) {
        dir.create(folder)
    }

    ogopList <- list()

    ## each interval receives a unique interval_id
    grl. <- GenomicRanges::GRangesList()
    for (i in 1:length(grl)) {
        a <- grl[[i]]
        a$interval_id <- 1:length(a)
        grl. <- c(grl., GenomicRanges::GRangesList(a))
        names(grl.)[i] <- names(grl)[i]
    }
    grl <- grl.
    rm(grl.)


    ## iterate through grl

    for (i in 1:length(grl)) {
        message(paste("processing", names(grl)[i]), "...")

        file1 <- paste0(folder, "/ints_", names(grl)[i], "_interacting.txt")
        file2 <- paste0(
            folder,
            "/",
            names(grl)[i],
            "_interacting_rmapFrag_annotated.txt"
        )
        file3 <- paste0(
            folder,
            "/oneGeneOnePeak_",
            names(grl)[i],
            "_interacting.txt"
        )
        file4 <- paste0(
            folder,
            "/",
            names(grl)[i],
            "_withID.txt"
        )

        ## Only generates OneGeneOnePeak files if they do not exist yet
        if (!file.exists(file3)) {
            ## Immediately test if the intervals overlap at all with the interactions.
            ## If not, construct artificial empty ogop table
            ov <- as.matrix(findOverlaps(
                grl[[i]],
                makeGRangesFromDataFrame(
                    ints,
                    seqnames.field = "seqnames_otherEnd",
                    start.field = "start_otherEnd",
                    end.field = "end_otherEnd",
                    ignore.strand = TRUE
                )
            ))

            if (nrow(ov) > 0) {
                ## for each interval: determine all restriction fragments it overlaps:
                ov <- as.matrix(GenomicRanges::findOverlaps(grl[[i]], rmapgr))

                ## only proceed if ov is not empty!
                if (nrow(ov) > 0) {
                    intervalToReFrags <- grl[[i]][ov[, 1], ]
                    intervalToReFrags$re_id <- rmapgr[ov[, 2]]$id
                    intervalToReFrags <- as.data.frame(intervalToReFrags)

                    ## aggregate the table by unique interval_id
                    reFragsPerInterval <-
                        aggregate(intervalToReFrags$re_id, by = list(intervalToReFrags$interval_id), FUN = concatenate)

                    ## resolve cases where two intervals overlap exactly the same restriction frags, they are then
                    ## assigned to be one single interval
                    intervalPerReFrags <-
                        aggregate(reFragsPerInterval$Group.1, by = list(reFragsPerInterval$x), FUN = concatenate)
                    names(intervalPerReFrags) <- c("re_id", "interval_id")
                    intervalPerReFrags <- intervalPerReFrags[order(intervalPerReFrags$interval_id), ]

                    ## If all IDs remain unique in their restriction frags overlap, we can proceed, if two IDs overlap,
                    ## then grl[[i]] needs adjustment -> we will merge the overlapping intervals and their ID will be the mixed ID
                    if ("FALSE" %in% names(summary(intervalPerReFrags$interval_id %in% grl[[i]]$interval_id))) {
                        tomerge <- intervalPerReFrags$interval_id[!intervalPerReFrags$interval_id %in% grl[[i]]$interval_id]
                        uniq <- grl[[i]][grl[[i]]$interval_id %in% intervalPerReFrags$interval_id]
                        for (x in tomerge) {
                            gr <- grl[[i]][grl[[i]]$interval_id %in% unlist(strsplit(x, split = ","))]
                            grr <- GenomicRanges::GRanges(unique(seqnames(gr)), IRanges::IRanges(min(start(gr)), max(end(gr))))
                            grr$interval_id <- x
                            uniq <- c(uniq, grr)
                        }
                        intervals <- uniq
                    } else {
                        intervals <- grl[[i]]
                    }

                    ## Now intervals, adjusted for those that overlap exactly the same re_ids, are assigned re_ids from intervalPerReFrags
                    intervals <- intervals[order(intervals$interval_id)]
                    intervalPerReFrags <- intervalPerReFrags[order(intervalPerReFrags$interval_id), ]

                    message("Interval overlaps resolved")

                    ## test if interval_id are in the same order as these ids in grl[[i]]
                    ## if so, then remove those that are different and issue a warning, retest,
                    ## if still different, stop teh function
                    if ("FALSE" %in% names(summary(intervals$interval_id == intervalPerReFrags$interval_id))) {
                        intervals <- intervals[intervals$interval_id %in% intervalPerReFrags$interval_id]
                        if ("FALSE" %in% names(summary(intervals$interval_id == intervalPerReFrags$interval_id))) {
                            stop("Intervals and reFrags overlapping intervals do not match...")
                        }
                    }

                    ## Annotate each interval with all restriction fragments that overlap with it
                    intervals$re_id <- intervalPerReFrags$re_id

                    ## for each interaction determine the name of the k27ac peak the gene interacts
                    # b <- grl[[i]]
                    # b <- b[as.matrix(findOverlaps(b, y))[, 1], ]

                    ## AF: 05.02.2026: Updated to an S4 function:
                    # ix = bed2gr(ints[, c("seqnames_otherEnd", "start_otherEnd", "end_otherEnd")])
                    ix <- makeGRangesFromDataFrame(
                        ints,
                        seqnames.field = "seqnames_otherEnd",
                        start.field = "start_otherEnd",
                        end.field = "end_otherEnd",
                        ignore.strand = TRUE
                    )

                    ov <- as.matrix(GenomicRanges::findOverlaps(ix, intervals, maxgap = maxgap))
                    ints.k <- ints[ov[, 1], ]
                    ints.k[, paste0("re_id_", names(grl)[i])] <- intervals[ov[, 2], ]$re_id
                    ints.k[, paste0("interval_id_", names(grl)[i])] <- intervals[ov[, 2], ]$interval_id
                    # ints.k = cbind(ints.k, gr2bed(grl[[i]][ov[, 2]]))  ##AF (20.11.2025): Commented out. grl[[i]] may not be the right one anymore here, the correct one is b

                    ## AF: 5.2.2026 Updated the redundant gr2bed function
                    # ints.k = cbind(ints.k, gr2bed(b[ov[, 2]]))
                    ints.k <- cbind(ints.k, as.data.frame(intervals[ov[, 2]])[, 1:3])
                    names(ints.k)[(ncol(ints.k) - 2):ncol(ints.k)] <- paste(
                        names(ints.k)[(ncol(ints.k) - 2):ncol(ints.k)],
                        names(grl)[i],
                        sep = "_"
                    )

                    message("Interactions overlapped with intervals")


                    write.table(
                        ints.k,
                        file1,
                        sep = "\t",
                        eol = "\n",
                        quote = FALSE,
                        row.names = FALSE
                    )
                    write.table(
                        as.data.frame(intervals),
                        file2,
                        sep = "\t",
                        eol = "\n",
                        quote = FALSE,
                        row.names = FALSE
                    )

                    ## make a table with one baitID one peak interactions.
                    ## Scores and reads won't be quantified properly yet!
                    ints.k$ogopID <- paste(ints.k$baitID, ints.k[, paste0("re_id_", names(grl)[i])], sep = ";")
                    ogop <- data.frame("intID" = unique(ints.k$ogopID))
                    ogop$baitID <- sapply(as.character(ogop$intID), FUN = function(x) unlist(strsplit(x, split = ";"))[1])
                    ogop$otherEndID <- sapply(as.character(ogop$intID), FUN = function(x) unlist(strsplit(x, split = ";"))[2])

                    baitdata <- do.call("rbind", lapply(as.character(ogop$baitID), FUN = function(x) {
                        as.data.frame(rmapgr[rmapgr$id == x])[, 1:3]
                    }))
                    names(baitdata) <- paste0(names(baitdata), "_bait")

                    oedata <- do.call("rbind", lapply(as.character(ogop$otherEndID), FUN = function(x) {
                        as.data.frame(reduce(rmapgr[rmapgr$id %in% unlist(strsplit(x, split = ","))]))[, 1:3]
                    }))
                    names(oedata) <- paste0(names(oedata), "_otherEnd")

                    intervalinfo <- do.call("rbind", lapply(ogop$intID, FUN = function(x) {
                        unique(ints.k[
                            ints.k$ogopID == x,
                            c(
                                paste0("interval_id_", names(grl)[i]),
                                paste0("seqnames_", names(grl)[i]),
                                paste0("start_", names(grl)[i]),
                                paste0("end_", names(grl)[i])
                            )
                        ])
                    }))

                    ogop <- cbind(ogop, baitdata, oedata, intervalinfo)

                    message("OneGeneOnePeak table created")


                    ## ogop is now a one-gene-one-peak file, where bait is the bait and otherEnd are all other ends that overlap
                    ## the entire interval. ogopID sumarizes this. interval_id gives the id of the intervals included

                    ## In the next step, we will annotate ogop with unique interval peak summits
                    ## that overlap the exact center of the interval fragments.
                    ## Note that this may result in ogop having two different summits.
                    ## We need to make sure to take care of it in the next step of data analysis
                    ## For instance for
                    ## a) the matrix, where we will have to reassign start_otherEnd and end_otherEnd
                    ## to refer to the summit
                    ## b) also we will have to adjust it for the quantification from the matrix

                    ## assign central re-frag to each interval
                    interval.summits <- resize(grl[[i]], 1, fix = "center")
                    ov <- as.matrix(findOverlaps(interval.summits, rmapgr))
                    intervals.withSummits <- grl[[i]][ov[, 1]]
                    intervals.withSummits$summit_re_id <- rmapgr[ov[, 2]]$id

                    write.table(
                        as.data.frame(intervals.withSummits),
                        file4,
                        sep = "\t",
                        eol = "\n",
                        quote = FALSE,
                        row.names = FALSE
                    )

                    ## annotate ogop with these central re_ids
                    ogop$summit_re_id <- unlist(lapply(as.character(ogop[, paste0("interval_id_", names(grl)[i])]), FUN = function(x) {
                        concatenate(unlist(intervals.withSummits[intervals.withSummits$interval_id %in% unlist(strsplit(x, split = ","))]$summit_re_id))
                    }))

                    message("Summits annotated")

                    ## save one gene one peak file:
                    ## it now contains the following columns:
                    "intID"
                    "baitID"
                    "otherEndID"
                    "seqnames_bait"
                    "start_bait"
                    "end_bait"
                    "seqnames_otherEnd"
                    "start_otherEnd"
                    "end_otherEnd"
                    "interval_id_testintervals"
                    "seqnames_testintervals"
                    "start_testintervals"
                    "end_testintervals"
                    "summit_re_id"
                    ## intID summarizes all otherEndIDs covering the interval peak
                    ## start_otherEnd and end_otherEnd are start of the first and end of the last fragment
                    ## This table can now be annotated using the new annotation function (quantifyWithinPeaks)
                    ## or the matrix can be created from it using getMatrix()

                    write.table(ogop, file3, sep = "\t", eol = "\n", quote = FALSE)

                    ## add ogop to the list:
                    d <- ogop
                }
            } else {
                message("Note: no overlap with ints!")

                ## create an empty table!
                d <- as.data.frame(matrix(1:14, nrow = 1))

                names(d) <- c(
                    "intID",
                    "baitID",
                    "otherEndID",
                    "seqnames_bait",
                    "start_bait",
                    "end_bait",
                    "seqnames_otherEnd",
                    "start_otherEnd",
                    "end_otherEnd",
                    "interval_id_testintervals",
                    "seqnames_testintervals",
                    "start_testintervals",
                    "end_testintervals",
                    "summit_re_id"
                )

                d <- d[-1, ]
            }
        } else {
            message(paste(
                "OneGeneOnePeak file exists for",
                names(grl)[i],
                ", skipping..."
            ))
            d <- read.delim(file3, stringsAsFactors = FALSE)
        }
        ogopList <- c(ogopList, list(d))
    }
    names(ogopList) <- names(grl)
    return(ogopList)
}


#############################################
#############################################
## Quantification in peaks:
#############################################
#############################################


## Main function:

quantifyWithinPeaks <- function(data, cdlist, cdlist.reps = NULL, rmapgr, minsize = NULL, maxsize = NULL, file = NULL) {
    #' @title Quantification of Reads and Scores within oneGeneOnePeak or AggregatedPeak Interactions.
    #' @description Annotates a table containing interactions between baits and distal sites where some interactions
    #' overlap with multiple restriction fragments, such as oneGenOnePeak or aggregatePeak tables.
    #' The function summarizes all read and score information for the interaction,
    #' including restriction fragments that are not detected as significant interactions.
    #' Output contains two types of quantifications:
    #' 1) total sum of reads (raw: _N and downsampled: _N_downsampled) and scores (_score) for each interaction, and
    #' 2) normalized quantifications (_N_downsampled_mean and _score_mean) that normalize for the total number of
    #' valid restriction fragments within an interaction peak.
    #' @param data Table with interactions derived either through oneGeneOnePeak() or aggregatePeaks().
    #' Usually contains the following columns:
    #' "intID","baitID","otherEndID","seqnames_bait","start_bait","end_bait","seqnames_otherEnd","start_otherEnd","end_otherEnd",
    #' minimal requirement: "intID","baitID","otherEndID".
    #' @param cdlist List of ChicagoData tables from which to make quantifications.
    #' @param cdlist.reps Optional: List of ChicagoData tables from which to make quantifications for individual replicates.
    #' @param rmapgr GenomicRanges object containing restriction fragments. IDs must be saved in a column called 'id'.
    #' @param minsize Optional, specifies the minimal size of valid restriction fragments and should correspond to the
    #' parameter chosen for creating the ChicagoData objects (mySettings$minFragLen) within the cdlist.
    #' @param maxsize Optional, specifies the maximal size of valid restriction fragments and should correspond to the
    #' parameter chosen for creating the ChicagoData objects (mySettings$maxFragLen) within the cdlist.
    #' @param file Optional; Path to the file in which the data should be saved.
    #' If specified, the final output will be saved under 'file', while up to four intermediate outputs will be saved under
    #' file1...4.txt. Note: Output will not be saved if file is not specified.
    #' @examples
    #' # example code
    #'
    #' ##Create an rmapgr object
    #' rmapgr <- GenomicRanges::GRanges(rep(1,5), IRanges::IRanges( 1:5,2:6  ))
    #' rmapgr$id=1:5
    #'
    #' ##create data with two interactions
    #' data <- data.frame(intID=c('1;2', '1;3,1;5'),
    #'                   baitID = as.character(c(1,1)),
    #'                   otherEndID=as.character(c(2, '3,5')))
    #'
    #' ##create a ChicagoData object with two samples
    #' cdlist <- list(
    #' sample1 = data.frame(baitID=rep(1,4),
    #'                       otherEndID=c(2,3:5),
    #'                       N=c(2,5,10,8),
    #'                       score=c(5, 1, 2, 4)),
    #'
    #' sample2 = data.frame(baitID=rep(1,4),
    #'                       otherEndID=c(2,3:5),
    #'                       N=c(20,3,1,6),
    #'                       score=c(3, 5, 3, 5))
    #'       )
    #'
    #' ##run quantifyWithinPeaks
    #' quantifyWithinPeaks(data, cdlist, rmapgr=rmapgr)
    #'
    #' @returns Interaction table annotated with the sum and mean of read and score counts for each interaction.
    #' @export

    ## total counts and scores:
    data <- annotateIntsWithCounts(data, cdlist)
    if (!is.null(file)) {
        write.table(data, paste0(file, 1, ".txt"), sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE)
    }
    gc()

    ## replicate counts and scores:
    if (!is.null(cdlist.reps)) {
        data <- annotateIntsWithCounts(data, cdlist.reps)
        if (!is.null(file)) {
            write.table(data, paste0(file, 2, ".txt"), sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE)
        }
        gc()
    }

    ## downsampled counts:
    data <- annotateIntsWithCounts(data, cdlist, downsampled = TRUE)
    if (!is.null(file)) {
        write.table(data, paste0(file, 3, ".txt"), sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE)
    }
    gc()

    ## downsampled rep counts:
    if (!is.null(cdlist.reps)) {
        data <- annotateIntsWithCounts(data, cdlist.reps, downsampled = TRUE)
        if (!is.null(file)) {
            write.table(data, paste0(file, 4, ".txt"), sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE)
        }
        gc()
    }

    ## Normalize to total fragment count (mean counts)
    ## while excluding all frags from the quantification that are below the assessed frag length:

    ## identify all validFrags in rmapgr:
    if (is.null(minsize)) {
        minsize <- 0
    }
    if (is.null(maxsize)) {
        maxsize <- max(width(rmapgr))
    }
    includedIDs <- rmapgr[width(rmapgr) >= minsize & width(rmapgr) <= maxsize]$id

    ## identify the total number of valid frags:
    data$totalValidFrags <- sapply(1:nrow(data), FUN = function(x) {
        length(as.vector(min(as.numeric(unlist(strsplit(as.character(data[x, ]$otherEndID), split = ",")))):
        max(as.numeric(unlist(strsplit(as.character(data[x, ]$otherEndID), split = ",")))))[
            as.vector(min(as.numeric(unlist(strsplit(as.character(data[x, ]$otherEndID), split = ",")))):
            max(as.numeric(unlist(strsplit(as.character(data[x, ]$otherEndID), split = ","))))) %in% includedIDs
        ])
    })

    ## normalize downsampled reads to this number:
    df <- data[, grep("_downsampled", names(data))]
    dff <- apply(df, 2, FUN = function(x) x / data$totalValidFrags)
    colnames(dff) <- paste0(colnames(dff), "_mean")
    data <- cbind(data, dff)

    ## normalize downsampled scores to this number:
    df <- data[, grep("_score", names(data))]
    dff <- apply(df, 2, FUN = function(x) x / data$totalValidFrags)
    colnames(dff) <- paste0(colnames(dff), "_mean")
    data <- cbind(data, dff)

    if (!is.null(file)) {
        write.table(data, file, sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE)
    }

    return(data)
}


loadGrl <- function(intDir = "intervals") {
    #' @title Load Intervals into a GenomicRangesList
    #' @description Reads in .rds and .bed files in intDir into a GRangesList object.
    #' @param intDir Directory containing .rds or .bed files with intervals.
    #' @examples
    #' # example code
    #'
    #' #requires precomputed matrices in the resfolder (see getMatrix())
    #' extdata <- system.file("extdata", package="PostChicago")
    #' intDir <- file.path(extdata,'intervals')
    #' dir(intDir)
    #'
    #' #Load intervals for which the boxplots should be quantified into a GRangesList
    #' grl <- suppressWarnings(loadGrl(intDir = intDir))
    #'
    #' @returns GenomicRanges object grl containing all intervals
    #' @export

    if (!file.exists(intDir)) {
        stop("enter a valid directory")
    }

    ## AF: 5.2.2026: Updated with full.names=TRUE, which also includes the path!:
    # ilist = list.files(paste0(intDir)
    ilist <- list.files(paste0(intDir), full.names = TRUE)
    grl <- GenomicRanges::GRangesList()
    names <- rep(NA, length(ilist))
    i <- 0
    for (f in ilist) {
        i <- i + 1
        if (length(grep(".rds$", f)) > 0) {
            ## AF: 5.2.2026: Updated to remove the redundant stripGR function:
            # gr = stripGR(readRDS(f))
            gr <- granges(readRDS(f))
            name <- strsplit_string(f, s2 = ".rds")
        } else if (length(grep(".bed$", f)) > 0) {
            # gr = bed2gr(read.delim(paste0(intDir, "/", f),header=FALSE))
            # ii need to test this still:
            ## Comment AF: Should be working correctly now! But needs to be tested with .rds as well!
            ## AF@ 5.2.2026@ Updated here, since I also had to change the old bed2gr function,
            ## which is now fully unnecessary, since we import directly to gr!
            # gr = bed2gr(rtracklayer::import(list.files(intDir, full.names = TRUE)[1]))
            gr <- rtracklayer::import(f, format = "BED")
            ## AF: 5.2.2026: Updated to extract the correct name from the full path:
            name <- unlist(strsplit(
                strsplit_string(f, s2 = ".bed"),
                split = "/"
            ))[length(
                unlist(strsplit(strsplit_string(f, s2 = ".bed"), split = "/"))
            )]
        } else {
            message(paste(f, "is an unclear file type, please add extension!"))
        }

        grl <- c(grl, GenomicRanges::GRangesList(gr))
        names[i] <- name
    }
    names(grl) <- names

    return(grl)
}


## NOTE: Changed 18.03.2026 AF to readRDS instead of load
## also simplified types and fls

loadCdList <- function(resultsDir = "results", baited_genes = NULL) {
    #' @title Loads List with Interaction Tables
    #' @description Reads in all ChicagoData objects in resultsDir and converts them to a list of ChicagoData tables.
    #' ChicagoData objects must be saved during the Chicago pipeline, which can be done by running Chicago
    #' via runChicagoForPostChicago()
    #' @param baited_genes Optional; Used to control for the compatibility between loaded ChicagoData and current
    #' experiment; GenomicRanges object with at least two metadata columns:
    #' 1) 're_id' (captured fragment ID (bait), corresponding to fragment numbers indicated in rmapgr),
    #' 2) genenames' (names of the corresponding baits, typically gene names).
    #' Can be created automatically from the Chicago baitmap table (see R Bioconductor package 'Chicago')
    #' using the function baitmap2baited_genes().
    #' @param resultsDir Default: 'results'; Directory containing ChicagoData output.
    #' @examples
    #' # example code:
    #' extdata <- system.file("extdata", package="PostChicago")
    #' chicagoOutputDir <- file.path(extdata,'results')
    #' designDir <- file.path(extdata,'designDir')
    #'
    #' #create baited_genes:
    #' baited_genes=baitmap2baited_genes(designDir,save=FALSE)
    #'
    #' #load ChicagoData object list:
    #' cdlist <- loadCdList(resultsDir = chicagoOutputDir, baited_genes = baited_genes)
    #'
    #' @returns List containing ChicagoData tables.
    #' @export

    if (!file.exists(resultsDir)) {
        stop("enter a valid directory")
    }

    fls <- list.files(pattern = "cd", path = resultsDir)

    if (length(grep(".rds", fls)) > 0) {
        types <- strsplit_string(fls, s1 = "cd_", s2 = ".rds")

        L <- list()
        for (f in fls) {
            cd <- readRDS(paste0(resultsDir, "/", f))
            L <- c(L, list(cd@x))
        }
    } else {
        types <- strsplit_string(fls, s1 = "cd_")

        L <- list()
        for (f in fls) {
            load(paste0(resultsDir, "/", f))
            L <- c(L, list(cd@x))
        }
    }

    names(L) <- types

    ids <- unique(L[[1]]$baitID)

    if (is.null(baited_genes)) {
        message(
            "No control for bait count... Add argument if you want to control output"
        )
    } else if (!length(ids) == length(baited_genes)) {
        message(
            "There is an unequal amount of baits in your list and bait files. Please check your Chicago Outputs..."
        )
    }

    return(L)
}


## NOTE: 18.03.2026 AF rewritten to save everything as .rds

runChicagoForPostChicago <- function(
      df,
      mySettings = NULL,
      outputDir = "results",
      DataPath = "dataPath",
      DesignDir = "designDir",
      overwrite = FALSE
) {
    #' @title Run Chicago and save ChicagoData objects.
    #' @description Runs Chicago such that the saved output files are compatible with our pipeline.
    #' Default settings are only executed if the files are not yet saved in outputDir.
    #' @param df data.frame containing at least two columns: 'filename' (.chinput filename) and 'type' (sample type)
    #' @param mySettings Chicago Settings for the pipeline. If not specified, default Chicago settings will be used.
    #' @param outputDir Default 'results'; Output directory in which ChicagoData objects are saved.
    #' @param DataPath Default 'dataPath'; Directory of input files (typically .chinput).
    #' @param DesignDir Default 'designDir'; Directory of Experiment Design Files, see Chicago documentation.
    #' @param overwrite Default: FALSE; Overwrites ChicagoData objects that are present in outputDir if set to TRUE.
    #' @examples
    #'
    #' ##Example code using package-supplied data
    #' extdata <- system.file("extdata", package = "PostChicago")
    #' chicagoOutputDir <- file.path(extdata,'results')
    #' designDir <- file.path(extdata,'designDir')
    #' dataPath <- file.path(extdata, "dataPath")
    #' chicagoOutputDir <- file.path(extdata, "results")
    #'
    #' ## checking which .chinput files we have stored:
    #' dir(dataPath)
    #'
    #' ##Adjust weights settings
    #'
    #' weightsPath <- file.path(
    #' system.file("extdata", package = "Chicago"),
    #' "weights"
    #' )
    #' weightSettings <- file.path(weightsPath, "mESC-2reps.settings")
    #' weightSettings <- read.delim(weightSettings, header = FALSE)
    #'
    #' mySettings <- Chicago::defaultSettings()
    #' mySettings[grep("weight", names(mySettings))] <- weightSettings[, 2]
    #' mySettings$minFragLen <- 100
    #'
    #' ##Create a table with sample info
    #'
    #' fls <- list.files(dataPath)
    #' fls
    #'
    #' df <- data.frame(
    #' filename = fls,
    #' RA = stringr::str_remove(fls, "_rep.*"),
    #' rep = stringr::str_remove(stringr::str_remove(fls, "^.*h_"), "\\.chinput$")
    #' )
    #' df$type <- df$RA
    #' head(df)
    #'
    #' ##Now you can run Chicago!
    #'
    #' runChicagoForPostChicago(
    #'  df = df, mySettings = mySettings, outputDir = chicagoOutputDir,
    #'  DataPath = dataPath, DesignDir = designDir
    #' )
    #'
    #' @returns Saves complete ChicagoData objects containing score and read counts for all ditags, including
    #' non-significant interactions and named after df$type. Also saves standard Chicago outputs in outputDir.
    #' @export

    if (is.null(mySettings)) {
        mySettings <- Chicago::defaultSettings()
    }

    types <- unique(df$type)

    for (type in types) {
        file <- paste0("cd", "_", type, ".rds") ## rewritten 18.3.2026

        ## To execute only if cd is not in outputDir or if overwrite is set to TRUE
        if (!file %in% dir(outputDir)) {
            files <- paste(DataPath, df[df$type == type, ]$filename, sep = "/")

            cd <- Chicago::setExperiment(designDir = DesignDir, settings = mySettings)
            cd <- Chicago::readAndMerge(files = files, cd = cd)
            cd <- Chicago::chicagoPipeline(
                cd,
                outprefix = paste(outputDir, "/", type, sep = "")
            )
            Chicago::exportResults(
                cd,
                file.path(outputDir, paste("vignetteOutput", type, sep = "_"))
            )
            saveRDS(cd, file = paste0(outputDir, "/", file)) ## changed 18.3.2026 using file
        } else if (overwrite) {
            files <- paste(DataPath, df[df$type == type, ]$filename, sep = "/")

            cd <- Chicago::setExperiment(designDir = DesignDir, settings = mySettings)
            cd <- Chicago::readAndMerge(files = files, cd = cd)
            cd <- Chicago::chicagoPipeline(
                cd,
                outprefix = paste(outputDir, "/", type, sep = "")
            )
            Chicago::exportResults(
                cd,
                file.path(outputDir, paste("vignetteOutput", type, sep = "_"))
            )
            # save(cd, file = paste(outputDir, "/cd", "_", type, sep = ""))
            saveRDS(cd, file = paste0(outputDir, "/", file)) ## changed 18.3.2026 using file
        } else {
            message(paste("cd file", file, "detected, skipping Chicago..."))
        }
    }
}


########################################
## Make tracks for
## Genome Browser visualization
########################################

## AF 06.02.2026: Updated function, using rtracklayer::export:
## Now also works with IGV, since the files are saved as .bdg

makeBedGraphs <- function(
      id,
      L_data,
      rmapgr,
      baited_genes,
      ints = NULL,
      L_norm = NULL,
      outputDir = "bedgraph",
      smoothingWindow = "off"
) {
    #' @title Create Bedgraphs From ChicagoData
    #' @description Takes ChicagoData objects supplied in the list L, normalizes them to read coverage mapped to baits.
    #' Creates two file types:
    #' 1) a bedgraph file containing normalized reads for each sample mapping to the bait
    #' defined by id,
    #' 2) a bed file of the view point (bait).
    #' @param id baitID, must correspond to the baitIDs defined in baited_genes and the supplied rmapgr object
    #' @param L_data List of ChicagoData (cd) objects.
    #' @param rmapgr GenomicRanges object containing restriction fragments. IDs must be saved in a column called 'id'.
    #' @param baited_genes GenomicRanges object with at least two metadata columns:
    #' 1) 're_id' (captured fragment ID (bait), corresponding to fragment numbers indicated in rmapgr),
    #' 2) genenames' (names of the corresponding baits, typically gene names).
    #' Can be created automatically from the Chicago baitmap table (see R Bioconductor package 'Chicago')
    #' using the function baitmap2baited_genes().
    #' @param ints optional; table containing all significant interactions across L_data, usually created using,
    #' makeIntsTable(). If supplied, a separate bed file containing all significant interactions mapping to the id
    #' is created
    #' @param L_norm Separate list of ChicagoData tables which should be used to normalize samples.
    #' Useful if bedgraphs should be created only for a subset of promoters or genomic locations, but the
    #' normalization performed based on the whole dataset. If not provided, L_data is used for normalization.
    #' @param outputDir Default: 'bedgraph'. Output directory for the .bedgraphs, automatically created if not currently present.
    #' @param smoothingWindow Default: 'off'. Over how many restriction fragments should the data be smoothed?
    #' Not recommended for UCSC visualization, but can improve visualization on IGV.
    #' @returns Saves bedgraphs and annotation bedfiles of interactions with id as view point for each sample in L_data.
    #' @examples
    #'
    #' # example code:
    #'
    #' extdata <- system.file("extdata", package="PostChicago")
    #' chicagoOutputDir <- file.path(extdata,'results')
    #' designDir <- file.path(extdata,'designDir')
    #'
    #' #create baited_genes and rmapgr files:
    #' baited_genes <- baitmap2baited_genes(designDir,save=FALSE)
    #' rmapgr <- rmap2rmapgr(designDir, save = FALSE)
    #'
    #' #load ChicagoData object list:
    #' cdlist <- loadCdList(resultsDir = chicagoOutputDir, baited_genes = baited_genes)
    #'
    #' ## extract all view points (baits):
    #' baits <- unique(unlist(lapply(cdlist,
    #'                              FUN = function(x) unique(x$baitID))))
    #'
    #' ## create bedgraphs for selected view points and store them in outdir:
    #' outdir <- "bedgraphs"
    #' dir.create(outdir)
    #' makeBedGraphs(baits[1], cdlist, rmapgr,
    #'                baited_genes = baited_genes,
    #'                outputDir = outdir
    #' )
    #'
    #' @export

    ## make sure genenames are supplied!
    if (is.null(names(baited_genes))) {
        if (!is.null(baited_genes$genename)) {
            names(baited_genes) <- baited_genes$genename
        } else if (!is.null(baited_genes$genenames)) {
            names(baited_genes) <- baited_genes$genenames
        } else {
            message("no genenames supplied with baited_genes!")
        }
    }

    ## make sure the outputdir exists, if not create one:
    if (!dir.exists(outputDir)) {
        dir.create(outputDir)
    }

    name <- baited_genes[baited_genes$re_id == id]$genename
    ## AF. 5.2.2026: updated to not be hardcoded chr names!
    ## Now they are extracted from baited_genes, since only these names need to be used (cis-chromosomal interactions!)
    # seqnames=paste('chr',c(1:19,'X','Y'),sep='')
    seqnames <- as.character(unique(seqnames(baited_genes)))

    ## get normalization parameters:
    L <- L_data
    if (is.null(L_norm)) {
        L_norm <- L
    }

    ## for each sample write a normalized bg
    samples <- names(L_data)

    bait <- rmapgr[rmapgr$id == id]

    for (i in 1:length(L)) {
        sample <- samples[i]

        ## only if bedgraphs do not already exist
        bgfile <- paste0(outputDir, "/", name, "_id", id, "_", sample, ".bdg")
        if (!file.exists(bgfile)) {
            print(i)

            ## AF: 6.2.2026: Change from here:

            ## Create bedgraph-like object, a GRanges:
            bg <- getnormbed(
                id,
                data = L[[i]],
                samplename = names(L)[i],
                L = L_norm,
                same.chromosome = TRUE,
                rmapgr = rmapgr
            )
            bg <- bg[order(start(bg), decreasing = FALSE)]
            bg <- bg[GenomicRanges::seqnames(bg) %in% seqnames]
            bg. <- granges(bg)
            bg.$score <- bg$N
            bg <- bg.
            rm(bg.)

            ## save bedgraph with a header (will be opening on UCSC):
            a <- unlist(strsplit(
                sapply(unlist(strsplit(bgfile, split = ".bdg$")), FUN = function(x) {
                    x[1]
                }),
                split = "/"
            ))
            trackname <- paste0(strsplit(a, split = "/")[[length(strsplit(
                a,
                split = "/"
            ))]])

            ## Trackline (header) of the bedgraph
            track <- paste0(
                "track type=bedGraph name=",
                trackname,
                " description=",
                trackname,
                " visibility=full maxHeightPixels=100 priority=10 autoScale=off alwaysZero=on gridDefault=off",
                " graphType=bar windowingFunction=mean viewLimits=0:100 smoothingWindow=",
                smoothingWindow
            )

            con <- file(bgfile, "w")
            writeLines(track, con)
            rtracklayer::export(bg, con, format = "bedGraph")
            close(con)
        }
    }

    ## make bedfiles of bait and of interactions, but only if the files do not exist yet:
    viewbedfile <- paste0(outputDir, "/", name, "_id", id, "_viewpoint.bed")

    if (!file.exists(viewbedfile)) {
        ## AF: 6.2.26: Updated the redundant changes to df etc
        ## # bedbait=gr2bed(bait)
        # bedbait=as.data.frame(bait)[,1:3]
        # saveBed(bedbait, viewbed)
        rtracklayer::export(bait, viewbedfile, format = "bed")
    }

    intsbedfile <- paste0(outputDir, "/", name, "_id", id, "_ints.bed")

    if (!is.null(ints) && !file.exists(intsbed)) {
        intsInId <- ints[ints$baitID == id, ]$otherEndID
        if (length(intsInId) > 0) {
            ## AF: 5.2.26: Updated the redundant gr2bed function
            # bedints=gr2bed(rmapgr[rmapgr$id %in% bedints])
            intsbed <- rmapgr[rmapgr$id %in% intsInId]
            rtracklayer::export(intsbed, con = intsbedfile, format = "bed")
        }
    }
}


######################################################
## Aggregate interaction peaks over a specific distance
######################################################

## Function to aggregate peaks for interactions from CaptureC analysis (version 16.08.2022):
# 08Feb25: SK added outdir
aggregatePeaks <- function(
      ints,
      dis,
      samples,
      fileprefix = "ints",
      outdir = NULL
) {
    #' @title Aggregate Interactions to Interaction Peaks
    #' @description Takes a table of interactions (ints) and aggregates significant interactions within the distance
    #' dis in restriction fragments for each bait. Works similarly to the reduce() function from GenomicRanges.
    #' Caution: All interactions, regardless the sample from which they originated, are aggregated, creating a 'consensus'
    #' interaction peak set. If interaction peaks should be aggregated for individual samples, we recommend a prior
    #' sample-specific filtering of ints.
    #' Only creates the table if it is not yet saved, otherwise loads the existing table.
    #' @param ints Table containing all interactions that should be aggregated.
    #' @param dis Numeric, distance in restriction fragments over which interactions should be aggregated.
    #' dis=1 will merge two adjacent fragments with significant interactions.
    #' @param samples Character vector of sample names, should correspond to the sample names used in ints.
    #' Typically corresponding to the names provided in the list of ChicagoData objects cd from which ints
    #' has been generated.
    #' @param fileprefix  Default: ints; Prefix of the filename for the aggregatePeaks table file, will be added the following extension:
    #' paste0('_aggregatePeaks_in_regions_',dis,'bp.txt').
    #' @param outdir Default: current directory; Directory where the final table is saved.
    #' Will not create a copy if already present.
    #' @examples
    #' # example code:
    #' extdata <- system.file("extdata", package="PostChicago")
    #' outputDir <- file.path(extdata,'postchicago')
    #'
    #' #read in the interactions table:
    #' ints=read.delim(paste0(outputDir,'/ints_all.txt'),stringsAsFactors=FALSE)
    #'
    #' #aggregate all peaks within a distance of 5 restriction fragments
    #' dis=5
    #' samples=names(ints)[grep('_N$',names(ints))]
    #' samples=strsplit_string(samples,s2='_N')
    #' aggregatePeaks(ints,dis,samples)
    #'
    #' @returns A table containing the positions of aggregated interactions. Redefines intIDs so that the ID contains all
    #' significant interactions that were aggregated.
    #' Data is automatically saved in a 'aggregatePeaks_dis' text file whose name is defined by fileprefix and dis.
    #' Also saved as a bedfile.
    #' @export

    if (is.null(outdir)) {
        outdir <- getwd()
    }

    file <- paste0(fileprefix, "_aggregatePeaks_dis", dis, ".txt")

    if (file %in% list.files(path = outdir)) {
        message(paste("dis", dis, ", aggregatePeaks file exists"))
        d <- read.delim(paste0(outdir, "/", file), stringsAsFactors = FALSE)
    } else {
        message(paste("dis", dis, ", aggregating peaks..."))
        d <- ints
        L_peaks <- list()

        ids <- unique(ints$baitID)

        for (id in ids) {
            a <- d[d$baitID == id, ]
            gr <- GenomicRanges::GRanges(
                a$seqnames_otherEnd,
                IRanges::IRanges(a$otherEndID, a$otherEndID)
            )
            values(gr) <- a
            gr <- GenomicRanges::reduce(gr, min.gapwidth = dis)

            L_ints <- list()
            for (i in 1:length(gr)) {
                oe <- start(gr[i]):end(gr[i])
                L_ints[[i]] <- a[a$otherEndID %in% oe, ]
            }

            ## aggregate to peaks

            ## 19.3.2026 (AF): Updated the following problem after talking with Ezo:
            ## Right now the function aggregates peaks and then takes the first
            ## aggregated interaction into the aggregated_peaks table, keeping reads, scores etc.
            ## This is misleading. Therefore, moving forward remove all unnecessary info and quantifications.
            ## Now we only keep columns that have nothing to do with quantification, quantification will be done later!

            for (i in 1:length(L_ints)) {
                b <- L_ints[[i]]
                c <- b[1, ]

                if (nrow(b) > 1) {
                    ## intID should contain all significant interactions that went into the aggregation
                    intid <- b[1, ]$intID
                    oe <- b[1, ]$otherEndID

                    for (j in 2:nrow(b)) {
                        intid <- paste(intid, b[j, ]$intID, sep = ",")
                    }
                    for (j in 2:nrow(b)) {
                        oe <- paste(oe, b[j, ]$otherEndID, sep = ",")
                    }

                    c$intID <- intid
                    c$otherEndID <- oe
                }

                c$start_otherEnd <- min(b$start_otherEnd)
                c$end_otherEnd <- max(b$end_otherEnd)

                ## remove all unnecessary columns from c, including scores, reads etc.
                ## only keep regional annotations
                keep <- c(
                    "intID",
                    "baitID",
                    "otherEndID",
                    "seqnames_bait",
                    "start_bait",
                    "end_bait",
                    "seqnames_otherEnd",
                    "start_otherEnd",
                    "end_otherEnd"
                )
                c <- c[, keep]

                L_peaks <- c(L_peaks, list(c))
            }
        }

        df <- do.call("rbind", L_peaks)
        d <- df
        d$clusters_refined <- d$clusters

        saveBed(
            d[, c("seqnames_otherEnd", "start_otherEnd", "end_otherEnd", "intID")],
            paste0(outdir, "/", file, ".bed")
        )
        write.table(
            d,
            paste0(outdir, "/", file),
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            row.names = FALSE
        )
    }

    return(d)
}


## aggregate by distance in bp:
# 08Feb25: SK added outdir

aggregatePeaks_regions <- function(
      ints,
      dis,
      samples,
      fileprefix = "ints",
      outdir = NULL
) {
    #' @title Aggregate Interactions to Interaction Peaks by Distance in bp
    #' @description Takes a table of interactions (ints) and aggregates significant interactions within the distance
    #' dis in base pairs (bp) for each bait. Works similarly to the reduce() function from GenomicRanges.
    #' Caution: All interactions, regardless the sample from which they originated, are aggregated, creating a 'consensus'
    #' interaction peak set. If interaction peaks should be aggregated for individual samples, we recommend a prior
    #' sample-specific filtering of ints.
    #' Only creates the table if it is not yet saved, otherwise loads the existing table.
    #' @param ints Table containing all interactions that should be aggregated.
    #' @param dis Numeric, distance in bp over which interactions should be aggregated.
    #' dis=1 will merge two adjacent fragments with significant interactions.
    #' @param samples Character vector of sample names, should correspond to the sample names used in ints.
    #' Typically corresponding to the names provided in the list of ChicagoData objects cd from which ints
    #' has been generated.
    #' @param fileprefix  Default: ints; Prefix of the filename for the aggregatePeaks table file, will be added the following extension:
    #' paste0('_aggregatePeaks_in_regions_',dis,'bp.txt').
    #' @param outdir Default: current directory; Directory where the final table is saved.
    #' Will not create a copy if already present.
    #' @examples
    #' # example code:
    #' extdata <- system.file("extdata", package="PostChicago")
    #' outputDir <- file.path(extdata,'postchicago')
    #'
    #' #read in the interactions table:
    #' ints=read.delim(paste0(outputDir,'/ints_all.txt'),stringsAsFactors=FALSE)
    #'
    #' #aggregate all peaks within a distance of 5 restriction fragments
    #' dis=5000
    #' samples=names(ints)[grep('_N$',names(ints))]
    #' samples=strsplit_string(samples,s2='_N')
    #' aggregatePeaks_regions(ints,dis,samples)
    #'
    #' @returns A table containing the position of aggregated interactions. Redefines intIDs so that the ID of the
    #' first interacting fragment in a peak stands in as otherEnd.
    #' Data is automatically saved in a 'aggregatePeaks_dis' text file whose name is defined by fileprefix and dis.
    #' Also saved as a bedfile.
    #' @export
    # Example usage:
    # ints=ints
    # dis=5000
    # samples=names(L_CapC10)
    # fileprefix="ints"

    if (is.null(outdir)) {
        outdir <- getwd()
    }

    file <- paste0(fileprefix, "_aggregatePeaks_dis", dis, "bp.txt")

    if (file %in% list.files(path = outdir)) {
        message(paste("distance of ", dis, ",aggregatePeaks file exists"))

        d <- read.delim(paste0(outdir, "/", file), stringsAsFactors = FALSE)
    } else {
        message(paste("distance of", dis, ", aggregating peaks..."))

        d <- ints
        L_peaks <- list()
        ids <- unique(ints$baitID)

        for (id in ids) {
            a <- d[d$baitID == id, ]
            a <- a[order(a$seqnames_otherEnd, a$start_otherEnd), ]
            gr <- GenomicRanges::GRanges(
                a$seqnames_otherEnd,
                IRanges::IRanges(a$start_otherEnd, a$end_otherEnd)
            )

            values(gr) <- a
            gr <- GenomicRanges::reduce(gr, min.gapwidth = dis, with.revmap = TRUE)

            revmap <- gr$revmap

            otherEndID <- a$otherEndID

            mapped <- lapply(revmap, function(idx) otherEndID[idx])

            L_ints <- list()

            for (i in 1:length(mapped)) {
                oe <- mapped[[i]]
                L_ints[[i]] <- a[a$otherEndID %in% oe, ]
            }

            ## aggregate to peaks

            for (i in 1:length(L_ints)) {
                b <- L_ints[[i]]
                c <- b[1, ]

                if (nrow(b) > 1) {
                    ## intID should contain all significant interactions that went into the aggregation
                    intid <- b[1, ]$intID
                    oe <- b[1, ]$otherEndID

                    for (j in 2:nrow(b)) {
                        intid <- paste(intid, b[j, ]$intID, sep = ",")
                    }
                    for (j in 2:nrow(b)) {
                        oe <- paste(oe, b[j, ]$otherEndID, sep = ",")
                    }

                    c$intID <- intid
                    c$otherEndID <- oe
                }

                c$start_otherEnd <- min(b$start_otherEnd)
                c$end_otherEnd <- max(b$end_otherEnd)

                ## AF: 19.03.2026:
                ## remove all unnecessary columns from c, including scores, reads etc.
                ## only keep regional annotations
                keep <- c(
                    "intID",
                    "baitID",
                    "otherEndID",
                    "seqnames_bait",
                    "start_bait",
                    "end_bait",
                    "seqnames_otherEnd",
                    "start_otherEnd",
                    "end_otherEnd"
                )
                c <- c[, keep]

                L_peaks <- c(L_peaks, list(c))
            }
        }

        df <- do.call("rbind", L_peaks)

        d <- df

        if (is.null(outdir)) {
            outdir <- getwd()
        }

        saveBed(
            d[, c("seqnames_otherEnd", "start_otherEnd", "end_otherEnd", "intID")],
            paste0(outdir, "/", file, ".bed")
        )
        write.table(
            d,
            paste0(outdir, "/", file),
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            row.names = FALSE
        )

        return(d)
    }
}


####################################################################################################
####################################################################################################
## main new getMatrix() function:
## replaced on 4.5.2026 (AF) due to the new makeOneGeneOnePeak()
####################################################################################################
####################################################################################################

## Newest function:

# 08Feb26: SK added outdir
# AF: 09.04.2026: Changed description, 4.5.2026: Added examples
getMatrix <- function(
      L,
      zoom = 100,
      d,
      resfolder = "matrices/",
      type = "reads",
      readnorm = "downsampling",
      name = "",
      norm = FALSE,
      normsam = 1,
      pseudo = 1,
      rmapgr = rmapgr
) {
    #' @title Create Interaction Matrices
    #' @description
    #' Function that creates and normalizes matrices containing reads or scores of fragment-to-fragment interactions
    #' surrounding the interactions that are provided in the interactions table `d`.
    #'
    #' Output is a list of matrices for the samples supplied in `L`, including distance-matched control matrices.
    #' The matrix list is automatically saved in the `resfolder`.
    #' Matrices are stranded, such that the promoter is to the right (+ direction) of the interaction, which is
    #' important because the read numbers are increasing with increased promoter proximity.
    #' A new list is only created and saved as `.rds` file if the matrix cannot be found in the specified `resfolder`,
    #' otherwise the preexisting matrix list is loaded.
    #' If the summit of an interaction peak or oneGeneOnePeak file is provided, getMatrix() extracts the summit info
    #' from the data, otherwise the function simply determines the central fragment of an interaction peak.
    #'
    #' If norm = TRUE, normalized matrices are calculated from pre-exsisting unnormalized matrices.
    #' All values are normalized to the interaction in the reference sample defined by `normsam`,
    #' which is recommended to be a positive control sample.
    #' @param L list of ChicagoData tables.
    #' @param zoom Distance around interaction in restriction fragments for which matrices should be generated.
    #' @param d A table with cis-chromosomal interactions. Should contain the columns `baitID`, `otherEndID`,
    #' defined as in Chicago and can have a column `intID`, which is expected to be in the format `baitID;otherEndID`.
    #' `otherEndID` should be the center of interaction, from which the matrices will be calculated.
    #' @param resfolder Default: 'matrices/'. Path to the directory, where matrices should be stored.
    #' Will not create a copy if one is already saved in the given directory. Creates the directory if nonexistant.
    #' @param type 'reads' (default)/'scores'. If reads, the function performs downsampling relative to the sample
    #' with the smallest amount of reads that are mapped to baits.
    #' @param readnorm Default: 'downsampling'. Read normalization, either 'downsampling' or 'scaling'.
    #' Scaling means the same type of normaliyation as for the lineplots in plotInteractions().
    #' @param name Name of the matrix.
    #' @param norm TRUE/FALSE (default). Should the data be normalized?
    #' @param normsam numeric; defines to which sample in the list L should the data be normalized.
    #' Default: 1 (the first sample!). Only used if norm=TRUE.
    #' @param pseudo Default: 1; pseudocount to be added during the normalization between matrices.
    #' @param rmapgr GenomicRanges object containing restriction fragments within the reference genome.
    #' IDs must be saved in a column called 'id'.
    #'
    #' @examples
    #' # example code
    #'
    #' ##create ints table with 4 interactions
    #' ints <- data.frame(intID=c('11;4', '11;7', '11;8','12;14'),
    #'                    baitID = as.character(c(11,11,11,12)),
    #'                    otherEndID=as.character(c(4,7,8,14)),
    #'                    seqnames_bait=as.character((paste0('chr',rep(1,4)))),
    #'                    start_bait=as.character(c(21,21,21,23)),
    #'                    end_bait=as.character(c(22,22,22,23)),
    #'                    seqnames_otherEnd=as.character((paste0('chr',rep(1,4)))),
    #'                    start_otherEnd=as.character(c(7,13,15,27)),
    #'                    end_otherEnd=as.character(c(8,14,16,28))
    #'                    )
    #'
    #' ##make rmapgr
    #' rmapgr <- GenomicRanges::GRanges(rep('chr1',20),
    #' IRanges::IRanges( c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39),
    #' c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40)  ) )
    #' rmapgr$id=1:length(rmapgr)
    #'
    #' cdlist <- list(
    #' sample1 = data.frame(baitID=c(rep(11,20),rep(12,20)),
    #'                       otherEndID=rep(c(1:20),2),
    #'                       N=sample(1:20,40,replace=TRUE),
    #'                       N.1=sample(1:20,40,replace=TRUE),
    #'                       N.2=sample(1:20,40,replace=TRUE),
    #'                       score=sample(1:5,40,replace=TRUE)),
    #'
    #' sample2 = data.frame(baitID=c(rep(11,20),rep(12,20)),
    #'                       otherEndID=rep(c(1:20),2),
    #'                       N=sample(15:45,40,replace=TRUE),
    #'                       N.1=sample(15:45,40,replace=TRUE),
    #'                       N.2=sample(15:45,40,replace=TRUE),
    #'                       score=sample(1:20,40,replace=TRUE))
    #'       )
    #'
    #' ##getMatrix
    #' getMatrix(L = cdlist, d = ints, zoom = 1, rmapgr=rmapgr)
    #'
    #' ##normalized matrices
    #' getMatrix(L = cdlist, d = ints, zoom = 1, rmapgr=rmapgr,norm=TRUE,normsam=2)
    #'
    #' ##Add summit info to interactions:
    #' ints$summit_re_id=c('4','7','8','13,14')
    #' getMatrix(L = cdlist, d = ints, zoom = 1, rmapgr=rmapgr, resfolder='matrices_summits/')
    #'
    #' @return List of interaction matrices.
    #' @export

    message(paste("Matrices of", type))

    ## if resfolder does not exist, create it!
    if (!file.exists(resfolder)) {
        dir.create(resfolder)
    }

    ## stop if type is not reads or scores:
    if (!type %in% c("reads", "scores")) {
        stop("type must be 'reads' or 'scores'")
    }

    ## define filename:
    filename <- paste0(resfolder, "matrices_", name, "_", type, ".rds") ## changed RDS 18.03.2026 (AF)

    ## include read normalization into the filename (downsampling or scaling!)
    if (type == "reads") {
        filename <- paste0(
            resfolder,
            "matrices_",
            name,
            "_",
            type,
            "_",
            readnorm,
            ".rds"
        ) ## changed RDS 18.03.2026 (AF)
    }

    if (!file.exists(filename)) {
        ## keep only interactions on the same chromosome!
        d <- d[d$seqnames_bait == d$seqnames_otherEnd, ]

        ## AF 4.5.2026: Adjusted to fit the new oneGeneOnePeak tables
        ## that contain as intID values like 1;1,2 indicating that the interaction
        ## occurs with an interval that overlaps re_id 1 and 2.
        ## Instead, we will redefine intID as bait;<central re_id ov interaction peak>
        ## These central IDs ('summits!' area lso saved in  the new ogop tables under
        ## summit_re_id)
        ## We remodel d, for cases where summit_re_id contains two entries
        ## (for instance: 3,4) to have two rows that contain both!
        ##
        ## We do this, however, only if the column summit_re_id exists in d
        if ("summit_re_id" %in% names(d)) {
            message("extracting summit info...")

            d$otherEndID <- d$summit_re_id
            d$intID <- paste(d$baitID, d$otherEndID, sep = ";")

            ## split up d into unique and multi re_id sites!
            ## only if there are summits covering ultiple restriction frags!
            if (length(grep(",", d$otherEndID)) > 0) {
                uniq <- d[-grep(",", d$otherEndID), ]
                multi <- d[grep(",", d$otherEndID), ]

                final <- uniq

                ## make two rows out of multi
                for (i in 1:nrow(multi)) {
                    v <- unlist(strsplit(as.character(multi$otherEndID), split = ","))
                    multi2 <- multi[i, ]
                    j <- 1
                    while (j < length(v)) {
                        multi2 <- rbind(multi2, multi[i, ])
                        j <- j + 1
                    }

                    multi2$otherEndID <- v
                    multi2$intID <- paste(multi2$baitID, multi2$otherEndID, sep = ";")

                    final <- rbind(final, multi2)
                }

                d <- final
            }
        }

        ## test d and d$intID for uniqueness
        ## In case of non-uniqueness: Keep just one intID
        if (length(d$intID) != length(unique(d$intID))) {
            intid.counts <- table(d$intID)
            final <- d[d$intID %in% names(intid.counts[intid.counts == 1]), ]

            for (id in names(intid.counts[intid.counts > 1])) {
                final <- rbind(final, d[d$intID %in% names(intid.counts[intid.counts > 1]), ][1, ])
            }

            d <- final
        }

        ## convert factors to characters or numbers:
        d$intID <- as.character(d$intID)
        d$baitID <- as.numeric(as.character(d$baitID))
        d$otherEndID <- as.numeric(as.character(d$otherEndID))

        ## step1: normalize reads in L
        message("normalizing reads...")

        if (readnorm == "downsampling") {
            message("downsampling reads...")

            nreads <- sapply(L, FUN = function(x) sum(x$N))
            downsamplingFacs <- min(nreads) / nreads * 0.9999
            Ldownsampled <- L
            for (i in 1:length(Ldownsampled)) {
                Ldownsampled[[i]]$N <- Ldownsampled[[i]]$N * downsamplingFacs[i]
            }
        } else {
            message(
                "scaling reads upon coverage normalisation using getsk() *100000 reads and *nproms..."
            )

            scalingFacs <- getsk(L)
            Ldownsampled <- L
            for (i in 1:length(Ldownsampled)) {
                Ldownsampled[[i]]$N <- Ldownsampled[[i]]$N * scalingFacs[i]
            }
        }

        ## step2: create a list of matrices Lm of length nrow(d) with all 0's; add intIDs as row names
        message("creating empty matrix lists...")

        m <- matrix(0, nrow = nrow(d), ncol = (2 * zoom + 1))
        row.names(m) <- d$intID

        k <- 1
        Lm <- list()
        while (k <= length(Ldownsampled)) {
            Lm <- c(Lm, list(m))
            k <- k + 1
        }

        names(Lm) <- names(Ldownsampled)

        ## step3: create a list of control matrices Lmc of length nrow(d) with all 0's; add intIDs as row names; define otherEnds

        Lmc <- Lm
        names(Lmc) <- paste0(names(Lm), "_con")

        ## step4: populate Lm and Lmc with downsampled values
        message("making matrices, this may take a while...")

        if (type == "reads") {
            ## annotate Lm, the list of normal matrices:
            for (i in 1:length(Ldownsampled)) {
                matr <- lapply(
                    d$intID,
                    FUN = annotateSampleMatrix,
                    downsampledReads = Ldownsampled[[i]],
                    zoom = zoom
                )
                matr <- do.call("rbind", matr)
                Lm[[i]] <- matr
                row.names(Lm[[i]]) <- d$intID
                gc()
            }

            ## annotate Lmc, the list of control matrices:
            message(
                "making distance-matched control matrices, this may take a while..."
            )
            for (i in 1:length(Ldownsampled)) {
                matr <- lapply(
                    d$intID,
                    FUN = annotateControlMatrix,
                    downsampledReads = Ldownsampled[[i]],
                    zoom = zoom
                )
                matr <- do.call("rbind", matr)
                Lmc[[i]] <- matr
                row.names(Lmc[[i]]) <- d$intID
                gc()
            }
        } else if (type == "scores") {
            ## annotate Lm, the list of normal matrices:
            for (i in 1:length(Ldownsampled)) {
                matr <- lapply(
                    d$intID,
                    FUN = annotateSampleMatrixScores,
                    downsampledReads = Ldownsampled[[i]],
                    zoom = zoom
                )
                matr <- do.call("rbind", matr)
                Lm[[i]] <- matr
                row.names(Lm[[i]]) <- d$intID
                gc()
            }

            ## annotate Lmc, the list of control matrices:
            message(
                "making distance-matched control matrices, this may take a while..."
            )
            for (i in 1:length(Ldownsampled)) {
                matr <- lapply(
                    d$intID,
                    FUN = annotateControlMatrixScores,
                    downsampledReads = Ldownsampled[[i]],
                    zoom = zoom
                )
                matr <- do.call("rbind", matr)
                Lmc[[i]] <- matr
                row.names(Lmc[[i]]) <- d$intID
                gc()
            }
        }

        ## step5: flip if the otherEndID is larger than baitID
        message("flipping the matrices...")

        biggerthan <- d$otherEndID > d$baitID
        smallerthan <- d$otherEndID < d$baitID

        Lmflipped <- Lm
        Lmcflipped <- Lmc

        for (i in 1:length(Lm)) {
            upstream <- matrix(Lm[[i]][smallerthan, ], ncol = ncol(Lm[[i]]))
            downstream <- t(apply(matrix(Lm[[i]][biggerthan, ], ncol = ncol(Lm[[i]])), 1, FUN = function(x) rev(x)))

            Lmflipped[[i]] <- rbind(upstream, downstream)

            row.names(Lmflipped[[i]]) <- c(
                row.names(Lm[[i]])[smallerthan],
                row.names(Lm[[i]])[biggerthan]
            )
        }

        for (i in 1:length(Lmc)) {
            downstream <- matrix(Lmc[[i]][biggerthan, ], ncol = ncol(Lm[[1]]))
            upstream <- t(apply(matrix(Lmc[[i]][smallerthan, ], ncol = ncol(Lm[[1]])), 1, FUN = function(x) rev(x)))

            Lmcflipped[[i]] <- rbind(upstream, downstream)

            row.names(Lmcflipped[[i]]) <- c(
                row.names(Lmc[[i]])[smallerthan],
                row.names(Lmc[[i]])[biggerthan]
            )
        }

        ## step6: combine matrices and save them in the resfolder
        message("saving the matrices...")
        Lm <- c(Lmflipped, Lmcflipped)
        Lm <- lapply(Lm, FUN = function(x) x[order(row.names(x)), ])
        saveRDS(Lm, file = filename)
    } else {
        ## if file exists already, load the data into an object called Lm
        message("raw matrix file exists, loading the matrices...")
        Lm <- readRDS(filename) ## changed 18.03.2026 (AF)
    }

    ## step7: normalize matrices to the central 3 frags of the normsam matrix, adding pseudocount of 1

    if (norm) {
        filename2 <- paste0(filename, "_normto.", names(L)[normsam], ".rds") ## included .rds 18.03.2026 (AF)

        ## only execute if file doesn't exist yet!
        if (!file.exists(filename2)) {
            message("making normalized matrices...")

            Lmnorm <- Lm

            normData <- Lm[[normsam]]
            normValues <- normData[, (2 * zoom + 1 - zoom)]

            for (i in 1:length(Lm)) {
                for (j in 1:nrow(Lm[[i]])) {
                    ## if peak summit equals zero, normalize to two adjacent frags:
                    if (normValues[j] == 0) {
                        normValues[j] <- mean(c(
                            normData[j, (2 * zoom + 1 - zoom - 1)],
                            normData[j, (2 * zoom + 1 - zoom + 1)]
                        ))
                        Lmnorm[[i]][j, ] <- (Lm[[i]][j, ] + pseudo) /
                            (normValues[j] + pseudo)
                    } else {
                        Lmnorm[[i]][j, ] <- (Lm[[i]][j, ] + pseudo) /
                            (normValues[j] + pseudo)
                    }
                }
            }

            # save(Lm, file=filename2)
            Lm <- Lmnorm
            saveRDS(Lm, file = filename2) ## changed 18.03.2026 (AF)
        } else {
            message("loading normalized matrices...")
            Lm <- readRDS(filename2) ## changed 18.03.2026 (AF)
        }
    }

    return(Lm)
}
