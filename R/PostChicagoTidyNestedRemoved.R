#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import GenomicRanges
#' @import S4Vectors
#' @import pheatmap
#' @import IRanges
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom Chicago defaultSettings
#' @importFrom Chicago setExperiment
#' @importFrom Chicago readAndMerge
#' @importFrom Chicago chicagoPipeline
#' @importFrom Chicago exportResults
#' @importFrom utils read.delim str write.table
#' @importFrom graphics abline barplot boxplot legend lines par rect text polygon
#' @importFrom grDevices dev.off palette.colors pdf png rgb
#' @importFrom stats density
## usethis namespace: end
NULL

##global variables:
utils::globalVariables(c("grl", "resfolder", "m", "d", "i", "a", "zoom",
                         "oeids", "rmapgr", "baited_genes", "L", "type",
                         "troubleshooting", "cd", "bg", "scalingL", "id", "bait", "Lm"))

testfunc <- function() {
  message("this totally worked")
}


colMedians <- function(df) {
  v = vector()
  for (i in 1:dim(df)[2]) {
    v = c(v, median(df[, i], na.rm = TRUE))
  }
  return(v)
}


colMeans <- function(m) {
  ncols <- ncol(m)
  v = vector()
  for (i in 1:ncols) {
    a = m[, i]
    num = mean(a)
    v = c(v, num)
  }
  return(v)
}


bed2gr <- function(bed) {
  #' @title Convert BED3 to Genomic Ranges
  #' @description Takes an object containing a table in the BED3 format and returns a GenomicRanges object
  #' @param bed R object containing a table in the BED3 format (chr, strat, end)
  #' @examples
  #' # basic usage of bed2gr:
  #' bed <- data.frame(chr='chr1',start=1000, end=2000)
  #' gr <- bed2gr(bed)
  #' gr
  #' @return GenomicRanges object
  #' @export

  gr = GenomicRanges::GRanges(bed[, 1], IRanges::IRanges(bed[, 2], bed[, 3]))
  if (dim(bed)[2] > 3) {
    names(gr) = bed[, 4]
  }
  if ("strand" %in% names(bed)) {
    strand(gr) = bed$strand
  } else if ("Strand" %in% names(bed)) {
    strand(gr) = bed$Strand
  }
  return(gr)
}

stripGR <- function(gr) {
  #' @title Remove Metadata From Genomic Ranges Objects
  #' @description Takes a GenomicRanges object with values (metadata) and returns a minimal GenomicRanges object containing only chr,start,end,strand info
  #' @param gr GenomicRanges object to be stripped
  #' @return minimal GenomicRanges object (chr,start,end,strand info)
  #' @examples
  #' # basic usage of stripGR:
  #' df <- data.frame(chr='chr1',start=1000, end=2000, id='xyz')
  #' gr <- bed2gr(df)
  #' gr$id=df$id
  #' gr
  #' stripGR(gr)
  #' @export
  a = GenomicRanges::GRanges(GenomicRanges::seqnames(gr), IRanges::IRanges(start(gr), end(gr)), BiocGenerics::strand(gr))
  names(a) = names(gr)
  return(a)
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

  a = x
  if (!is.null(s1)) {
    a = unlist(strsplit(a, split = s1))
    a = a[1:length(a)%%2 == 0]
  }
  if (!is.null(s2)) {
    a = unlist(strsplit(a, split = s2))
  }
  return(a)
}


getnormbed <- function(id, data, samplename, L, same.chromosome = FALSE, p = 0, rmapgr, sk=NULL) {
  #' title Get Normalised Bed From ChicagoData Objects
  #' description Normalises the read counts between samples in one Capture-C experiment by bait-mapped library size. Returns a
  #' GenomicRanges object containing normalized read counts for one bait.
  #' param id baitID
  #' param data chicago data, 'cd(at)x', for a specific sample
  #' param samplename name of the sample, has to be the same as the name of the samples in L
  #' param L list of chicago data tables, 'cd(at)x', for which to calculate scaling factors sk
  #' param same.chromosome should only reads from the same chromosome be taken?; Default: FALSE
  #' param p pseudocount to add before normalization; Default:0
  #' param rmapgr GenomicRanges object containing restriction fragments. IDs must be saved in a column called 'id'.
  #' param sk optional, scaling factors can be provided separately.
  #' return normalised read counts in Genomic Ranges format
  #' @noRd

  if(is.null(sk)){
    sk = getsk(L)
  }
  names(sk) = names(L)
  s = sk[samplename]

  a = data[data$baitID == id, ]
  a = a[order(a$otherEndID), ]
  bed = rmapgr[rmapgr$id %in% a$otherEndID]
  S4Vectors::values(bed) = cbind(as.data.frame(S4Vectors::values(bed)), as.data.frame(a)[, grep("N", names(a))])
  summary(bed$id == a$otherEndID)
  bait = IRanges::resize(rmapgr[rmapgr$id %in% a$baitID], 1, fix = "center")

  if (same.chromosome == TRUE) {
    bed = bed[GenomicRanges::seqnames(bed) == GenomicRanges::seqnames(bait)]
  }
  bed$intID = paste(bait$id, bed$id, sep = ";")

  ## normalize N using the scaling factors:
  bed$N = (bed$N + p) * s

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

  bed = bed[order(GenomicRanges::start(bed), decreasing = FALSE)]

  ## cis interactions: for runmean include also all frags that have 0 interactions:
  b = rmapgr[GenomicRanges::seqnames(rmapgr) == GenomicRanges::seqnames(bait)]
  min = min(bed[GenomicRanges::seqnames(bed) == GenomicRanges::seqnames(bait)]$id)
  max = max(bed[GenomicRanges::seqnames(bed) == GenomicRanges::seqnames(bait)]$id)
  b = b[b$id %in% c(min:max)]
  b$N = 0
  b[b$id %in% bed$id]$N = bed$N

  ## trans interactions separately

  # trans=bed[seqnames(bed) %in% seqnames(bait)==FALSE]
  trans = bed[!as.character(GenomicRanges::seqnames(bed)) %in% as.character(GenomicRanges::seqnames(bait))]
  S4Vectors::values(trans) = S4Vectors::values(trans)[names(S4Vectors::values(b))]

  ## combine cis and trans interactions:
  bed = c(b, trans)

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

  x = (S4Vectors::runmean(S4Vectors::Rle(bed$N), k = k, endrule = "constant"))

  ylab = "% reads per promoter"
  xlab = paste(as.character(GenomicRanges::seqnames(bait)), ":", xlim[1], "-", xlim[2], sep = "")

  bed$x = x
  bed = bed[start(bed) > (xlim[1] - (0.1 * abs(xlim[1] - xlim[2]))) & end(bed) < (xlim[2]) + (0.1 * abs(xlim[1] -
                                                                                                          xlim[2]))]
  return(bed)
}


getbait <- function(id, rmapgr) {
  #' title Get Bait
  #' description Returns a GenomicRanges object with the position of the bait
  #' param id baitID
  #' param rmapgr  GenomicRanges object containing restriction fragments. IDs must be saved in a column called 'id'.
  #' return GenomicRanges object with bait position
  #' @noRd

  bait = IRanges::resize(rmapgr[rmapgr$id %in% id], 1, fix = "center")
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

  ##only executed if baited_genes.rds doesn't exist yet, otherwise loads the existing file:
  if ('baited_genes.rds' %in% list.files(designDir)){
    message('file exists, loading...')
    bg=readRDS(f)
  } else {
    baitmap = paste(designDir, list.files(designDir)[grep("baitmap$", list.files(designDir))], sep = "/")
    baitmap = read.delim(baitmap, header = FALSE, stringsAsFactors = FALSE)
    bg = bed2gr(baitmap)
    bg$re_id = baitmap$V4
    bg$genename = baitmap$V5
    if (save) {
        saveRDS(bg, f)
    }
  }
  return(bg)
}


rmap2rmapgr <- function(designDir = "designDir", save = FALSE) {
  #' @title Convert Rmap to the GenomicRanges object 'rmapgr'
  #' @description Searches for rmap file in designDir and converts to a GenomicRanges object rmapgr.
  #' Only executed if the GenomicRanges object doesn't exist yet, otherwise loads the existing file
  #' @param designDir  Path to the directory with Chicago design files. Default: 'designDir'
  #' @param save If TRUE, saves the converted object as .rds in the designDir,
  #' using the same name as the .rmap file Default: FALSE
  #' @examples
  #' # Example running rmap2rmapgr
  #' extdata <- system.file("extdata", package="PostChicago")
  #' designDir <- file.path(extdata,'designDir')
  #' rmap2rmapgr(designDir,save=FALSE)
  #' @return Returns a GenomicRanges object with the positions of the baits, annotated with fragment IDs (re_id) and gene names
  #' @export

  f <- paste0(designDir, "/", list.files(designDir)[grep("rmap$", list.files(designDir))], ".rds")

  ##only executed if baited_genes.rds doesn't exist yet, otherwise loads the existing file:
  if (paste0(list.files(designDir)[grep("rmap$", list.files(designDir))], ".rds") %in% list.files(designDir)){
    message('file exists, loading...')
    rmapgr=readRDS(f)
  } else {
    rmap = paste(designDir, list.files(designDir)[grep("rmap$", list.files(designDir))], sep = "/")
    rmap = read.delim(rmap, header = FALSE, stringsAsFactors = FALSE)
    rmapgr = bed2gr(rmap)
    rmapgr$id = rmap$V4
    if (save) {
      saveRDS(rmapgr, f)
    }
  }
  return(rmapgr)
}


getbaitplus = function(id, rmapgr) {
  #' title Get bait +/- 1
  #' description returns a gr object with the positons of the bait +/- 1 fragment
  #' param id baitID
  #' param rmapgr GenomicRanges object containing restriction fragments. IDs must be saved in a column called 'id'
  #' return GRanges Object of bait +/- 1 fragment
  #' @noRd


  min = as.numeric(id) - 1
  max = as.numeric(id) + 1
  bait = GenomicRanges::reduce(rmapgr[rmapgr$id %in% min:max])
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

  cov = vector()
  nprom = vector()
  for (i in 1:length(L)) {
    a = L[[i]]
    b = sum(a$N)
    cov = c(cov, b)
    nprom = c(nprom, length(unique(a$baitID)))
  }

  names(cov) = names(L)
  sk = 1/cov * nprom * 100 * 1000
  names(sk) = names(L)

  return(sk)
}


gr2bed <- function(gr) {
  #' @title Genomic Ranges To Bed
  #' @description Converts a Genomic Ranges Object to a BED3-formatted table
  #' @param gr Genomic Ranges Object
  #' @examples
  #' # example running gr2bed
  #' gr <- GenomicRanges::GRanges(1,IRanges::IRanges(1,100))
  #' gr
  #' gr2bed(gr)
  #' @return Bed files with columns names, start, end
  #' @export

  return(data.frame(chr = as.character(GenomicRanges::seqnames(gr)), start = GenomicRanges::start(gr), end = GenomicRanges::end(gr)))
}


getnormscores <- function(id, data, samplename, L, same.chromosome = FALSE, bait, p = 0, rmapgr = rmapgr) {
  ## returns a gr object with chicago scores for a specific bait arguments: id = baitID data=cd@x for a
  ## specific sample samplename = name of the sample, has to be the same as the name of the scaling factors
  ## sk L=list of cd@x objects for which to calculate sk same.chromosome: should only reads from teh same
  ## chromosome be taken? Default:FALSE p:pseudocount to add before normalization, Default:0

  a = data[data$baitID == id, ]
  a = a[order(a$otherEndID), ]
  bed = rmapgr[rmapgr$id %in% a$otherEndID]
  S4Vectors::values(bed) = cbind(S4Vectors::values(bed), data.frame(N = as.data.frame(a)[, grep("score", names(a))]))
  summary(bed$id == a$otherEndID)
  bait = resize(rmapgr[rmapgr$id %in% a$baitID], 1, fix = "center")
  if (same.chromosome == TRUE) {
    bed = bed[GenomicRanges::seqnames(bed) == GenomicRanges::seqnames(bait)]
  }
  bed$intID = paste(bait$id, bed$id, sep = ";")

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
    if ("data.frame" %in% class(df1) & data.table::is.data.table(df1) == FALSE) {
      gr1 = GenomicRanges::GRanges(rep(1, dim(df1)[1]), IRanges::IRanges(df1[, col1], df1[, col1]))
    } else if (data.table::is.data.table(df1)) {
      gr1 = GenomicRanges::GRanges(rep(1, dim(df1)[1]), IRanges::IRanges(df1[[col1]], df1[[col1]]))
    } else {
      gr1 = GenomicRanges::GRanges(rep(1, length(df1)[1]), IRanges::IRanges(S4Vectors::values(df1)[, col1],
                                                                            S4Vectors::values(df1)[, col1]))
    }

    ## get gr2 for overlaps
    if ("data.frame" %in% class(df2) & data.table::is.data.table(df2) == FALSE) {
      gr2 = GenomicRanges::GRanges(rep(1, dim(df2)[1]), IRanges::IRanges(df2[, col2], df2[, col2]))
    } else if (data.table::is.data.table(df2)) {
    } else {
      gr2 = GenomicRanges::GRanges(rep(1, length(df2)[1]), IRanges::IRanges(S4Vectors::values(df2)[, col2],
                                                                            S4Vectors::values(df2)[, col2]))
    }

  }

  if (length(col1) == 3) {

    ## get gr1 for overlaps:
    if ("data.frame" %in% class(df1)) {
      gr1 = GenomicRanges::GRanges(df1[, col1[1]], IRanges::IRanges(df1[, col1[2]], df1[, col1[3]]))
    } else {
      gr1 = GenomicRanges::GRanges(S4Vectors::values(df1)[, col1[1]], IRanges::IRanges(S4Vectors::values(df1)[,
                                                                                                              col1[2]], S4Vectors::values(df1)[, col1[3]]))
    }

    ## get gr2 for overlaps
    if ("data.frame" %in% class(df2)) {
      gr2 = GenomicRanges::GRanges(df2[, col2[1]], IRanges::IRanges(df2[, col2[2]], df2[, col2[3]]))
    } else {
      gr2 = GenomicRanges::GRanges(S4Vectors::values(df2)[, col2[1]], IRanges::IRanges(S4Vectors::values(df2)[,
                                                                                                              col2[2]], S4Vectors::values(df2)[, col2[3]]))
    }

  }

  ov = as.matrix(GenomicRanges::findOverlaps(gr1, gr2))

  return(ov)

}


## annotate interactions with peaks from your grl:

annotateInts <- function(ints, grl, minov = NULL, maxgp = NULL) {
  #' @title Annotate Interactions
  #' @description annotates bait and otherEnd from a table with interactions (ints) with features within a GRangesList (grl)
  #' @param ints a table with interactions, for example derived from Chicago Data tables, for instance by using getInts().
  #' assumes six named columns:
  #' 'seqnames_bait','start_bait','end_bait' (coordinates of the baited restriction fragment) and
  #' 'seqnames_otherEnd','start_otherEnd','end_otherEnd' (coordinates of the interacting restriction fragment)
  #' @param grl GRangesList containing features with which ints should be annotated.
  #' Names of the grl correspond to the names that will be given to table columns
  #' @param minov Default:NULL,  minimal overlap (minoverlap in findOverlaps) between otherEnd or bait fragment and the feature
  #' @param maxgp Default:NULL, maximum gap (maxgap in findOverlaps) between otherEnd or bait fragment and the feature
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

  b = bed2gr(ints[, c("seqnames_bait", "start_bait", "end_bait")])

  for (i in 1:length(grl)) {
    iv = grl[[i]]
    n = names(grl)[i]
    # initialise variable values before
    S4Vectors::values(b)[, n] = FALSE
    if (!is.null(minov)) {
      S4Vectors::values(b[as.matrix(GenomicRanges::findOverlaps(b, iv, minoverlap = minov))[, 1]])[, n] = TRUE
    } else if (!is.null(maxgp)) {
      S4Vectors::values(b[as.matrix(GenomicRanges::findOverlaps(b, iv, maxgap = maxgp))[, 1]])[, n] = TRUE
    } else {
      S4Vectors::values(b[as.matrix(GenomicRanges::findOverlaps(b, iv))[, 1]])[, n] = TRUE
    }
  }

  names(S4Vectors::values(b)) = paste0(names(S4Vectors::values(b)), "_bait")
  ## ensure there is no duplication of the same column names. If ints already contains info on interval
  ## overlaps, they will now be replaced:
  ints = ints[, names(ints) %in% names(S4Vectors::values(b)) == FALSE]
  ints = cbind(ints, as.data.frame(S4Vectors::values(b)))

  b = bed2gr(ints[, c("seqnames_otherEnd", "start_otherEnd", "end_otherEnd")])

  for (i in 1:length(grl)) {
    iv = grl[[i]]
    n = names(grl)[i]
    S4Vectors::values(b)[, n] = FALSE
    if (!is.null(minov)) {
      S4Vectors::values(b[as.matrix(GenomicRanges::findOverlaps(b, iv, minoverlap = minov))[, 1]])[, n] = TRUE
    } else if (!is.null(maxgp)) {
      S4Vectors::values(b[as.matrix(GenomicRanges::findOverlaps(b, iv, maxgap = maxgp))[, 1]])[, n] = TRUE
    } else {
      S4Vectors::values(b[as.matrix(GenomicRanges::findOverlaps(b, iv))[, 1]])[, n] = TRUE
    }
  }

  names(S4Vectors::values(b)) = paste0(names(S4Vectors::values(b)), "_otherEnd")
  ## ensure there is no duplication of the same column names. If ints already contains info on interval
  ## overlaps, they will now be replaced:
  ints = ints[, names(ints) %in% names(S4Vectors::values(b)) == FALSE]
  ints = cbind(ints, as.data.frame(S4Vectors::values(b)))

  return(ints)
}


annotateIntsRepscores <- function(L, ints) {
  #' @title Annotate Interactions with Replicate Scores
  #' @description annotates ints table with replicate scores; only works if scores were determined for each of the replicates
  #' @param L List of replicate cd objects
  #' @param ints interactions table from getInts
  #' @return Table of interactions with annotations
  #' @noRd

  ids = ints$intID
  ints = ints[order(ints$intID), ]

  for (i in 1:length(L)) {
    L[[i]] = L[[i]][L[[i]]$baitID %in% ints$baitID, ]
    L[[i]]$intID = paste(L[[i]]$baitID, L[[i]]$otherEndID, sep = ";")
    L[[i]] = L[[i]][order(L[[i]]$intID), ]
    L[[i]] = L[[i]][L[[i]]$intID %in% ids, ]
    ints[, paste(names(L)[i], "score", sep = "_")] = 0
    ints[ints$intID %in% L[[i]]$intID, paste(names(L)[i], "score", sep = "_")] = L[[i]]$score

  }

  return(ints)
}


## I need to take out troubleshooting later? I actually think it is kind of nice
plotInteractions <- function(L, id, k, zoom=200000, rmapgr, ylim = NULL, show.legend = TRUE, d = NULL, preselected = FALSE,
                             name = NULL, intervals = NULL, show.legend.intervals = TRUE, xlim = NULL, col = NULL, colintervals = NULL,
                             troubleshooting = FALSE, lwd = 2, lty=NULL, show.legend.outside=FALSE, cex.intervals=0.7,
                             cex.legend=1) {
  #' @title Plot Lineplots of Interactions
  #' @description Plots line plots of normalized Capture-C read counts around one view point (bait) from one or many samples.
  #' Optional annotation with genomic intervals of choice and interactions.
  #' @param L List with the chicagoData tables (cd(at)x)
  #' @param id baitID: corresponds to a specific view point.
  #' @param k smoothing factor of the line plot (running average)
  #' @param zoom area around the bait in +/-bp to plot, Default: 200000
  #' @param rmapgr GenomicRanges object containing restriction fragments. IDs must be saved in a column called 'id'
  #' @param ylim optional: limit of the y axis, usually c(0,1), Default c(0, max(all samples))
  #' @param show.legend logical, should the sample legend be shown? Default: TRUE
  #' @param d optional: interactions table to highlight significant interactions as transparent rectangles. Default highlights any interaction.
  #' If different categories of interactions should be highlighted, they should be specified in a column called 'clusters_refined'. Example usage would be lost, gained and stable interactions between
  #' two samples. Up to 8 categories can be defined which then receive a different color each. Default: NULL
  #' @param preselected optional: are interactions in d preselected for the correct baitID?
  #' Useful in cases when the ids in the interactions table do not match the ids in L, Default: FALSE
  #' @param name name of the promoter to plot, optional: Default: name = baitID
  #' @param intervals optional: intervals to be plotted as a named GenomicRangesList grl, useful for annotation of the plot with ChIP-Seq data
  #' @param show.legend.intervals logical input, should the intervals legend be shown? Default: TRUE
  #' @param xlim optional: x limit of the plot as xlim=c(30000,200000), overrides zoom, Default bait +/- zoom
  #' @param col colors of lines in plot, optional: Default: 'Okabe-Ito' with a max of 9 different samples
  #' @param colintervals colors for interval plots, optional: Default: c('grey','pink','purple','green','red','blue','orange','black','darkred','darkgreen','cornflowerblue')
  #' @param troubleshooting if TRUE will return messages showing how far the process is, Default: FALSE
  #' @param lwd line width, Default 2pt
  #' @param lty optional: line type, Default: 1
  #' @param show.legend.outside If TRUE, shows the legend outside of the plot. Default: FALSE
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
  #' #load list with chicago data objects:
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

  #name to use in the plot title:

  if (!is.null(name)) {
    if (length(grep(',',name))>0 || length(grep('\n',name))>0){ ##only add comma before ID if name doesn't contain comma or a line break!
      main = paste0(name, 'id=',id)
    } else {
      main = paste0(name, ',id=',id)
    }
  } else {
    main = paste0('id=',id)
  }

  #optional zoom use in the plot title:

  if(!is.na(zoom)){
    main = paste0(main, ",TSS+/-", zoom/1e+06, "Mb,k=", k)
  } else {
    main = paste0(main, ",k=", k)
  }


  if (!is.null(d)) {
    if ("clusters_refined" %in% names(d) == FALSE) {
      d$clusters_refined = "interaction"
    }
    clust = unique(unlist(strsplit(d$clusters_refined, split = ",")))
  }


  if (is.null(col)) {
    colors = palette.colors(length(L), recycle = TRUE)
  } else {
    colors = col
  }

  if (is.null(lty)) {
    lty = rep(1,length(L))
  }

  if (troubleshooting) {
    message(paste("getting bait..."))
  }
  bait = getbait(id, rmapgr = rmapgr)

  if (troubleshooting) {
    message(paste("bait retrieved..."))
  }

  if (is.null(xlim)) {
    xlim = c(GenomicRanges::start(bait) - zoom, GenomicRanges::end(bait) + zoom)
  }

  if (troubleshooting) {
    message(paste("Bait: \n", print(bait)))
  }



  LL = list()  ##make a list with bed files for all samples
  legend = vector()
  all = vector()

  for (i in 1:length(L)) {
    # for every interaction

    if (troubleshooting) {
      message(paste("sample", i))
    }
    samplename = names(L)[i]
    legend = c(legend, samplename)
    if (troubleshooting) {
      message(paste("getting normbed..."))
    }

    bed = getnormbed(id, data = L[[i]], samplename = names(L)[i], L = L, same.chromosome = TRUE,
                     rmapgr = rmapgr)

    if (troubleshooting) {
      message(paste("normbed retrieved, coverage..."))
      message(str(bed))  # Check structure of bed
    }

    bed = bed2cov(bed, bait, rmapgr = rmapgr)

    if (troubleshooting) {
      message(paste("cov retrieved, getrunmean..."))
      message(str(bed))  # Check structure of bed
      message(paste0("Bed is between ", min(GenomicRanges::start(bed)), " - ", max(GenomicRanges::start(bed))))
      message(paste0("Bait is ", GenomicRanges::start(bait)))
    }


    bed = getrunmeanbed(bed, bait, xlim, k)
    if (troubleshooting) {
      message(paste("runmean retrieved"))
    }

    LL = c(LL, list(bed))
    all = c(all, as.vector(bed$x))

  }
  if (troubleshooting) {
    message(paste("list made"))
  }

  ## make plots:

  # check for graph specifications
  if (is.null(ylim)) {
    ylim = max(LL[[1]]$x, na.rm = TRUE)

    if (length(LL) > 1) {
      for (j in 2:length(LL)) {
        ylim = c(ylim, max(LL[[j]]$x, na.rm = TRUE))
      }
    }

    ylim = as.numeric(c(0, max(ylim, na.rm = TRUE) * 1.2))
  }


  if (!is.null(intervals)) {
    len = length(intervals)
    ylim = c(-ylim[2]/20 * len, ylim[2])
  }

  if (troubleshooting) {
    message(paste("xlab doing, Bait: \n", print(bait)))
  }

  ylab = "normalised read counts (PRPP)"
  xlab = paste(as.character(GenomicRanges::seqnames(bait)), ":", xlim[1], "-", xlim[2], sep = "")

  if (troubleshooting) {
    message(paste("xlab done, Bait: \n", print(bait)))
  }


  bed = LL[[1]]


  plot((GenomicRanges::end(bed) + GenomicRanges::start(bed))/2, bed$x, type = "l", lwd = lwd, xlim = xlim, col = colors[1],
       ylab = ylab, xlab = xlab, main = main, ylim = ylim,lty=lty[1])
  for (j in 2:length(LL)) {
    bed = LL[[j]]
    lines((GenomicRanges::end(bed) + GenomicRanges::start(bed))/2, bed$x, type = "l", xlim = xlim, col = colors[j],
          lwd = lwd,lty=lty[j])
  }

  ##plot bait:
  b = getbaitplus(id, rmapgr = rmapgr)
  # rect(xleft = GenomicRanges::start(b), xright = GenomicRanges::end(b), ybottom = 0, ytop = ylim[2], col = "green",
       # border = NA)
  baitposition=mean(c(  GenomicRanges::start(b), GenomicRanges::end(b) ))
  x=c(baitposition-0.02*(xlim[2]-xlim[1]), baitposition+0.02*(xlim[2]-xlim[1]), baitposition )
  y=c(0, 0, 0.075*max(ylim))
  polygon(x,y, col = "darkgrey", border = NA)

  ## plot signficant interactions

  if (!is.null(d)) {
    if (preselected == FALSE) {
      ib = d[d$baitID == id, ]
    } else {
      ib = d
    }
    ibb = GenomicRanges::GRanges(ib$seqnames_otherEnd, IRanges::IRanges(ib$start_otherEnd, ib$end_otherEnd))
    values(ibb) = ib
    ib = ibb

    cols = list(rgb(0, 0, 0, alpha = 0.1), rgb(1, 0, 0, alpha = 0.1), rgb(0, 1, 0, alpha = 0.1),
                rgb(0, 0, 1, alpha = 0.1),  rgb(1, 1, 0, alpha = 0.1), rgb(1, 0, 1, alpha = 0.1), rgb(0, 1, 1, alpha = 0.1), rgb(0, 0.5, 0.5, alpha = 0.1))
    clust = clust[order(clust)]
    #message('plotting interaction types:', paste(clust))

    for (j in 1:length(clust)) {

      cl = clust[j]
      col = cols[[j]]

      lost = ib[grep(cl, ib$clusters_refined)]

      if (length(lost) > 0) {

        rect(xleft = GenomicRanges::start(lost), xright = GenomicRanges::end(lost), ybottom = -0, ytop = (max(all) +
                                                                                                            1), col = col, border = NA)

      }


    }

    legend("topleft", legend = clust, bty = "n", cex = 0.8, fill = unlist(cols)[1:length(clust)], ncol = 1)

  }

  ## Plot intervals

  if (!is.null(intervals)) {

    abline(h = 0)

    if (is.null(colintervals)) {
      # does this work???  colintervals = palette.colors(length(intervals), recycle = TRUE)
      colintervals = c("grey", "pink", "purple", "green", "red", "blue", "orange", "black", "darkred", "darkgreen",
                       "cornflowerblue")
    }

    grl = intervals
    col = colintervals[1:length(grl)]

    for (i in 1:length(grl)) {
      gr = grl[[i]]
      gr = gr[as.character(GenomicRanges::seqnames(gr)) == as.character(GenomicRanges::seqnames(bait))]
      if (length(gr) > 0) {
        st = 0.05 * ylim[2]
        step = st * i
        if (length(grep("TAD", names(grl)[i])) > 0) {
          rect(xleft = start(gr), xright = end(gr), ybottom = (-st - step), ytop = (-st - step) + step/i, col = col[i])
        } else {
          rect(xleft = start(gr), xright = end(gr), ybottom = (-st - step), ytop = (-st - step) + step/i, col = col[i],
               border = NA)
        }
      }
    }
    if (show.legend.intervals){
      legend("bottomright", fill = col, legend = names(grl), bty = "n", cex = cex.intervals)
    }
  }

  if (show.legend == TRUE) {
    legend("topright", legend = legend, col = colors, lwd = lwd, bty = "n", ncol = 1, cex = cex.legend,lty=lty)
  }

  ##plot the legend outside of the plot, for both: intervals and samples
  if (show.legend.outside){
    plot(0,0,col='white',axes=FALSE,xlab='',ylab='')
    legend('topleft',ncol = 1,col=colors,lwd=lwd,legend=names(L),bty='n', cex = cex.legend,lty=lty)
    if(!is.null(intervals)){
      legend('bottomleft',fill=colintervals,legend=names(grl),bty='n',cex = cex.intervals)
    }
  }
}


makeIntsTable <- function(L, baited_genes, repscores = FALSE, LL = NULL, ngroups = 1, scorecut = 5, readcut = 1,
                           mode = "score", intervals = NULL, outfolder = NULL, printer = TRUE, overwrite = FALSE) {
  #' @title Create interactions table
  #' @description
  #' makeIntsTable() integrates significant interactions from multiple data sets supplied in the list L.
  #' The table is populated by scores, raw and downsampled reads and optionally with replicate data.
  #' @param L list of chicago data objects (cd(at)x) of summarized interactions
  #' @param baited_genes GenomicRanges object, with at least two metadata (or Values) columns: 're_id' (fragment ID covering the bait, corresponding
  #' to fragment numbers indicated in rmapgr), genenames' (names of the baited genes). Can be created automatically from the baitmap table (see Chicago)
  #' using the function baitmap2baited_genes()
  #' @param repscores logical input, should replicate scores be added? if TRUE, LL must be supplied,
  #' Default: FALSE
  #' @param LL optional: list of chicago data objects (cd(at)x) of the replicates of summarized interactions
  #' @param ngroups amount of groups into which data should be split, is ca. 200 entries per group, optional: Default: 1
  #' @param printer Should correlation plots and ints be printed or saved? Default: FALSE, correlation plots are saved.
  #' @param scorecut Chicago score cutoff, Default: 5 (as recommended by Chicago)
  #' @param readcut reads cutoff, optional: Default: 1
  #' @param mode used for cutoff when getting interactions, by 'score', 'read' or 'both', Default: 'score'
  #' @param intervals optional: A GenomicRanges object with intervals with which the otherEnds should overlap.
  #' Useful if mode='reads', to reduce the number of interactions in ints
  #' @param outfolder optional: Path to the directory in which the output is saved, Default: current folder
  #' @param overwrite If TRUE runs the function even if an ints_all.txt is present in the outputdir, otherwise loads ints.txt file if present.
  #' Default: FALSE
  #' @examples
  #' # example code
  #' #define directories
  #' extdata <- system.file("extdata", package="PostChicago")
  #' chicagoOutputDir <- file.path(extdata,'results')
  #' designDir <- file.path(extdata,'designDir')
  #' outputDir <- file.path(extdata,'postchicago')
  #' 
  #' #create baited_genes and rmapgr:
  #' baited_genes=baitmap2baited_genes(designDir,save=FALSE)
  #' rmapgr=rmap2rmapgr(designDir,save=FALSE)
  #' 
  #' #load list with chicago data objects:
  #' cdlist <- loadCdList(resultsDir = chicagoOutputDir, baited_genes = baited_genes)
  #' L=cdlist[-grep('rep',names(cdlist))] ##pooled samples
  #' LL=cdlist[grep('rep',names(cdlist))] ##replicates
  #' 
  #' #run function:
  #' ints = makeIntsTable(L,baited_genes,repscores=TRUE,LL=LL, ngroups = 4, outfolder = outputDir)
  #' head(ints)
  #' 
  #' @returns table of all interactions, which can be saved as .txt and saves a QC plot via correlation heatmap between samples
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

  if (printer == FALSE && is.null(outfolder)) {
    outfolder = getwd()
    message("Output saved here: ", outfolder)
  }

  ## Checking if ints_all.txt is already in the outfolder
  file = "ints_all.txt"

  if (!file %in% dir(outfolder) || overwrite == TRUE) {

    # why is there a nested function??
    f = function(y, scorecut = scorecut, readcut = readcut, mode = "score", intervals = NULL) {

      filtered_L <- lapply(L, function(df) {
        # Check if baitID is in any of the arrays in y
        match_found <- sapply(y, function(arr) any(df$baitID %in% arr))

        # If any match is found in y, filter the rows
        if (any(match_found)) {
          return(df[df$baitID %in% unlist(y[match_found]), ])
        } else {
          return(NULL)  # If no match found, return NULL (or empty data frame)
        }
      })

      # Clean up the list: remove NULL entries (if any)
      filtered_L <- filtered_L[!sapply(filtered_L, is.null)]
      return(getInts(filtered_L, baited_genes[baited_genes$re_id %in% unlist(y)], repscores, LL, scorecut = scorecut,
                     readcut = readcut, mode = mode, intervals = intervals))
    }

    ## split baitIDs into groups of ~100:
    v = vector()
    for (i in 1:length(L)) {
      v = c(v, unique(L[[i]]$baitID))
    }

    # check so that the samples measured are a multiple of the ngroups
    a = split(baited_genes$re_id[baited_genes$re_id %in% v], 1:ngroups)

    if (!is.null(outfolder)) {
      pdf(paste0(outfolder, "/Library_sizes_capturedFrags.pdf"))
      par(mfrow = c(2, length(L)))
    }

    ints = do.call("rbind", BiocGenerics::lapply(a, FUN = f, scorecut = scorecut, readcut = readcut, mode = mode,
                                                 intervals = intervals))
    row.names(ints) = ints$intID

    if (!is.null(outfolder)) {
      dev.off()
    }

    ## Heatmaps, correlations:

    m = ints[, grep("N[.]", names(ints))]
    # m=m[,-grep('downsampled',names(m))]
    m = m[, grep("downsampled", names(m))]
    m = log2(m + 1)
    if (!is.null(outfolder)) {
      pheatmap::pheatmap(cor(m, use = "pairwise.complete.obs"), main = "frags log2(reads+1), cor pearson",
                         filename = paste0(outfolder, "/Frags_Reads_CorPearson_ints.png"))
      pheatmap::pheatmap(cor(m, method = "spearman", use = "pairwise.complete.obs"), main = "frags log2(reads+1), cor spearman",
                         filename = paste0(outfolder, "/Frags_Reads_CorSpearman_ints.png"))
    }
    if (printer == TRUE) {
      pheatmap::pheatmap(cor(m, use = "pairwise.complete.obs"), main = "frags log2(reads+1), cor pearson")
      pheatmap::pheatmap(cor(m, method = "spearman", use = "pairwise.complete.obs"), main = "frags log2(reads+1), cor spearman")
    }

    if (repscores) {
      m = ints[, grep("score", names(ints))]
      m = m[, grep("rep", names(m))]
      m = asinh(m)
      pheatmap::pheatmap(cor(m, use = "pairwise.complete.obs"), main = "frags asinh(scores), cor pearson",
                         filename = paste0(outfolder, "/Frags_Scores_CorPearson_ints.png"))
      pheatmap::pheatmap(cor(m, method = "spearman", use = "pairwise.complete.obs"), main = "frags asinh(scores), cor spearman",
                         filename = paste0(outfolder, "/Frags_Scores_CorSpearman_ints.png"))
    }



    if (!is.null(outfolder)) {
      write.table(ints, paste0(outfolder, "/ints_all.txt"), sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE)
    }

    dev.off()

  } else {
    message("ints_all.txt found in outfolder, loading...")
    ints = read.delim(paste(outfolder, file, sep = "/"), stringsAsFactors = FALSE)
  }

  return(ints)

}


getInts <- function(L, baited_genes, repscores = FALSE, LL = NULL, scorecut = 5, readcut = 1, mode = "score",
                    intervals = NULL) {
  #' @title Create interactions data table
  #' @description
  #' Normalises, sorts, downsamples and saves chicago data of interactions
  #' @param L List with chicago data (cd(at)x)
  #' @param baited_genes Genomic Ranges object or data.frame, with columns 're_id' (bait restrictions fragments),
  #' 'genenames' (names of baited genes)
  #' @param repscores logical input, should replicate scores be added, if TRUE, LL (list of replicate chicago data) is required,
  #' optional: Default: FALSE
  #' @param LL list of chicago data objects (cd(at)x) of replicates of summarized interactions, optional: Default: FALSE
  #' @param scorecut score cutoff, optional: Default: 5
  #' @param readcut reads cutoff, optional: Default: 1
  #' @param mode used for cutoff when getting interactions, by 'score', by 'read' or 'both', optional: Default: 'score'
  #' @param intervals A genomic ranges object with intervals with which the otherEnds should overlap.
  #' Useful if mode='reads', to reduce the number of interactions in ints, optional: Default: NULL
  #' @returns 'int' variable, table of interactions data, saved .png files of quality control: library size both by reps and total
  #' @noRd

  ids = vector()

  for (i in 1:length(L)) {
    # L[[i]]=L[[i]][L[[i]]$baitID %in% baited_genes$re_id,]
    L[[i]]$intID = paste(L[[i]]$baitID, L[[i]]$otherEndID, sep = ";")

    if (mode == "score") {
      ids = unique(c(ids, L[[i]][L[[i]]$intID %in% ids == FALSE & L[[i]]$score >= scorecut, ]$intID))
    } else if (mode == "read") {
      if (is.null(intervals)) {
        ids = unique(c(ids, L[[i]][L[[i]]$intID %in% ids == FALSE & L[[i]]$N >= readcut, ]$intID))
      } else {
        ids = unique(c(ids, L[[i]][L[[i]]$intID %in% ids == FALSE & L[[i]]$N >= readcut & L[[i]]$otherEndID %in%
                                     oeids, ]$intID))
      }
    } else {
      ids = unique(c(ids, L[[i]][L[[i]]$intID %in% ids == FALSE & L[[i]]$N >= readcut & L[[i]]$score >= scorecut,
      ]$intID))
    }

  }

  ids = unique(ids)

  if (length(ids) > 0) {

    ints = data.frame(intID = ids[order(ids)])

    for (i in 1:length(L)) {
      a = L[[i]][L[[i]]$intID %in% ids, ]
      a = a[order(a$intID), ]
      ints[, paste(names(L)[i], "score", sep = "_")] = 0
      ints[, paste(names(L)[i], "N", sep = "_")] = 0

      ## get pooled read numbers and scores
      ints[ints$intID %in% a$intID, paste(names(L)[i], "score", sep = "_")] = a$score
      ints[ints$intID %in% a$intID, paste(names(L)[i], "N", sep = "_")] = a$N

      ## get replicate read numbers
      x = names(a)
      x = x[grep("N", x)]
      x = x[grep("[.]", x)]
      if (length(x) > 0) {
        n = as.numeric(unlist(strsplit(x, split = "[.]"))[length(unlist(strsplit(x, split = "[.]")))])
        for (j in 1:n) {
          ints[, paste0(names(L)[i], "_N.", j)] = 0
          ints[ints$intID %in% a$intID, paste0(names(L)[i], "_N.", j)] = as.data.frame(a)[, paste0("N.",
                                                                                                   j)]

        }
      }

    }


    ints$intID = as.character(ints$intID)
    ints$baitID = as.numeric(unlist(strsplit(ints$intID, split = ";"))[1:length(unlist(strsplit(ints$intID,
                                                                                                split = ";")))%%2 == 1])
    ints$otherEndID = as.numeric(unlist(strsplit(ints$intID, split = ";"))[1:length(unlist(strsplit(ints$intID,
                                                                                                    split = ";")))%%2 == 0])

    oe = gr2bed(rmapgr[findOverlapsDf(ints, rmapgr, "otherEndID", "id")[, 2]])
    b = gr2bed(rmapgr[findOverlapsDf(ints, rmapgr, "baitID", "id")[, 2]])
    names(oe) = paste(names(oe), "otherEnd", sep = "_")
    names(b) = paste(names(b), "bait", sep = "_")
    names(oe)[1] = "seqnames_otherEnd"
    names(b)[1] = "seqnames_bait"
    ints = cbind(ints, oe, b)
    row.names(ints) = ints$intID

    ints = ints[order(ints$otherEndID), ]
    ints = ints[order(ints$baitID), ]

    ## normalize read counts of replicates by the read numbers in the respective libraries (written to
    ## vector v):
    m = ints[, grep("[.]", names(ints))]

    v = vector()
    for (i in 1:length(L)) {
      x = names(L[[i]])
      x = x[grep("N", x)]
      x = x[grep("[.]", x)]
      if (length(x) > 0) {
        n = as.numeric(unlist(strsplit(x, split = "[.]"))[length(unlist(strsplit(x, split = "[.]")))])
        for (j in 1:n) {
          v = c(v, sum(as.data.frame(L[[i]])[, paste0("N.", j)]))
          names(v)[length(v)] = paste0(names(L)[i], "_N.", j)
        }
      }
    }
    libsize = v

    ## TO DO: spearate as function png(paste0(outfolder, '/Library_sizes_capturedFrags.png')) #maybe I
    ## would like this to be external??

    barplot(libsize, las = 3, main = "library sizes captured frags reps")

    ## downsample to smallest:
    if (length(summary(names(m) == names(libsize))[names(summary(names(m) == names(libsize))) == FALSE]) ==
        0) {
      mat = m
      for (i in 1:ncol(m)) {
        m[, i] = m[, i]/libsize[i] * min(libsize)
      }

      names(m) = paste(names(m), "downsampled", sep = "_")
      ints = cbind(ints, m)
    }

    ## normalize total read counts by the read numbers in the respective libraries (written to vecor v):

    m = ints[, grep("_N$", names(ints))]

    v = vector()
    for (i in 1:length(L)) {

      v = c(v, sum(as.data.frame(L[[i]])[, "N"]))
      names(v)[length(v)] = paste0(names(L)[i], "_N")

    }
    libsize = v

    barplot(libsize, main = "library sizes captured frags weighted average read numbers", las = 3)


    ## TO DO: Allow for separate normalization of grouped samples (samplegroups=NULL), such as Med13fl and
    ## Pcgf2fl samples norm to each othe ronly downsample to smallest:
    names(m) == names(libsize)
    mat = m
    for (i in 1:ncol(m)) {
      m[, i] = m[, i]/libsize[i] * min(libsize)
    }

    names(m) = paste(names(m), "downsampled", sep = "_")
    ints = cbind(ints, m)

    ## annotate ints with genenames:
    ints$genename = "NA"
    i = 0
    for (b in unique(ints$baitID)) {
      ## check that baited_genes is valid:
      nonames = c("seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular", "start", "end",
                  "width", "element")
      if (TRUE %in% (names(S4Vectors::values(baited_genes)) %in% nonames)) {
        names(S4Vectors::values(baited_genes)) -> v
        x = as.data.frame(S4Vectors::values(baited_genes))
        x = x[, v %in% nonames == FALSE]
        S4Vectors::values(baited_genes) = x
      }
      n = unique(as.character(as.data.frame(baited_genes)[as.data.frame(baited_genes)$re_id == b, grep("genename",
                                                                                                        names(as.data.frame(baited_genes)))]))

      if (length(n) > 1) {
        nl = n[1]
        for (j in 2:length(n)) {
          nl = paste0(nl, ", ", n[j])
        }
        n = nl
      }
      ints[ints$baitID == b, ]$genename = n
    }

    if (repscores == TRUE) {
      ints = annotateIntsRepscores(LL, ints)
    }

    return(ints)
  }
}


makeOneGeneOnePeak <- function(grl, rmapgr, ints, maxgap = 1000, folder = 'oneGeneOnePeak') {
  #' @title Make one-Gene-one-Peak Interaction Table
  #' @description This function takes intervals in grl, gives them a unique ID and identifies all interaction with a given interval.
  #' @param grl GRangesList containing features with which ints should be annotated.
  #' Names of the grl correspond to the names that will be given to table columns
  #' @param rmapgr GenomicRanges object containing restriction fragments. IDs must be saved in a column called 'id'
  #' @param ints chicago data over interactions (as returned by getInts)
  #' @param maxgap minimum gap for two different peaks/maximum separation for two maxima to be binned together, optional: Default 1000
  #' @param folder Output folder
  #' @examples
  #' # example code
  #' 
  #' #define directories:
  #' extdata <- system.file("extdata", package="PostChicago")
  #' outputDir <- file.path(extdata,'postchicago')
  #' designDir <- file.path(extdata,'designDir')
  #' intDir <- file.path(extdata,'intervals')
  #' 
  #' #read in the interactions file
  #' ints=read.delim(paste0(outputDir,'/ints_all.txt'))
  #' 
  #' #construct rmapgr
  #' rmapgr=rmap2rmapgr(designDir,save=FALSE)
  #' 
  #' #construct grl
  #' grl <- suppressWarnings(loadGrl(intDir = intDir))
  #' 
  #' #create oneGeneOnePeak tables
  #' quantfolder=paste0(outputDir,'/oneGeneOnePeak')
  #' makeOneGeneOnePeak(grl = grl,rmapgr = rmapgr, ints = ints, folder = quantfolder, maxgap=1000)
  #' 
  #' @returns annotated OneGeneOnePeak interaction table as .txt.
  #' @export

  if (!file.exists(folder)) {
      dir.create(folder)
  }

  for (i in 1:length(grl)) {

    file1=paste0(folder, "/ints_", names(grl)[i], "_interacting.txt")
    file2=paste0(folder, "/", names(grl)[i], "_interacting_rmapFrag_annotated.txt")
    file3=paste0(folder, "/oneGeneOnePeak_", names(grl)[i], "_interacting.txt")

    ##Only generates OneGeneOnePeak files if they do not exist yet
    if (!file.exists(file3)){

        ## for each k27ac peak: determine the restriction fragment its center overlaps:
        b = IRanges::resize(grl[[i]], 1, fix = "center")
        ov = as.matrix(GenomicRanges::findOverlaps(b, rmapgr))
        y = b[ov[, 1], ]
        y$re_id = rmapgr[ov[, 2]]$id

        ## for each interaction determine the name of the k27ac peak the gene interacts
        b = grl[[i]]
        b = b[as.matrix(findOverlaps(b,y))[,1],]
        ix = bed2gr(ints[, c("seqnames_otherEnd", "start_otherEnd", "end_otherEnd")])

        ov = as.matrix(GenomicRanges::findOverlaps(ix, b, maxgap = maxgap))
        ints.k = ints[ov[, 1], ]
        ints.k[, paste0("re_id_", names(grl)[i])] = y[ov[, 2], ]$re_id
        #ints.k = cbind(ints.k, gr2bed(grl[[i]][ov[, 2]]))  ##AF (20.11.2025): Commented out. grl[[i]] may not be the right one anymore here, the correct one is b
        ints.k = cbind(ints.k, gr2bed(b[ov[, 2]]))
        names(ints.k)[(ncol(ints.k) - 2):ncol(ints.k)] = paste(names(ints.k)[(ncol(ints.k) - 2):ncol(ints.k)], names(grl)[i],
                                                               sep = "_")

        write.table(ints.k, file1, sep = "\t", eol = "\n",
                    quote = FALSE, row.names = FALSE)
        write.table(y, file2, sep = "\t", eol = "\n",
                    quote = FALSE, row.names = FALSE)

        ## make a table with one baitID one peak interactions. Scores and reads won't be quantified properly
        ## yet!
        peaks = read.delim(file1, stringsAsFactors = FALSE)

        ## make data frame d with clusters_refined and other end data encompassing the peak summit!
        d = peaks
        d$otherEndID = d[, paste0("re_id_", names(grl)[i])]
        d[, c("seqnames_otherEnd", "start_otherEnd", "end_otherEnd")] = d[, c(paste0("chr_", names(grl)[i]), paste0("start_",
                                                                                                                    names(grl)[i]), paste0("end_", names(grl)[i]))]
        dim(d)
        d$intID = paste(d$baitID, d$otherEndID, sep = ";")

        id = unique(d$intID)[1]
        a = d[d$intID == id, ]
        df = a[1, ]
        for (id in unique(d$intID)[2:length(unique(d$intID))]) {
          df = rbind(df, d[d$intID == id, ][1, ])
        }

        dim(df)
        d = df
        row.names(d) = d$intID

        gr = bed2gr(d[, paste(c("chr", "start", "end"), names(grl)[i], sep = "_")])
        for (n in names(grl)) {
          ov = as.matrix(GenomicRanges::findOverlaps(gr, grl[[n]]))
          d[, paste(n, "otherEnd", sep = "_")] = FALSE
          d[ov[, 1], paste(n, "otherEnd", sep = "_")] = TRUE
        }

        # write.table(d,paste0(folder,
        # '/oneGeneOnePeak_',names(grl)[i],'_maxgap.',maxgap,'bp_interacting.txt'),sep='\t',eol='\n',quote=FALSE)
        write.table(d, file3, sep = "\t", eol = "\n",
                    quote = FALSE)

    } else {

       message(paste('OneGeneOnePeak file exists for', names(grl)[i], ', skipping...'))

        }

    }

}


###########################
##Quantification in peaks:
###########################

## a)Function to identify fragments flanking peak summit that also overlap the peaks


getPeakFrags <- function(i, infolder = '', grl, rmapgr) {
  ##Helper function for getDataInPeaks()

  d = read.delim(paste0(infolder, "/oneGeneOnePeak_", names(grl)[i], "_interacting.txt"), stringsAsFactors = FALSE)

  ## only take peaks that are on the same chromosome:
  d = d[d$seqnames_bait == d$seqnames_otherEnd, ]

  a = read.delim(paste0(infolder, "/", paste(names(grl)[i], "interacting_rmapFrag_annotated.txt", sep = "_")),
                 stringsAsFactors = FALSE)
  aa = bed2gr(a)
  aa$id = a$re_id
  a = aa
  rm(aa)
  gc()

  ov = as.matrix(findOverlaps(grl[[i]], a))
  aa = grl[[i]][ov[, 1]]
  aa$id = a[ov[, 2]]$id
  a = aa
  rm(aa)
  gc()

  a = a[a$id %in% d[, paste("re_id", names(grl)[i], sep = "_")], ]

  a = unlist(GenomicRanges::GRangesList(lapply(split(a, f = 1:length(a)), FUN = identifyFragSteps, rmapgr = rmapgr)))
  saveRDS(a, paste0(infolder, "/", paste(names(grl)[i], "intInData_rmapFrag_annotated.rds", sep = "_")))  ##TO DO: This saves sometimes two entries per row!
}



identifyFragSteps <- function(x, rmapgr) {
  ##Helper function for getPeakFrags()
  ##Formerly nested as f and used in

  y = rmapgr[as.matrix(findOverlaps(x, rmapgr))[, 2]]$id
  x$min = min(y) - x$id
  x$max = max(y) - x$id
  return(x)
}


## b)Function to calculate reads/scores within all peaks

calculateUnderAllPeaks <- function(i, what = "reads", infolder, resfolder, grl, zoom) {
  ##Helper function for getDataInPeaks
  ##Formerly nested

  message(paste("calculating", what, "in", names(grl)[i], "interactions..."))

  ## load interaction data:
  d = read.delim(paste0(infolder, "/oneGeneOnePeak_", names(grl)[i], "_interacting.txt"), stringsAsFactors = FALSE)

  ## only take peaks that are on the same chromosome:
  d = d[d$seqnames_bait == d$seqnames_otherEnd, ]
  d = d[order(d$intID), ]

  ## load info about where the fragments overlapping the peak are relative to the summit
  a = readRDS(paste0(infolder, "/", paste(names(grl)[i], "intInData_rmapFrag_annotated.rds", sep = "_")))
  load(paste0(resfolder, "matrices_", names(grl)[i], '_',what))
  m=Lm

  ## bring matrices (m) and d in the same order
  m = lapply(m, FUN = function(x) x = x[order(row.names(x)), ])

  df = do.call(rbind, lapply(1:nrow(d), FUN = calculateUnderOnePeak, what = what, m = m, grl = grl, i = i, zoom = zoom, d = d, a = a))
  dat = cbind(d, df)

  message("finished calculations...")

  return(dat)
}


## c)Function to calculate reads/scores within one peak

## TO DO: Check, here is the problem with x=85 (K27ac peaks) - turns out that there are two such
## peaks!
calculateUnderOnePeak <- function(x, what, m, grl, i, zoom, d, a) {
  ##Helper function for calculateUnderAllPeaks
  ##formerly nested

  df = data.frame(matrix(nrow = 1, ncol = length(m) * 3))
  colnames(df) = paste(rep(names(m), each = 3), c(paste0("sum_bednorm_", what), paste0("mean_bednorm_",
                                                                                       what), paste0("summit_bednorm_", what)), sep = ".")
  id = d[x, paste("re_id", names(grl)[i], sep = "_")]
  # ax=a[a$id==id ] ##TO DO: Here is the problem: in x=85 (K27ac peaks) we get two entries, should
  # be one!  ##AF 22.07.2025: Fix: Select the one ax object whose start site overlaps the interval
  # peak site ax=a[a$id==id & start(a)==d[x,'start_otherEnd']] Fix: changed AF 20.08.2025, id can be
  # misleading if we have two peaks that overlap the same rmapgr_id: ax=a[a$id==id]
  ax = a[seqnames(a) == d[x, paste0("chr_", names(grl[i]))] & start(a) == d[x, paste0("start_", names(grl[i]))] &
           end(a) == d[x, paste0("end_", names(grl[i]))]]

  ## Matrices are flipped if oE>bait id! That means that now we also need to flip a for interactions
  ## with baits with a smaller re_id
  ##potentially change: reverse flipping of controls? -> probably not necessary, since distance-matched
  if (d[x, paste("re_id", names(grl)[i], sep = "_")] > d[x, "baitID"]) {
    names(values(ax)) = c("id", "max", "min")
    ax$max = -ax$max
    ax$min = -ax$min
  }
  for (j in 1:length(m)) {
    if ((zoom + 1 + as.data.frame(ax)$min) > 0 & (zoom + 1 + as.data.frame(ax)$max) <= (zoom * 2 + 1)) {
      df[, paste(names(m)[j], paste0("sum_bednorm_", what), sep = ".")] = sum(m[[j]][x, ][(zoom + 1 +
                                                                                             as.data.frame(ax)$min):(zoom + 1 + as.data.frame(ax)$max)])
      df[, paste(names(m)[j], paste0("mean_bednorm_", what), sep = ".")] = mean(m[[j]][x, ][(zoom +
                                                                                               1 + as.data.frame(ax)$min):(zoom + 1 + as.data.frame(ax)$max)])
      df[, paste(names(m)[j], paste0("summit_bednorm_", what), sep = ".")] = m[[j]][x, ][(zoom + 1)]
    } else {

      if ((zoom + 1 + as.data.frame(ax)$min) < 0) {
        mn = 1
      }
      if ((zoom + 1 + as.data.frame(ax)$max) > (zoom * 2 + 1)) {
        ma = (zoom * 2 + 1)
      }

      df[, paste(names(m)[j], paste0("sum_bednorm_", what), sep = ".")] = sum(m[[j]][x, ])
      df[, paste(names(m)[j], paste0("mean_bednorm_", what), sep = ".")] = mean(m[[j]][x, ])
      df[, paste(names(m)[j], paste0("summit_bednorm_", what), sep = ".")] = m[[j]][x, ][(zoom + 1)]

    }
  }
  return(df)

  ##-> turns out some peaks are larger than the assessed matrix size!

}


##AF 27.08.2025: Changed this function to only overwrite ReadsScoresInPeaks.txt if at least one is not existing!

getDataInPeaks <- function(grl, resfolder = "matrices/", overwrite = FALSE, infolder = '', zoom) {
  #' @title Quantify reads and scores in interaction peaks
  #' @description
  #' Quantifies reads and scores under interaction peaks, which can be defined
  #' either as interactions with some features (i.e. Ring1b peaks) or as interaction peak empirically determined.
  #' Reads and scores are quantified in 3 ways: 1) sum value of all restriction fragments under the peak, 2) mean
  #' value of all restriction fragments under a peak and 3) value in the peak summit.
  #' @param grl Genomic Ranges List containing intervals for which interaction frequencies should be quantified. Typically a list of ChIP-Seq peaks.
  #' @param resfolder Path to the directory in which matrices are stored. Default: 'matrices/'. Creates the directory if necessary. Note that the names of the matrix lists are extracted from the names of intervals in grl. Therefore, if creating matrices from Capture-C replicates, it is recommended to use a separate resfolder, such as 'matrices_reps/'
  #' @param overwrite Should helper tables be overwritten? Refers to the helper tables that detail under which restriction fragments the quantification should occur. Default: FALSE.
  #' @param infolder Path to the directory in which oneGeneOnePeak files are stored. Default: Current working directory.
  #' @param zoom Zoom in restriction fragments, for which the matrices were quantified by getMatrix()
  #' @examples
  #' # example code
  #' 
  #' #requires precomputed matrices in the resfolder (see getMatrix())
  #' extdata <- system.file("extdata", package="PostChicago")
  #' outputDir <- file.path(extdata,'postchicago')
  #' quantfolder=paste0(outputDir,'/oneGeneOnePeak')
  #' resfolder=paste0(outputDir,'/matrices/')
  #' intDir <- file.path(extdata,'intervals')
  #' 
  #' #construct grl
  #' grl <- suppressWarnings(loadGrl(intDir = intDir))
  #' 
  #' ##run function
  #' dip=getDataInPeaks(grl[grep('K27ac',names(grl))],resfolder=resfolder, 
  #'     infolder = quantfolder, overwrite=FALSE, zoom = 100)
  #' 
  #' @returns Saved annotated OneGeneOnePeak tables. Extension from resfolder is used to save the table (e.g. if resfolder='matrices_reps/', the extension '_reps' is used).
  #' @export
  ## Requires: 1) premade oneGeneOnePeak interaction tables and rmapgr frag annotated interval tables (both
  ## are created when running makeOneGeneOnePeak(grl,rmapgr,ints)) 2) interaction matrices 3) rmapgr,
  ## an rds containing ids and coordinates of restriction fragments The function uses and annotates
  ## oneGeneOnePeak tables and matrices previously generated using the makeOneGeneOnePeak() and the
  ## getMatrix() functions.  arguments: grl: GenomicRangesList object containing intervals for which the
  ## calculations should be made zoom: distance in restriction fragments from the summit of intervals for
  ## which the matrices have been calculated resfolder: folder containing the interaction matrices, Default:
  ## 'matrices/' overwrite: should existing tables detailing which intervals interact with promoters in teh
  ## dataset be overwritten? Default: FALSE maxgap: What is the maximum distance to an interval from grl
  ## that is acceptable? Must be the same as defined in makeOneGeneOnePeak. Default: 1000 infolder:
  ## typical workflow: makeOneGeneOnePeak(grl,rmapgr,ints)
  ## getMatrix(L,zoom,d,resfolder,type='reads',norm=FALSE,name=names(grl)[i])
  ## getMatrix(L,zoom,d,resfolder,type='scores',norm=FALSE,name=names(grl)[i]) getDataInPeaks(grl,zoom)

  # if (is.null(infolder)) {
  #   infolder = getwd()
  # }

  fileprefix=paste0(infolder, "/quantReadsScoresInPeak")

  ## a)Identify all fragments that overlap the peak summit and the rest of the peak

  ## Check if rmapgrFrag_annotated tables are present in the infolder, load them if overwrite==FALSE, otherwise create de novo:

  if (overwrite == FALSE) {
    v = vector()
    for (i in 1:length(grl)) {
      v = c(v, paste(names(grl)[i], "intInData_rmapFrag_annotated.rds", sep = "_") %in% list.files(infolder))
    }
    vv = data.frame(1:length(v), v)
  } else {
    v = rep(FALSE, length(grl))
    vv = data.frame(1:length(v), v)
  }

  if (FALSE %in% v) {
    lapply(vv[v == FALSE, 1], FUN = getPeakFrags, infolder = infolder, grl = grl, rmapgr = rmapgr)
  }

  message("getPeakFrags finished...")

  ##Addition: AF 27.08.2025
  ## optional: only create the readsScoresInPeak.txt tables if they have not already been created:
  if (overwrite == FALSE) {
    v = vector()
    for (i in 1:length(grl)) {
      v = c(v, paste0(fileprefix, '_', names(grl)[i], ".txt") %in% list.files(infolder))
    }
    vv = data.frame(1:length(v), v)
  } else {
    v = rep(FALSE, length(grl))
    vv = data.frame(1:length(v), v)
  }

  if (FALSE %in% v) {
      ##Previous code was within the brackets above (before 27.08.2025 AF):

      newlist = lapply(1:length(grl), FUN = calculateUnderAllPeaks, what = "reads", infolder = infolder, resfolder = resfolder, grl = grl, zoom = zoom)
      newlist2 = lapply(1:length(grl), FUN = calculateUnderAllPeaks, what = "scores", infolder = infolder, resfolder = resfolder, grl = grl, zoom = zoom)

      message(print("combining reads and scores..."))

      newlist3 = newlist
      for (k in 1:length(newlist)) {
        d = read.delim(paste0(infolder, "/oneGeneOnePeak_", names(grl)[k], "_interacting.txt"), stringsAsFactors = FALSE)
        newlist3[[k]] = cbind(newlist[[k]], newlist2[[k]][, (ncol(d) + 1):ncol(newlist2[[k]])])
      }
      names(newlist3) = names(grl)

      message(print("writing tables..."))

      for (i in 1:length(newlist3)) {
        newlist3[[i]] = annotateInts(newlist3[[i]], grl)
        write.table(newlist3[[i]], paste0(fileprefix, '_', names(newlist3)[i], ".txt"), sep = "\t", eol = "\n", quote = FALSE,
                    row.names = FALSE)
      }

    } else {

      newlist3=list()
      for (i in 1:length(grl)) {
        newlist3=c(newlist3,list(
          read.delim(paste0(fileprefix, '_', names(newlist3)[i], ".txt"),
                     stringsAsFactors=FALSE)
        ))

        # newlist3[[i]] = annotateInts(newlist3[[i]], grl)
        # write.table(newlist3[[i]], paste0(infolder, "/oneGeneOnePeak_", names(newlist3)[i], strsplit_string(resfolder,
        #                                                                                                     s1 = "matrices", s2 = "/"), "_interacting_readsScoresInPeak.txt"), sep = "\t", eol = "\n", quote = FALSE,
        #             row.names = FALSE)
      }

    }

    return(newlist3)

}




## function to plot aggregate peaks from CaptureC analysis matrices:

plotAggregatePeaks <- function(Lm, mainprefix = NULL, col = NULL, ylim = NULL, lty = NULL, xlim = NULL, k = 3,
                               ylab = "normalised reads",plotwhat = "all") {
  #' @title Plot Aggregate Interaction Profiles
  #' @description Takes a list of matrices and plots aggregate profiles
  #' @param Lm List of matrices containing interaction reads or scores
  #' @param mainprefix Optional; text added to the plot title.
  #' @param col Optional; colour for the aggregate plots. Default contains four colours.
  #' @param lty Optional; line type for the aggregate plots
  #' @param ylim Optional; Default data-driven
  #' @param xlim Optional; should be added in restriction fragments surrounding the interval/peak center. Example: c(-50,50).
  #' Default: All restriction fragments represented in the matrix.
  #' @param k Numeric; over how many fargments should the running mean be calculated? Default: 3
  #' @param ylab Optional. Default: 'normalised reads'
  #' @param plotwhat What plots to plot? Options are 'all','mean','median'
  #' 'mean' and 'median' plot running averages over k fragments.
  #' Default: 'all': Four plots representing mean and median read counts both smoothed and unsmoothed.
  #' @examples
  #' # example code
  #' 
  #' #requires precomputed matrices in the resfolder (see getMatrix())
  #' extdata <- system.file("extdata", package="PostChicago")
  #' outputDir <- file.path(extdata,'postchicago')
  #' resfolder=paste0(outputDir,'/matrices/')
  #' 
  #' #load list of matrices supplied in the PostChicago package (Lm):
  #' load(paste0(resfolder,'matrices_K27ac_peaks_reads'))
  #' 
  #' #plot aggregate profiles:
  #' plotAggregatePeaks(Lm,plotwhat='median',k=3,xlim=c(-50,50),mainprefix= 'K27ac\n', ylab='reads', 
  #'                     col=rep(c('red','blue'),2), lty=c(1,1,3,3))
  #' 
  #' @return Plot containing aggregate matric profiles
  #' @export
  ## arguments: optional: xlim: Region (in restriction fragments) surrounding the peak summit to plot.
  ## Default: All restriction fragments represented in the matrix. plotwhat: What plots to plot? Options are
  ## 'all','mean','median' Default: four plots representing mean and median read counts both smoothed and
  ## unsmoothed.  Mean and median plot running averages over k fragments.  which: which of the plots to
  ## plot?

  ## removed: AF 25.08.2025 if (is.null(ylim)){ ylim=c(0,1.5) }
  if (is.null(col)) {
    col = rep(c("blue", "red", "green", "grey"), 2)
  }
  if (is.null(lty)) {
    lty = rep(c(1, 2), each = 4)
  }

  ## columns for the legend (up to 6 samples per column)
  ncol = 1
  if (length(Lm) > 8) {
    ncol = ceiling(length(Lm)/8)
  }

  if (plotwhat == "all") {

    for (type in c("mean", "median")) {

      if (type == "mean") {
        main = paste(mainprefix, "mean enrichment", sep = "")
      } else {
        main = paste(mainprefix, "median enrichment", sep = "")
      }

      ## aggregate plots:
      Lx = list()
      for (i in 1:length(Lm)) {
        if (type == "mean") {
          Lx = c(Lx, list(colMeans(Lm[[i]])))
        } else {
          Lx = c(Lx, list(colMedians(Lm[[i]])))
        }
      }
      names(Lx) = names(Lm)

      Lxx = Lx
      for (i in 1:length(Lxx)) {
        Lxx[[i]] = runmean(Rle(Lxx[[i]]), k = k, endrule = "constant")
      }

      ## plot:

      dis = (length(Lx[[1]]) - 1)/2
      x = -dis:dis

      if (!is.null(xlim)) {

        for (i in 1:length(Lx)) {
          a = Lx[[i]]
          df = data.frame(a, x)
          Lx[[i]] = df[df$x %in% xlim[1]:xlim[2], ]$a

          a = Lxx[[i]]
          df = data.frame(a, x)
          Lxx[[i]] = df[df$x %in% xlim[1]:xlim[2], ]$a
        }
        x = xlim[1]:xlim[2]
      }

      ## Added AF 26.08.2025: Data-driven ylim:
      ylim = c(0, max(unlist(Lx)) * 1.1)


      plot(x, Lx[[1]], type = "l", ylim = ylim, lwd = 2, ylab = ylab, xlab = "distance from summit [RE frags]",
           col = col[1], main = main)
      for (i in 1:length(Lx)) {
        lines(x, Lx[[i]], lwd = 2, col = col[i], lty = lty[i])
      }
      legend("topright", cex = 0.8, bty = "n", lwd = 2, col = col, legend = names(Lm), lty = lty, ncol = ncol)
      legend("topleft", bty = "n", cex = 0.8,
             legend = c(paste("n_ints=", nrow(Lm[[1]]), sep = ""),
                        paste0("n_prom=", length(unique(strsplit_string(row.names(Lm[[1]]), s2 = ";")[1:length(strsplit_string(row.names(Lm[[1]]),
                              s2 = ";"))%%2 == 1])))))



      plot(x, Lxx[[1]], type = "l", ylim = ylim, lwd = 2, ylab = ylab, xlab = "distance from summit [RE frags]",
           col = col[1], main = paste0(main, " k=", k))
      for (i in 1:length(Lxx)) {
        lines(x, Lxx[[i]], lwd = 2, col = col[i], lty = lty[i])
      }
      legend("topright", cex = 0.7, bty = "n", lwd = 2, col = col, legend = names(Lm), lty = lty, ncol = ncol)
      legend("topleft", bty = "n", cex = 0.8,
             legend = c(paste("n_ints=", nrow(Lm[[1]]), sep = ""),
                        paste0("n_prom=", length(unique(strsplit_string(row.names(Lm[[1]]), s2 = ";")[1:length(strsplit_string(row.names(Lm[[1]]),
                         s2 = ";"))%%2 == 1])))))
    }

  } else {


    type = plotwhat

    if (type == "mean") {
      main = paste(mainprefix, "mean enrichment", sep = "")
    } else {
      main = paste(mainprefix, "median enrichment", sep = "")
    }

    ## aggregate plots:
    Lx = list()
    for (i in 1:length(Lm)) {
      if (type == "mean") {
        Lx = c(Lx, list(colMeans(Lm[[i]])))
      } else {
        Lx = c(Lx, list(colMedians(Lm[[i]])))
      }
    }
    names(Lx) = names(Lm)

    Lxx = Lx
    for (i in 1:length(Lxx)) {
      Lxx[[i]] = runmean(Rle(Lxx[[i]]), k = k, endrule = "constant")
    }

    ## plot:

    dis = (length(Lx[[1]]) - 1)/2
    x = -dis:dis

    if (!is.null(xlim)) {

      for (i in 1:length(Lx)) {
        a = Lx[[i]]
        df = data.frame(a, x)
        Lx[[i]] = df[df$x %in% xlim[1]:xlim[2], ]$a

        a = Lxx[[i]]
        df = data.frame(a, x)
        Lxx[[i]] = df[df$x %in% xlim[1]:xlim[2], ]$a
      }
      x = xlim[1]:xlim[2]
    }

    ## Added AF 26.08.2025: Data-driven ylim:
    ylim = c(0, max(unlist(Lx)) * 1.1)

    plot(x, Lxx[[1]], type = "l", ylim = ylim, lwd = 2, ylab = ylab, xlab = "distance from summit [RE frags]",
         col = col[1], main = paste0(main, " k=", k))
    for (i in 1:length(Lxx)) {
      lines(x, Lxx[[i]], lwd = 2, col = col[i], lty = lty[i])
    }
    legend("topright", cex = 0.7, bty = "n", lwd = 2, col = col, legend = names(Lm), lty = lty, ncol = ncol)
    legend("topleft", bty = "n", cex = 0.8,
           legend = c(paste("n_ints=", nrow(Lm[[1]]), sep = ""),
                      paste0("n_prom=", length(unique(strsplit_string(row.names(Lm[[1]]), s2 = ";")[1:length(strsplit_string(row.names(Lm[[1]]),
                                                                                                s2 = ";"))%%2 == 1])))))


  }
}


##AF 23.11.2025: Changed grl to NULL
boxplotsCC <- function(grl, resfolder, infolder = '', l = NULL, col=NULL) {
  #' @title Plot Quantification of OneGeneOnePeak Interactions
  #' @description
  #' Plots scores and reads of each experiment
  #' @param grl GRangesList containing features that were quantified for the boxplot.
  #' @param resfolder Folder where interaction matrices are saved
  #' @param infolder Folder where the tables with feature interaction quantifications are stored.
  #' @param l Optional, list with readsScoresInPeak tables from which to make plots. Default: 
  #' No list supplied, then the information from grl and saved tables generated by getDataInPeaks() is used for the plots.
  #' @param col Optional color
  #' @examples
  #' # example code
  #' 
  #' #requires precomputed matrices in the resfolder (see getMatrix())
  #' extdata <- system.file("extdata", package="PostChicago")
  #' outputDir <- file.path(extdata,'postchicago')
  #' quantfolder=paste0(outputDir,'/oneGeneOnePeak')
  #' resfolder=paste0(outputDir,'/matrices/')
  #' intDir <- file.path(extdata,'intervals')
  #' 
  #' #Load intervals for which the boxplots should be quantified into a GRangesList 
  #' grl <- suppressWarnings(loadGrl(intDir = intDir))
  #' 
  #' #Quantify interactions under the peaks in grl to generate quantification tables:
  #' #dip=getDataInPeaks(grl[grep('K27ac',names(grl))],resfolder=resfolder, 
  #' #infolder = quantfolder, overwrite=FALSE, zoom = 100)
  #' 
  #' #Read in interaction quantification tables into a list where each table represents the quantification
  #' #over one set of intercals in grl:
  #' df=read.delim(paste0(quantfolder, '/quantReadsScoresInPeak_K27ac_peaks.txt'),stringsAsFactors=FALSE)
  #' quantlist=list(K27ac_peaks=df)
  #' 
  #' #Plot boxplots of these quantifications
  #' boxplotsCC(grl[grep('K27ac',names(grl))],resfolder=resfolder, infolder=outputDir, l=quantlist)
  #' 
  #' @returns boxplots of scores and reads of each experiment of a specific feature
  #' @export

  if (is.null(col)){
    col = c(RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))
  }

  if (is.null(l)) {
    l = list()
    for (i in 1:length(grl)) {
      l = c(l, list(read.delim(paste0(infolder, '/quantReadsScoresInPeak_',names(grl)[i],'.txt'), stringsAsFactors = FALSE)))
    }
    names(l) = names(grl)
  }

  for (i in 1:length(l)) {
    par(mfrow = c(2, 3))

    nints = nrow(l[[i]])
    nproms = length(unique(l[[i]]$baitID))

    for (what in c("summit_bednorm_scores", "sum_bednorm_scores", "mean_bednorm_scores")) {

      a = l[[i]][, grep(what, names(l[[i]]))]
      names(a) = strsplit_string(names(a), s2 = paste0("[.]", what))

      ylim = c(0, max(asinh(a),na.rm=TRUE))
      boxplot(asinh(a[, c(1:(ncol(a)/2), ((ncol(a)/2 + 1):ncol(a)))] + 1), las = 3, outline = FALSE,
              ylab = paste("asinh", what), main = names(l)[i], ylim = ylim, col = rep(col[1:(ncol(a)/2)], 2), notch = TRUE,
              border = rep(c("black", "lightgrey"), each = ncol(a)/2))

      legend("topright", legend = c(paste("n_ints =", nints), paste("n_proms =", nproms)), bty = "n")
    }

    for (what in c("summit_bednorm_reads", "sum_bednorm_reads", "mean_bednorm_reads")) {
      a = l[[i]][, grep(what, names(l[[i]]))]
      names(a) = strsplit_string(names(a), s2 = paste0("[.]", what))

      ylim = c(0, max(log2(a + 1),na.rm=TRUE))
      boxplot(log2(a[, c(1:(ncol(a)/2), ((ncol(a)/2 + 1):ncol(a)))] + 1), las = 3, outline = FALSE, ylab = paste("log2", what, "+1"),
              main = names(l)[i], ylim = ylim, col = rep(col[1:(ncol(a)/2)], 2), notch = TRUE, border = rep(c("black", "lightgrey"), each = ncol(a)/2))


      legend("topright", legend = c(paste("n_ints =", nints), paste("n_proms =", nproms)), bty = "n")

    }

  }

}


makeQCplot <- function(ints, folder = NULL, plot = FALSE) {
  #' @title Quality Control Plot
  #' @description
  #' Will create and save heatmaps of downsampled and raw data from the interactions matrix by read
  #' and the scores of each replicate
  #' @param ints Interaction Matrix
  #' @param folder Output folder; optional: Default current folder
  #' @param plot Should the plots be printed? If 'TRUE', then the QC plots are saved and printed, otherwise only saved. Default: FALSE
  #' @examples 
  #' # example code:
  #' extdata <- system.file("extdata", package="PostChicago")
  #' outputDir <- file.path(extdata,'postchicago')
  #' 
  #' #read in interactions:
  #' ints=read.delim(paste0(outputDir,'/ints_all.txt'),stringsAsFactors=FALSE)
  #' 
  #' #assess interactions:
  #' makeQCplot(ints, outputDir, plot = TRUE)
  #' 
  #' @returns png figures of association heatmaps of raw and downsampled data as well as the scores by replicate
  #' @export

  if (is.null(folder)) {
    folder = getwd()
  }

  m = ints[, grep("N[.]", names(ints))]
  mraw = log2(m[, -grep("downsampled", names(m))] + 1)
  mnorm = log2(m[, grep("downsampled", names(m))] + 1)
  if (plot) {
    plt <- pheatmap::pheatmap(mraw, show_rownames = F, main = "log2(raw_reads+1)",border_color=NA)
    show(plt)
    ggplot2::ggsave(filename = paste0(folder, "/pheatmap_rawReads_byInt.png"), plot = plt)
    plt <- pheatmap::pheatmap(mnorm, show_rownames = F, main = "log2(downsampled_reads+1)",border_color=NA)
    show(plt)
    ggplot2::ggsave(filename = paste0(folder, "/pheatmap_normReads_byInt.png"), plot = plt)
  } else {
    pheatmap::pheatmap(mraw, show_rownames = F, filename = paste0(folder, "/pheatmap_rawReads_byInt.png"), main = "log2(raw_reads+1)"
                      ,border_color=NA)
    pheatmap::pheatmap(mnorm, show_rownames = F, filename = paste0(folder, "/pheatmap_normReads_byInt.png"),
                       main = "log2(downsampled_reads+1)",border_color=NA)
  }
  ## correlate quantification between datasets: scores
  m = ints[, grep("score", names(ints))]
  msc = asinh(m[, grep("rep", names(m))])

  if (plot) {
    plt <- pheatmap::pheatmap(msc, show_rownames = F, main = "asinh(scores)",border_color=NA)
    show(plt)
    ggplot2::ggsave(filename = paste0(folder, "/pheatmap_Scores_byInt.png"), plot = plt)
  } else {
    pheatmap::pheatmap(msc, show_rownames = F, filename = paste0(folder, "/pheatmap_Scores_byInt.png"), 
              main = "asinh(scores)",border_color=NA)

  }

}


##TO TEST

## heatmaps to look at replicate correlations:
heatmapsCC <- function(grl, resfolder = "matrices_reps/", num) {
  #' @title Plot Replicate Correlation Heatmaps of Feature-overlapping Interactions
  #' @description
  #' Plots and saves correlation of replicates of each experiment and condition given as quality control
  #' @param grl GRangesList containing features by which ints should be plotted
  #' @param resfolder folder where matrices created by getMatrix() are saved for separate reps, optional: Default: 'matrices_reps/'
  #' @param num indices of features to be plotted
  #' @returns Correlation heatmaps of replicates of intervals by features
  #' @export


  if (!(dir.exists("pheatmaps"))) {
    dir.create("pheatmaps")
  }
  # make this an if?
  l = list()
  for (i in num) {
    l = c(l, list(read.delim(paste0("oneGeneOnePeak_", names(grl)[i], strsplit_string(resfolder, s1 = "matrices",
                                                                                      s2 = "/"), "_interacting_readsScoresInPeak.txt"), stringsAsFactors = FALSE)))
  }
  names(l) = names(grl)[num]

  for (i in 1:length(l)) {

    for (what in c("summit_bednorm_scores", "sum_bednorm_scores", "mean_bednorm_scores")) {

      a = l[[i]][, grep(what, names(l[[i]]))]
      names(a) = strsplit_string(names(a), s2 = paste0("[.]", what))

      typ = sapply(strsplit(names(a), split = "_rep"), FUN = function(x) x[1])

      annotation_col = data.frame(type = typ)
      row.names(annotation_col) = colnames(a)

      pheatmap::pheatmap(asinh(a[, 1:(ncol(a)/2)]), annotation_col = annotation_col, show_rownames = FALSE,
                         main = paste(names(l)[i], "asinh", what, " interactions"), file = paste0("pheatmaps/", names(l)[i],
                                                                                                  "_", what, "_interactions.png"))
      pheatmap::pheatmap(asinh(a[, (ncol(a)/2 + 1):ncol(a)]), annotation_col = annotation_col, show_rownames = FALSE,
                         main = paste(names(l)[i], "asinh", what, "controls"), file = paste0("pheatmaps/", names(l)[i], "_",
                                                                                             what, "_controls.png"))

    }

    for (what in c("summit_bednorm_reads", "sum_bednorm_reads", "mean_bednorm_reads")) {

      a = l[[i]][, grep(what, names(l[[i]]))]
      names(a) = strsplit_string(names(a), s2 = paste0("[.]", what))

      typ = sapply(strsplit(names(a), split = "_rep"), FUN = function(x) x[1])

      annotation_col = data.frame(type = typ)
      row.names(annotation_col) = colnames(a)

      pheatmap::pheatmap(log2(a[, 1:(ncol(a)/2)] + 1), annotation_col = annotation_col, show_rownames = FALSE,
                         main = paste(names(l)[i], "log2+1", what, " interactions"), file = paste0("pheatmaps/", names(l)[i],
                                                                                                   "_", what, "_interactions.png"))
      pheatmap::pheatmap(log2(a[, (ncol(a)/2 + 1):ncol(a)] + 1), annotation_col = annotation_col, show_rownames = FALSE,
                         main = paste(names(l)[i], "log2+1", what, "controls"), file = paste0("pheatmaps/", names(l)[i],
                                                                                              "_", what, "_controls.png"))

    }

  }
  dev.off()
}


loadGrl <- function(intDir = "intervals") {
  #' @title Load Annotation into a GenomicRangesList
  #' @description reads all intervals and reads them into a GRanges List, converts bed files into GRanges
  #' @param intDir folder containing all interval .rds files
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
  #' @returns GenomicRanges object Grl containing all intervals
  #' @export

  if (!file.exists(intDir)) {
    stop("enter a valid directory")
  }

  ilist = list.files(paste0(intDir))
  grl = GenomicRanges::GRangesList()
  names=vector()
  for (f in ilist) {

    if (length(grep('.rds$',f))>0){
      gr = stripGR(readRDS(paste0(intDir, "/", f)))
      name=strsplit_string(f, s2 = ".rds")
    } else if (length(grep('.bed$',f))>0){
      gr = bed2gr(read.delim(paste0(intDir, "/", f),header=FALSE))
      name=strsplit_string(f, s2 = ".bed")
    } else {
      message(paste(f,'is an unclear file type, please add extension!'))
    }

    grl = c(grl, GenomicRanges::GRangesList(gr))
    names=c(names,name)
  }
  names(grl) = names

  return(grl)
}


loadCdList <- function(resultsDir = "results",baited_genes = NULL) {
  #' @title Load Interactions List
  #' @description Reads all relevant Chicago files, (concatenated, not singular reps) and put them into a list
  #' @param baited_genes List of bait to control whether the baits in Chicago data is equal to the bait input of the
  #' Post Chicago pipeline, optional
  #' @param resultsDir Directory with Chicago data output. Default: 'results'
  #' @examples 
  #' # example code:
  #' extdata <- system.file("extdata", package="PostChicago")
  #' chicagoOutputDir <- file.path(extdata,'results')
  #' designDir <- file.path(extdata,'designDir')
  #' 
  #' #create baited_genes:
  #' baited_genes=baitmap2baited_genes(designDir,save=FALSE)
  #' 
  #' #load chicago data object list:
  #' cdlist <- loadCdList(resultsDir = chicagoOutputDir, baited_genes = baited_genes)
  #' 
  #' @returns List L of Chicago data
  #' @export

  if (!file.exists(resultsDir)) {
    stop("enter a valid directory")
  }

  types = list.files(pattern = "cd", path = resultsDir)
  types = strsplit_string(types, s1 = "cd_")

  L = list()
  for (type in types) {
    load(paste(paste0(resultsDir, "/cd"), type, sep = "_"))
    L = c(L, list(cd@x))
  }
  names(L) = types

  ids = unique(L[[1]]$baitID)

  if (is.null(baited_genes)) {
    message("No control for bait count... Add argument if you want to control output")
  } else if (!length(ids) == length(baited_genes)) {
    message("There is an unequal amount of baits in your list and bait files. Please check your Chicago Outputs...")
  }

  return(L)
}




runChicagoForPostChicago <- function(df, mySettings = NULL, outputDir = "results", DataPath = "dataPath", DesignDir = "designDir",
                       overwrite = FALSE) {
  #' @title Run Chicago
  #' @description to run chicago such that the output is compatible with our pipeline. Default settings are only executed if the files are not yet saved in outputDir.
  #' It is also possible to run chicago as is, but this output respects naming conventions for easier post chicago
  #' @param df data.frame containing at least two columns: 'filename' (.chinput filename) and 'type' (sample type)
  #' @param mySettings Chicago Settings for the pipeline, if not specified, the Default settings will be used
  #' @param outputDir Output directory, Default 'results'
  #' @param DataPath Directory of input files, Default 'dataPath'
  #' @param DesignDir Directory of Experiment Design Files, Default 'designDir'
  #' @param overwrite Overwrites cd that are present in outputDir if set to TRUE. Default: FALSE
  #' @returns Saves Chicago output files named after df$type, including oeNorm, techNoise, distFun as well as vignette and complete chicago data with scores as cd in outputDir
  #' @export

  if (is.null(mySettings)) {
    mySettings <- Chicago::defaultSettings()
  }

  types=unique(df$type)

  for (type in types) {

    file = paste("cd", "_", type, sep = "")

    ## To execute only if cd is not in outputDir or if overwrite is set to TRUE
    if (!file %in% dir(outputDir)) {

      files <- paste(DataPath, df[df$type == type, ]$filename, sep = "/")

      cd <- Chicago::setExperiment(designDir = DesignDir, settings = mySettings)
      cd <- Chicago::readAndMerge(files = files, cd = cd)
      cd <- Chicago::chicagoPipeline(cd, outprefix = paste(outputDir, "/", type, sep = ""))
      Chicago::exportResults(cd, file.path(outputDir, paste("vignetteOutput", type, sep = "_")))
      save(cd, file = paste(outputDir, "/cd", "_", type, sep = ""))

    } else if (overwrite) {
      files <- paste(DataPath, df[df$type == type, ]$filename, sep = "/")

      cd <- Chicago::setExperiment(designDir = DesignDir, settings = mySettings)
      cd <- Chicago::readAndMerge(files = files, cd = cd)
      cd <- Chicago::chicagoPipeline(cd, outprefix = paste(outputDir, "/", type, sep = ""))
      Chicago::exportResults(cd, file.path(outputDir, paste("vignetteOutput", type, sep = "_")))
      save(cd, file = paste(outputDir, "/cd", "_", type, sep = ""))

    } else {
      message(paste("cd file", file, "detected, skipping Chicago..."))
    }

  }

}





########################################
##Genome Browser visualization
########################################

###################################################
##normalization between cell lines using getnormbed
###################################################

##make normalized bedGraphs (between cell lines, conditions etc.)
##normalize to an average number of reads between cell lines/conditions etc.
##make bedGraphs for duplicates or triplicates
##make BedGraphs which are uploadable on UCSC
##renamed, previously: getBedGraphsNormGetbednorm

makeBedGraphs <- function(id, L_data,rmapgr,baited_genes,ints=NULL,L_norm=NULL, outputDir="bedgraph"){
  #' @title Create Bedgraphs From ChicagoData
  #' @description Takes ChicagoData objects supplied in the list L, normalizes them to read coverage mapped to baits.
  #' Creates two file types. 1) A bedgraph file containing normalized reads for each sample mapping to the bait defined by id.
  #' 2) Creates a bed file of the view point (bait).
  #' @param id baitID, must correspond to the baitIDs defined in baited_genes and the supplied rmapgr object
  #' @param L_data list of ChicagoData (cd) objects
  #' @param rmapgr GenomicRanges object containing restriction fragments. IDs must be saved in a column called 'id'
  #' @param baited_genes GenomicRanges object containing the information on view points (baited genes);
  #' usually generated from the .baitmap (see Chicago package), using for instance baitmap2baited_genes()
  #' @param ints optional, table containing all interactions across L_data, if supplied, a separate bed file containing all
  #' interactions mapping to the bait defined by id is created
  #' @param L_norm list of cd objects which should be used to normalize samples, useful if bedgraphs are to be created only
  #' for a subset of promoters or genomic locations, as it keeps the normalization constant to the whole dataset. Default: L_data
  #' @param outputDir output directory for bedgraphs, automatically created if not surrently present; Default: 'bedgraph'
  #' @returns Saves bedgraphs for the bait defined by id for each sample and a bedfile containing the view point.
  #' @export

  ##make sure genenames are supplied!
  if(is.null(names(baited_genes))){

    if(!is.null(baited_genes$genename)){
      names(baited_genes)=baited_genes$genename
    } else if (!is.null(baited_genes$genenames)) {
      names(baited_genes)=baited_genes$genenames
    } else {
      message('no genenames supplied with baited_genes!')
    }

  }

  ##make sure the outputdir exists, if not create one:
  if(!dir.exists(outputDir)){
      dir.create(outputDir)
  }

  name=names(baited_genes[baited_genes$re_id==id])
  seqnames=paste('chr',c(1:19,'X','Y'),sep='')

  ##get normalization parameters:
  L=L_data
  if (is.null(L_norm)){
    L_norm=L
  }

  ##for each sample write a normalized bg
  samples=names(L_data)

  bait=rmapgr[rmapgr$id==id]

  for (i in 1:length(L)){

    ##only if bedgraphs do not already exist
    if (!file.exists(paste0(outputDir,'/',name,'_id',id,'_',sample,'_bg.bed'))){

        print(i)

        sample=samples[i]

        bed=getnormbed(id,data=L[[i]],samplename = names(L)[i],L=L_norm,same.chromosome=TRUE,rmapgr=rmapgr)
        bed=bed[order(start(bed),decreasing=FALSE)]
        bed=bed[seqnames(bed) %in% seqnames]

        ##bg
        bg=data.frame(as.character(seqnames(bed)),start(bed),end(bed),bed$N)
        names(bg)=c('chromosome','start','end','score')

        ##save all bed files:
        saveBedGraph(bg, paste0(outputDir,'/',name,'_id',id,'_',sample,'_bg.bed'))

      }

  }


  ##make bedfiles of bait and of interactions, but only if the files do not exist yet:

  if (!file.exists(paste0(outputDir,'/',name,'_id',id,'_viewpoint.bed'))){
    bedbait=gr2bed(bait)
    saveBed(bedbait, paste0(outputDir,'/',name,'_id',id,'_viewpoint.bed'))
  }

  if (!is.null(ints) && !file.exists(paste0(outputDir,'/',name,'_id',id,'_ints.bed'))){

    bedints=ints[ints$baitID==id,]$otherEndID
    if (length(bedints)>0){
      bedints=gr2bed(rmapgr[rmapgr$id %in% bedints])
      saveBed(bedints, paste0(outputDir,'/',name,'_id',id,'_ints.bed'))
    }

  }


}




saveBedGraph <- function(bg, bgfile, filetype='bg',viewlimit=100,autoscale='on',color='0,0,0'){
  ##Helper function of getBedgraphsNormGedBedNorm()

  ##define the name of the header:
  if ( filetype=='bg' ){
    a=sapply(unlist(strsplit(bgfile,split='_bg.bed')),FUN=function(x) x[1])
    header=paste0("track type=bedGraph name=",a," description=",a, " visibility=full track height=100 color=",color," priority=10 autoScale=",autoscale," alwaysZero=on gridDefault=off maxHeightPixels=1000:100:20 graphType=bar windowingFunction=mean viewLimits=0:",viewlimit,sep='')
  } else if (filetype=='bed' ){
    a=sapply(unlist(strsplit(bgfile,split='.bed')),FUN=function(x) x[1])
    header=paste("track name=",a,"_bedid description=",a,"_baitID color=150,0,0 visibility=full",sep='')
  }

  f <- file(bgfile, "w")

  writeLines(header,f)
  saveBed(bg, f)

  close(f)

}


saveBed <- function(bed,bedfile){
  #' @title Save Bed File
  #' @description Saves a table as bedfile
  #' @param bed table in a bed format or regular table that should be saved without headers
  #' @param bedfile name of the file to which to save the table
  #' @examples
  #' # basic usage of saveBed:
  #' bed <- data.frame(chr='chr1',start=1000, end=2000)
  #' saveBed(bed, 'myfile.bed')
  #' @return saves a table as bed file (no header)
  #' @export


  write.table(bed,bedfile,sep='\t', eol='\n', quote=FALSE,row.names=FALSE,col.names=FALSE)
}


######################################################
##Aggregate interaction peaks over a specific distance
######################################################

##Function to aggregate peaks for interactions from CaptureC analysis (version 16.08.2022):

aggregatePeaks <- function(ints, dis, samples, fileprefix='ints'){
  #' @title Aggregate Interactions to Interaction Peaks
  #' @description Takes a table of interactions (ints) and aggregates interactions with the same bait over the distance dis in
  #' restriction fragments. Only creates the table if it is not yet saved. Otherwise loads the existing table.
  #' Caution: All interactions, regardless the sample from which they originated, are aggregated,
  #' similarly to the reduce() function from GenomicRanges.
  #' @param ints table containing all interactions that should be aggregated
  #' @param dis numeric, distance in restriction fragments over which interactions should be aggregated.
  #' dis=1 will merge two adjacent fragments with significant interactions.
  #' @param samples Character vector of sample names, should correspond to the sample names used in ints.
  #' Typically corresponding to the names provided in the list of chicagoData objects cd
  #' @param fileprefix start of the filename for the aggregatePeaks table file, will be added the following extension:
  #' paste0('_aggregatePeaks_in_regions_',dis,'bp.txt'). Should contain path if required to be saved in a different folder. Default: ints.
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
  #' @returns A table containing the position of aggregated interactions. Redefines intIDs so that the ID of the
  #' first interacting fragment in a peak stands in as otherEnd.
  #' Data is automatically saved in a 'aggregatePeaks_dis' text file whose name is defined by fileprefix and dis.
  #' Also saved as a bedfile.
  #' @export

  file=paste0(fileprefix,'_aggregatePeaks_dis',dis,'.txt')

  if (file %in% list.files()){

    message(paste('dis',dis,', aggregatePeaks file exists'))
    d=read.delim(file,stringsAsFactors = FALSE)

  } else {

    message(paste('dis',dis,', aggregating peaks...'))
    d=ints
    L_peaks=list()

    ids=unique(ints$baitID)

    for (id in ids){

      a=d[d$baitID==id,]
      gr=GenomicRanges::GRanges(a$seqnames_otherEnd, IRanges::IRanges(a$otherEndID,a$otherEndID))
      values(gr)=a
      gr=GenomicRanges::reduce(gr,min.gapwidth=dis)

      L_ints=list()
      for (i in 1:length(gr)){
        oe=start(gr[i]):end(gr[i])
        L_ints[[i]]=a[a$otherEndID %in% oe,]
      }

      ##aggregate to peaks, averaging everything

      for (i in 1:length(L_ints)){

        b=L_ints[[i]]
        c=b[1,]

        if (nrow(b)>1){
          oe=b[1,]$otherEndID
          for (j in 2:nrow(b)){oe=paste(oe,b[j,]$otherEndID,sep=',')}
          c$otherEndID=oe

          oe=b[1,]$intID
          for (j in 2:nrow(b)){oe=paste(oe,b[j,]$intID,sep=',')}
          c$intID=oe
        }

        c$start_otherEnd=min(b$start_otherEnd)
        c$end_otherEnd=max(b$end_otherEnd)

        L_peaks=c(L_peaks,list(c))

      }

    }

    df=do.call('rbind',L_peaks)
    d=df
    d$clusters_refined=d$clusters

    write.table(d,file,sep='\t',eol='\n',quote=FALSE,row.names=FALSE)
    saveBed(d[,c('seqnames_otherEnd', 'start_otherEnd', 'end_otherEnd','intID')], paste0(file,'.bed'))

  }

  return(d)

}



##aggregate by distance in bp:

aggregatePeaks_regions <- function(ints,dis,samples,fileprefix='ints'){

  #' @title Aggregate Interactions to Interaction Peaks by Distance in bp
  #' @description Takes a table of interactions (ints) and aggregates interactions with the same bait over the distance dis in bp.
  #' Only creates the table if it is not yet saved. Otherwise loads the existing table.
  #' Caution: All interactions, regardless the sample from which they originated, are aggregated, similarly to the reduce() function
  #' from GenomicRanges.
  #' @param ints Interactions table containing all interactions that should be aggregated
  #' @param dis numeric in bp, fragments which lie within this distance from each other will be aggregated into regions.
  #' @param samples Character vector of sample names, should correspond to the sample names used in ints.
  #' Typically corresponding to the names provided in the list of chicagoData objects cd
  #' @param fileprefix start of the filename for the aggregatePeaks table file, will be added the following extension:
  #' paste0('_aggregatePeaks_in_regions_',dis,'bp.txt'). Should contain path if required to be saved in a different folder. Default: ints.
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

  file=paste0(fileprefix,'_aggregatePeaks_in_regions_',dis,'bp.txt')

  if (file %in% list.files()){
    message(paste('distance of ',dis,',aggregatePeaks file exists'))

    d=read.delim(file,stringsAsFactors = FALSE)

  } else {

    message(paste('distance of',dis,', aggregating peaks...'))

    d=ints
    L_peaks=list()
    ids=unique(ints$baitID)

    for (id in ids){

      a=d[d$baitID==id,]
      a <- a[order(a$seqnames_otherEnd, a$start_otherEnd),]
      gr=GenomicRanges::GRanges(a$seqnames_otherEnd, IRanges::IRanges(a$start_otherEnd,a$end_otherEnd))

      values(gr)=a
      gr=GenomicRanges::reduce(gr, min.gapwidth = dis, with.revmap = TRUE)

      revmap <- gr$revmap

      otherEndID <- a$otherEndID

      mapped <- lapply(revmap, function(idx) otherEndID[idx])

      L_ints=list()

      for (i in 1:length(mapped)){

        oe=mapped[[i]]
        L_ints[[i]]=a[a$otherEndID %in% oe,]

      }

      ##aggregate to peaks, averaging everything

      for (i in 1:length(L_ints)){

        b=L_ints[[i]]
        c=b[1,]

        if (nrow(b)>1){

          oe=b[1,]$otherEndID

          for (j in 2:nrow(b)){oe=paste(oe,b[j,]$otherEndID,sep=',')}

          c$otherEndID=oe
          oe=b[1,]$intID

          for (j in 2:nrow(b)){oe=paste(oe,b[j,]$intID,sep=',')}

          c$intID=oe

        }

        c$start_otherEnd=min(b$start_otherEnd)
        c$end_otherEnd=max(b$end_otherEnd)

        L_peaks=c(L_peaks,list(c))

      }

    }

    df=do.call('rbind',L_peaks)

    d=df

    write.table(d,file,sep='\t',eol='\n',quote=FALSE,row.names=FALSE)
    saveBed(d[,c('seqnames_otherEnd','start_otherEnd','end_otherEnd','intID')], paste0(file,'.bed'))

  }

  return(d)

}



###########################################################
##Function to plot statistics from the chicagoData lists
###########################################################

##helper functions:
nGenesSig<-function(x){return(length(table(x[x$score>=5,]$baitID)))}
nreads <- function(x){return(sum(x$N))}
fracReadsInSig <- function(x){return(  round(sum(x[x$score>=5,]$N)/sum(x$N)*100,2) )}

plotSigIntsStats <- function(L,col=NULL,lty=NULL,filename=NULL,plotDistribution=TRUE,plotExamples=FALSE) {
  #' @title Plot Interaction Statistics
  #' @description Takes a list of ChicagoData objects and plots a per sample statistial summary.
  #' @param L list of ChicagoData tables for which to extract statistics
  #' @param col optional; vector, colours for the plots
  #' @param lty optional; vector, line type
  #' @param filename optional; if provided, the plots will be saved under filename; must have the extension '.pdf'
  #' @param plotDistribution if TRUE (default), read distribution per bait is plotted.
  #' @param plotExamples if TRUE, example distributions will be plotted based on data from Mahara et al, default: FALSE
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
  #' @returns Plots with numbers of interactions etc. for each sample and bait.
  #' @export


  if(is.null(col)){
    col=c(RColorBrewer::brewer.pal(8,'Set2'),RColorBrewer::brewer.pal(9,'Set1'))
  }

  if(is.null(lty)){
    lty=rep(1,length(L))
  }

  if(!is.null(filename)){
    pdf(filename,width=14,height=10)
  }

  par(mfrow=c(3,3))

  ##a) total reads

  bp=barplot(sapply(L,FUN=nreads),las=3,main='total reads mapped to baits',col=col,ylim=c(0,1.4*max(sapply(L,FUN=nreads))))
  text(bp,sapply(L,FUN=nreads)*1.2,labels=sapply(L,FUN=nreads),cex=0.75)

  ##b) total genes with sig ints

  bp=barplot(sapply(L,FUN=nGenesSig),las=3,main='total genes with significant ints',ylim=c(0,1.4*max(sapply(L,FUN=nGenesSig))),col=col)
  text(bp,sapply(L,FUN=nGenesSig)*1.2,labels=sapply(L,FUN=nGenesSig),cex=0.75)

  ##c) total unique ints

  bp=barplot(sapply(L,FUN=function(x) nrow(x)),las=3,main='total unique ints',ylim=c(0,1.4*max(sapply(L,FUN=function(x) nrow(x)))),col=col)
  text(bp,sapply(L,FUN=function(x) nrow(x))*1.3,labels=sapply(L,FUN=function(x) nrow(x)),cex=0.75)

  ##d) total sig ints

  bp=barplot(sapply(L,FUN=function(x) nrow(x[x$score>=5,])),las=3,main='total significant ints',ylim=c(0,1.4*max(sapply(L,FUN=function(x) nrow(x[x$score>=5,])))),col=col)
  text(bp,sapply(L,FUN=function(x) nrow(x[x$score>=5,]))*1.3,labels=sapply(L,FUN=function(x) nrow(x[x$score>=5,])),cex=0.75)

  ##e) % reads within sig ints

  bp=barplot(sapply(L,FUN=fracReadsInSig),las=3,main='% reads within sig ints',ylim=c(0,1.4*max(sapply(L,FUN=fracReadsInSig))),col=col)
  text(bp,sapply(L,FUN=fracReadsInSig)*1.3,labels=sapply(L,FUN=fracReadsInSig),cex=0.75)

  ##if distribution should not be plotted, plot examples

  if (plotDistribution){

    ##f) distribution total reads per bait

    nreadsPerBait=function(x){
      return(sapply(unique(x$baitID), FUN=function(y) sum(x[x$baitID==y,]$N) ))
    }
    dat=lapply(L,FUN=nreadsPerBait)

    plot(density(dat[[1]]),lwd=2,xlab='total reads per bait', main='total reads per bait distribution',col=col[1],lty=lty[1],xlim=c(0,1.25*max( unlist(dat)  )))
    for (i in 2:length(L)){
      lines(density(dat[[i]]),col=col[i],lwd=2,lty=lty[i])
    }
    legend('topright',cex=0.75,bty='n',col=col[1:length(L)],legend = paste0(names(L),',median=',round(sapply(dat, FUN='median'),0)),lwd=2,lty=lty)

    ##g) distribution total ints per bait

    dat=lapply(L,FUN=function(x) table(x$baitID))
    plot(density(dat[[1]]),lwd=2,xlab='unique ints per bait', main='unique ints per bait distribution',col=col[1],lty=lty[1],xlim=c(0,1.25*max( unlist(dat)  )))
    for (i in 2:length(L)){
      lines(density(dat[[i]]),col=col[i],lwd=2,lty=lty[i])
    }
    legend('topright',cex=0.75,bty='n',col=col[1:length(L)],legend = paste0(names(L),',median=',round(sapply(dat, FUN='median'),0)),lwd=2,lty=lty)

    ##h) distribution sig ints per bait

    dat=lapply(L,FUN=function(x) table(x[x$score>=5,]$baitID))
    plot(density(dat[[1]]),lwd=2,ylim=c(0,0.05),xlab='significant ints per bait', main='sig ints per bait distribution',col=col[1],lty=lty[1])
    for (i in 2:length(L)){
      lines(density(dat[[i]]),col=col[i],lwd=2,lty=lty[i])
    }
    legend('topright',cex=0.75,bty='n',col=col[1:length(L)],legend = paste0(names(L),',median=',round(sapply(dat, FUN='median'),0)),lwd=2,lty=lty)

    ##i) distribution reads per sig int per bait (log2)

    nSigReadsPerIntPerBait=function(x){
      return(sapply(unique(x$baitID), FUN=function(y) median(x[x$baitID==y & x$score>=5,]$N,na.rm = TRUE) ))
    }
    dat=lapply(L,FUN=nSigReadsPerIntPerBait)
    dat=lapply(dat,FUN=function(x) replace(x,is.na(x),0))

    plot(density(dat[[1]]),lwd=2,xlab='Median reads in sig ints per bait', main='Median reads in sig ints per bait distribution',col=col[1],lty=lty[1])
    for (i in 2:length(L)){
      lines(density(dat[[i]]),col=col[i],lwd=2,lty=lty[i])
    }
    legend('topright',cex=0.75,bty='n',col=col[1:length(L)],legend = paste0(names(L),',median=',round(sapply(dat, FUN='median'),0)),lwd=2,lty=lty)

  }  else if (plotExamples) {

    ##f) distribution total reads per bait

    load(paste0(chicagoOutputDir,'/fnreadsPerBait'))

    plot(density(dat[[1]]),lwd=2,xlab='total reads per bait', main='Example: \n total reads per bait distribution',col=col[1],lty=lty[1],xlim=c(0,1.25*max( unlist(dat)  )))
    for (i in 2:length(L)){
      lines(density(dat[[i]]),col=col[i],lwd=2,lty=lty[i])
    }
    legend('topright',cex=0.75,bty='n',col=col[1:length(L)],legend = paste0(names(L),',median=',round(sapply(dat, FUN='median'),0)),lwd=2,lty=lty)

    ##g) distribution total ints per bait

    load(paste0(chicagoOutputDir,'/gDistrTotIntsPerBait'))

    plot(density(dat[[1]]),lwd=2,xlab='unique ints per bait', main='Example: \n unique ints per bait distribution',col=col[1],lty=lty[1],xlim=c(0,1.25*max( unlist(dat)  )))
    for (i in 2:length(L)){
      lines(density(dat[[i]]),col=col[i],lwd=2,lty=lty[i])
    }
    legend('topright',cex=0.75,bty='n',col=col[1:length(L)],legend = paste0(names(L),',median=',round(sapply(dat, FUN='median'),0)),lwd=2,lty=lty)

    ##h) distribution sig ints per bait

    load(paste0(chicagoOutputDir,'/hSigIntsPerBait'))

    plot(density(dat[[1]]),lwd=2,ylim=c(0,0.05),xlab='significant ints per bait', main='Example: \n sig ints per bait distribution',col=col[1],lty=lty[1])
    for (i in 2:length(L)){
      lines(density(dat[[i]]),col=col[i],lwd=2,lty=lty[i])
    }
    legend('topright',cex=0.75,bty='n',col=col[1:length(L)],legend = paste0(names(L),',median=',round(sapply(dat, FUN='median'),0)),lwd=2,lty=lty)

    ##i) distribution reads per sig int per bait (log2)

    load(paste0(chicagoOutputDir,'/iDistReadsPerSigIntPerBait'))

    plot(density(dat[[1]]),lwd=2,xlab='Median reads in sig ints per bait', main='Example: \n Median reads in sig ints per bait distribution',col=col[1],lty=lty[1])
    for (i in 2:length(L)){
      lines(density(dat[[i]]),col=col[i],lwd=2,lty=lty[i])
    }
    legend('topright',cex=0.75,bty='n',col=col[1:length(L)],legend = paste0(names(L),',median=',round(sapply(dat, FUN='median'),0)),lwd=2,lty=lty)

  }
  if(!is.null(filename)){
    dev.off()
  }

}


##################################################
##################################################
##getMatrix() replaces the old getMatrix() and all
##of its helper functions
##
##doesn't make plots
##################################################
##################################################

####################################################################################################
##Helpers for getMatrix reads
####################################################################################################

findOthers <- function(intID,zoom){
  oe=as.numeric(unlist(strsplit(intID,split=';'))[2])
  minoe=oe-zoom
  maxoe=oe+zoom
  return(minoe:maxoe)
}

extractValues <- function(downsampledReads, bait,others){
  df=data.frame(others=others,Ndownsampled=0)
  reads=downsampledReads[downsampledReads$baitID==bait & downsampledReads$otherEndID %in% others,c('otherEndID','N')][order(
    downsampledReads[downsampledReads$baitID==bait & downsampledReads$otherEndID %in% others,c('otherEndID','N')]$otherEndID),]
  df[df$others %in% reads$otherEndID,]$Ndownsampled=reads$N
  return(df$Ndownsampled)
}

annotateSampleMatrix <- function(intID,downsampledReads,zoom){
  bait=as.numeric(unlist(strsplit(intID,split=';'))[1])
  others=findOthers(intID,zoom)
  return(extractValues(downsampledReads,bait,others))
}

annotateControlMatrix <- function(intID,downsampledReads,zoom){
  bait=as.numeric(unlist(strsplit(intID,split=';'))[1])
  otherEnd=as.numeric(unlist(strsplit(intID,split=';'))[2])
  ##remake intID to be distance-matched on the other side of the bait:
  intID=paste(bait,(bait+bait-otherEnd),sep=';')
  others=findOthers(intID,zoom)   ##no need to check chromosome, values will be 0 if different chromosome!
  return(extractValues(downsampledReads,bait,others))
}

####################################################################################################
##Helpers for getMatrix scores
####################################################################################################

extractValuesScores <- function(downsampledReads, bait,others){
  df=data.frame(others=others,Ndownsampled=0)
  reads=downsampledReads[downsampledReads$baitID==bait & downsampledReads$otherEndID %in% others,c('otherEndID','score')][order(
    downsampledReads[downsampledReads$baitID==bait & downsampledReads$otherEndID %in% others,c('otherEndID','score')]$otherEndID),]
  df[df$others %in% reads$otherEndID,]$Ndownsampled=reads$score
  return(df$Ndownsampled)
}

annotateSampleMatrixScores <- function(intID,downsampledReads,zoom){
  bait=as.numeric(unlist(strsplit(intID,split=';'))[1])
  others=findOthers(intID,zoom)
  return(extractValuesScores(downsampledReads,bait,others))
}

annotateControlMatrixScores <- function(intID,downsampledReads,zoom){
  bait=as.numeric(unlist(strsplit(intID,split=';'))[1])
  otherEnd=as.numeric(unlist(strsplit(intID,split=';'))[2])
  ##remake intID to be distance-matched on the other side of the bait:
  intID=paste(bait,(bait+bait-otherEnd),sep=';')
  others=findOthers(intID,zoom)   ##no need to check chromosome, values will be 0 if different chromosome!
  return(extractValuesScores(downsampledReads,bait,others))
}

####################################################################################################
####################################################################################################
##main new getMatrix() function:
####################################################################################################
####################################################################################################


getMatrix <- function(L, zoom=100, d, resfolder = "matrices/", type = 'reads', readnorm='downsampling',name = '', norm = FALSE,
                       normsam = 1, pseudo = 1, rmapgr = rmapgr){
  #' @title Create Interaction Matrices
  #' @description
  #' Function that creates, normalizes and plots matrices with reads or scores surrounding an interaction as determined in d. Interaction is always
  #' determined as a restriction fragment to restriction fragment interaction. If interaction peaks or interval peaks are used, it is advised
  #' to either determine the peak summit or the central fragment that overlaps an interval.
  #' CAUTION: Only creates a matrix if the matrix doesn't exist yet in the specified resfolder
  #' The function outputs a list with interaction matrices and distance-matched control matrices (con). For each interaction the control is considered to be
  #' an interaction between the same promoter and a fragment at the same distance (in restriction fragments) on the opposite side of this promoter.
  #' Output matrices are stranded, such that the promoter is to the right (+ direction) of the interaction. The matrices are saved as rds in the resfolder and can be loaded later.
  #' Figures (png) are saved in the same folder. Before recalculating the matrices the function checks the resfolder if they exist. In that case the matrices are loaded and only the figures are recreated.
  #' Normalised matrices are calculated from unnormalised matrices if they exist.
  #' Normalisation happens to the center of interaction in the reference sample normsam. Normsam is recommended to be a positive control, expected to have the highest values in the samples.
  #' Plots can be interpreted as fraction of interaction strength compared to the center of interaction in sample normsam.
  #' @param L list with chicago data objects (cd(at)x)
  #' @param zoom distance around interaction in restriction fragments
  #' @param d A table with cis-chromosomal interactions. Should contain the columns baitID, otherEndID, defined as in Chicago and can have a column intID, which is expected to be in the format baitID;otherEndID.
  #' otherEndID should be the center of interaction, from which the matrices will be calculated.
  #' @param resfolder path to the directory, where the matrices and plots should be stored,
  #' optional, Default: 'matrices/'. Creates the directory if necessary.
  #' @param type 'reads' or 'scores'. If reads, then downsampling to the sample with the smallest amount of reads mapped to baits is
  #' performed. Default: 'reads'
  #' @param readnorm Read normalization, either 'downsampling' or 'scaling'. Scaling is the same as for Lineplots. Default: 'downsampling'
  #' @param name Name of the matrix.
  #' @param norm Logical input. Should the data be normalised?, optional: Default: TRUE
  #' @param normsam To which sample in the list L should the data be normalised? optional: Default: 1 (the first sample!).
  #' @param pseudo pseudocount to be added for normalisation between matrices. Default: 1
  #' @param rmapgr GenomicRanges object containing restriction fragments. IDs must be saved in a column called 'id'
  #' @return List of interaction matrices (Lm), also saved in the resfolder. Loadable as Lm.
  #' @export


  message(paste('Matrices of', type))

  ##if resfolder does not exist, create it!
  if(!file.exists(resfolder)){dir.create(resfolder)}

  ##define filename:
  filename=paste0(resfolder,'matrices_',name,'_',type)

  if(type=='reads' && readnorm=='scaling'){
    filename=paste0(resfolder,'matrices_',name,'_',type,'_',readnorm)
  }

  if(!file.exists(filename)){

    ##keep only interactions on the same chromosome!
    d=d[d$seqnames_bait==d$seqnames_otherEnd,]

    ##step1: normalize reads in L
    message('normalizing reads...')

    if(readnorm=='downsampling'){
      message('downsampling reads...')

      nreads=sapply(L,FUN=function(x) sum(x$N))
      downsamplingFacs=min(nreads)/nreads*0.9999
      Ldownsampled=L
      for (i in 1:length(Ldownsampled)){
        Ldownsampled[[i]]$N=Ldownsampled[[i]]$N*downsamplingFacs[i]
      }

    } else {
      message('scaling reads upon coverage normalisation using getsk() *100000 reads and *nproms...')

      scalingFacs=getsk(L)
      Ldownsampled=L
      for (i in 1:length(Ldownsampled)){
        Ldownsampled[[i]]$N=Ldownsampled[[i]]$N*scalingFacs[i]
      }

    }

    ##step2: create a list of matrices Lm of length nrow(d) with all 0's; add intIDs as row names
    message('creating empty matrix lists...')

    m=matrix(0,nrow=nrow(d),ncol=(2*zoom+1))
    row.names(m)=d$intID

    k=1
    Lm=list()
    while(k<=length(Ldownsampled)){
      Lm=c(Lm,list(m))
      k=k+1
    }

    names(Lm) = names(Ldownsampled)

    ##step3: create a list of control matrices Lmc of length nrow(d) with all 0's; add intIDs as row names; define otherEnds

    Lmc=Lm
    names(Lmc)=paste0(names(Lm),'_con')

    ##step4: populate Lm and Lmc with downsampled values
    message('making matrices, this may take a while...')

    if (type == 'reads'){

      ##annotate Lm, the list of normal matrices:
      for (i in 1:length(Ldownsampled)){
        matr=lapply(d$intID, FUN=annotateSampleMatrix, downsampledReads=Ldownsampled[[i]], zoom=zoom)
        matr=do.call('rbind',matr)
        Lm[[i]]=matr
        row.names(Lm[[i]])=d$intID
        gc()
      }

      ##annotate Lmc, the list of control matrices:
      message('making distance-matched control matrices, this may take a while...')
      for (i in 1:length(Ldownsampled)){
        matr=lapply(d$intID, FUN=annotateControlMatrix, downsampledReads=Ldownsampled[[i]], zoom=zoom)
        matr=do.call('rbind',matr)
        Lmc[[i]]=matr
        row.names(Lmc[[i]])=d$intID
        gc()
      }

    } else if (type == 'scores'){

      ##annotate Lm, the list of normal matrices:
      for (i in 1:length(Ldownsampled)){
        matr=lapply(d$intID, FUN=annotateSampleMatrixScores, downsampledReads=Ldownsampled[[i]],zoom=zoom)
        matr=do.call('rbind',matr)
        Lm[[i]]=matr
        row.names(Lm[[i]])=d$intID
        gc()
      }

      ##annotate Lmc, the list of control matrices:
      message('making distance-matched control matrices, this may take a while...')
      for (i in 1:length(Ldownsampled)){
        matr=lapply(d$intID, FUN=annotateControlMatrixScores, downsampledReads=Ldownsampled[[i]],zoom=zoom)
        matr=do.call('rbind',matr)
        Lmc[[i]]=matr
        row.names(Lmc[[i]])=d$intID
        gc()
      }

    }

    ##step5: flip if the otherEndID is larger than baitID
    message('flipping the matrices...')

    biggerthan=d$otherEndID>d$baitID
    smallerthan=d$otherEndID<d$baitID

    Lmflipped=Lm
    Lmcflipped=Lmc

    for (i in 1:length(Lm)){
      Lmflipped[[i]]=rbind(Lm[[i]][smallerthan,],
                           t(apply(Lm[[i]][biggerthan,],1,FUN=function(x) rev(x))))
    }

    for (i in 1:length(Lmc)){
      Lmcflipped[[i]]=rbind(Lmc[[i]][biggerthan,],
                            t(apply(Lmc[[i]][smallerthan,],1,FUN=function(x) rev(x))))
    }

    ##step6: combine matrices and save them in the resfolder
    message('saving the matrices...')
    Lm=c(Lmflipped,Lmcflipped)
    save(Lm,file=filename)

  } else {

    ##if file exists already, load the data into an object called Lm
    message('file exists, loading the matrices...')
    load(filename)

  }

  ##step7: normalize matrices to the central 3 frags of the normsam matrix, adding pseudocount of 1

  if (norm){

    filename2=paste0(filename,'_normto.',names(L)[normsam])

    ##only execute if file doesn't exist yet!
    if(!file.exists(filename2)) {

      message('making normalized matrices...')

      Lmnorm=Lm

      normData=Lm[[normsam]]
      pseudo=1
      normValues=normData[,(2*zoom+1-zoom)]

      for (i in 1:length(Lm)){
        for (j in 1:nrow(Lm[[i]])){
          Lmnorm[[i]][j,]=(Lm[[i]][j,])/(normValues[j])
        }
      }

      Lm=Lmnorm

      save(Lm, file=filename2)

    } else {

      message('loading normalized matrices...')
      load(filename2)

    }

  }

  return(Lm)

}


