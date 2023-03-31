library(RUnit)
library(CytofkitNormalization)
library(ggplot2)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_markerNames()
    test_normalizeMarker()
    test_normalizeMarker_dashInName()
    test_getCluster()
    test_createTableForViolinPlot()
    test_createTableForViolinPlot_matrixSupplied()

    test_calculateColorBoundaries()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    f <- system.file(package="CytofkitNormalization", "extdata", "cytofkit-leukemia.RData")
    checkTrue(file.exists(f))
    x <- CytofkitNormalization$new(f)
    checkTrue(all(c("R6", "CytofkitNormalization") %in% class(x)))

    mtx <- x$getMatrix()
    checkEquals(dim(mtx), c(50000, 58))

    tbl.tsne <- x$getTsne()
    checkEquals(dim(tbl.tsne), c(50000, 2))

    clusters <- x$getClusters()
    checkEquals(length(clusters), 50000)
    checkEquals(clusters[["CBday4_2022-09-16_ery-MLL_Histone_LIVE_46719"]], 1)

} # test_ctor
#----------------------------------------------------------------------------------------------------
# give each intricate marker name a standard
test_markerNames <- function()
{
    message(sprintf("--- test_markerNames"))

    f <- system.file(package="CytofkitNormalization", "extdata", "cytofkit-leukemia.RData")
    checkTrue(file.exists(f))
    x <- CytofkitNormalization$new(f)

    x$createSimpleMarkerNames()

    markers <- x$getMarkers()
    h3.markers <- grep("H3", names(markers), v=TRUE)
    checkEquals(length(h3.markers), 20)
    h4.markers <- grep("H4", names(markers), value=TRUE)
    checkEquals(length(h4.markers), 4)

    checkTrue("H3" %in% h3.markers)
    checkTrue("H4" %in% h4.markers)

    checkEquals(markers[["H3"]], "Yb176Di<176Yb_H3>")
    checkEquals(markers[["H4"]], "Nd142Di<142Nd_H4>")

    checkTrue(length(markers) > 35)
    checkTrue(length(markers) < 45)

       # we later encountered long marker names with two underscores.
       # make sure we handle them properly

    f <- system.file(package="CytofkitNormalization", "extdata",
                     "cytofkit-dashNameExample.RData")
    checkTrue(file.exists(f))
    x <- CytofkitNormalization$new(f)

    x$createSimpleMarkerNames()

    markers <- x$getMarkers()
    h3.markers <- grep("H3", names(markers), v=TRUE)
    checkEquals(length(h3.markers), 12)
    h4.markers <- grep("H4", names(markers), value=TRUE)
    checkEquals(length(h4.markers), 2)

    checkTrue("H3" %in% h3.markers)
    checkTrue("H4" %in% h4.markers)

    checkEquals(markers[["H3"]], "Yb176Di<176Yb_H3>")
    checkEquals(markers[["H4"]], "Nd142Di<142Nd_H4>")


    checkTrue("c_Myc" %in% names(markers))
    checkTrue("MLL_Nterm" %in% names(markers))
    checkEquals(markers[["c_Myc"]], "Dy163Di<163Dy_c_Myc>")
    checkEquals(markers[["MLL_Nterm"]], "Gd160Di<160Gd_MLL_Nterm>")

    mtx <- x$getMatrix()
    moi <- c("H3", "H4", "c_Myc", "MLL_Nterm")
    longNames <- unlist(lapply(moi, function(shortName)  markers[[shortName]]))
    checkTrue(all(longNames %in% colnames(mtx)))

    moi <- h3.markers
    longNames <- unlist(lapply(moi, function(shortName)  markers[[shortName]]))
    checkTrue(all(longNames %in% colnames(mtx)))

    moi <- h4.markers
    longNames <- unlist(lapply(moi, function(shortName)  markers[[shortName]]))
    checkTrue(all(longNames %in% colnames(mtx)))

    checkTrue(length(markers) > 35)
    checkTrue(length(markers) < 45)

} # test_markerNames
#----------------------------------------------------------------------------------------------------
test_normalizeMarker <- function()
{
    message(sprintf("--- test_normalizeMarker"))

    f <- system.file(package="CytofkitNormalization", "extdata", "cytofkit-leukemia.RData")
    checkTrue(file.exists(f))
    x <- CytofkitNormalization$new(f)

    x$createSimpleMarkerNames()
    markers <- x$getMarkers()

    h3.markers <- grep("H3", markers, value=TRUE)
    h4.markers <- grep("H4", markers, value=TRUE)

    target <- "H3K4me3"
    h3.reference <- "H3"
    h4.reference <- "H4"

        #-----------------------------------
        # first just against H3
        #-----------------------------------

    new.col.name <- x$normalizeMarker(target, h3.reference)
    checkEquals(new.col.name, "H3K4me3.regress.H3")
        # make sure the new column is also in the markers data structure
    markers <- x$getMarkers()
    checkTrue(new.col.name %in% names(markers))
    checkTrue(markers[[new.col.name]] == new.col.name)

      # a hasty un-nuanced check is to get the quartiles of
      # of the distribution before and after normaliaztion,
      # check their sums

    h3.normalized.vec <- x$getMatrix()[, new.col.name]
    quick.assay.normalized <- sum(as.numeric(fivenum(h3.normalized.vec)))
       # -4.32963408 -0.25168910  0.01600153  0.29289136  1.89312134
    checkTrue(quick.assay.normalized < -2)
    raw.vec <- x$getMatrix()[, markers[["H3K4me3"]]]
    quick.assay.raw <- sum(as.numeric(fivenum(raw.vec)))
       #  -0.006413926  3.504356736  3.804160211  4.115765684  5.413913164
    checkTrue(quick.assay.raw > 16)

        #-----------------------------------
        # now against H4
        #-----------------------------------

    new.col.name <- x$normalizeMarker(target, h4.reference)
    checkEquals(new.col.name, "H3K4me3.regress.H4")

      # a hasty un-nuanced check is to get the quartiles of
      # of the distribution before and after normaliaztion,
      # check their sums

    h4.normalized.vec <- x$getMatrix()[, new.col.name]
    quick.assay.normalized <- sum(as.numeric(fivenum(h4.normalized.vec)))
       # -4.08416812 -0.23665835  0.01581525  0.28050035  1.60220508
    checkTrue(quick.assay.normalized < -2)
    quick.assay.raw <- sum(as.numeric(fivenum(raw.vec)))
       #  -0.006413926  3.504356736  3.804160211  4.115765684  5.413913164
    checkTrue(quick.assay.raw > 16)

        #-----------------------------------
        # H3 and  H4
        #-----------------------------------

    new.col.name <- x$normalizeMarker(target, c(h3.reference, h4.reference))
    checkEquals(new.col.name, "H3K4me3.regress.H3+H4")

      # a hasty un-nuanced check is to get the quartiles of
      # of the distribution before and after normaliaztion,
      # check their sums

    h3h4.normalized.vec <- x$getMatrix()[, new.col.name]
    quick.assay.normalized <- sum(as.numeric(fivenum(h3h4.normalized.vec)))
       # -4.13746033 -0.23606436  0.01371299  0.27937228  1.59271114
    checkTrue(quick.assay.normalized < -2)
    quick.assay.raw <- sum(as.numeric(fivenum(raw.vec)))
       #  -0.006413926  3.504356736  3.804160211  4.115765684  5.413913164
    checkTrue(quick.assay.raw > 16)

} # test_normalizeMarker
#----------------------------------------------------------------------------------------------------
test_normalizeMarker_dashInName <- function()
{
    message(sprintf("--- test_normalizeMarker_dashInName"))

        # this file has MLL-2 protein, a problem for lm formula
    f <- system.file(package="CytofkitNormalization", "extdata",
                     "cytofkit-dashNameExample.RData")
    checkTrue(file.exists(f))
    x <- CytofkitNormalization$new(f)

    x$createSimpleMarkerNames()
    markers <- x$getMarkers()

    h3.markers <- grep("H3", markers, value=TRUE)
    h4.markers <- grep("H4", markers, value=TRUE)

    target <- "MLL-2"
    checkTrue(target %in% names(markers))
    h3.reference <- "H3"
    h4.reference <- "H4"

        #-----------------------------------
        # H3 and  H4
        #-----------------------------------

    new.col.name <- x$normalizeMarker(target, c(h3.reference, h4.reference))
    checkEquals(new.col.name, "MLL-2.regress.H3+H4")

      # a hasty un-nuanced check is to get the quartiles of
      # of the distribution before and after normaliaztion,
      # check their sums.  normalized should approach zero.
      # raw (unnormalized) should be a largish positive number

    h3h4.normalized.vec <- x$getMatrix()[, new.col.name]
    quick.assay.normalized <- abs(sum(as.numeric(fivenum(h3h4.normalized.vec))))
       # -2.19692586 -0.36953319  0.04917722  0.41320936  2.00188338
    checkTrue(quick.assay.normalized < 0.2)  # 0.1021
    raw.vec <- x$getMatrix()[, markers[["MLL-2"]]]
    quick.assay.raw <- sum(as.numeric(fivenum(raw.vec)))
       #  -0.006413926  3.504356736  3.804160211  4.115765684  5.413913164
    checkTrue(quick.assay.raw > 9)

       #------------------------------------------------------------
       # audrey reports problems with  MLL_Nterm and c_myc as well
       #------------------------------------------------------------

    target <- "MLL-Cterm"
    checkTrue(target %in% names(markers))
    h3.reference <- "H3"
    h4.reference <- "H4"

        #-----------------------------------
        # H3 and  H4
        #-----------------------------------

    new.col.name <- x$normalizeMarker(target, c(h3.reference, h4.reference))
    checkEquals(new.col.name, "MLL-Cterm.regress.H3+H4")

      # a hasty un-nuanced check is to get the quartiles of
      # of the distribution before and after normaliaztion,
      # check their sums.  normalized should approach zero.
      # raw (unnormalized) should be a largish positive number

} # test_normalizeMarker_dashInName
#----------------------------------------------------------------------------------------------------
test_getCluster <- function()
{
    message(sprintf("--- test_getCluster"))

    f <- system.file(package="CytofkitNormalization", "extdata", "cytofkit-leukemia.RData")
    checkTrue(file.exists(f))
    x <- CytofkitNormalization$new(f)

    x$createSimpleMarkerNames()
    markers <- x$getMarkers()

    clusters <- x$getClusters()
    checkEquals(length(clusters), 50000)
    checkEquals(clusters[["CBday4_2022-09-16_ery-MLL_Histone_LIVE_46719"]], 1)

    cells.in.cluster <- x$getCluster(17)
    head(cells.in.cluster)
    small.cluster.size <- table(clusters)[[17]]
    checkEquals(length(cells.in.cluster), small.cluster.size)

    checkEquals(unique(as.numeric(clusters[cells.in.cluster])), 17)

} # test_getCluster
#----------------------------------------------------------------------------------------------------
test_createTableForViolinPlot <- function()
{
    message(sprintf("--- test_createTableForViolinPlot"))

        #--------------------------------------
        # first test and plot: H3 in cluster 3
        #--------------------------------------

    f <- system.file(package="CytofkitNormalization", "extdata",
                     "cytofkit-leukemia.RData")
    checkTrue(file.exists(f))
    x <- CytofkitNormalization$new(f)

    x$createSimpleMarkerNames()
    markers <- x$getMarkers()
    clusters <- x$getClusters()

      # create a tbl from clusters 5, 3, 1 in that order
      # then use ordered factors on those cluster names to
      # plot the violins in 1,3,5 order
    tbl.violin <- x$createTableForViolinPlot(clusters=c(5,3,1), marker="H3")
    tbl.violin$name <- factor(tbl.violin$name,
                              levels = c("H3.c1", "H3.c3", "H3.c5"), ordered = TRUE)

    checkEquals(dim(tbl.violin), c(8081, 2))
    if(interactive()){
        p <- ggplot(tbl.violin, aes(x=name, y=value, fill=name)) +
             geom_violin() +
             geom_boxplot(width=.1) +
             theme(axis.text = element_text(size = 14)) +
             ggtitle("H3 - clusters 5,3,1")
        p <- p + theme(legend.position="none")
        p <- p + scale_fill_manual(values=rep("beige",  13))
        p
           # by visual inspection: are these ranges displayed, in order,
           # in the 3 violins?
        round(range(subset(tbl.violin, name=="H3.c1")$value), digits=2) # 0.32 6.54
        round(range(subset(tbl.violin, name=="H3.c3")$value), digits=2) # 0.00 6.67
        round(range(subset(tbl.violin, name=="H3.c5")$value), digits=2) # 0.00 4.71
        } # interactive




      #--------------------------------------------------
      # second test and plot: H3 in all clusters, 1:20
      #--------------------------------------------------

    tbl.violin <- x$createTableForViolinPlot(clusters=1:20, marker="H3")

    if(interactive()){
        ggplot(tbl.violin, aes(x=name, y=value, fill=name)) + geom_violin() +
            #coord_cartesian(ylim=c(-5,5)) +
            theme(axis.text = element_text(size = 14)) +
            ggtitle("H3 - all clusters")
        }


        #-----------------------------------------------------------------------
        # third test and plot: H3K4me3 in all clusters, 1:20, no normalization
        #-----------------------------------------------------------------------

    tbl.violin <- x$createTableForViolinPlot(clusters=1:20, marker="H3K4me3")
    checkTrue(sum(fivenum(tbl.violin$value)) > 16)

    if(interactive()){
        ggplot(tbl.violin, aes(x=name, y=value, fill=name)) + geom_violin() +
            #coord_cartesian(ylim=c(-5,5)) +
            theme(axis.text = element_text(size = 14)) +
            ggtitle("H3K4me3 - all clusters")
        }


        #-----------------------------------------------------------------------
        # fourth test and plot: H3K4me3 in all clusters, 1:20, regression
        # normalized against H3 and H4
        #-----------------------------------------------------------------------

    target <- "H3K4me3"
    h3.reference <- "H3"
    h4.reference <- "H4"

    new.col.name <- x$normalizeMarker(target, c(h3.reference, h4.reference))
    markers <- x$getMarkers()
    checkTrue(new.col.name %in% names(markers))
    checkTrue(markers[[new.col.name]] == new.col.name)

    checkEquals(new.col.name, "H3K4me3.regress.H3+H4")
    mtx.extended <- x$getMatrix()
    checkTrue(new.col.name %in% colnames(mtx.extended))

    tbl.violin <- x$createTableForViolinPlot(clusters=1:20, marker=new.col.name)
    checkTrue(sum(fivenum(tbl.violin$value)) < -2.0)

    if(interactive()){
        ggplot(tbl.violin, aes(x=name, y=value, fill=name)) + geom_violin() +
            #coord_cartesian(ylim=c(-5,5)) +
            theme(axis.text = element_text(size = 14)) +
            ggtitle(sprintf("%s - all clusters", new.col.name))
        }

} # test_createTableForViolinPlot
#----------------------------------------------------------------------------------------------------
# woratree requests separate violin based on her introduced ("N" vs "Thal") labelled cells
# to minimize code changes and clutter, my current solution is to simply add an optional
# matrix argument to x$createTableForViolinPlot(clusters=3, marker="H3", matrix=mtx.low)
# the returned tbl.violin, from successive calls, rbind together, but DO take care
# to append, e.g., "N", "Thal", to the name column of the tbl.violin before doing the rbind.
# all this is demonstrated below.
test_createTableForViolinPlot_matrixSupplied <- function()
{
    message(sprintf("--- test_createTableForViolinPlot_matrixSupplied"))

        #--------------------------------------
        # first test and plot: H3 in cluster 3
        #--------------------------------------

    f <- system.file(package="CytofkitNormalization", "extdata", "cytofkit-leukemia.RData")
    checkTrue(file.exists(f))
    x <- CytofkitNormalization$new(f)

    x$createSimpleMarkerNames()
    markers <- x$getMarkers()
    clusters <- x$getClusters()
    mtx <- x$getMatrix()
    low.h3 <- which(mtx[, markers["H3"]] < 2)
    high.h3 <- which(mtx[, markers["H3"]] >= 2)
    length(low.h3)   # 34813
    length(high.h3)  # 63406
    mtx.low <- mtx[low.h3,]
    mtx.high <- mtx[high.h3,]

    tbl.violin.low <- x$createTableForViolinPlot(clusters=3, marker="H3", mtx.low)
    checkEquals(length(which(is.na(tbl.violin.low$value))), 0)
    checkEquals(dim(tbl.violin.low), c(453, 2))
    if(interactive()){
        ggplot(tbl.violin.low, aes(x=name, y=value, fill=name)) + geom_violin() +
               theme(axis.text = element_text(size = 14)) +
               ggtitle("H3 - cluster 3")
        } # interactive

    tbl.violin.high <- x$createTableForViolinPlot(clusters=3, marker="H3", mtx.high)
    checkEquals(dim(tbl.violin.high), c(2162, 2))
    checkEquals(length(which(is.na(tbl.violin.high$value))), 0)
    if(interactive()){
        ggplot(tbl.violin.high, aes(x=name, y=value, fill=name)) + geom_violin() +
               theme(axis.text = element_text(size = 14)) +
               ggtitle("H3 - cluster 3")
        } # interactive

    tbl.violin.low$name <- paste0(tbl.violin.low$name, ".low")
    tbl.violin.high$name <- paste0(tbl.violin.high$name, ".high")
    tbl.violin.both <- rbind(tbl.violin.low, tbl.violin.high)
    checkEquals(sort(unique(tbl.violin.both$name)), c("H3.c3.high", "H3.c3.low"))

    if(interactive()){
        ggplot(tbl.violin.both, aes(x=name, y=value, fill=name)) + geom_violin() +
               theme(axis.text = element_text(size = 14)) +
               ggtitle("H3 - cluster 3")
        } # interactive

} # test_createTableForViolinPlot_matrixSupplied
#----------------------------------------------------------------------------------------------------
test_calculateColorBoundaries <- function()
{
    message(sprintf("--- test_calculateColorBoundaries"))

    f <- system.file(package="CytofkitNormalization", "extdata", "cytofkit-leukemia.RData")
    checkTrue(file.exists(f))
    x <- CytofkitNormalization$new(f)

       #------------------------------------------------------------
       # the simplest case: 10 integers, uniformly distributed,
       # 2 colors, mode="quantile"
       #
       #   start  end color
       # 1   1.0  5.5 black
       # 2   5.5 10.0 white
       #------------------------------------------------------------

    colors <- c("black", "white")
    vec <- 1:10
    tbl <- x$calculateColorBoundaries(vec, colors, mode="quantile")
    checkEquals(dim(tbl), c(2,3))
    checkEquals(as.list(tbl[1,]), list(start=1, end=5.5, color="black"))
    checkEquals(as.list(tbl[2,]), list(start=5.5, end=10, color="white"))
    checkEquals(min(tbl$start), min(vec))
    checkEquals(max(tbl$end), max(vec))

       #------------------------------------------------------------
       # normal distribution, 10 numbers, clumped at the mean
       # 2 colors, mode="quantile"
       #------------------------------------------------------------

    set.seed(17)
    vec <- sort(rnorm(n=10, mean=5, sd=2))
    tbl <- x$calculateColorBoundaries(vec, colors, mode="quantile")
    checkEquals(dim(tbl), c(2,3))
      # low half of vector should be black, top half white,
      # no overlaps, nothing missing
     black.values <-  vec[vec>=tbl$start[1] & vec < tbl$end[1]]
     white.values <-  vec[vec>=tbl$start[2] & vec <= tbl$end[2]]
     checkTrue(all(black.values < 5))
     checkTrue(all(white.values >= 5))

       #------------------------------------------------------------
       # normal distribution, 100 numbers, clumped at the mean
       # 4 colors, mode="quantile", should be +/1 equal number of
       # numbers for each color
       #------------------------------------------------------------

    set.seed(31)
    vec <- sort(rnorm(n=100, mean=50, sd=15))
    colors <- c("black", "darkgray", "lightgray", "white")
    tbl <- x$calculateColorBoundaries(vec, colors, mode="quantile")

    black.values <-  vec[vec>=tbl$start[1] & vec < tbl$end[1]]
    darkGray.values <-  vec[vec>=tbl$start[2] & vec <= tbl$end[2]]
    lightGray.values <-  vec[vec>=tbl$start[3] & vec <= tbl$end[3]]
    white.values <-  vec[vec>=tbl$start[4] & vec <= tbl$end[4]]

    checkEqualsNumeric(length(black.values), 25, tol=1)
    checkEqualsNumeric(length(darkGray.values), 25, tol=1)
    checkEqualsNumeric(length(lightGray.values), 25, tol=1)
    checkEqualsNumeric(length(white.values), 25, tol=1)
         # mode=="quantile" puts an even number of elements into each bin
         # but bins will be of very different size when the distribution
         # of the elements is non-uniform.  test that
    checkEquals(round(with(tbl, end-start)), c(28, 9, 10, 25))

       #------------------------------------------------------------
       # normal distribution, 100 numbers, clumped at the mean
       # 4 colors, mode="interval".
       # with evenly sized bins, the elements-per-bin will vary
       #------------------------------------------------------------
    tbl <- x$calculateColorBoundaries(vec, colors, mode="interval")
    black.values <-  vec[vec>=tbl$start[1] & vec < tbl$end[1]]
    darkGray.values <-  vec[vec>=tbl$start[2] & vec <= tbl$end[2]]
    lightGray.values <-  vec[vec>=tbl$start[3] & vec <= tbl$end[3]]
    white.values <-  vec[vec>=tbl$start[4] & vec <= tbl$end[4]]

       #--------------------------------------------
       # note: very uneven distribution across bins
       #--------------------------------------------

    checkEquals(length(black.values), 8)
    checkEquals(length(darkGray.values), 39)
    checkEquals(length(lightGray.values), 42)
    checkEquals(length(white.values), 11)

    checkEquals(min(vec), tbl$start[1])
    checkEquals(max(vec), tbl$end[nrow(tbl)])

} # test_calculateColorBoundaries
#----------------------------------------------------------------------------------------------------
test_createTableForTsnePlot <- function()
{
    message(sprintf("--- test_createTableForTsnePlot"))

    #cytofkit.results.file <- "cytofkit-leukemia.RData"
    cytofkit.results.file <- "exp54-NvsThal.RData"
    f <- system.file(package="CytofkitNormalization", "extdata", cytofkit.results.file)
    checkTrue(file.exists(f))
    x <- CytofkitNormalization$new(f)

    x$createSimpleMarkerNames()
    markers <- x$getMarkers()
    marker <- "H3K9cr"
    checkTrue(marker %in% names(markers))
    vec <- x$getMatrix()[, markers[[marker]]]
    round(fivenum(as.numeric(vec)), digits=3)
      # -0.007  0.948  1.415  1.846  4.237
    require(RColorBrewer)
    colors <- rev(brewer.pal(11, "Spectral"))
     # colors <- c("black", "darkgray", "lightgray", "white")
    boundary.mode <- "quantile"
    boundary.mode <- "interval"

    tbl.colors <- x$calculateColorBoundaries(vec, colors, mode=boundary.mode)
    tbl.tsneWithColors <- x$createTableForTsnePlot("H3K9cr", tbl.colors, matrix.sub=NA)

    if(interactive()){
       quartz(width=10, height=10)
       par(mar=c(5,5,5,5)) # generous margins
       #tbl.show <- subset(tbl.tsneWithColors, value > 2.0)
       tbl.show <- tbl.tsneWithColors
       with(tbl.show, plot(tsne_1, tsne_2, col=color,
                           main=sprintf("%s (%s colors) %s", marker, boundary.mode, cytofkit.results.file),
                           cex=0.4, pch=16,
                       ylim=c(-50, 50), xlim=c(-50, 50)))

       # starts <- tbl.colors$start[seq_len(nrow(tbl.colors))[-1]]
       # starts <- c(starts, tail(tbl.colors$end, n=1))
       # boundaries <- round(tbl.colors$end, digits=2)
       boundaries <- round(tbl.colors$end, digits=2)
       legend(35, 40, round(rev(boundaries), digits=3), rev(tbl.colors$color))
       } # interactive

} # test_createTableForTsnePlot
#----------------------------------------------------------------------------------------------------
test_test_reproduce.known.plot <- function()
{
   message(sprintf("--- test_reproduce.known.plot"))

   f <- system.file(package="CytofkitNormalization", "extdata", "cytofkit-leukemia.RData")
   f <- "../extdata/cytofkit-all-abs.RData"
   checkTrue(file.exists(f))
   x <- CytofkitNormalization$new(f)

   x$createSimpleMarkerNames()
   target <- "H3K9Ac"
   h3.reference <- "H3"
   new.col.name <- x$normalizeMarker(target, c(h3.reference)) #, h4.reference))
   checkEquals(new.col.name, "H3K9Ac.regress.H3")

   markers <- x$getMarkers()
   clusters <- x$getClusters()

       #-----------------------------------------------------------
       # show clusters 1 and 3, nominated by marjorie, for H3K9Ac
       # first with raw data, then regression normalized against H3
       #-----------------------------------------------------------

   tbl.violin.01 <- x$createTableForViolinPlot(clusters=c(1,3), marker=target)
   tbl.violin.02 <- x$createTableForViolinPlot(clusters=c(1,3), marker=new.col.name)
   tbl.violin <- rbind(tbl.violin.01, tbl.violin.02)

   if(interactive()){
       ggplot(tbl.violin, aes(x=name, y=value, fill=name)) + geom_violin() +
            theme(axis.text = element_text(size = 14)) +
            ggtitle(sprintf("%s", new.col.name))
        }



} # test_reproduce.known.plot
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
