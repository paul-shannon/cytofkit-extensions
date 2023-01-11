library(RUnit)
library(CytofkitNormalization)
library(ggplot2)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_markerNames()
    test_normalizeMarker()
    test_getCluster()
    test_createTableForViolinPlot()
    test_createTableForViolinPlot_matrixSupplied()

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

    f <- system.file(package="CytofkitNormalization", "extdata", "cytofkit-leukemia.RData")
    checkTrue(file.exists(f))
    x <- CytofkitNormalization$new(f)

    x$createSimpleMarkerNames()
    markers <- x$getMarkers()
    clusters <- x$getClusters()

    tbl.violin <- x$createTableForViolinPlot(clusters=3, marker="H3")
    checkEquals(dim(tbl.violin), c(2615, 2))
    if(interactive()){

        ggplot(tbl.violin, aes(x=name, y=value, fill=name)) + geom_violin() +
            #coord_cartesian(ylim=c(-5,5)) +
            theme(axis.text = element_text(size = 14)) +
            ggtitle("H3 - cluster 3")
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
    checkEquals(dim(tbl.violin.low), c(2615, 2))
    if(interactive()){
        ggplot(tbl.violin.low, aes(x=name, y=value, fill=name)) + geom_violin() +
               theme(axis.text = element_text(size = 14)) +
               ggtitle("H3 - cluster 3")
        } # interactive

    tbl.violin.high <- x$createTableForViolinPlot(clusters=3, marker="H3", mtx.high)
    checkEquals(dim(tbl.violin.high), c(2615, 2))
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
test_reproduce.known.plot <- function()
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
