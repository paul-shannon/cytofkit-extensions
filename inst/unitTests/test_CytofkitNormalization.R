library(RUnit)
library(CytofkitNormalization)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_markerNames()
    test_normalizeMarker()

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

    target <- markers[["H3K4me3"]]
    h3.reference <- markers[["H3"]]
    h4.reference <- markers[["H4"]]

    target <- "H3K4me3"
    h3.reference <- "H3"
    h4.reference <- "H4"

    new.col.name <- x$normalizeMarker(target, h3.reference)
    checkEquals(new.col.name, "H3K4me3.regress.H3")

      # a hasty un-nuanced check is to get the quartiles of
      # of the distribution before and after normaliaztion,
      # check their sums

    normalized.vec <- x$getMatrix()[, new.col.name]
    quick.assay.normalized <- sum(as.numeric(fivenum(normalized.vec)))
       # -4.32963408 -0.25168910  0.01600153  0.29289136  1.89312134
    checkTrue(quick.assay.normalized < -2)
    raw.vec <- x$getMatrix()[, markers[["H3K4me3"]]]
    quick.assay.raw <- sum(as.numeric(fivenum(raw.vec)))
       #  -0.006413926  3.504356736  3.804160211  4.115765684  5.413913164
    checkTrue(quick.assay.raw > 16)

} # test_normalizeMarker
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
