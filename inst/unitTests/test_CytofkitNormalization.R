library(RUnit)
library(CytofkitNormalization)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()

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
if(!interactive())
    runTests()
