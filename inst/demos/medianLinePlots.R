library(RUnit)
library(CytofkitNormalization)


#-----------------------------------------------------------------------------------------------------------------
getNormalizedMarkerMediansAcrossClusters <- function (cytofkitNormalizer,
                                                      targetMarker,
                                                      referenceMarkers,
                                                      clustersOfInterest)
{
   markers <- cytofkitNormalizer$getMarkers()
   clusters <- cytofkitNormalizer$getClusters()

   stopifnot(targetMarker %in% names(markers))
   stopifnot(all(referenceMarkers %in% names(markers)))
   normalizedTargetName <- cytofkitNormalizer$normalizeMarker(targetMarker, referenceMarkers)

   mtx <- cytofkitNormalizer$getMatrix()[, normalizedTargetName, drop=FALSE]
   dim(mtx)
   tbl.violin <- cytofkitNormalizer$createTableForViolinPlot(clustersOfInterest,
                                                             marker=normalizedTargetName,
                                                             matrix=mtx)

   medians <- list()
   for(cluster in clustersOfInterest){
      cluster.sig <- sprintf("\\.c%d$", cluster)
      tbl.cluster <- subset(tbl.violin, grepl(cluster.sig, tbl.violin$name, fixed=FALSE))
      new.value <- median(tbl.cluster$value)
      title <- sprintf("c.%d", cluster)
      medians[[title]] <- new.value
      } # for cluster

   names(medians) <- clustersOfInterest

   medians

} # getNormalizedMarkerMediansAcrossClusters
#-----------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_getNormalizedMarkerMediansAcrossClusters()

} # runTests
#-----------------------------------------------------------------------------------------------------------------
test_getNormalizedMarkerMediansAcrossClusters <- function()
{
   clustersOfInterest <- c(1,7,3,8,5,6,12,15,14,16,13,17,18) # missing 2, 4, 9, 10 , 11
   f <- system.file(package="CytofkitNormalization", "extdata", "exp54-NvsThal.RData")
   checkTrue(file.exists(f))
   x <- CytofkitNormalization$new(f)
   x$createSimpleMarkerNames()
   markers <- x$getMarkers()

        #-----------------
        # first H3K9me3
        #-----------------

   target <-   "H3K9me3"
   h3.reference <- "H3"
   h4.reference <- "H4"
   ref.markers <- c(h3.reference, h4.reference)

   checkTrue(all(c(target, h3.reference, h4.reference) %in% names(markers)))

   vec.h3k9me3 <- getNormalizedMarkerMediansAcrossClusters(x, target, ref.markers,
                                                           clustersOfInterest)
   checkEquals(as.numeric(names(vec.h3k9me3)), clustersOfInterest)

        #-----------------
        # now H3K9ac
        #-----------------

   target <-   "H3K27ac"
   h3.reference <- "H3"
   h4.reference <- "H4"
   ref.markers <- c(h3.reference, h4.reference)

   checkTrue(all(c(target, h3.reference, h4.reference) %in% names(markers)))

   vec.h3k27ac <- getNormalizedMarkerMediansAcrossClusters(x, target, ref.markers,
                                                           clustersOfInterest)
   checkEquals(as.numeric(names(vec.h3k9me3)), clustersOfInterest)

        #-----------------
        # now plot
        #-----------------

   title <- sprintf("median marker expression, normalized against %s",
                    paste(ref.markers, collapse=" & "))
   plot(as.numeric(vec.h3k9me3), type="b", col="blue", ylim=c(-2, 2), pch=16,
        main=title,
        xlab="Cluster", ylab="median expression",
        xaxt="n")
   lines(as.numeric(vec.h3k27ac), type="b", col="darkgreen", pch=16)
   axis(1,
        at=seq_len(length(clustersOfInterest)),
        labels=sprintf("%s", clustersOfInterest),
        col.axis="black", las=0)
   legend(10, 2, c("H3K9me3", "H3K27ac"), c("blue", "darkgreen"))

} # test_getNormalizedMarkerMediansAcrossClusters
#-----------------------------------------------------------------------------------------------------------------

