library(RUnit)
library(CytofkitNormalization)


#-----------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  #test_getNormalizedMarkerMediansAcrossClusters()
  test_getRatioOfNormalizedMarkerMediansAcrossClusters()

} # runTests
#-----------------------------------------------------------------------------------------------------------------
# matrix.filter.rownames.string, with our current data, should be either "-N_" or "-Thal_"
getNormalizedMarkerMediansAcrossClusters <- function (cytofkitNormalizer,
                                                      targetMarker,
                                                      referenceMarkers,
                                                      clustersOfInterest,
                                                      matrix.rowname.filter=NA)
{
   stopifnot(matrix.rowname.filter %in% c(NA, "-N_", "-Thal_"))
   markers <- cytofkitNormalizer$getMarkers()
   clusters <- cytofkitNormalizer$getClusters()

   stopifnot(targetMarker %in% names(markers))
   stopifnot(all(referenceMarkers %in% names(markers)))
   normalizedTargetName <- cytofkitNormalizer$normalizeMarker(targetMarker, referenceMarkers)

   mtx <- cytofkitNormalizer$getMatrix()[, normalizedTargetName, drop=FALSE]
   if(!is.na(matrix.rowname.filter)){
      filtered.row.names <- grep(matrix.rowname.filter, rownames(mtx), fixed=TRUE)
      mtx <- mtx[filtered.row.names,,drop=FALSE]
      }

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
test_getNormalizedMarkerMediansAcrossClusters <- function()
{
   message(sprintf("--- test_getNormalizedMarkerMediansAcrossClusters"))

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
   vec.h3k9me3 <- getNormalizedMarkerMediansAcrossClusters(x, target, ref.markers,clustersOfInterest)
   checkEquals(as.numeric(names(vec.h3k9me3)), clustersOfInterest)

        #---------------------
        # now H3K9me3 Normal
        #---------------------

   target <-   "H3K9me3"
   h3.reference <- "H3"
   h4.reference <- "H4"
   ref.markers <- c(h3.reference, h4.reference)
   checkTrue(all(c(target, h3.reference, h4.reference) %in% names(markers)))
   vec.h3k9me3.norm <- getNormalizedMarkerMediansAcrossClusters(x, target, ref.markers,clustersOfInterest,
                                                                matrix.rowname.filter="-N_")

        #---------------------
        # now H3K9me3 Thal
        #---------------------

   vec.h3k9me3.thal <- getNormalizedMarkerMediansAcrossClusters(x, target, ref.markers,clustersOfInterest,
                                                                matrix.rowname.filter="-Thal_")

        #-----------------
        # now H3K9ac
        #-----------------

   target <-   "H3K27ac"
   h3.reference <- "H3"
   h4.reference <- "H4"
   ref.markers <- c(h3.reference, h4.reference)
   checkTrue(all(c(target, h3.reference, h4.reference) %in% names(markers)))
   vec.h3k27ac <- getNormalizedMarkerMediansAcrossClusters(x, target, ref.markers,clustersOfInterest)
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
   lines(as.numeric(vec.h3k9me3.norm), type="b", col="black", pch=16)
   lines(as.numeric(vec.h3k9me3.thal), type="b", col="red", pch=16)
   lines(as.numeric(vec.h3k27ac), type="b", col="darkgreen", pch=16)

   axis(1,
        at=seq_len(length(clustersOfInterest)),
        labels=sprintf("%s", clustersOfInterest),
        col.axis="black", las=0)
   legend(9, 2,
          c("H3K9me3 all", "H3K9me3.norm", "H3K9me3.thal", "H3K27ac all"),
          c("blue", "black", "red", "darkgreen"))

} # test_getNormalizedMarkerMediansAcrossClusters
#-----------------------------------------------------------------------------------------------------------------
# matrix.filter.rownames.string, with our current data, should be
# either "-N_" or "-Thal_"
getRatioOfNormalizedMarkerMediansAcrossClusters <- function (cytofkitNormalizer,
                                                           targetMarker.1,
                                                           targetMarker.2,
                                                           referenceMarkers,
                                                           clustersOfInterest,
                                                           matrix.rowname.filter=NA)
{
   stopifnot(matrix.rowname.filter %in% c(NA, "-N_", "-Thal_"))
   markers <- cytofkitNormalizer$getMarkers()
   clusters <- cytofkitNormalizer$getClusters()

   stopifnot(targetMarker.1 %in% names(markers))
   stopifnot(targetMarker.2 %in% names(markers))
   stopifnot(all(referenceMarkers %in% names(markers)))
   normalizedTargetName.1 <- cytofkitNormalizer$normalizeMarker(targetMarker.1, referenceMarkers)
   normalizedTargetName.2 <- cytofkitNormalizer$normalizeMarker(targetMarker.2, referenceMarkers)

   colnames <- c(normalizedTargetName.1, normalizedTargetName.2)
   mtx <- cytofkitNormalizer$getMatrix()[, colnames, drop=FALSE]
   ratio <- mtx[, colnames[1]] / mtx[, colnames[2]]
   mtx <- cbind(mtx, ratio)
      # mtx access is via colnames, where the original names are long and unwieldy,
      # for which we have provided shorter, human-friendly names.
      # colnames for newly calculated columns are often (always) short, and shared
   cytofkitNormalizer$addMarker("ratio", "ratio")

   if(!is.na(matrix.rowname.filter)){
      filtered.row.names <- grep(matrix.rowname.filter, rownames(mtx), fixed=TRUE)
      mtx <- mtx[filtered.row.names,,drop=FALSE]
      }

   tbl.violin <- cytofkitNormalizer$createTableForViolinPlot(clustersOfInterest,
                                                             marker="ratio", #normalizedTargetName,
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

} # getRatioOfNormalizedMarkerMediansAcrossClusters
#-----------------------------------------------------------------------------------------------------------------
# email from woratree (8 feb 2023)
#   The next task would be generating a violin plot and also this
#   median line plot but the data is the ratio of 2 histone marks of
#   interest (after normalization). For example the ratio of
#   H3K9cr/H3K9ac, H3K27cr/H3K27ac.
test_getRatioOfNormalizedMarkerMediansAcrossClusters <- function()
{
   message(sprintf("--- test_getRatioOfNormalizedMarkerMediansAcrossClusters"))

   clustersOfInterest <- c(1,7,3,8,5,6,12,15,14,16,13,17,18)  # missing 2, 4, 9, 10 , 11
   f <- system.file(package="CytofkitNormalization", "extdata", "exp54-NvsThal.RData")
   checkTrue(file.exists(f))
   x <- CytofkitNormalization$new(f)
   x$createSimpleMarkerNames()
   markers <- x$getMarkers()

        #-----------------
        # first H3K9me3
        #-----------------

   target.1 <-   "H3K9cr"
   target.2 <-   "H3K9Ac"
   h3.reference <- "H3"
   h4.reference <- "H4"
   ref.markers <- c(h3.reference, h4.reference)
   checkTrue(all(c(target.1, target.2, h3.reference, h4.reference) %in% names(markers)))

   vec.h3k9cr <- getNormalizedMarkerMediansAcrossClusters(x, target.1, ref.markers,
                                                          clustersOfInterest)
   vec.h3k9ac <- getNormalizedMarkerMediansAcrossClusters(x, target.2, ref.markers,
                                                          clustersOfInterest)

   vec.ratios <- as.numeric(vec.h3k9cr)/as.numeric(vec.h3k9ac)

      # getRatioOfNormalizedMarkerMediansAcrossClusters(x, target.1, target.2,
      # ref.markers, clustersOfInterest)

      # for easy QC
   tbl <- data.frame(h3k9cr=as.numeric(vec.h3k9cr),
                     h3k9ac=as.numeric(vec.h3k9ac),
                     ratio=as.numeric(vec.ratios))
   tbl <- round(tbl, digits=3)
   tbl$cluster <- clustersOfInterest

   plot(as.numeric(vec.ratios), type="b", ylim=c(-3,3), pch=16,
        xlab="Cluster", ylab="median expression ratios",
        xaxt="n", main="Normalized Median Expression by Cluster")
   axis(1,
        at=seq_len(length(clustersOfInterest)),
        labels=sprintf("%s", clustersOfInterest),
        col.axis="black", las=0)

   lines(as.numeric(vec.h3k9cr), type="b", col="blue", pch=16)
   lines(as.numeric(vec.h3k9ac), type="b", col="darkgreen", pch=16)

   legend(9, -2,
          c("H3K9cr", "H3K9Ac", "H3K9cr/H3K9Ac"),
          c("blue", "darkgreen", "black"))


} # test_getRatioOfNormalizedMarkerMediansAcrossClusters
#-----------------------------------------------------------------------------------------------------------------
