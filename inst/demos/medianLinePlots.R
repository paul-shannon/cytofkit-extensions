library(RUnit)
library(CytofkitNormalization)
#-----------------------------------------------------------------------------------------------------------------
runDemos <- function()
{
  demo_getNormalizedMarkerMediansCollapsedAcrossClusters()
  demo_getMedianRatiosByCluster.h3k9cr.vs.h3k9ac()
  demo_getMedianRatiosByCluster.h3k9cr.thal.vs.norm()
  demo_getNormalizedMarkerRatiosByClusterBySample()

} # runDemos
#-----------------------------------------------------------------------------------------------------------------
# matrix.filter.rownames.string, with our current data, should be either "-N_" or "-Thal_"
getNormalizedMarkerMediansCollapsedAcrossClusters <- function (cytofkitNormalizer,
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

} # getNormalizedMarkerMediansCollapsedAcrossClusters
#-----------------------------------------------------------------------------------------------------------------
demo_getNormalizedMarkerMediansCollapsedAcrossClusters <- function()
{
   message(sprintf("--- demo_getNormalizedMarkerMediansCollapsedAcrossClusters"))

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
   vec.h3k9me3 <- getNormalizedMarkerMediansCollapsedAcrossClusters(x, target, ref.markers,clustersOfInterest)
   checkEquals(as.numeric(names(vec.h3k9me3)), clustersOfInterest)

        #---------------------
        # now H3K9me3 Normal
        #---------------------

   target <-   "H3K9me3"
   h3.reference <- "H3"
   h4.reference <- "H4"
   ref.markers <- c(h3.reference, h4.reference)
   checkTrue(all(c(target, h3.reference, h4.reference) %in% names(markers)))
   vec.h3k9me3.norm <- getNormalizedMarkerMediansCollapsedAcrossClusters(x, target, ref.markers,
                                                                clustersOfInterest,
                                                                matrix.rowname.filter="-N_")

        #---------------------
        # now H3K9me3 Thal
        #---------------------

   vec.h3k9me3.thal <- getNormalizedMarkerMediansCollapsedAcrossClusters(x, target, ref.markers,
                                                                clustersOfInterest,
                                                                matrix.rowname.filter="-Thal_")

        #-----------------
        # now H3K9ac
        #-----------------

   target <-   "H3K27ac"
   h3.reference <- "H3"
   h4.reference <- "H4"
   ref.markers <- c(h3.reference, h4.reference)
   checkTrue(all(c(target, h3.reference, h4.reference) %in% names(markers)))
   vec.h3k27ac <- getNormalizedMarkerMediansCollapsedAcrossClusters(x, target, ref.markers,
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

} # demo_getNormalizedMarkerMediansCollapsedAcrossClusters
#-----------------------------------------------------------------------------------------------------------------
# email from woratree (8 feb 2023)
#   The next task would be generating a violin plot and also this
#   median line plot but the data is the ratio of 2 histone marks of
#   interest (after normalization). For example the ratio of
#   H3K9cr/H3K9ac, H3K27cr/H3K27ac.
demo_getMedianRatiosByCluster.h3k9cr.vs.h3k9ac <- function()
{
   message(sprintf("--- demo_getMedianRatiosByCluster.h3k9cr.vs.h3k9ac"))

   clustersOfInterest <- c(1,7,3,8,5,6,12,15,14,16,13,17,18)  # missing 2, 4, 9, 10 , 11
   f <- system.file(package="CytofkitNormalization", "extdata", "exp54-NvsThal.RData")
   checkTrue(file.exists(f))
   x <- CytofkitNormalization$new(f)
   x$createSimpleMarkerNames()
   markers <- x$getMarkers()

        #-----------------
        # H3K9cr/H3K9Ac
        #-----------------

   target.1 <-   "H3K9cr"
   target.2 <-   "H3K9Ac"
   h3.reference <- "H3"
   h4.reference <- "H4"
   ref.markers <- c(h3.reference, h4.reference)
   checkTrue(all(c(target.1, target.2, h3.reference, h4.reference) %in% names(markers)))

   vec.h3k9cr <- getNormalizedMarkerMediansCollapsedAcrossClusters(x, target.1, ref.markers,
                                                          clustersOfInterest)
   vec.h3k9ac <- getNormalizedMarkerMediansCollapsedAcrossClusters(x, target.2, ref.markers,
                                                          clustersOfInterest)

   vec.ratios <- as.numeric(vec.h3k9cr)/as.numeric(vec.h3k9ac)

      #---------------------------------------------
      # this data.frame enables easy inspection, QC
      #---------------------------------------------

   tbl <- data.frame(h3k9cr=as.numeric(vec.h3k9cr),
                     h3k9ac=as.numeric(vec.h3k9ac),
                     ratio=as.numeric(vec.ratios))
   tbl <- round(tbl, digits=3)
   tbl$cluster <- clustersOfInterest

      #---------------------------------------------------------
      # now plot.  note that the very large spikes, + and -
      # at clusters 6 and 7 are artifacts of small denominators
      # one negative, one positive, and of doubtful biological
      # import
      #---------------------------------------------------------

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


} # demo_getMedianRatiosByCluster.h3k9cr.vs.h3k9ac
#-----------------------------------------------------------------------------------------------------------------
demo_getMedianRatiosByCluster.h3k9cr.thal.vs.norm <- function()
{
   message(sprintf("--- demo_getMedianRatiosByCluster.h3k9cr.thal.vs.norm"))

   clustersOfInterest <- c(1,7,3,8,5,6,12,15,14,16,13,17,18)  # missing 2, 4, 9, 10 , 11
   f <- system.file(package="CytofkitNormalization", "extdata", "exp54-NvsThal.RData")
   checkTrue(file.exists(f))
   x <- CytofkitNormalization$new(f)
   x$createSimpleMarkerNames()
   markers <- x$getMarkers()

        #-----------------
        # normal cells
        #-----------------

   target <-   "H3K9Ac"
   h3.reference <- "H3"
   h4.reference <- "H4"
   ref.markers <- c(h3.reference, h4.reference)
   checkTrue(all(c(target, h3.reference, h4.reference) %in% names(markers)))
   vec.h3k9ac.norm <- getNormalizedMarkerMediansCollapsedAcrossClusters(x, target, ref.markers,
                                                               clustersOfInterest,
                                                               matrix.rowname.filter="-N_")

        #---------------------
        # thalassemic cells
        #---------------------

   vec.h3k9ac.thal <- getNormalizedMarkerMediansCollapsedAcrossClusters(x, target, ref.markers,
                                                               clustersOfInterest,
                                                               matrix.rowname.filter="-Thal_")
   vec.ratios <- as.numeric(vec.h3k9ac.thal)/as.numeric(vec.h3k9ac.norm)

      #---------------------------------------------
      # a data.frame for easy inspection, QC
      #---------------------------------------------

   tbl <- data.frame(h3k9ac.norm=as.numeric(vec.h3k9ac.norm),
                     h3k9ac.thal=as.numeric(vec.h3k9ac.thal),
                     ratio=as.numeric(vec.ratios))
   tbl <- round(tbl, digits=3)
   fivenum(as.matrix(tbl))  # -23.6600  -0.0575   0.1120   0.3630   7.5250

   tbl$cluster <- clustersOfInterest

      #---------------------------------------------------------
      # now plot.  note that the very large spikes, + and -
      # at clusters 6 and 7 are artifacts of small denominators
      # one negative, one positive, and of doubtful biological
      # import
      #---------------------------------------------------------

   plot(as.numeric(vec.ratios), type="b", ylim=c(-30,10), pch=16,
        xlab="Cluster", ylab="median expression ratios",
        xaxt="n", main="H3K9Ac THAL/NORM Normalized Median Expression by Cluster")
   axis(1,
        at=seq_len(length(clustersOfInterest)),
        labels=sprintf("%s", clustersOfInterest),
        col.axis="black", las=0)

   lines(as.numeric(vec.h3k9ac.thal), type="b", col="red", pch=16)
   lines(as.numeric(vec.h3k9ac.norm), type="b", col="darkgreen", pch=16)

   legend(9, -15,
          c("H3K9Ac.norm", "H3K9Ac.thal", "H3K9Ac.thal/H3K9Ac.norm"),
          c("red", "darkgreen", "black"))

} # demo_getMedianRatiosByCluster.h3k9cr.thal.vs.norm
#--------------------------------------------------------------------------------------------------------------
# matrix.filter.rownames.string, with our current data, should be either "-N_" or "-Thal_"
getNormalizedMarkerValuesByCluster <- function (cytofkitNormalizer,
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
   invisible(tbl.violin)

} # getNormalizedMarkerValuesByClusters
#-----------------------------------------------------------------------------------------------------------------
demo_getNormalizedMarkerRatiosByClusterBySample <- function()
{
   clustersOfInterest <- c(1,7,3,8,5,6,12,15,14,16,13,17,18) # missing 2, 4, 9, 10 , 11
   f <- system.file(package="CytofkitNormalization", "extdata", "exp54-NvsThal.RData")
   checkTrue(file.exists(f))
   x <- CytofkitNormalization$new(f)
   x$createSimpleMarkerNames()
   markers <- x$getMarkers()

   target.1 <-   "H3K9cr"
   target.2 <-   "H3K9Ac"
   h3.reference <- "H3"
   h4.reference <- "H4"
   ref.markers <- c(h3.reference, h4.reference)
   checkTrue(all(c(target.1, target.2, h3.reference, h4.reference) %in% names(markers)))


        #-----------------
        # first target.1
        #-----------------

   tbl.1 <- getNormalizedMarkerValuesByCluster(x, target.1, c(h3.reference, h4.reference),
                                               clustersOfInterest,
                                               matrix.rowname.filter=NA)

   tbl.2 <- getNormalizedMarkerValuesByCluster(x, target.2, c(h3.reference, h4.reference),
                                               clustersOfInterest,
                                               matrix.rowname.filter=NA)

        #--------------------------------------------------------
        # guarantee that the tables have compatible rows
        #--------------------------------------------------------
   checkEquals(sub(target.1, "", tbl.1$name),
               sub(target.2, "", tbl.2$name))
   clusters <- sub("H3K9cr.regress.H3+H4.","", tbl.1$name, fixed=TRUE)
   tbl.ratios <- data.frame(target.1=tbl.1$value, target.2=tbl.2$value,
                            cluster=clusters)
   tbl.ratios$ratio <- with(tbl.ratios, target.1/target.2)
   dim(tbl.ratios)     # 77218 4

   tbls.byCluster <- list()

   for(cluster in clustersOfInterest){
      cluster.name <- paste0("c", cluster)
      ratios <- subset(tbl.ratios, cluster==cluster.name)$ratio
      values.1 <- subset(tbl.ratios, cluster==cluster.name)$target.1
      values.2 <- subset(tbl.ratios, cluster==cluster.name)$target.2
      tbls.byCluster[[cluster.name]] <- data.frame(cluster=cluster.name,
                                                   value.1=median(values.1),
                                                   value.2=median(values.2),
                                                   ratioBySample=median(ratios))
      } # for cluster
   tbl.byCluster <- do.call(rbind, tbls.byCluster)
   tbl.byCluster$ratioByCluster <- with(tbl.byCluster, value.1/value.2)

   tbl.byCluster$value.1 <- round(tbl.byCluster$value.1, digits=3)
   tbl.byCluster$value.2 <- round(tbl.byCluster$value.2, digits=3)
   tbl.byCluster$ratioBySample <- round(tbl.byCluster$ratioBySample, digits=3)
   tbl.byCluster$ratioByCluster <- round(tbl.byCluster$ratioByCluster, digits=3)

   tbl.byCluster$cluster <- as.numeric(sub("c", "", tbl.byCluster$cluster))
   dim(tbl.byCluster)

   plot.title <- sprintf("%s/%s median expression by cluster", target.1, target.2)

   vec <- tbl.byCluster$ratioByCluster
   plot(vec, type="b", ylim=c(-3,3), pch=16,
        xlab="Cluster", ylab="median expression ratios",
        xaxt="n", main=plot.title)
   axis(1,
        at=seq_len(length(clustersOfInterest)),
        labels=sprintf("%s", clustersOfInterest),
        col.axis="black", las=0)

   vec <- tbl.byCluster$ratioBySample
   lines(vec, type="b", col="red", lwd=5, pch=16)
   vec <- rep(0, length(clustersOfInterest))
   lines(vec, col="gray")

   vec <- tbl.byCluster$value.1
   lines(vec, type="b", col="darkgreen", pch=16)
   vec <- tbl.byCluster$value.2
   lines(vec, type="b", col="darkblue", pch=16)

   legend(10, 2.8,
          c("ratio by sample", "ratio by cluster",
            sprintf("medians of %s", target.1),
            sprintf("medians of %s", target.2)),
          c("red", "black", "darkgreen", "darkblue"))


} # demo_getNormalizedMarkerRatiosByClusterBySample
#-----------------------------------------------------------------------------------------------------------------
demo_getRawMarkerRatiosByClusterBySample <- function()
{
   clustersOfInterest <- c(1,7,3,8,5,6,12,15,14,16,13,17,18) # missing 2, 4, 9, 10 , 11
   f <- system.file(package="CytofkitNormalization", "extdata", "exp54-NvsThal.RData")
   checkTrue(file.exists(f))
   x <- CytofkitNormalization$new(f)
   x$createSimpleMarkerNames()
   markers <- x$getMarkers()

   target.1 <-   "H3K9cr"
   target.2 <-   "H3K9Ac"
   h3.reference <- "H3"
   h4.reference <- "H4"
   ref.markers <- c(h3.reference, h4.reference)
   checkTrue(all(c(target.1, target.2, h3.reference, h4.reference) %in% names(markers)))


        #-----------------
        # first target.1
        #-----------------

   tbl.1 <- x$createTableForViolinPlot(clustersOfInterest, marker=target.1)
   tbl.2 <- x$createTableForViolinPlot(clustersOfInterest, marker=target.2)

        #--------------------------------------------------------
        # guarantee that the tables have compatible rows
        #--------------------------------------------------------

   checkEquals(sub(target.1, "", tbl.1$name),
               sub(target.2, "", tbl.2$name))
   clusters <- sub("H3K9cr.","", tbl.1$name, fixed=TRUE)
   tbl.ratios <- data.frame(target.1=tbl.1$value, target.2=tbl.2$value, cluster=clusters)
   tbl.ratios$ratio <- with(tbl.ratios, target.1/target.2)
   dim(tbl.ratios)     # 77218 4

   tbls.byCluster <- list()

   for(cluster in clustersOfInterest){
      cluster.name <- paste0("c", cluster)
      ratios <- subset(tbl.ratios, cluster==cluster.name)$ratio
      values.1 <- subset(tbl.ratios, cluster==cluster.name)$target.1
      values.2 <- subset(tbl.ratios, cluster==cluster.name)$target.2
      tbls.byCluster[[cluster.name]] <- data.frame(cluster=cluster.name,
                                                   value.1=median(values.1),
                                                   value.2=median(values.2),
                                                   ratioBySample=median(ratios))
      } # for cluster

   tbl.byCluster <- do.call(rbind, tbls.byCluster)
   tbl.byCluster$ratioByCluster <- with(tbl.byCluster, value.1/value.2)

   tbl.byCluster$value.1 <- round(tbl.byCluster$value.1, digits=3)
   tbl.byCluster$value.2 <- round(tbl.byCluster$value.2, digits=3)
   tbl.byCluster$ratioBySample <- round(tbl.byCluster$ratioBySample, digits=3)
   tbl.byCluster$ratioByCluster <- round(tbl.byCluster$ratioByCluster, digits=3)

   tbl.byCluster$cluster <- as.numeric(sub("c", "", tbl.byCluster$cluster))
   checkEquals(dim(tbl.byCluster), c(13, 5))

      #----------------------------------------------------------
      # spot check the two strategies for calculating medians
      # first: reducing all samples to their median, then divide
      #----------------------------------------------------------
   target1.c1.median <- median(subset(tbl.ratios, cluster=="c1")$target.1)
   target2.c1.median <- median(subset(tbl.ratios, cluster=="c1")$target.2)
   c1.median.by.cluster <- target1.c1.median/target2.c1.median
   checkEqualsNumeric(c1.by.cluster, 0.4303003, tolerance=1e-6)

      #----------------------------------------------------------
      # second, get ratio of every sample pair, then find median
      #----------------------------------------------------------
   ratios <- with(subset(tbl.ratios, cluster=="c1"), target.1/target.2)
   c1.median.by.sample <- median(ratios)
   checkEqualsNumeric(c1.median.by.sample, 0.4358723, tolerance=1e-6)

   plot.title <- sprintf("%s/%s median expression by cluster", target.1, target.2)

   vec <- tbl.byCluster$ratioByCluster
   plot(vec, type="b", ylim=c(-3,4), pch=16,
        xlab="Cluster", ylab="median expression ratios",
        xaxt="n", main=plot.title)
   axis(1,
        at=seq_len(length(clustersOfInterest)),
        labels=sprintf("%s", clustersOfInterest),
        col.axis="black", las=0)

   vec <- tbl.byCluster$ratioBySample
   lines(vec, type="b", col="red", lwd=1, pch=16)
   vec <- rep(0, length(clustersOfInterest))
   lines(vec, col="gray")

   vec <- tbl.byCluster$value.1
   lines(vec, type="b", col="darkgreen", pch=16)
   vec <- tbl.byCluster$value.2
   lines(vec, type="b", col="darkblue", pch=16)

   legend(10, -1.0,
          c("ratio by sample", "ratio by cluster",
            sprintf("medians of %s", target.1),
            sprintf("medians of %s", target.2)),
          c("red", "black", "darkgreen", "darkblue"))

} # demo_getRawMarkerRatiosByClusterBySample
#-----------------------------------------------------------------------------------------------------------------
