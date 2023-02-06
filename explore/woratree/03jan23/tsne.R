library(CytofkitNormalization)
library(ggplot2)
library(RUnit)

f <- "exp54-NvsThal.RData"
checkTrue(file.exists(f))
x <- CytofkitNormalization$new(f)
x$createSimpleMarkerNames()
markers <- x$getMarkers()


mtx <- x$getMatrix()
dim(mtx)  # 98219    58

normal.rows <- grep("-N_", rownames(mtx), fixed=TRUE)
thalassemia.rows <- grep("-Thal_", rownames(mtx), fixed=TRUE)
checkEquals(nrow(mtx), sum(length(normal.rows), length(thalassemia.rows)))

mtx.normal <- mtx[normal.rows,]
mtx.thal   <- mtx[thalassemia.rows,]

   # drawing code from
   # ~~/github/TrenaProjectErythropoiesis/explore/cytof/cytofit.runs/cdMarkersOnly-cell5k-noNormalization-tsne/normalizeByCluster.R
   # see Cd114Di<114Cd_Hb_tot> Expression Level Plot in
   # Histone normalization test-RpackagePaul-WK_3Jan2023.pptx, slide 15, this directory
   # use H3K9Cr, Thal|Normal|Thal+Normal where, by default, colored scales differ

tbl.tsne <- x$getTsne()
dim(tbl.tsne)  # [1] 98219     2
marker <- "H3K9cr"
stopifnot(marker %in% names(markers))

full.marker.name <- markers[[marker]]
plot.marker.on.tsne.coordinates(mtx, tbl.tsne, full.marker.name)
#----------------------------------------------------------------------------------------------------
create.tsne.plotting.table <- function(mtx, tbl.tsne, full.marker.name,
                                       colors, color.boundaries)
{
   stopifnot(nrow(mtx) == nrow(tbl.tsne))
   stopifnot(full.marker.name %in% colnames(mtx))
   stopifnot(length(colors) == length(color.boundaries))
   stopifnot(length(colors) >= 3 && length(colors) <= 11)

   tbl.plot <- cbind(tbl.tsne[, c("tsne_1", "tsne_2")],
                     mtx[, full.marker.name, drop=FALSE])
   require(RColorBrewer)
   # display.brewer.all()
   colors <- rev(brewer.pal(11, "Spectral"))
   bands <- 11

   #colors <- colorRampPalette(c("blue", "blue", "pink", "red"))(9)
   color.bins <- seq(-5, 5, length.out=length(colors))

   marker.vec <- tbl.plot[, full.marker.name]
   #qs <- fivenum(marker.vec)
   qs <- quantile(marker.vec, probs=seq(0, 1, length.out=10))
   marker.color.indices <- rep(NA_real_, length(marker.vec))
   for(i in 1:9){
     largerThanThreshold <- which(marker.vec >= qs[i])
     marker.color.indices[largerThanThreshold] <- i
     }
   table(marker.color.indices)

   tbl.plot$color <- colors[marker.color.indices]
   par(mar=c(5,5,5,5))#  + 0.1)

   with(tbl.plot, plot(tsne_1, tsne_2, col=color, main=full.marker.name,
                       ylim=c(-50, 50), xlim=c(-50, 50)))

   legend(40, 15, rev(as.character(round(as.numeric(qs), digits=2))),
              rev(colors))


   quartz()
   marker.vec <- tbl.plot$cr
   qs <- fivenum(marker.vec)
   marker.color.indices <- rep(1, length(marker.vec))
   marker.color.indices[marker.vec > qs[2]] <- 3
   marker.color.indices[marker.vec > qs[3]] <- 5
   marker.color.indices[marker.vec > qs[4]] <- 9
   table(marker.color.indices)

   tbl.plot$color.cr <- colors[marker.color.indices]
   with(tbl.plot, plot(tsne_1, tsne_2, col=color.cr, main="H3K9Cr"))

   quartz()
   marker.vec <- tbl.plot$h3
   qs <- fivenum(marker.vec)
   marker.color.indices <- rep(1, length(marker.vec))
   marker.color.indices[marker.vec > qs[2]] <- 3
   marker.color.indices[marker.vec > qs[3]] <- 5
   marker.color.indices[marker.vec > qs[4]] <- 9
   table(marker.color.indices)

   tbl.plot$color.h3 <- colors[marker.color.indices]
   with(tbl.plot, plot(tsne_1, tsne_2, col=color.h3, main="H3"))

   color.indices <- tbl.plot$cr %/% round(max(tbl.plot$cr)/bands) + 1
   table(color.indices)
   tbl.plot$color.cr <- colors[color.indices]
   with(tbl.plot, plot(tsne_1, tsne_2, col=color.cr, main="H3K9Cr"))

   color.indices <- tbl.plot$h3 %/% round(max(tbl.plot$h3)/bands) + 1
   table(color.indices)
   tbl.plot$color.h3 <- colors[color.indices]
   with(tbl.plot, plot(tsne_1, tsne_2, col=color.h3, main="H3"))

} # plot.marker.on.tsne.coordinates
#----------------------------------------------------------------------------------------------------
