library(CytofkitNormalization)
library(RColorBrewer)
library(RUnit)

cytofkit.results.file <- "exp54-NvsThal.RData"
f <- system.file(package="CytofkitNormalization", "extdata", cytofkit.results.file)
checkTrue(file.exists(f))
x <- CytofkitNormalization$new(f)

x$createSimpleMarkerNames()
markers <- x$getMarkers()
marker <- "H3K18cr"
checkTrue(marker %in% names(markers))

#boundary.mode <- "quantile"
boundary.mode <- "interval"

target <- "H3K18cr"
h3.reference <- "H3"
h4.reference <- "H4"

normalized.marker.name <- x$normalizeMarker(target, c(h3.reference, h4.reference))
markers <- x$getMarkers()
checkTrue(normalized.marker.name %in% names(markers))
checkTrue(markers[[normalized.marker.name]] == normalized.marker.name)

mtx <- x$getMatrix()
dim(mtx) # 98219    59
mtx.full <- mtx[, normalized.marker.name, drop=FALSE]
dim(mtx.full)

mtx.norm <- mtx.full[grep("-N_", rownames(mtx.full)),,drop=FALSE]
mtx.thal <- mtx.full[grep("-Thal_", rownames(mtx.full)),,drop=FALSE]

checkEquals(nrow(mtx.full), nrow(mtx.norm) + nrow(mtx.thal))

fivenum(mtx.full)  # -3.91267858 -0.42954024 -0.01679341  0.45182397  3.80952729
fivenum(mtx.thal)  # -3.91267858 -0.43966255 -0.02107368  0.51595033  3.80952729
fivenum(mtx.norm)  # -3.74391611 -0.42035791 -0.01402064  0.39867349  2.77348580

    #--------------------------------------------------
    # copied from cytofkit/R/cytofkit_postProcess.R
    #--------------------------------------------------

spectral1 <- c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4",
               "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61",
               "#F46D43", "#D53E4F", "#9E0142")
spectral2 <- rev(c("#7F0000","red","#FF7F00","yellow","white",
                   "cyan", "#007FFF", "blue","#00007F"))


#colors <- rev(brewer.pal(11, "Spectral"))
colors <- spectral1

tbl.colors <- x$calculateColorBoundaries(mtx.full, colors, mode=boundary.mode)

tbl.tsneNorm <- x$createTableForTsnePlot(normalized.marker.name, tbl.colors, matrix.sub=mtx.norm)


    #--------------------------------------------------
    # plot marker tsne Normal
    #--------------------------------------------------

quartz(width=10, height=10)
par(mar=c(5,5,5,5)) # generous margins
title <- sprintf("%s NORM (%s colors) %s", normalized.marker.name, boundary.mode, cytofkit.results.file)
with(tbl.tsneNorm, plot(tsne_1, tsne_2, col=color,
                        main=title,
                        cex=0.4, pch=16,
                        ylim=c(-50, 50), xlim=c(-50, 50)))

boundaries <- round(tbl.colors$end, digits=2)
legend(35, 40, round(rev(boundaries), digits=3), rev(tbl.colors$color))

    #--------------------------------------------------
    # plot marker tsne Thalassemia
    #--------------------------------------------------

tbl.tsneThal <- x$createTableForTsnePlot(normalized.marker.name, tbl.colors, matrix.sub=mtx.thal)

quartz(width=10, height=10)
par(mar=c(5,5,5,5)) # generous margins
title <- sprintf("%s THAL (%s colors) %s", normalized.marker.name, boundary.mode, cytofkit.results.file)
with(tbl.tsneThal, plot(tsne_1, tsne_2, col=color,
                        main=title,
                        cex=0.4, pch=16,
                        ylim=c(-50, 50), xlim=c(-50, 50)))

boundaries <- round(tbl.colors$end, digits=2)
legend(35, 40, round(rev(boundaries), digits=3), rev(tbl.colors$color))



