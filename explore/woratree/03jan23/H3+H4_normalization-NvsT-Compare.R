library(RUnit)
library(CytofkitNormalization)
library(ggplot2)
f <- "/Users/woratreekaewsakulthong/Dropbox/Family\ Room/CYTOF-Ottawa/EXPERIMENTS/54-Thal\ vs\ normal\ histone\ panel/cytofkit/cytofkit.RData"
f <- "exp54-NvsThal.RData"

x <- CytofkitNormalization$new(f)
message(sprintf("--- test_createTableForViolinPlot"))

myClusters <- c(1,7,3,8,5,6,12,15,14,16,13,17,18) # missing 2, 4, 9, 10 , 11

x$createSimpleMarkerNames()
markers <- x$getMarkers()
clusters <- x$getClusters()

#-----------------------------------------------------------------------
# normalized against H3 and H4
#-----------------------------------------------------------------------
markers <- x$getMarkers()

  #--------------------------------------------------
  # the expression matrix uses long obscure names
  # like "Er167Di<167Er_H3K9me3>"
  # which we shorten to "H3K9me3" the "markers"
  # variable maps short names to the long names
  #--------------------------------------------------
target <-   "H3K9me3"
h3.reference <- "H3"
h4.reference <- "H4"

stopifnot(target %in% names(markers))
stopifnot(h3.reference %in% names(markers))
stopifnot(h4.reference %in% names(markers))
new.col.name <- x$normalizeMarker(target, c(h3.reference, h4.reference))

mtx <- x$getMatrix()[, new.col.name, drop=FALSE]
dim(mtx)
normal.rows <- grep("-N_", rownames(mtx), fixed=TRUE)
thalassemia.rows <- grep("-Thal_", rownames(mtx), fixed=TRUE)
checkEquals(nrow(mtx), sum(length(normal.rows), length(thalassemia.rows)))

mtx.normal <- mtx[normal.rows,, drop=FALSE]
mtx.thal   <- mtx[thalassemia.rows,, drop=FALSE]
dim(mtx.normal); dim(mtx.thal)

tbl.violin.normal <-
  x$createTableForViolinPlot(myClusters, marker=new.col.name, matrix=mtx.normal)
tbl.violin.thal <-
  x$createTableForViolinPlot(myClusters, marker=new.col.name, matrix=mtx.thal)


tbl.violin.normal$name <- paste0(tbl.violin.normal$name, ".N")
tbl.violin.thal$name <- paste0(tbl.violin.thal$name, ".T")

tbl.violin.both <- rbind(tbl.violin.normal, tbl.violin.thal)
long.names <- tbl.violin.both$name
shorter.names <- sub("H3K9me3.regress.H3+H4.", "", long.names, fixed=TRUE)

tbl.violin.both$name <- factor(shorter.names)

ggplot(tbl.violin.both, aes(x=name, y=value, fill=name)) +
  geom_violin() +
  theme(axis.text = element_text(size = 14)) +
  stat_summary(fun=median, geom="point", size=1, color="white") +
  theme_bw() +
  ggtitle(sprintf("%s - erythroid trajectory", new.col.name)) +
  theme_grey(base_size = 18)

# export::graph2ppt(file = "compare groups-H3+H4-normalization-2023-01-10_Woratree.pptx")

