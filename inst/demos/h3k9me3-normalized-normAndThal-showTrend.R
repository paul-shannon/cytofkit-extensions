library(RUnit)
library(CytofkitNormalization)
library(ggplot2)

f <- system.file(package="CytofkitNormalization", "extdata", "exp54-NvsThal.RData")
checkTrue(file.exists(f))

x <- CytofkitNormalization$new(f)
message(sprintf("--- test_createTableForViolinPlot"))

ordered.clusters.of.interest <- c(1,7,3,8,5,6,12,15,14,16,13,17,18) # omit 2, 4, 9, 10, 11

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
tbl.tsne <- x$getTsne()
clusters <- x$getClusters()
cells.in.cluster.7 <- names(clusters[clusters==7])
length(cells.in.cluster.7)  # 12659
vec.7 <- as.numeric(mtx[cells.in.cluster.7,1])
round(fivenum(vec.7), digits=2) # -1.82 -0.18  0.23  0.60  2.70

tbl.violin <- x$createTableForViolinPlot(ordered.clusters.of.interest,
                                         marker=new.col.name, matrix=mtx)
    #--------------------------------------------------------------------------
    # simplify the cluster names, which come back as factors, with long names
    #--------------------------------------------------------------------------

lapply(tbl.violin, class)
tbl.violin$name <- as.character(tbl.violin$name)
coi <- "H3K9me3.regress.H3+H4.c7"
fivenum(subset(tbl.violin, name==coi)$value)
   #  -1.8202801 -0.1788207  0.2256116  0.5964198  2.6982170

tbl.violin$name <- sub("H3K9me3.regress.H3+H4.", "", fixed=TRUE, tbl.violin$name)
head(tbl.violin)
coi.short <- "c7"
fivenum(subset(tbl.violin, name==coi.short)$value)
   # -1.8202801 -0.1788207  0.2256116  0.5964198  2.6982170
tbl.violin$name <- as.factor(tbl.violin$name)
fivenum(subset(tbl.violin, name==coi.short)$value)
   # -1.8202801 -0.1788207  0.2256116  0.5964198  2.6982170

tbl.violin$name <- factor(tbl.violin$name,
                          levels=paste("c", ordered.clusters.of.interest, sep=""))
head(tbl.violin)

for(cluster in ordered.clusters.of.interest){
   cluster.sig <- sprintf("c%d$", cluster)
   tbl.cluster <- subset(tbl.violin, grepl(cluster.sig, as.character(tbl.violin$name),
                                           fixed=FALSE))
   stats <- fivenum(tbl.cluster$value)
   quartiles <- sprintf("%6.2f %6.2f %6.2f %6.2f %6.2f", stats[1], stats[2], stats[3],
                        stats[4], stats[5])
   printf("cluster %2d, %5d rows, median: %s", cluster, nrow(tbl.cluster), quartiles)
   }


lapply(tbl.violin, class)
levels(tbl.violin$name)

p <- ggplot(tbl.violin,
            aes(x=name, y=value, fill=name)) +
            geom_violin() +
            geom_boxplot(width=.1) +
            theme(axis.text = element_text(size = 14)) +
            theme_bw() +
            ggtitle(sprintf("%s - erythroid trajectory", new.col.name)) +
            theme_grey(base_size = 18)
#p <- p + scale_fill_manual(values=rep("beige", 13))
p <- p + theme(legend.position="none")
p <- p + scale_fill_manual(values=rep("beige",  13))
p
