library(RUnit)
library(CytofkitNormalization)
library(ggplot2)

f <- system.file(package="CytofkitNormalization", "extdata", "exp54-NvsThal.RData")
checkTrue(file.exists(f))

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
tbl.violin <- x$createTableForViolinPlot(myClusters, marker=new.col.name, matrix=mtx)

   #----------------------------------------------------------------
   # explore woratree's claim that median expression drops from c15
   #---------------------------------------------------------------
for(cluster in myClusters){
   cluster.sig <- sprintf("\\.c%d$", cluster)
   tbl.cluster <- subset(tbl.violin, grepl(cluster.sig, tbl.violin$name, fixed=FALSE))
   quartiles <- paste(as.character(round(fivenum(tbl.cluster$value), digits=2)), collapse=" ")
   printf("cluster %d, %d rows, median: %s", cluster, nrow(tbl.cluster), quartiles)
   }

levels(tbl.violin$name)
new.levels <- as.character(lapply(myClusters,
                                  function(cluster)
                                      grep(sprintf("\\.c%d$", cluster),  levels(tbl.violin$name),
                                           fixed=FALSE, value=TRUE)
                                  )
                           )
levels(tbl.violin$name) <- new.levels

normal.rows <- grep("-N_", rownames(mtx), fixed=TRUE)
thalassemia.rows <- grep("-Thal_", rownames(mtx), fixed=TRUE)
checkEquals(nrow(mtx), sum(length(normal.rows), length(thalassemia.rows)))

mtx.normal <- mtx[normal.rows,, drop=FALSE]
mtx.thal   <- mtx[thalassemia.rows,, drop=FALSE]
dim(mtx.normal); dim(mtx.thal)

tbl.violin.normal <-
  x$createTableForViolinPlot(myClusters, marker=new.col.name, matrix=mtx.normal)
checkEquals(length(which(is.na(tbl.violin.normal$values))), 0)
tbl.violin.normal$status <- "NORM"

tbl.violin.thal <-
    x$createTableForViolinPlot(myClusters, marker=new.col.name, matrix=mtx.thal)
checkEquals(length(which(is.na(tbl.violin.thal$values))), 0)

tbl.violin.normal$name <- paste0(tbl.violin.normal$name, ".N")
tbl.violin.thal$name <- paste0(tbl.violin.thal$name, ".T")
tbl.violin.thal$status <- "THAL"

for(cluster in myClusters){
   cluster.sig <- sprintf(".c%d", cluster)
   tbl.cluster <- subset(tbl.violin.thal, grepl(cluster.sig, tbl.violin.thal$name, fixed=TRUE))
   printf("thal cluster %d, %d rows, median: %f", cluster, nrow(tbl.cluster), median(tbl.cluster$value))
   }

for(cluster in myClusters){
   cluster.sig <- sprintf(".c%d", cluster)
   tbl.cluster <- subset(tbl.violin.normal, grepl(cluster.sig, tbl.violin.normal$name, fixed=TRUE))
   printf("normal cluster %d, %d rows, median: %f", cluster, nrow(tbl.cluster), median(tbl.cluster$value))
   }


tbl.violin.both <- rbind(tbl.violin.normal, tbl.violin.thal)
long.names <- tbl.violin.both$name
shorter.names <- sub("H3K9me3.regress.H3+H4.", "", long.names, fixed=TRUE)

tbl.violin.both$name <- factor(shorter.names)
tbl.violin.both$status <- as.factor(tbl.violin.both$status)

for(cluster in myClusters){
   cluster.sig <- sprintf("c%d.", cluster)
   tbl.cluster <- subset(tbl.violin.both, grepl(cluster.sig, tbl.violin.both$name, fixed=TRUE))
   printf("both cluster %d, %d rows, median: %f", cluster, nrow(tbl.cluster), median(tbl.cluster$value))
   }



  # the name column is now a factor.  the order of its "levels" controls their order in the plot
levels(tbl.violin.both$name) <- c("c1.N", "c1.T", "c7.N", "c7.T", "c3.N", "c3.T", "c8.N", "c8.T",
                                  "c5.N", "c5.T", "c6.N", "c6.T", "c12.N", "c12.T", "c15.N", "c15.T",
                                  "c14.N", "c14.T", "c16.N", "c16.T", "c13.N", "c13.T", "c17.N", "c17.T",
                                  "c18.N", "c18.T")
tbl.violin.x <- tbl.violin.both
tbl.violin.x <- subset(tbl.violin.both, name=="c1.N")

p <- ggplot(tbl.violin.x, aes(x=name, y=value, fill=name)) +
            geom_violin() +
            theme(axis.text = element_text(size = 14)) +
            stat_summary(fun=median, geom="point", size=3, color="black") +
            theme_bw() +
            ggtitle(sprintf("%s - erythroid trajectory", new.col.name)) +
            theme_grey(base_size = 18)

p <- p + scale_fill_manual(values=rep(c("gray", "red"), 13)) + theme(legend.position="none")
p

# export::graph2ppt(file = "compare groups-H3+H4-normalization-2023-01-10_Woratree.pptx")

