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
tbl.violin.shortClusterNames <- tbl.violin
tbl.violin.shortClusterNames$name <- sub("H3K9me3.regress.H3+H4.", "",
                                         tbl.violin.shortClusterNames$name,
                                         fixed=TRUE)
tbl.violin.shortClusterNames$name <- as.factor(tbl.violin.shortClusterNames$name)
levels(tbl.violin.shortClusterNames$name)
new.levels <- as.character(lapply(myClusters,
                                  function(cluster)
                                      grep(sprintf("c%d$", cluster),
                                           levels(tbl.violin.shortClusterNames$name),
                                           fixed=FALSE, value=TRUE)
                                  )
                           )

levels(tbl.violin.shortClusterNames$name) <- new.levels
p <- ggplot(tbl.violin.shortClusterNames,
            aes(x=name, y=value, fill=name)) +
            geom_violin() + # draw_quantiles = c(0.25, 0.5, 0.75)) +
            geom_boxplot(width=.1) +
            theme(axis.text = element_text(size = 14)) +
            #stat_summary(fun.y=median, geom="point", size=3, color="black") +
            theme_bw() +
            ggtitle(sprintf("%s - erythroid trajectory", new.col.name)) +
            theme_grey(base_size = 18)
p <- p + scale_fill_manual(values=rep("beige", 13)) + theme(legend.position="none")
p
