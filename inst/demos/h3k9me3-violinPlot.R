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

# levels(tbl.violin$name)
# new.levels <- as.character(lapply(myClusters,
#                                   function(cluster)
#                                       grep(sprintf("\\.c%d$", cluster),  levels(tbl.violin$name),
#                                            fixed=FALSE, value=TRUE)
#                                   )
#                            )
# levels(tbl.violin$name) <- new.levels

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
checkEquals(as.character(lapply(tbl.violin.normal, class)), c("factor", "numeric", "character"))
checkEquals(levels(tbl.violin.normal$name),
   c("H3K9me3.regress.H3+H4.c1", "H3K9me3.regress.H3+H4.c7", "H3K9me3.regress.H3+H4.c3",
     "H3K9me3.regress.H3+H4.c8", "H3K9me3.regress.H3+H4.c5", "H3K9me3.regress.H3+H4.c6",
     "H3K9me3.regress.H3+H4.c12", "H3K9me3.regress.H3+H4.c15", "H3K9me3.regress.H3+H4.c14",
     "H3K9me3.regress.H3+H4.c16", "H3K9me3.regress.H3+H4.c13","H3K9me3.regress.H3+H4.c17",
     "H3K9me3.regress.H3+H4.c18"))

tbl.violin.thal <-
    x$createTableForViolinPlot(myClusters, marker=new.col.name, matrix=mtx.thal)
checkEquals(length(which(is.na(tbl.violin.thal$values))), 0)

tbl.violin.normal$name <- paste0(tbl.violin.normal$name, ".N")
tbl.violin.thal$name <- paste0(tbl.violin.thal$name, ".T")
tbl.violin.thal$status <- "THAL"

tbl.violin.both <- rbind(tbl.violin.normal, tbl.violin.thal)
long.names <- tbl.violin.both$name

   #--------------------------------------------------------
   # remove the long & redundant prefix: same for all cells
   #--------------------------------------------------------

shorter.names <- sub("H3K9me3.regress.H3+H4.", "", long.names, fixed=TRUE)

   #--------------------------------------------------------
   # remove Normal Thal N & T: this is now encoded in status
   #--------------------------------------------------------
shorter.names <- sub("\\.[NT]$", "", shorter.names)

   #--------------------------------------------------------
   # use factors with ordered levels so ggplot will order the
   # clusters in x properly
   #--------------------------------------------------------

preferred.order <- c("c1", "c3", "c5", "c6", "c7",
                     "c8", "c12", "c13", "c14", "c15",
                     "c16", "c17", "c18")
tbl.violin.both$name <- factor(shorter.names, levels = preferred.order, ordered = TRUE)

p <- ggplot(tbl.violin.both,
            aes(x=name, y=value, fill=status)) +
            geom_violin() +
            theme(axis.text = element_text(size = 14)) +
            #stat_summary(fun=median, geom="point", size=3, color="black") +
            theme_bw() +
            ggtitle(sprintf("%s - erythroid trajectory", new.col.name)) +
            theme_grey(base_size = 18)

p <- p + scale_fill_manual(values=c("darkgreen", "red"))
print(p)



