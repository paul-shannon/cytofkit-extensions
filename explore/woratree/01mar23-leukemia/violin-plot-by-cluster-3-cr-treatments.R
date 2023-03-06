library(RUnit)
library(CytofkitNormalization)
library(ggplot2)

f <- "CrotonateTreatment.RData"
checkTrue(file.exists(f))
x <- CytofkitNormalization$new(f)

myClusters <- c(1,2,10,6) # CFUe,ProEB,BasoEB1,BasoEB2

x$createSimpleMarkerNames()
markers <- x$getMarkers()
clusters <- x$getClusters()

   #------------------------------------------------
   # normalize against H3 and H4
   #
   # the expression matrix uses long obscure names
   # e.g. "Er167Di<167Er_H3K9cr>"
   # which we shorten to "H3K9cr" the "markers"
   # variable maps short names to the long names
   #------------------------------------------------

markers <- x$getMarkers()

target <-   "H3K9cr"
h3.reference <- "H3"
h4.reference <- "H4"

stopifnot(target %in% names(markers))
stopifnot(h3.reference %in% names(markers))
stopifnot(h4.reference %in% names(markers))

    #-------------------------------------------------
    # normalize and get the resulting 1-column matrix
    #-------------------------------------------------
new.col.name <- x$normalizeMarker(target, c(h3.reference, h4.reference))

mtx <- x$getMatrix()[, new.col.name, drop=FALSE]
dim(mtx)

   #-------------------------------------------------------------
   # create violin plots, 1 cluster at a time, 1 violin
   # for each of 3 conditions:  untreated, 2.5mM, 10mM crotonate.
   #-------------------------------------------------------------

untreated.rows <- grep("-unt_", rownames(mtx), fixed=TRUE)
length(untreated.rows)  # 40000

cr2.5mM.rows <- grep("-2.5mM_", rownames(mtx), fixed=TRUE)
length(cr2.5mM.rows)    # 40000

cr10mM.rows <- grep("-10mM_", rownames(mtx), fixed=TRUE)
length(cr10mM.rows)     # 40000
nrow(mtx)
checkEquals(nrow(mtx),
            sum(length(untreated.rows), length(cr2.5mM.rows), length(cr10mM.rows)))

mtx.0 <- mtx[untreated.rows,, drop=FALSE]
mtx.25 <- mtx[cr2.5mM.rows,, drop=FALSE]
mtx.10  <- mtx[cr10mM.rows,, drop=FALSE]
dim(mtx.0); dim(mtx.25); dim(mtx.10)

cluster <- 1

tbl.violin.0 <-  x$createTableForViolinPlot(cluster, marker=new.col.name, matrix=mtx.0)
tbl.violin.25 <- x$createTableForViolinPlot(cluster, marker=new.col.name, matrix=mtx.25)
tbl.violin.10 <- x$createTableForViolinPlot(cluster, marker=new.col.name, matrix=mtx.10)

tbl.violin.0$treatment <- "cr 0"
tbl.violin.25$treatment <- "cr 2.5"
tbl.violin.10$treatment <- "cr 10"

dim(tbl.violin.0)    # 2521 3
dim(tbl.violin.25)   # 2795 3
dim(tbl.violin.10)   # 3420 3

     #------------------------------------------------------
     # to sanity check the violin plots, to make sure
     # they are labeled properly, get the quartiles on each
     #------------------------------------------------------

round(fivenum(tbl.violin.0$value), digits=2)   #  -1.98 -1.07 -0.65 -0.24  1.15
round(fivenum(tbl.violin.25$value), digits=2)  #  -1.95 -0.64 -0.26  0.10  1.48
round(fivenum(tbl.violin.10$value), digits=2)  #  -1.78  0.08  0.40  0.70  2.07

tbl.violin.combined <- rbind(tbl.violin.0, tbl.violin.25, tbl.violin.10)
tbl.violin.combined$treatment <- factor(tbl.violin.combined$treatment,
                                        levels=c("cr 0", "cr 2.5", "cr 10"))
tbl.violin.combined$name <- as.character(tbl.violin.combined$name)

levels(tbl.violin.combined$treatment)
table(tbl.violin.combined$treatment)
lapply(tbl.violin.combined, class)

p <- ggplot(tbl.violin.combined, aes(x=treatment, y=value,  color=treatment)) +
            geom_violin() +
            theme(axis.text = element_text(size = 14)) +
            theme_bw() +
            ggtitle(sprintf("crotynate: %s, cluster  %s", new.col.name, cluster)) +
            theme_grey(base_size = 18)

p <- p + scale_color_brewer(palette="Dark2")
p <- p +  theme(legend.position="none")
p
