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

mtx.untreated <- mtx[untreated.rows,, drop=FALSE]
mtx.cr2.5mM <- mtx[cr2.5mM.rows,, drop=FALSE]
mtx.cr10mM  <- mtx[cr10mM.rows,, drop=FALSE]
dim(mtx.untreated); dim(mtx.cr2.5mM); dim(mtx.cr10mM)

tbl.violin.0 <-  x$createTableForViolinPlot(1, marker=new.col.name, matrix=mtx.untreated)
tbl.violin.25 <- x$createTableForViolinPlot(1, marker=new.col.name, matrix=mtx.cr2.5mM)
tbl.violin.10 <- x$createTableForViolinPlot(1, marker=new.col.name, matrix=mtx.cr10mM)

tbl.violin.0$treatment <- "cr 0"
tbl.violin.25$treatment <- "cr 2.5"
tbl.violin.10$treatment <- "cr 10"


dim(tbl.violin.0)    # 2521 3
dim(tbl.violin.25)   # 2795 3
dim(tbl.violin.10)   # 3420 3

tbl.violin.combined <- rbind(tbl.violin.0, tbl.violin.25, tbl.violin.10)
tbl.violin.combined$treatment <- factor(tbl.violin.combined$treatment)

levels(tbl.violin.combined$treatment) <- c("cr 0", "cr 2.5", "cr 10")

p <- ggplot(tbl.violin.combined, aes(x=treatment, y=value, fill=name)) +
            geom_violin() +
            theme(axis.text = element_text(size = 14)) +
            stat_summary(fun=median, geom="point", size=1, color="white") +
            theme_bw() +
            # ggtitle(sprintf("%s - erythroid trajaectory", new.col.name)) +
            theme_grey(base_size = 18)

p <- p + scale_fill_manual(values=rep(c("gray", "red"), 13)) + theme(legend.position="none")
p


# levels(tbl.violin$name)
# new.levels <- as.character(lapply(myClusters,
#                                   function(cluster)
#                                       grep(sprintf("\\.c%d$", cluster),  levels(tbl.violin$name),
#                                            fixed=FALSE, value=TRUE)
#                                   )
#                            )
# levels(tbl.violin$name) <- new.levels

unt.rows <- grep("-unt_", rownames(mtx), fixed=TRUE)
cr2.5mM.rows <- grep("-2.5mM_", rownames(mtx), fixed=TRUE)
cr10mM.rows <- grep("-10mM_", rownames(mtx), fixed=TRUE)
checkEquals(nrow(mtx), sum(length(unt.rows), length(cr10mM.rows)))

mtx.unt <- mtx[unt.rows,, drop=FALSE]
mtx.cr2.5mM <- mtx[cr2.5mM.rows,, drop=FALSE]
mtx.cr10mM  <- mtx[cr10mM.rows,, drop=FALSE]
dim(mtx.unt); dim(mtx.cr2.5mM); dim(mtx.cr10mM)

tbl.violin.unt <-
  x$createTableForViolinPlot(myClusters, marker=new.col.name, matrix=mtx.unt)
checkEquals(length(which(is.na(tbl.violin.unt$values))), 0)
tbl.violin.unt$status <- "UNT"
checkEquals(as.character(lapply(tbl.violin.unt, class)), c("factor", "numeric", "character"))
checkEquals(levels(tbl.violin.unt$name),
            c("H3K9cr.regress.H3+H4.c1", "H3K9cr.regress.H3+H4.c2", "H3K9cr.regress.H3+H4.c10",
              "H3K9cr.regress.H3+H4.c6"))

tbl.violin.cr2.5mM <-
  x$createTableForViolinPlot(myClusters, marker=new.col.name, matrix=mtx.cr2.5mM)
checkEquals(length(which(is.na(tbl.violin.cr2.5mM$values))), 0)

tbl.violin.cr10mM <-
  x$createTableForViolinPlot(myClusters, marker=new.col.name, matrix=mtx.cr10mM)
checkEquals(length(which(is.na(tbl.violin.cr10mM$values))), 0)

tbl.violin.unt$name <- paste0(tbl.violin.unt$name, ".unt")
tbl.violin.cr2.5mM$name <- paste0(tbl.violin.cr2.5mM$name, ".cr2.5mM")
tbl.violin.cr10mM$name <- paste0(tbl.violin.cr10mM$name, ".cr10mM")
tbl.violin.cr2.5mM$status <- "cr2.5mM"
tbl.violin.cr10mM$status <- "cr10mM"

tbl.violin.both <- rbind(tbl.violin.unt, tbl.violin.cr2.5mM, tbl.violin.cr10mM)
long.names <- tbl.violin.both$name

#--------------------------------------------------------
# remove the long & redundant prefix: same for all cells
#--------------------------------------------------------

shorter.names <- sub("H3K9cr.regress.H3+H4.", "", long.names, fixed=TRUE)

#--------------------------------------------------------
# remove unt cr10mM N & T: this is now encoded in status
#--------------------------------------------------------
shorter.names <- sub("\\.[U2.510]$", "", shorter.names)

#--------------------------------------------------------
# use factors with ordered levels so ggplot will order the
# clusters in x properly
#--------------------------------------------------------

preferred.order <- c("c1")
tbl.violin.both$name <- factor(shorter.names, levels = preferred.order, ordered = TRUE)

p <- ggplot(tbl.violin.both,
            aes(x=name, y=value, fill=status)) +
  geom_violin() +
  theme(axis.text = element_text(size = 14)) +
  #stat_summary(fun=median, geom="point", size=3, color="black") +
  theme_bw() +
  ggtitle(sprintf("%s - erythroid trajectory", new.col.name)) +
  theme_grey(base_size = 18)

p <- p + scale_fill_manual(values=c("#7ACB8D", "#E09AAD", "blue"))
print(p)

