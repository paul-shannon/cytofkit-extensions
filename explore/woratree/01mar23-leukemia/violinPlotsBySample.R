library(RUnit)
library(CytofkitNormalization)
library(ggplot2)
f <- "2023-01-31_cytofkit.RData"

x <- CytofkitNormalization$new(f)
myClusters <- c(1,7,3,8,5,6,12,15,14,16,13,17,18) # missing 2, 4, 9, 10 , 11

x$createSimpleMarkerNames()
markers <- x$getMarkers()
clusters <- x$getClusters()
mtx <- x$getMatrix()
dim(mtx)   # 30000 58
head(rownames(mtx))

length(grep("c12_2023-01-31_Helios4-4_Sample_LIVE_", rownames(mtx)))  # 5000

sample.names <- sub("_[0-9]+$", "", rownames(mtx))
length(table(sample.names))  # 6
as.data.frame(table(sample.names))
  #                            sample.names Freq
  # 1 c12_2023-01-31_Helios4-4_Sample_LIVE 5000
  # 2 c13_2023-01-31_Helios4-4_Sample_LIVE 5000
  # 3 c14_2023-01-31_Helios4-4_Sample_LIVE 5000
  # 4 c15_2023-01-31_Helios4-4_Sample_LIVE 5000
  # 5 c16_2023-01-31_Helios4-4_Sample_LIVE 5000
  # 6 c17_2023-01-31_Helios4-4_Sample_LIVE 5000


   # c12_2023-01-31_Helios4-4_Sample_LIVE ->ALL-SIL c12 
   # c13_2023-01-31_Helios4-4_Sample_LIVE ->Patient13
   # c14_2023-01-31_Helios4-4_Sample_LIVE ->Patient14
   # c15_2023-01-31_Helios4-4_Sample_LIVE ->EOL1 
   # c16_2023-01-31_Helios4-4_Sample_LIVE ->MV4 
   # c17_2023-01-31_Helios4-4_Sample_LIVE ->Jurkat 

# woratree says: we need a violin plot of histone marks after normalization
# showing by sample (not by cluster).

#-----------------------------------------------------------------------
# normalized against H3 and H4
#-----------------------------------------------------------------------

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
checkEquals(length(which(is.na(tbl.violin.normal$values))), 0)
tbl.violin.thal <-
    x$createTableForViolinPlot(myClusters, marker=new.col.name, matrix=mtx.thal)
checkEquals(length(which(is.na(tbl.violin.thal$values))), 0)


tbl.violin.normal$name <- paste0(tbl.violin.normal$name, ".N")
tbl.violin.thal$name <- paste0(tbl.violin.thal$name, ".T")

tbl.violin.both <- rbind(tbl.violin.normal, tbl.violin.thal)
long.names <- tbl.violin.both$name
shorter.names <- sub("H3K9me3.regress.H3+H4.", "", long.names, fixed=TRUE)

tbl.violin.both$name <- factor(shorter.names)

  # the name column is now a factor.  the order of its "levels" controls their order in the plot
levels(tbl.violin.both$name) <- c("c1.N", "c1.T", "c7.N", "c7.T", "c3.N", "c3.T", "c8.N", "c8.T",
                                  "c5.N", "c5.T", "c6.N", "c6.T", "c12.N", "c12.T", "c15.N", "c15.T",
                                  "c14.N", "c14.T", "c16.N", "c16.T", "c13.N", "c13.T", "c17.N", "c17.T",
                                  "c18.N", "c18.T")

p <- ggplot(tbl.violin.both, aes(x=name, y=value, fill=name)) +
            geom_violin() +
            theme(axis.text = element_text(size = 14)) +
            stat_summary(fun=median, geom="point", size=1, color="white") +
            theme_bw() +
            ggtitle(sprintf("%s - erythroid trajectory", new.col.name)) +
            theme_grey(base_size = 18)

p <- p + scale_fill_manual(values=rep(c("gray", "red"), 13)) + theme(legend.position="none")
p

# export::graph2ppt(file = "compare groups-H3+H4-normalization-2023-01-10_Woratree.pptx")

