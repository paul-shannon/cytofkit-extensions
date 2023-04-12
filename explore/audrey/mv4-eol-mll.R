library(RUnit)
library(CytofkitNormalization)
library(ggplot2)

# f <- "C:/Users/koppa/Desktop/PhD/MLL-PTD/Protocols-techniques/CyTOF-single cells/2023-01-30 CyToF MLL myeloid/R file only cell lines+MLL-PTD cd markers/R file all cd markers with 45RA+Cter-Nter/2023-01-31_cytofkit.RData"
f <- "2023-01-31_cytofkit.RData"
checkTrue(file.exists(f))

x <- CytofkitNormalization$new(f)


x$createSimpleMarkerNames()
markers <- x$getMarkers()
clusters <- x$getClusters()



#### H3+H4 -normalization ####
target <- "MLL-Cterm"
checkTrue(target %in% names(markers))

h3.reference <- "H3"
h4.reference <- "H4"

new.col.name <- x$normalizeMarker(target, c(h3.reference, h4.reference))
markers <- x$getMarkers()

mtx <- x$getMatrix()

# from audrey:
#   c12 is ALL-SIL
#   c13 is MLL-PTD BM
#   c14 is MLL-PTD spleen
#   c15 is EOL1 cells
#   c16 is MV4
#   c17 Jurkat cells

mtx.eol1 <- mtx[grep("c15_", rownames(mtx)),]
dim(mtx.eol1)  # 5000 59

mtx.mv4 <-  mtx[grep("c16_", rownames(mtx)),]
dim(mtx.mv4)   # 5000 59

   # not sure if you want "BM" or "spleen" so chosing  MLL-PTD BM
mtx.mll <- mtx[grep("c13", rownames(mtx)),]
dim(mtx.mll)   # 5000 59

   # from audrey, the celltypes/clusters she is intereested in
   #   list(`MV4-11`= c(9,12,5),
   #        `EOL1` = c(5,8),
   #        `MLL-PTD` = c(6,7,10,11))

clusters.for.mv4 <- c(9, 12, 5)
clusters.for.eol1 <- c(5, 8)
clusters.for.mll <- c(6,7,10,11)

tbl.violin.eol1 <-  x$createTableForViolinPlot(clusters.for.eol1, new.col.name, mtx.eol1)
dim(tbl.violin.eol1)  # 3870

tbl.violin.mv4 <-  x$createTableForViolinPlot(clusters.for.mv4, new.col.name, mtx.mv4)
dim(tbl.violin.mv4)   # 4385

tbl.violin.mll <-  x$createTableForViolinPlot(clusters.for.mll, new.col.name, mtx.mll)
dim(tbl.violin.mll)   # 4599

tbl.violin.eol1$name <- paste0(tbl.violin.eol1$name, ".eol1")
tbl.violin.mv4$name <- paste0(tbl.violin.mv4$name, ".mv4")
tbl.violin.mll$name <- paste0(tbl.violin.mll$name, ".mll")

tbl.violin.combined <- rbind(tbl.violin.eol1, tbl.violin.mv4, tbl.violin.mll)
dim(tbl.violin.combined)  # 12854

   # remove the redundant parts of the violin names
tbl.violin.combined$name <- sub("MLL-Cterm.regress.H3+H4.", "",
                                tbl.violin.combined$name,
                                fixed=TRUE)

   # check the distribution of names (celltype, cluster, count)
   # the cNN numbers should match those defined above, for each celltype
table(tbl.violin.combined$name);
   #  c6.mll: 2523
   #  c7.mll: 1633
   # c10.mll: 409
   # c11.mll: 34
   #
   #  c9.mv4: 2719
   # c12.mv4: 833
   #  c5.mv4: 833
   #
   # c5.eol1: 894
   # c8.eol1: 2976

   # ggplot2 requires this fancy footwork to get the violins
   # presented in the order we want

preferred.order <- c("c9.mv4",
                     "c12.mv4",
                     "c5.mv4",
                     "c5.eol1",
                     "c8.eol1",
                     "c6.mll",
                     "c7.mll",
                     "c10.mll",
                     "c11.mll")

setdiff(preferred.order, tbl.violin.combined$name)  # should be character(0)
checkTrue(all(preferred.order %in% tbl.violin.combined$name))

tbl.violin.combined$name <- factor(tbl.violin.combined$name,
                                   levels = preferred.order,
                                   ordered = TRUE)

gg <- ggplot(tbl.violin.combined, aes(x=name, y=value, fill=name)) + geom_violin() +
         theme(axis.text = element_text(size = 14)) +
         ggtitle(sprintf("%s", new.col.name))
