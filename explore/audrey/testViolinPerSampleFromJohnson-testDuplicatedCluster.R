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

mtx.eol1 <- mtx[grep("c15", rownames(mtx)),]
dim(mtx.eol1)  # 5000 59

mtx.mv4 <-  mtx[grep("c16", rownames(mtx)),]
dim(mtx.mv4)   # 5000 59

cell.names <- sample(rownames(mtx.extended), size=5)
cell.names <- rownames(mtx.extended)

clusters <- c(5,12)
tbl.violin.eol1 <-  x$createTableForViolinPlot(clusters, new.col.name, mtx.eol1)
tbl.violin.mv4 <-  x$createTableForViolinPlot(clusters, new.col.name, mtx.mv4)

tbl.violin.eol1$name <- paste0(tbl.violin.eol1$name, ".eol1")
tbl.violin.mv4$name <- paste0(tbl.violin.mv4$name, ".mv4")

tbl.violin.both <- rbind(tbl.violin.eol1, tbl.violin.mv4)

ggplot(tbl.violin.both, aes(x=name, y=value, fill=name)) + geom_violin() +
    theme(axis.text = element_text(size = 14)) +
    ggtitle(sprintf("cluster 5: %s", new.col.name))


# tbl.violin$name <-  factor(tbl.violin$name, levels = paste0(target, ".regress.H3+H4.c", clusters))
# ggplot(tbl.violin, aes(x=name, y=value, fill=name)) + geom_violin() +
#   #coord_cartesian(ylim=c(-5,5)) +
#   theme(axis.text = element_text(size = 14)) +
#   ggtitle(sprintf("%s - AML histones", new.col.name))+
#   geom_boxplot(width=0.1, outlier.colour = NA)
# export::graph2ppt(file = "New-working-MLLsamples-H3+H4_normalization_n1_noALL-Jurkat.pptx", append = TRUE)
#
#
# #### Plot cells instead of clusters ####
#
# cells <-
#   list(`MV4-11`= c(9,12,5),
#        `EOL1` = c(5,8),
#        `MLL-PTD` = c(6,7,10,11))
# samples <-
#   unlist(cells)
#
# tbl.violin <-
#   x$createTableForViolinPlot(samples, new.col.name)
#
# new.tbl.violin <-
#   do.call(rbind, lapply(names(cells), function(cellname){
#     dat <-
#       tbl.violin[tbl.violin$name %in% paste0(target, ".regress.H3+H4.c", cells[[cellname]]), ]
#     dat[, "cell"] <-
#       cellname
#     return(dat)
#   })
#   )
# new.tbl.violin$cell <-
#   factor(new.tbl.violin$cell, levels = names(cells))
#
# ggplot(new.tbl.violin, aes(x=cell, y=value, fill=cell)) + geom_violin() +
#   #coord_cartesian(ylim=c(-5,5)) +
#   theme(axis.text = element_text(size = 14)) +
#   ggtitle(sprintf("%s - AML histones", new.col.name))+
#   geom_boxplot(width=0.1, outlier.colour = NA)
# export::graph2ppt(file = "New-working-MLLsamples-H3+H4_normalization_cell_n1_noALL-Jurkat.pptx", append = TRUE)
#
