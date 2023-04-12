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

mtx.mv4 <-  mtx[grep("c17", rownames(mtx)),]
dim(mtx.mv4)   # 5000 59

cell.names <- sample(rownames(mtx.extended), size=5)
cell.names <- rownames(mtx.extended)

clusters <- c(5,12)
tbl.violin.eol1 <-  x$createTableForViolinPlot(clusters, new.col.name, mtx.eol1)
tbl.violin.mv4 <-  x$createTableForViolinPlot(clusters, new.col.name, mtx.mv4)

tbl.violin.eol1$name <- paste0(tbl.violin.eol1$name, ".eol1")
tbl.violin.mv4$name <- paste0(tbl.violin.mv4$name, ".mv4")

tbl.violin.both <- rbind(tbl.violin.eol1, tbl.violin.mv4)
dim(tbl.violin.both)  # 2638 2

   # learn the distribution of names
table(tbl.violin.both$name)
   # MLL-Cterm.regress.H3+H4.c12.eol1     78
   # MLL-Cterm.regress.H3+H4.c12.mv4     833
   # MLL-Cterm.regress.H3+H4.c5.eol1     894
   # MLL-Cterm.regress.H3+H4.c5.mv4      833

   # ggplot2 requires this fancy footwork to get the violins
   # presented in the order we want

preferred.order <- c("MLL-Cterm.regress.H3+H4.c5.eol1",
                     "MLL-Cterm.regress.H3+H4.c5.mv4",
                     "MLL-Cterm.regress.H3+H4.c12.eol1",
                     "MLL-Cterm.regress.H3+H4.c12.mv4")

setdiff(preferred.order, tbl.violin.both$name)

checkTrue(all(preferred.order %in% tbl.violin.both$name))

tbl.violin.both$name <- factor(tbl.violin.both$name,
                               levels = preferred.order,
                               ordered = TRUE)

ggplot(tbl.violin.both, aes(x=name, y=value, fill=name)) + geom_violin() +
    theme(axis.text = element_text(size = 14)) +
    ggtitle(sprintf("%s", new.col.name))
