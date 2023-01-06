library(CytofkitNormalization)
library(ggplot2)
library(RUnit)

f <- "exp54-NvsThal.RData"
checkTrue(file.exists(f))
x <- CytofkitNormalization$new(f)

mtx <- x$getMatrix()
dim(mtx)  # 98219    58

head(rownames(mtx))
   # [1] "day6-N_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_25907"
   # [2] "day6-N_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_4214"
   # [3] "day6-N_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_46609"
   # [4] "day6-N_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_57935"
   # [5] "day6-N_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_21447"
   # [6] "day6-N_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_43474"
tail(rownames(mtx))
   # [1] "day14-Thal_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_3631"
   # [2] "day14-Thal_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_2859"
   # [3] "day14-Thal_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_11140"
   # [4] "day14-Thal_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_4078"
   # [5] "day14-Thal_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_13376"
   # [6] "day14-Thal_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_3661"

tbl.tsne <- x$getTsne()
dim(tbl.tsne)  # 98219     2
head(rownames(tbl.tsne))
   # [1] "day6-N_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_25907"
   # [2] "day6-N_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_4214"
   # [3] "day6-N_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_46609"
   # [4] "day6-N_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_57935"
   # [5] "day6-N_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_21447"
   # [6] "day6-N_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_43474"

tail(rownames(tbl.tsne))
   # [1] "day14-Thal_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_3631"
   # [2] "day14-Thal_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_2859"
   # [3] "day14-Thal_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_11140"
   # [4] "day14-Thal_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_4078"
   # [5] "day14-Thal_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_13376"
   # [6] "day14-Thal_2022-11-25_Herlios3_3_Normal-vs-Thalassaemia_LIVE_3661"

x$createSimpleMarkerNames()
markers <- x$getMarkers()

h3.markers <- grep("H3", markers, value=TRUE)
h4.markers <- grep("H4", markers, value=TRUE)

target <- markers[["H3K9Ac"]]
h3.reference <- markers[["H3"]]
h4.reference <- markers[["H4"]]

mtx <- x$getMatrix()
normal.rows <- grep("-N_", rownames(mtx), fixed=TRUE)
thalassemia.rows <- grep("-Thal_", rownames(mtx), fixed=TRUE)
checkEquals(nrow(mtx), sum(length(normal.rows), length(thalassemia.rows)))

mtx.normal <- mtx[normal.rows,]
mtx.thal   <- mtx[thalassemia.rows,]

clusters.requested <- 1:10
tbl.violin.normal <- x$createTableForViolinPlot(clusters=clusters.requested, marker="H3", matrix=mtx.normal)
tbl.violin.thal <- x$createTableForViolinPlot(clusters=clusters.requested, marker="H3", matrix=mtx.thal)

tbl.violin.normal$name <- paste0(tbl.violin.normal$name, ".N")
tbl.violin.thal$name <- paste0(tbl.violin.thal$name, ".T")

tbl.violin.both <- rbind(tbl.violin.normal, tbl.violin.thal)

ggplot(tbl.violin.both, aes(x=name, y=value, fill=name)) + geom_violin() +
    theme(axis.text = element_text(size = 14)) +
    ggtitle("H3 - cluster 3")

