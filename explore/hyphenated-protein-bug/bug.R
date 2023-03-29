library(RUnit)
library(CytofkitNormalization)

f <- system.file(package="CytofkitNormalization",
                 "extdata", "cytofkit-dashNameExample.RData")

checkTrue(file.exists(f))
x <- CytofkitNormalization$new(f)
x$createSimpleMarkerNames()
markers <- x$getMarkers()
head(markers)
target <- "MLL-2"
grep("MLL-2", names(markers), v=TRUE)


h3.reference <- "H3"
h4.reference <- "H4"
new.col.name <- x$normalizeMarker(target, c(h3.reference, h4.reference))
checkEquals(new.col.name, "MLL-2.regress.H3+H4")

