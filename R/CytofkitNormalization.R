#' @title CytofkitNormalization
#' @description A template for building documented, tested R6 classes
#' @name CytofkitNormalization

#' @field id identifier for a class object
#'
#' @examples
#'   rt <- R6Template$new(id="abc")
#' @export

CytofkitNormalization = R6Class("CytofkitNormalization",

    #--------------------------------------------------------------------------------
    private = list(filename=NULL,
                   mtx=matrix(),
                   tbl.tsne=data.frame(),
                   clusters=list(),
                   markers=list()
                   ),

    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #' @param filename character, the full path to cytofkit results
         #' @return a new instance of CytofkitNormalization
        initialize = function(filename){
            private$filename <- filename
            x <- get(load(filename))
            private$mtx <- x$expressionData
            private$tbl.tsne <- as.data.frame(x$dimReducedRes[[1]])
            private$clusters <- x$clusterRes[[1]]
            },
        #------------------------------------------------------------
        #' @description accessor for the object's filename
        #' @return the current value of the filename
        getMatrix = function(){
           invisible(private$mtx)
           },

        #------------------------------------------------------------
        #' @description accessor for the object's tbl.tsne
        #' @return the data.frame
        getTsne = function(){
            invisible(private$tbl.tsne)
            },

        #------------------------------------------------------------
        #' @description accessor for the object's cluster list
        #' @return the matrix
        getClusters = function(){
            invisible(private$clusters)
        },
        #------------------------------------------------------------
        #' @description accessor for markers used in the cytof
        #' @return nothing
        createSimpleMarkerNames = function(){
            candidates <- colnames(private$mtx)

            for(notUseful in c("barcode", "center", "offset", "width", "residual", "environ",
                               "dna", "live-dead", "bckg", "bi209di", "beads")){
               deleters <- grep(notUseful, candidates, ignore.case=TRUE)
               if(length(deleters) > 0)
                   candidates <- candidates[-deleters]
               }
            markers <- candidates
            tokens <- strsplit(markers, "_")
            tokens.2 <- lapply(tokens, "[", 2)
            short.names <- sub(">", "", tokens.2)
            names(markers) <- short.names
            private$markers <- markers
            }, # getMarkers

        #------------------------------------------------------------
        #' @description accessor for named list, e.g. list("H3"="Yb176Di<176Yb_H3>")
        #' @return character list
        getMarkers = function(){
            private$markers
            }, # getMarkerShortNames


        #------------------------------------------------------------
        #' @description accessor for named list, e.g. list("H3"="Yb176Di<176Yb_H3>")
        #' @param clusterNumber numeric - must be in range
        #' @return character list
        getCluster = function(clusterNumber){
            stopifnot(clusterNumber %in% unique(as.numeric(private$clusters)))
            names(private$clusters[private$clusters==clusterNumber])
            }, # getMarkerShortNames


        #------------------------------------------------------------
        #' @description calculate and append a normalized vector column for the specified
        #'   histone marker, regressed against the base ("total") H3 and/or H4 vector
        #' @param target.marker character, a column name from the matrix, e.g., "Sm154Di<154Sm_H3K27me2>"
        #' @param reference.markers character, the marker names for reference ("total") H3 and/or H4
        #' @return character the newly created column name
        normalizeMarker = function(target.marker, reference.markers){
            if(!target.marker %in% names(private$markers))
                stop(sprintf("target marker '%s' is not recognized", target.marker))
            for(reference.marker in reference.markers)
                if(!reference.marker %in% names(private$markers))
                    stop(sprintf("reference marker '%s' is not recognized", target.marker))
            columns.requested <- c(target.marker, reference.markers)
            columns.actual <- as.character(lapply(columns.requested,
                                                  function(short.name) private$markers[[short.name]]))
            mtx.2 <- private$mtx[, columns.actual]
            colnames(mtx.2) <- columns.requested
            model <- sprintf("%s ~ 1 + %s", target.marker, paste(reference.markers, collapse=" + "))
            fit <- lm(model, data=as.data.frame(mtx.2))
            marker.resid <- as.numeric(residuals(fit))
            new.col.name <- sprintf("%s.regress.%s", target.marker, paste(reference.markers, collapse="+"))
            private$mtx <- cbind(private$mtx, new.col.name=marker.resid)
            colnames(private$mtx)[ncol(private$mtx)] <- new.col.name
            private$markers[[new.col.name]] <- new.col.name
            new.col.name
            }, # normalizeMarker

        #------------------------------------------------------------
        #' @description create a ggplot2 violin plot data.frame
        #' @param clusters numeric a vector of cluster numbers, in the order you prefer
        #' @param marker character, the marker name
        #' @param matrix matrix, a subset of the full matrix, used only if you wish to make violins
        #'        out of different classes of cells in each cluster see unitTests
        #' @return data.frame ready for call to ggplot
        #' @examples
        #' tbl.violin <- createTableForViolinPlot(1:20, "H3")
        #' ggplot(tbl.violin, aes(x=name, y=value, fill=name)) +
        #'     geom_violin() +
        #'     coord_cartesian(ylim=c(-5,5)) +
        #'     theme(axis.text = element_text(size = 14)) +
        #      ggtitle("example")
        createTableForViolinPlot = function(clusters, marker, matrix=NA){
           stopifnot(all(clusters %in% as.numeric(private$clusters)))
           if(all(is.na(matrix)))
              matrix <- self$getMatrix()
           tbls <- list()
           for(c in clusters){
              vec <- matrix[, private$markers[[marker]]]
              cells.in.cluster <- self$getCluster(c)
              cells.in.cluster.in.vector <- intersect(names(vec), cells.in.cluster)
              values <- as.numeric(vec[cells.in.cluster.in.vector])
              cluster.name <- sprintf("%s.c%d", marker, c)
              tbl.violin <- data.frame(name=cluster.name, value=values, stringsAsFactors=FALSE)
              tbls[[cluster.name]] <- tbl.violin
              } # for c
           tbl.violin <- do.call(rbind, tbls)
           rownames(tbl.violin) <- NULL
              # set the levels, to keep ggplot2 honest.  explicit level ordering
              # here, as a parameter to factor(), is apparently required
           tbl.violin$name <- factor(tbl.violin$name,
                                     levels=sprintf("%s.c%d", marker, clusters))
           tbl.violin
           }, # createTableForViolinPlot

        #------------------------------------------------------------
        #' @description calculates color boundaries for a distribution of numeric values
        #' @param vector numeric
        #' @param colors character, a vector of color names, e.g. "#5E4FA2" "#3288BD" ...
        #' @param mode character either "quantile" or "interval"
        #' @return data.frame with min, max, color columns

        calculateColorBoundaries = function(vector, colors, mode){
           stopifnot(mode %in% c("quantile", "interval"))
           if(mode == "quantile"){
               # we need one more boundary than there are colors:
              color.quantiles <- quantile(vector, probs=seq(0, 1, length.out=1+length(colors)))
              tbl <- data.frame(start=color.quantiles[1:length(colors)],
                                end=color.quantiles[2:(1+length(colors))],
                                color=colors)
              rownames(tbl) <- NULL
              return(tbl)
              } # quantile mode
           if(mode == "interval"){
              boundaries <- seq(min(vector), max(vector), length.out=(1+length(colors)))
              tbl <- data.frame(start=boundaries[1:length(colors)],
                                end=boundaries[2:(1+length(colors))],
                                color=colors)
              rownames(tbl) <- NULL
              return(tbl)
              } # interval mode
           }, # calculateColorBoundaries

        #------------------------------------------------------------
        #' @description create a tsne (x,y,color) data.frame with specified color mapping,
        #'    for the specified marker
        #' @param marker character, the (short) marker name
        #' @param tbl.colors data.frame with colnames "start", "end", "color"
        #' @param matrix.sub matrix, an optional subset of the full matrix,
        #' @return data.frame (for tsne plotting)

        createTableForTsnePlot = function(marker, tbl.colors, matrix.sub=NA){
           stopifnot(marker %in% names(private$markers))
           marker.longName <- private$markers[[marker]]
           mtx <- self$getMatrix()
           if(!all(is.na(matrix.sub)))
               mtx <- matrix.sub
           stopifnot(marker.longName %in% colnames(mtx))
           marker.vec <- mtx[, marker.longName, drop=TRUE]
           tbl.tsne <- self$getTsne()
           stopifnot(all(names(marker.vec) %in% rownames(tbl.tsne)))
           tbl.tsne <- tbl.tsne[names(marker.vec),]
           tbl.tsne$value <- as.numeric(marker.vec[rownames(tbl.tsne)])
           color.lookup <- function(value){
              which(tbl.colors$end >= value)[1]
              }
           colors.nums <- unlist(lapply(tbl.tsne$value, function(value) color.lookup(value)))
           tbl.tsne$colorNum <- colors.nums
           tbl.tsne$color <- tbl.colors$color[tbl.tsne$colorNum]
           # marker.color.indices <- rep(NA_real_, length(marker.vec.full))
           # color.quantiles <- quantile(marker.vec.full, probs=seq(0, 1, length.out=length(colors)))
           # marker.color.indices <- rep(NA_real_, length(marker.vec))
           # for(i in seq_len(length(colors))){
           #     largerThanThreshold <- which(marker.vec >= color.quantiles[i])
           #     marker.color.indices[largerThanThreshold] <- i
           #    }
           # table(marker.color.indices)
           # tbl.tsne <- self$getTSNE()
           # tbl.tsne$color <- colors[marker.color.indices]
           # list(tbl.tsne=tbl.tsne, color.quantiles=color.quantiles)
           invisible(tbl.tsne)
           } # createTableForTsnePlot
       ) # public

    ) # class
#--------------------------------------------------------------------------------
