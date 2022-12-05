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
        #' @description calculate and append a normalized vector column for the specified
        #'   histone marker, regressed against the base ("total") H3 and/or H4 vector
        #' @param target.marker character, a column name from the matrix, e.g., "Sm154Di<154Sm_H3K27me2>"
        #' @param reference.markers character, the marker names for reference ("total") H3 and/or H4
        #' @return character the newly created column name
        normalizeMarker = function(target.marker, reference.markers){
            printf("--- entering normalizeMarker")
            stopifnot(all(c(target.marker, reference.markers)%in% names(private$markers)))
            columns.requested <- c(target.marker, reference.markers)
            columns.actual <- as.character(lapply(columns.requested,
                                                   function(shortName) private$markers[[shortName]]))

            mtx.2 <- private$mtx[, columns.actual]
            colnames(mtx.2) <- columns.requested
            model <- sprintf("%s ~ 1 + %s", target.marker, paste(reference.markers, collapse=" + "))
            fit <- lm(model, data=as.data.frame(mtx.2))
            marker.resid <- as.numeric(residuals(fit))
            new.col.name <- sprintf("%s.regress.%s", target.marker, paste(paste(reference.markers, collapse="+")))
            private$mtx <- cbind(private$mtx, new.col.name=marker.resid)
            colnames(private$mtx)[ncol(private$mtx)] <- new.col.name
            new.col.name
            } # normalizeMarker

       ) # public

    ) # class
#--------------------------------------------------------------------------------
