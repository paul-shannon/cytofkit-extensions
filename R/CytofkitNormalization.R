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
                   clusters=list()
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
            }
       ) # public

    ) # class
#--------------------------------------------------------------------------------
