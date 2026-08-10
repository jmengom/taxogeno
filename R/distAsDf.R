#' Convert Distance Object to Data Frame
#'
#' Converts an R \code{dist} object to a data frame with rows, columns, and values.
#'
#' @param inDist A \code{dist} object.
#'
#' @return A data frame with three columns:
#'   \item{row}{Row label}
#'   \item{col}{Column label}
#'   \item{value}{Distance value}
#'
#' @examples
#' d <- dist(matrix(rnorm(100), nrow = 10))
#' df <- distAsDf(d)
#' head(df)
#'
#' @export
distAsDf <- function(inDist) {
    ## https://stackoverflow.com/questions/23474729/convert-object-of-class-dist-into-data-frame-in-r
    if (class(inDist) != "dist") stop("wrong input type")
    A <- attr(inDist, "Size")
    B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
    if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
    if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
    data.frame(
        row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
        col = rep(B[-length(B)], (length(B)-1):1),
        value = as.vector(inDist)
    ) ## end data.frame
}
