#' Split Vector into Equal-Sized Fragments
#'
#' Splits a vector into fragments of approximately equal size. This is
#' used for memory-efficient processing when comparing proteomes in chunks.
#'
#' @param x A vector to split.
#' @param n Integer. The maximum size of each fragment.
#'
#' @return A list of vector fragments.
#'
#' @examples
#' chunkVectorInEqualSizeFragments(1:100, n = 20)
#'
#' @export
chunkVectorInEqualSizeFragments <- function(x,n) {
        split(x, ceiling(seq_along(x)/n))
}
