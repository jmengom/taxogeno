#' Calculate Euclidean Distances Between Proteomes
#'
#' Computes Euclidean distances between proteomes based on their normalized
#' annotation vectors (GO Slim or Keywords). Uses sparse matrix representation
#' for memory efficiency.
#'
#' @param normalizedAnnotationDf Data frame with columns:
#'   \item{proteomeid}{Proteome identifier}
#'   \item{annotkwid}{Annotation term ID}
#'   \item{annotkwrelval}{Normalized annotation value (count/genecount)}
#'
#' @return A data frame with distance information:
#'   \item{proteomeid_greatest}{First proteome ID (greater numeric value)}
#'   \item{proteomeid_least}{Second proteome ID (lesser numeric value)}
#'   \item{distance}{Euclidean distance between the two proteomes}
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Builds a sparse matrix: proteomeid × annotkwid
#'   \item Calculates Euclidean distance between rows (proteomes)
#'   \item Converts the result to a data frame
#' }
#'
#' @examples
#' \dontrun{
#' norm_df <- data.frame(
#'   proteomeid = c(1, 1, 2, 2),
#'   annotkwid = c("GO:0003674", "GO:0005575", "GO:0003674", "GO:0008150"),
#'   annotkwrelval = c(0.1, 0.05, 0.2, 0.15)
#' )
#' distances <- calculateEuclideanDistances(norm_df)
#' }
#'
#' @export
calculateEuclideanDistances<-function (normalizedAnnotationDf){   
    euclideanDistancesDf<-distAsDf(
        dist(
            ## ####################
            ##            kw
            ## proteomeid a b c d
            ##          1 1 . . .
            ##          2 . 3 . .
            ##          3 . . 5 .
            ##          4 . . . 7
            xtabs( annotkwrelval~proteomeid+annotkwid, normalizedAnnotationDf, sparse=TRUE),
            ## De esa matriz dispersa la distancia euclídea entre las filas
            method="euclidean",
            diag=FALSE,
            upper=FALSE
        ) ## end dist
    ) ## end distAsDf
    colnames(euclideanDistancesDf)<-c("proteomeid_greatest","proteomeid_least","distance")
    euclideanDistancesDf
}
