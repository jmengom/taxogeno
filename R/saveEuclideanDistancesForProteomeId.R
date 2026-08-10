#' Save Euclidean Distances for a Proteome
#'
#' Calculates and caches Euclidean distances between a specific proteome
#' and all other proteomes in the database for the specified annotation types.
#' Uses chunking for memory-efficient processing.
#'
#' @param dbConn A DBI database connection object.
#' @param proteomeId Integer. The proteome ID to compare.
#' @param annotationTypeVec Character vector. Types of annotations to use:
#'   "generated_goslim" and/or "keyword".
#'
#' @return NULL (invisible). The function writes distances to database tables.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Validates the input parameters
#'   \item Identifies proteomes not yet compared with the target
#'   \item Processes comparisons in chunks to manage memory
#'   \item Writes results to \code{euclidean_distances_generated_goslim}
#'     and/or \code{euclidean_distances_keyword}
#' }
#'
#' @importFrom DBI dbGetQuery dbWriteTable SQL
#'
#' @examples
#' \dontrun{
#' conn <- dbConnect(RPostgres::Postgres(), dbname = "taxogeno")
#' saveEuclideanDistancesForProteomeId(conn, 1, c("generated_goslim", "keyword"))
#' }
#'
#' @export
saveEuclideanDistancesForProteomeId <- function(dbConn, proteomeId, annotationTypeVec) {
    ## Integrity tests
    allowedAnnotationTypes <- c("generated_goslim", "keyword")
    if (!(all(annotationTypeVec %in% allowedAnnotationTypes))) {
        stop("!(annotationTypeVec %in% allowedAnnotationTypes)")
    }
    
    ## Only one proteome in proteomeId vector
    if (length(proteomeId) > 1) {
        stop("length(proteomeId) > 1")
    }

    ## Only proceed if proteomeId exists in database
    existsProteomeId <- as.logical(dbGetQuery(
        dbConn,
        "SELECT EXISTS( SELECT proteomeid FROM taxogeno.proteome WHERE proteomeid = $1 )",
        params = list(proteomeId)
    )[1, 1])
    
    if (!existsProteomeId) {
        stop("!existsProteomeId")
    }
    
    ## Process each annotation type
    for (annotationType in annotationTypeVec) {
        ## Get all potential proteomes for comparison
        otherProteomeIdVec <- dbGetQuery(
            dbConn,
            "SELECT proteomeid FROM taxogeno.proteome"
        )[, "proteomeid", drop = TRUE]
        
        ## Get proteomes that have already been compared with the target
        ## It doesn't make sense to compare against those already stored
        alreadyExistsVec <- dbGetQuery(
            dbConn,
            sprintf(
                "SELECT proteomeid_least AS proteomeid
                 FROM taxogeno.euclidean_distances_%1$s
                 WHERE proteomeid_least = $1
                 UNION
                 SELECT proteomeid_greatest AS proteomeid
                 FROM taxogeno.euclidean_distances_%1$s
                 WHERE proteomeid_greatest = $1",
                annotationType
            ),
            params = list(proteomeId)
        )[, "proteomeid"]

        ## Keep only those not yet compared with the target
        otherProteomeIdVecFiltered <- setdiff(otherProteomeIdVec, alreadyExistsVec)
        
        ## Process in chunks to avoid memory issues when generating
        ## the sparse matrix and dist object (chunk size = 20)
        for (chunkVec in chunkVectorInEqualSizeFragments(otherProteomeIdVecFiltered, n = 20)) {
            ## Add the target proteome to the chunk for comparison
            comparingProteomeIdVec <- c(proteomeId, chunkVec)

            ## Get normalized annotations for all proteomes in this chunk
            normalizedKeywordAnnotationDf <- getNormalizedAnnotationForProteomeIdVec(
                dbConn,
                annotationType,
                comparingProteomeIdVec
            )

            ## Calculate Euclidean distances from the sparse matrix
            euclideanDistancesDf <- calculateEuclideanDistances(normalizedKeywordAnnotationDf)

            ## Only insert distances involving the target proteome
            ## to avoid violating the primary key constraint
            dbWriteTable(
                dbConn,
                SQL(sprintf("taxogeno.euclidean_distances_%s", annotationType)),
                euclideanDistancesDf[
                    euclideanDistancesDf[, "proteomeid_greatest"] == proteomeId |
                    euclideanDistancesDf[, "proteomeid_least"] == proteomeId,
                ],
                append = TRUE,
                row.names = FALSE
            )
        } ## end for(chunkVec in chunkVectorInEqualSizeFragments(...))
    } ## end for(annotationType in annotationTypeVec)
}

