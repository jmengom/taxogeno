#' Save GO Slim-Based Euclidean Distances for a Proteome
#'
#' Calculates and caches Euclidean distances based on GO Slim annotations
#' between a specific proteome and all other proteomes in the database.
#'
#' @param dbConn A DBI database connection object.
#' @param proteomeId Integer. The proteome ID to compare.
#'
#' @return NULL (invisible). The function writes distances to
#'   \code{taxogeno.euclidean_distances_generated_goslim}.
#'
#' @details
#' This function is a wrapper around the more general
#' \code{\link{saveEuclideanDistancesForProteomeId}} specifically for
#' GO Slim annotations. It performs chunked processing to manage memory.
#'
#' @importFrom DBI dbGetQuery dbWriteTable SQL
#'
#' @examples
#' \dontrun{
#' conn <- dbConnect(RPostgres::Postgres(), dbname = "taxogeno")
#' saveGeneratedGoslimEuclideanDistancesForProteomeId(conn, 1)
#' }
#'
#' @export
saveGeneratedGoslimEuclideanDistancesForProteomeId<-function(dbConn,proteomeId){
    ## integrity tests
    ## only one proteome in proteomeId vector
    if(length(proteomeId)>1){
        stop("length(proteomeId)>1")
    }

    ## only proceed if proteomeId exists in databse
    existsProteomeId<-as.logical( dbGetQuery(dbConn,"SELECT EXISTS( SELECT proteomeid FROM taxogeno.proteome WHERE proteomeid=$1 )",
                                             params=list(proteomeId))[1,1] )
    if(!existsProteomeId){
        stop("!existsProteomeId")
    }


    ## examining candidates for euclidean distance comparision
    otherProteomeIdVec<-dbGetQuery(dbConn,"SELECT proteomeid FROM taxogeno.proteome WHERE proteomeid<>$1",
                                   params=list(proteomeId))[,"proteomeid",drop=TRUE]
    ## Porque no tiene sentido comparar entre sí los que ya están guardados
    ## Hay que sacar un vector con los que ya existen
    alreadyExistsVec<-dbGetQuery(dbConn,
                                "SELECT proteomeid_least as proteomeid
                                   FROM taxogeno.taxogeno.euclidean_distances_generated_goslim
                                  WHERE proteomeid_least=$1
                                  UNION
                                 SELECT proteomeid_greatest as proteomeid
                                   FROM taxogeno.taxogeno.euclidean_distances_generated_goslim
                                  WHERE proteomeid_greatest=$1",
                                params=list(proteomeId))[,"proteomeid"]

    ## Y dejar únicamente los que no estén ya comparados entre sí
    otherProteomeIdVecFiltered<-setdiff(otherProteomeIdVec, alreadyExistsVec)
    
    ## Como puede petar la memoria al generar la matriz dispersa y el objeto dist, ir de 20 en 20
    for(chunkVec in chunkVectorInEqualSizeFragments(otherProteomeIdVecFiltered, n=20)){
        ## Se añade el proteoma que se está introduciendo al trozo de ids de proteoma contra los que se va a comparar
        comparingProteomeIdVec<-c(proteomeId,chunkVec)
        ## Obtener las anotaciones normalizadas para cada uno de ellos en un dataframe en formato mapa
        normalizedGeneratedGoslimAnnotationDf<-getGeneratedGoslimNormalizedAnnotationForProteomeIdVec(dbConn, comparingProteomeIdVec)
        ## Calcular la matriz dispersa reconvertida a data frame que tiene las distancias euclidianas
        euclideanDistancesDf<-calculateEuclideanDistances(normalizedGeneratedGoslimAnnotationDf)
        ## Sólo añadir las distancias calculadas para el proteoma que se está introduciendo, porque si no se viola la clave primaria
        dbWriteTable( dbConn,
                     SQL("taxogeno.euclidean_distances_generated_goslim"),
                     euclideanDistancesDf [euclideanDistancesDf[,"proteomeid_greatest"]==proteomeId | euclideanDistancesDf[,"proteomeid_least"]==proteomeId , ],
                     append=TRUE, row.names=FALSE )
    }
}
