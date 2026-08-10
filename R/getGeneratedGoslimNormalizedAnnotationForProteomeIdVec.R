#' Get Generated GO Slim Normalized Annotation for Proteomes
#'
#' Retrieves normalized GO Slim annotation vectors for a set of proteomes
#' from the generated_goslim_summary table.
#'
#' @param dbConn A DBI database connection object.
#' @param proteomeIdVec Integer vector. Proteome IDs to retrieve.
#'
#' @return A data frame with columns:
#'   \item{proteomeid}{Proteome identifier}
#'   \item{annotkwid}{GO Slim term ID}
#'   \item{annotkwcount}{Raw annotation count}
#'   \item{annotkwrelval}{Normalized value (count/genecount)}
#'
#' @importFrom DBI dbGetQuery
#'
#' @examples
#' \dontrun{
#' norm_df <- getGeneratedGoslimNormalizedAnnotationForProteomeIdVec(conn, c(1, 2, 3))
#' }
#'
#' @export
getGeneratedGoslimNormalizedAnnotationForProteomeIdVec<-function(dbConn,proteomeIdVec){
    dbGetQuery(dbConn,
               "SELECT
                  taxogeno.proteome.proteomeid as proteomeid,
                  annotkwid                     as annotkwid,
                  annotkwcount                  as annotkwcount,
                  annotkwcount::float/genecount as annotkwrelval

                FROM taxogeno.generated_goslim_summary
                INNER JOIN taxogeno.proteome
                  ON taxogeno.generated_goslim_summary.proteomeid=taxogeno.proteome.proteomeid
                  AND taxogeno.generated_goslim_summary.proteomeid IN ($1)

                ORDER BY taxogeno.proteome.proteomeid, annotkwid", params=list(proteomeIdVec))
}
