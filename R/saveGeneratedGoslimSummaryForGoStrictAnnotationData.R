#' Save GO Slim Summary for Strict GO Annotations
#'
#' Generates a GO Slim summary from strict GO annotations (C, F, P aspects)
#' and saves it to the database. This is the main function for creating
#' the simplified GO Slim representation of a proteome.
#'
#' @param dbConn A DBI database connection object.
#' @param proteomeId Integer. The proteome ID.
#' @param goStrictDf Data frame with strict GO annotations:
#'   \item{geneid}{Gene identifier}
#'   \item{annotkwid}{GO term ID}
#'   \item{aspect}{GO aspect (C, F, or P)}
#'
#' @return NULL (invisible). The function writes to the database.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Generates a GAF data frame from the strict GO annotations
#'   \item Maps the GAF to GO Slim using OWLTools
#'   \item Converts the mapped GAF to a summary data frame
#'   \item Writes the summary to \code{taxogeno.generated_goslim_summary}
#' }
#'
#' @importFrom DBI dbWriteTable SQL
#'
#' @examples
#' \dontrun{
#' conn <- dbConnect(RPostgres::Postgres(), dbname = "taxogeno")
#' saveGeneratedGoslimSummaryForGoStrictAnnotationData(conn, 1, go_df)
#' }
#'
#' @export
saveGeneratedGoslimSummaryForGoStrictAnnotationData<-function(dbConn,proteomeId,goStrictDf){
        gafDf<-generateGafDf(dbConn,goStrictDf)
        mappedGafDf<-owltoolsMap2SlimGafDf(gafDf,owl="goslim_generic.owl",subset="goslim_generic")
        summaryDf<-convertGafToSummaryDf(mappedGafDf, proteomeId)
        dbWriteTable(dbConn,SQL("taxogeno.generated_goslim_summary"), summaryDf,append=TRUE,row.names=FALSE)
}
