#' Delete Proteomes from the Database
#'
#' Deletes one or more proteomes and all associated data using the
#' \code{taxogeno.delete_proteome} stored procedure.
#'
#' @param dbConn A DBI database connection object.
#' @param proteomeIdVec Integer vector. Proteome IDs to delete.
#'
#' @return NULL (invisible). The function deletes data from the database.
#'
#' @details
#' The deletion is performed via a PostgreSQL stored procedure that
#' handles cascade deletion of all related records including:
#' \itemize{
#'   \item GO Slim summaries
#'   \item Distance matrices
#'   \item Gene records
#'   \item All annotation relationships (GO, Keywords, EC, Pathways)
#' }
#'
#' @importFrom DBI dbExecute
#'
#' @examples
#' \dontrun{
#' conn <- dbConnect(RPostgres::Postgres(), dbname = "taxogeno")
#' deleteProteome(conn, c(1, 2, 3))
#' }
#'
#' @export
deleteProteome<-function(dbConn,proteomeIdVec){
    for(proteomeId in proteomeIdVec){
        dbExecute(dbConn,"CALL taxogeno.delete_proteome($1)",params=list(proteomeId))
    }
}
