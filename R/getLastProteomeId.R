#' Get the Last (Maximum) Proteome ID
#'
#' Retrieves the highest proteome ID from the proteome table.
#'
#' @param dbConn A DBI database connection object.
#'
#' @return Integer. The maximum proteome ID, or NA if no proteomes exist.
#'
#' @importFrom DBI dbGetQuery
#'
#' @examples
#' \dontrun{
#' conn <- dbConnect(RPostgres::Postgres(), dbname = "taxogeno")
#' last_id <- getLastProteomeId(conn)
#' }
#'
#' @export
getLastProteomeId<-function(dbConn){
    dbGetQuery(dbConn,"SELECT MAX(proteomeid) as max_proteomeid FROM taxogeno.proteome")[1,"max_proteomeid"]
}

