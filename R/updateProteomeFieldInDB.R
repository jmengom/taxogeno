#' Update a Proteome Field in the Database
#'
#' Updates a specific field in the proteome table for a single proteome.
#'
#' @param dbConn A DBI database connection object.
#' @param proteomeId Integer. The proteome ID to update.
#' @param fieldName Character string. The field name to update.
#'   Allowed values: "is_userproteome", "gcaid", "ncbitaxid", "scientific_name".
#' @param newValue The new value for the field.
#'
#' @return NULL (invisible). The function updates the database directly.
#'
#' @importFrom DBI dbExecute
#'
#' @examples
#' \dontrun{
#' conn <- dbConnect(RPostgres::Postgres(), dbname = "taxogeno")
#' updateProteomeFieldInDB(conn, 1, "is_userproteome", TRUE)
#' }
#'
#' @export
updateProteomeFieldInDB<-function(dbConn,proteomeId,fieldName, newValue){
    if(length(proteomeId)>1){
        stop("length(proteomeId)>1")
    }
    allowedFieldNames<-c("is_userproteome", "gcaid", "ncbitaxid", "scientific_name")
    if(!(all(fieldName %in% allowedFieldNames))){
        error("!(fieldName %in% allowedFieldNames)")
    }
    dbExecute(dbConn, sprintf("UPDATE taxogeno.proteome SET %s = $1 WHERE proteomeid = $2", fieldName), params(list(newValue,proteomeId)))
    
}
