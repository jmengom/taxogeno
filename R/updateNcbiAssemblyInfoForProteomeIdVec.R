#' Update NCBI Assembly Information for Proteomes
#'
#' Updates NCBI assembly information (GCA, taxonomy ID, source URL) for
#' a set of proteomes by querying the NCBI assembly summary database.
#'
#' @param dbConn A DBI database connection object.
#' @param proteomeIdVec Integer vector. Proteome IDs to update.
#'
#' @return NULL (invisible). The function updates the database directly.
#'
#' @details
#' For each proteome:
#' \enumerate{
#'   \item Retrieves existing GCA ID or uses basename as fallback
#'   \item Queries NCBI assembly summary for updated information
#'   \item Updates taxonomy if taxon has changed
#'   \item Inserts new taxon if it doesn't exist
#'   \item Updates proteome record with new assembly info
#' }
#'
#' @importFrom DBI dbGetQuery dbExecute
#'
#' @examples
#' \dontrun{
#' conn <- dbConnect(RPostgres::Postgres(), dbname = "taxogeno")
#' updateNcbiAssemblyInfoForProteomeIdVec(conn, c(1, 2, 3))
#' }
#'
#' @export
updateNcbiAssemblyInfoForProteomeIdVec<-function(dbConn,proteomeIdVec){
    ## sacar la información del/de los proteomas de la base de datos
    proteomeInfoDf<-dbGetQuery(dbConn,"SELECT * FROM taxogeno.proteome WHERE proteomeid IN ($1)", params=list(proteomeIdVec))

    for(rowNum in seq_along(nrow(proteomeInfoDf))){
        gcaIdContainingString=NA
        ## si se proporcionó gcaid usarlo primero
        if(!is.na(proteomeInfoDf[rowNum,][["gcaid"]])){
            gcaIdContainingString<-proteomeInfoDf[rowNum,"gcaid"]
        }
        ## si no hay gcaid, tirar del nombre del fichero
        else if(!is.na(proteomeInfoDf[rowNum,"basename"])){
            gcaIdContainingString<-proteomeInfoDf[rowNum,"basename"]
        }
        
        ## tomar la información nueva
        ncbiAssemblyInfoList<-ncbiAssemblyInfoFromGcaIdContainingString(dbConn, gcaIdContainingString)
        str(ncbiAssemblyInfoList)
        if(length(ncbiAssemblyInfoList)>0){
            print("Entra en lista de información ncbi")
            ## ^############
            ## comprobar si había un taxón antiguo almacenado
            oldNcbiTaxId<-dbGetQuery(dbConn,
                                     "SELECT ncbitaxid FROM taxogeno.taxonomy WHERE ncbitaxid = $1 AND is_ancestor=FALSE",
                                     params=list(proteomeInfoDf[rowNum,"ncbitaxid"]) ) [1,"ncbitaxid"]

            if(!is.na(oldNcbiTaxId) && oldNcbiTaxId != ncbiAssemblyInfoList[["taxid"]]){
                dbExecute(dbConn, "CALL taxogeno.delete_taxon($1)", params=list(oldNcbiTaxId))
            }
            ## $##########

            ## ^############
            ## comprobar si había un taxón antiguo almacenado igual que el que se quiere insertar
            oldNcbiTaxId<-dbGetQuery(dbConn,
                                     "SELECT ncbitaxid FROM taxogeno.taxonomy WHERE ncbitaxid = $1 AND is_ancestor=FALSE",
                                     params=list(ncbiAssemblyInfoList[["taxid"]]) ) [1,"ncbitaxid"]
            
            if(is.na(oldNcbiTaxId)){
                dbExecute(dbConn, "CALL taxogeno.insert_taxon($1)", params=list(ncbiAssemblyInfoList[["taxid"]]))
            }
            ## $############

            dbExecute(dbConn,
                      "UPDATE taxogeno.proteome SET (gcaid,ncbitaxid,sourceurl,dbname) = ($2,$3,$4,$5) WHERE proteomeid=$1",
                      params=list(proteomeInfoDf[rowNum,"proteomeid"],
                                  ncbiAssemblyInfoList[["assembly_accession"]],
                                  ncbiAssemblyInfoList[["taxid"]],
                                  ncbiAssemblyInfoList[["ftp_path"]],
                                  "NCBI:GenBank") )
        }
    }
}
