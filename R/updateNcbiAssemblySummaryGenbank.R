#' Update NCBI Assembly Summary GenBank Table
#'
#' Downloads the latest NCBI assembly summary from GenBank and updates
#' the database table \code{ncbi.assembly_summary_genbank}.
#'
#' @param dbConn A DBI database connection object.
#'
#' @return NULL (invisible). The function updates the database directly.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Downloads \code{assembly_summary_genbank.txt} from NCBI FTP
#'   \item Reads the tab-separated file
#'   \item Deletes the existing table contents
#'   \item Inserts the new data
#' }
#'
#' @importFrom DBI dbExecute dbWriteTable SQL
#'
#' @examples
#' \dontrun{
#' conn <- dbConnect(RPostgres::Postgres(), dbname = "taxogeno")
#' updateNcbiAssemblySummaryGenbank(conn)
#' }
#'
#' @export
updateNcbiAssemblySummaryGenbank<-function(dbConn){
    ## semiconstants
    assemblySummaryGenbankUrl<-"ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
    assemblySummaryGenbankFilePath<-"assembly_summary_genbank.txt"
    assemblySummaryGenbankColNames<-c("assembly_accession", "bioproject", "biosample", "wgs_master",
                                      "refseq_category", "taxid", "species_taxid", "organism_name",
                                      "infraspecific_name", "isolate", "version_status", "assembly_level",
                                      "release_type", "genome_rep", "seq_rel_date", "asm_name", "submitter",
                                      "gbrs_paired_asm", "paired_asm_comp", "ftp_path", "excluded_from_refseq",
                                      "relation_to_type_material")
    
    write("Downloading assembly summary",stdout())
    system2("wget", args=c(assemblySummaryGenbankUrl,"-O",assemblySummaryGenbankFilePath), wait=TRUE)
    write("Loading assembly summary table",stdout())
    assemblySummaryGenbankDf<-read.table(assemblySummaryGenbankFilePath,
                                         comment.char="#",
                                         header=FALSE,
                                         row.names=NULL,
                                         col.names = assemblySummaryGenbankColNames,
                                         colClasses = rep("character",22),
                                         fill=TRUE,
                                         quote=NULL,
                                         sep="\t",
                                         stringsAsFactors=FALSE)
    unlink(assemblySummaryGenbankFilePath)
    write("Inserting assembly summary table into database",stdout())
    dbExecute(dbConn,"DELETE FROM ncbi.assembly_summary_genbank")
    dbWriteTable(dbConn,SQL("ncbi.assembly_summary_genbank"), assemblySummaryGenbankDf, append=TRUE, row.names=FALSE)
}

