#' Update BioSQL NCBI Taxonomy
#'
#' Downloads the latest NCBI taxonomy dump and loads it into the biosql
#' schema using the \code{load_ncbi_taxonomy.pl} Perl script.
#'
#' @param dbConn A DBI database connection object.
#'
#' @return NULL (invisible). The function updates the database directly.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Downloads \code{taxdump.tar.gz} from NCBI FTP
#'   \item Extracts the archive
#'   \item Executes the Perl script \code{load_ncbi_taxonomy.pl}
#'   \item Cleans up temporary files
#' }
#'
#' @examples
#' \dontrun{
#' conn <- dbConnect(RPostgres::Postgres(), dbname = "taxogeno")
#' updateBiosqlNcbiTaxonomy(conn)
#' }
#'
#' @export
updateBiosqlNcbiTaxonomy<-function(dbConn){
    ncbiTaxonomyDumpUrl="ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
    ncbiTaxonomyDumpFile="taxdump.tar.gz"

    write("Downloading taxonomy dump",stdout())
    system2("wget", args=c(ncbiTaxonomyDumpUrl,"-O",ncbiTaxonomyDumpFile), wait=TRUE)

    write("Uncompressing taxonomy dump",stdout())
    taxonomyDumpFileList<-untar(ncbiTaxonomyDumpFile,list=TRUE)
    untar(ncbiTaxonomyDumpFile, exdir=dirname(ncbiTaxonomyDumpFile))
    unlink(ncbiTaxonomyDumpFile)
    
    write("Loading taxonomy dump into biosql schema of database",stdout())
    system2("perl",
            args=c("./bin/load_ncbi_taxonomy.pl","--driver","Pg","--schema","biosql","--dbname" ,"taxogeno","--directory",dirname(ncbiTaxonomyDumpFile)),
            wait=TRUE)
}
