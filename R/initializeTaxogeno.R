#' Initialize Taxogeno Database with External Reference Data
#'
#' Performs a complete initialization of the taxogeno database by downloading
#' and installing all external reference datasets required for the system to
#' function properly. This is a one-time setup step that must be completed
#' before using the taxogeno package for proteome analysis.
#'
#' @param dbConn A DBI database connection object to the taxogeno database.
#'
#' @return NULL (invisible). The function updates the database directly.
#'
#' @details
#' This function sequentially updates four essential external data sources:
#'
#' \enumerate{
#'   \item \strong{NCBI Assembly Summary GenBank}:
#'     Downloads the complete assembly summary from NCBI GenBank containing
#'     metadata for all publicly available genome assemblies. This is stored
#'     in the \code{ncbi.assembly_summary_genbank} table and is used for
#'     mapping GCA accessions to taxonomic and assembly information.
#'
#'   \item \strong{NCBI Taxonomy (BioSQL)}:
#'     Downloads the complete NCBI taxonomy dump and loads it into the
#'     \code{biosql} schema using the BioSQL taxonomy loading script.
#'     This provides the full taxonomic hierarchy used for phylogenetic
#'     analysis and ancestor/descendant queries.
#'
#'   \item \strong{UniProt Keywords}:
#'     Downloads both OBO and TSV formats of UniProt keywords, which are
#'     used for functional annotation of proteins. These are stored in the
#'     \code{uniprot.uniprot_keyword} and \code{uniprot.uniprot_keyword_category}
#'     tables.
#'
#'   \item \strong{Gene Ontology (GO)}:
#'     Downloads the complete GO ontology in OWL format along with the
#'     GO Slim generic subset. These are parsed using PostgreSQL's XML
#'     capabilities and stored in the \code{gene_ontology} schema for
#'     annotation mapping and enrichment analysis.
#' }
#'
#' @section Dependencies:
#' \strong{System Requirements:}
#' \itemize{
#'   \item \code{wget}: Required for downloading files from FTP/HTTP servers
#'   \item \code{Perl}: Required for running the BioSQL taxonomy loading script
#'   \item \code{PostgreSQL} with \code{XML} extension enabled
#' }
#'
#' \strong{R Packages (already imported):}
#' \itemize{
#'   \item \code{DBI}: For database operations
#'   \item \code{RPostgres}: For PostgreSQL connectivity
#' }
#'
#' @section Time and Resources:
#' The initialization process is resource-intensive:
#' \itemize{
#'   \item \strong{Download size}: Approximately 350-400 MB total
#'   \item \strong{Database size}: Approximately 500-700 MB after insertion
#'   \item \strong{Time}: 10-30 minutes depending on network speed and system performance
#' }
#'
#' @section Error Handling:
#' If any of the four updates fails:
#' \itemize{
#'   \item The function will stop execution with an error message
#'   \item The database may be left in a partially updated state
#'   \item You can rerun the function to complete the failed update
#'   \item Check network connectivity and database permissions if errors occur
#' }
#'
#' @section Usage Recommendations:
#' \itemize{
#'   \item Run this function once during initial system setup
#'   \item Run monthly to keep reference data current with NCBI/GO releases
#'   \item Use with \code{\link{withProgress}} for progress tracking in Shiny apps
#' }
#'
#' @references
#' \itemize{
#'   \item NCBI Assembly: \url{ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/}
#'   \item NCBI Taxonomy: \url{ftp://ftp.ncbi.nih.gov/pub/taxonomy/}
#'   \item UniProt Keywords: \url{http://www.uniprot.org/keywords/}
#'   \item Gene Ontology: \url{http://current.geneontology.org/ontology/}
#' }
#'
#' @seealso
#' \code{\link{updateNcbiAssemblySummaryGenbank}} for NCBI assembly updates
#' \code{\link{updateBiosqlNcbiTaxonomy}} for taxonomy updates
#' \code{\link{updateUniprotKeyword}} for UniProt keyword updates
#' \code{\link{updateGeneOntology}} for GO updates
#' \code{\link{saveSma3sFileSet}} for loading proteome data after initialization
#'
#' @examples
#' \dontrun{
#' # Connect to the database
#' library(DBI)
#' library(RPostgres)
#' 
#' conn <- dbConnect(RPostgres::Postgres(), 
#'                    dbname = "taxogeno",
#'                    host = "localhost",
#'                    port = 5432,
#'                    user = "postgres",
#'                    password = "password")
#'
#' # Initialize the database (takes 10-30 minutes)
#' initializeTaxogeno(conn)
#'
#' # Verify the installation
#' dbGetQuery(conn, "SELECT COUNT(*) FROM ncbi.assembly_summary_genbank")
#' dbGetQuery(conn, "SELECT COUNT(*) FROM biosql.taxon")
#' dbGetQuery(conn, "SELECT COUNT(*) FROM uniprot.uniprot_keyword")
#' dbGetQuery(conn, "SELECT COUNT(*) FROM gene_ontology.gene_ontology")
#'
#' # Disconnect when done
#' dbDisconnect(conn)
#' }
#'
#' @export
initializeTaxogeno <- function(dbConn) {
    # Log start of initialization
    message("========================================")
    message("Starting Taxogeno Database Initialization")
    message("========================================")
    message("This process will take 10-30 minutes.")
    message("Please ensure you have a stable internet connection.")
    message("")
    
    # 1. Update NCBI Assembly Summary
    message("[1/4] Updating NCBI Assembly Summary GenBank...")
    message("Downloading assembly_summary_genbank.txt (~50 MB)...")
    start_time <- Sys.time()
    updateNcbiAssemblySummaryGenbank(dbConn)
    elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
    message(sprintf("[1/4] Complete! (%.2f minutes)", elapsed))
    message("")
    
    # 2. Update NCBI Taxonomy
    message("[2/4] Updating NCBI Taxonomy (BioSQL)...")
    message("Downloading taxdump.tar.gz (~100 MB)...")
    start_time <- Sys.time()
    updateBiosqlNcbiTaxonomy(dbConn)
    elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
    message(sprintf("[2/4] Complete! (%.2f minutes)", elapsed))
    message("")
    
    # 3. Update UniProt Keywords
    message("[3/4] Updating UniProt Keywords...")
    message("Downloading keyword.obo and keyword.tsv (~5 MB)...")
    start_time <- Sys.time()
    updateUniprotKeyword(dbConn)
    elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
    message(sprintf("[3/4] Complete! (%.2f minutes)", elapsed))
    message("")
    
    # 4. Update Gene Ontology
    message("[4/4] Updating Gene Ontology (GO)...")
    message("Downloading go.owl and goslim_generic.owl (~200 MB)...")
    message("Parsing OWL with PostgreSQL XMLTABLE...")
    start_time <- Sys.time()
    updateGeneOntology(dbConn)
    elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
    message(sprintf("[4/4] Complete! (%.2f minutes)", elapsed))
    message("")
    
    # Log completion
    total_time <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
    message("========================================")
    message("Taxogeno Database Initialization Complete!")
    message(sprintf("Total time: %.2f minutes", total_time))
    message("========================================")
    message("")
    message("The following external databases have been updated:")
    message("  - NCBI Assembly Summary GenBank")
    message("  - NCBI Taxonomy (BioSQL)")
    message("  - UniProt Keywords")
    message("  - Gene Ontology (GO and GO Slim)")
    message("")
    message("You can now load proteomes with saveSma3sFileSet()")
    
    invisible(NULL)
}
