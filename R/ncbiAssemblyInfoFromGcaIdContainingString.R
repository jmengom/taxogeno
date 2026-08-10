#' Retrieve NCBI Assembly Information from GCA ID
#'
#' Queries the NCBI assembly summary database for information about an assembly
#' identified by a GCA accession number.
#'
#' @param dbConn A DBI database connection object.
#' @param gcaIdContainingString Character string containing a GCA ID
#'   (e.g., "GCA_000001405.28_GRCh38.p14"). The function extracts the GCA ID
#'   using regex.
#'
#' @return A list containing assembly information from the NCBI summary table:
#'   \item{assembly_accession}{GCA accession}
#'   \item{bioproject}{BioProject ID}
#'   \item{biosample}{BioSample ID}
#'   \item{wgs_master}{WGS master accession}
#'   \item{refseq_category}{RefSeq category}
#'   \item{taxid}{NCBI taxonomy ID}
#'   \item{species_taxid}{Species taxonomy ID}
#'   \item{organism_name}{Scientific name}
#'   \item{infraspecific_name}{Strain/isolate name}
#'   \item{isolate}{Isolate identifier}
#'   \item{version_status}{Assembly version status}
#'   \item{assembly_level}{Assembly level (Complete, Chromosome, Scaffold, Contig)}
#'   \item{release_type}{Release type}
#'   \item{genome_rep}{Genome representation}
#'   \item{date}{Release date}
#'   \item{asm_name}{Assembly name}
#'   \item{submitter}{Submitting organization}
#'   \item{gbrs_paired_asm}{GenBank/RefSeq paired assembly}
#'   \item{paired_asm_comp}{Paired assembly comparison}
#'   \item{ftp_path}{FTP path for data download}
#'   \item{excluded_from_refseq}{Excluded from RefSeq flag}
#'   \item{relation_to_type_material}{Type material relationship}
#'
#' @importFrom DBI dbGetQuery
#'
#' @examples
#' \dontrun{
#' conn <- dbConnect(RPostgres::Postgres(), dbname = "taxogeno")
#' info <- ncbiAssemblyInfoFromGcaIdContainingString(conn, "GCA_000001405.28")
#' print(info$organism_name)
#' }
#'
#' @export
ncbiAssemblyInfoFromGcaIdContainingString<-function(dbConn,gcaIdContainingString){
    ## extractedGcaId<-gsub('^.*(GCA_[^.]*)[.].*$','\\1',gcaIdContainingString)
    ## assemblyInfoDf<-dbGetQuery(
    ##     dbConn,
    ##     "SELECT
    ##       *
    ##      FROM ncbi.assembly_summary_genbank
    ##      WHERE assembly_accession LIKE $1||'%'",
    ##     params=list(extractedGcaId))

    extractedGcaId<-gsub('^.*(GCA_[0123456789]+[.][0123456789]+).*$','\\1',gcaIdContainingString)
    print("Extracted gcaid")
    print(extractedGcaId)
    assemblyInfoList<-dbGetQuery(
        dbConn,
        "SELECT
           assembly_accession,
           bioproject,
           biosample,
           wgs_master,
           refseq_category,
           taxid,
           species_taxid,
           organism_name,
           infraspecific_name,
           isolate,
           version_status,
           assembly_level,
           release_type,
           genome_rep,
           seq_rel_date date,
           asm_name,
           submitter,
           gbrs_paired_asm,
           paired_asm_comp,
           ftp_path,
           excluded_from_refseq,
           relation_to_type_material
     FROM ncbi.assembly_summary_genbank
     WHERE assembly_accession IN($1)",
     params=list(extractedGcaId))[1,]
    print("AssemblyInfoList")
    print(assemblyInfoList)
    assemblyInfoList
}
