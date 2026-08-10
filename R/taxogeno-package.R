#' taxogeno: Comparative Phylogenetic Analysis Based on Functional Annotations
#'
#' @description
#' Taxogeno is an ETL (Extract-Transform-Load) system for processing, storing,
#' and analyzing complete proteomes functionally annotated by the Sma3s program.
#' It builds a normalized PostgreSQL database with:
#'
#' \itemize{
#'   \item Functional annotations (GO, GO Slim, Keywords, EC, Pathways)
#'   \item Protein sequences
#'   \item Taxonomic relationships
#'   \item Euclidean distance matrices between proteomes
#' }
#'
#' @section Key Functions:
#'
#' **Data Import:**
#' \itemize{
#'   \item \code{\link{read_sma3s_annotation}}: Read Sma3s TSV annotation files
#'   \item \code{\link{read_multifasta}}: Read large FASTA files line by line
#' }
#'
#' **Data Transformation:**
#' \itemize{
#'   \item \code{\link{extract_column_df}}: Convert semicolon-separated fields to 1:N relations
#'   \item \code{\link{generate_gaf_df}}: Generate Gene Association Format data
#'   \item \code{\link{owltools_map2slim_gaf_df}}: Map GO terms to GO Slim using OWLTools
#' }
#'
#' **Database Operations:**
#' \itemize{
#'   \item \code{\link{save_sma3s_fileset}}: Main ETL pipeline function
#'   \item \code{\link{update_ncbi_assembly_info}}: Update NCBI assembly metadata
#'   \item \code{\link{update_ncbi_assembly_summary_genbank}}: Download NCBI assembly summaries
#'   \item \code{\link{update_uniprot_keyword}}: Download UniProt keywords
#'   \item \code{\link{update_gene_ontology}}: Download and parse GO ontologies
#' }
#'
#' **Distance Calculations:**
#' \itemize{
#'   \item \code{\link{calculate_euclidean_distances}}: Compute Euclidean distances between proteomes
#'   \item \code{\link{get_normalized_annotation}}: Retrieve normalized annotation vectors
#' }
#'
#' **Proteome Management:**
#' \itemize{
#'   \item \code{\link{update_proteome_field}}: Update proteome metadata
#'   \item \code{\link{delete_proteome}}: Delete a proteome and all associated data
#' }
#'
#' @section Database Schema:
#' The package uses the following main schemas:
#' \itemize{
#'   \item \code{taxogeno}: Core data (proteomes, genes, annotations, distances)
#'   \item \code{biosql}: NCBI taxonomy (BioSQL schema with recursive functions)
#'   \item \code{ncbi}: Assembly summary data
#'   \item \code{uniprot}: UniProt keywords and proteome catalogs
#'   \item \code{gene_ontology}: GO and GO Slim ontologies (parsed from OWL)
#' }
#'
#' @docType package
#' @name taxogeno
#' @keywords internal
"_PACKAGE"
