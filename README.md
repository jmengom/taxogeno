# taxogeno

**Comparative Phylogenetic Analysis System Based on Functional Annotations of Complete Proteomes**

---

## Description

**taxogeno** is a comprehensive R package implementing an ETL (Extract-Transform-Load) system for processing, storing, and analyzing complete proteomes functionally annotated by [Sma3s](https://github.com/UPOBioinfo/sma3s) — developed by Dr. Antonio J. Pérez Pulido at the UPO Department of Bioinformatics. The package builds a normalized PostgreSQL database with:

- Functional annotations (GO, GO Slim, Keywords, EC, Pathways)
- Protein sequences
- Taxonomic relationships
- Euclidean distance matrices between proteomes

The system was developed during the **Higher Degree in Web Application Development** training period at the **Department of Bioinformatics, Universidad Pablo de Olavide**.

---

## System Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                          TAXOGENO PIPELINE                             │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  ┌───────────────┐     ┌───────────────┐     ┌────────────────────┐    │
│  │  Sma3s TSV    │────▶│  readSma3s    │────▶│  PostgreSQL        │    │
│  │  (annotation) │     │  Annotation   │     │  (taxogeno db)     │    │
│  └───────────────┘     └───────────────┘     └────────────────────┘    │
│                                                                         │
│  ┌───────────────┐     ┌───────────────┐     ┌────────────────────┐    │
│  │  FASTA        │────▶│  readMultif-  │────▶│  Tables:           │    │
│  │  (proteome)   │     │  asta         │     │  - proteome        │    │
│  └───────────────┘     └───────────────┘     │  - gene            │    │
│                                              │  - gene_goc_rel    │    │
│  ┌───────────────┐     ┌───────────────┐     │  - gene_gof_rel    │    │
│  │  OWLTools     │────▶│  map2slim     │────▶│  - gene_gop_rel    │    │
│  │  (GO Slim)    │     │               │     │  - gene_keyword_rel│    │
│  └───────────────┘     └───────────────┘     │  - gene_enzyme_rel │    │
│                                              │  - gene_pathway_rel│    │
│  ┌───────────────┐     ┌───────────────┐     │  - generated_goslim│    │
│  │  Euclidean    │────▶│  Distance     │────▶│    _summary        │    │
│  │  Distance     │     │  Matrix       │     │  - euclidean_      │    │
│  │  Calculation  │     │  (cached)     │     │    distances_*     │    │
│  └───────────────┘     └───────────────┘     └────────────────────┘    │
│                                              │                         │
│                                              ▼                         │
│                              ┌───────────────────────────────┐         │
│                              │  Shiny Application (R/Shiny)  │         │
│                              │  - Interactive queries        │         │
│                              │  - Tree visualization         │         │
│                              │  - Proteome comparison        │         │
│                              └───────────────────────────────┘         │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Package Structure

```
taxogeno/
├── R/                              # Core R functions
│   ├── taxogeno-package.R          # Package documentation
│   ├── readSma3sAnnotation.R       # Read Sma3s TSV files
│   ├── readMultifasta.R            # Read FASTA files (streaming)
│   ├── extractColumnDfFromInsertedGeneAnnotationDf.R
│   ├── ncbiAssemblyInfoFromGcaIdContainingString.R
│   ├── updateNcbiAssemblyInfoForProteomeIdVec.R
│   ├── generateGafDf.R             # Gene Association Format
│   ├── owltoolsMap2SlimGafDf.R     # GO Slim mapping
│   ├── convertGafToSummaryDf.R
│   ├── saveGeneratedGoslimSummaryForGoStrictAnnotationData.R
│   ├── chunkVectorInEqualSizeFragments.R
│   ├── getNormalizedAnnotationForProteomeIdVec.R
│   ├── getKeywordNormalizedAnnotationForProteomeIdVec.R
│   ├── getGeneratedGoslimNormalizedAnnotationForProteomeIdVec.R
│   ├── distAsDf.R                  # Convert dist objects
│   ├── calculateEuclideanDistances.R
│   ├── saveEuclideanDistancesForProteomeId.R
│   ├── saveKeywordEuclideanDistancesForProteomeId.R
│   ├── saveGeneratedGoslimEuclideanDistancesForProteomeId.R
│   ├── saveSma3sFileSet.R          # Main ETL pipeline
│   ├── updateProteomeFieldInDB.R
│   ├── getLastProteomeId.R
│   ├── deleteProteome.R
│   ├── updateNcbiAssemblySummaryGenbank.R
│   ├── updateBiosqlNcbiTaxonomy.R
│   ├── updateUniprotKeyword.R
│   ├── updateGeneOntology.R
│   ├── initializeTaxogeno.R        # One-time database setup
│   └── webApp/
│       ├── webAppFunctions.R
│       ├── webAppServer.R
│       └── webAppUi.R
├── inst/
│   ├── sql/
│   │   ├── taxogeno_full.sql       # Complete database schema
│   │   └── biosql_pgsql.sql        # BioSQL taxonomy schema
│   └── webApp/
│       └── webApp.R                # Shiny application entry point
├── DESCRIPTION
├── NAMESPACE
├── README.md
└── TODO.org
```

---

## Data Processing: ETL Pipeline

### 1. Extraction: Reading Unstructured Data

#### `readSma3sAnnotation()`

Converts the Sma3s-generated TSV into a structured data frame. The TSV contains tab-separated fields including:

```r
# Columns read:
fastaheader | genename | genedescription | enzyme | goc | gocname | gof | gofname | gop | gopname | keyword | pathway | goslim
```

The function extracts the short FASTA header (`fastashortheader`) for use as a unique identifier within the proteome.

#### `readMultifasta()`

Processes large FASTA files **line by line** to avoid loading the entire contents into RAM:

```r
# Incremental reading with file() and readLines()
fileConn <- file(multifastaPath, "r")
while (TRUE) {
    oneLine <- readLines(fileConn, n = 1)
    # Process headers (>id) and sequences
    # ...
}
close(fileConn)
```

**Memory Management Strategy:** Each proteome may contain thousands of sequences. Line-by-line processing allows handling of multi-GB files without exceeding system memory limits.

---

### 2. Transformation: From Unstructured to Structured

#### `extractColumnDfFromInsertedGeneAnnotationDf()`

Converts semicolon-separated fields into 1:N relationships:

```r
# Input:  geneid | enzyme
#         1      | EC1;EC2;EC3
# Output: geneid | ec
#         1      | EC1
#         1      | EC2
#         1      | EC3
colVecList <- strsplit(insertedGeneAnnotationDf[,colName], ";")
geneidVec <- rep(insertedGeneAnnotationDf[,"geneid"], lengths(colVecList))
colVec <- unlist(colVecList)
outputDf <- data.frame(geneidVec, colVec)
```

This pattern is applied to:
- GO terms (goc, gof, gop) → `gene_goc_rel`, `gene_gof_rel`, `gene_gop_rel`
- GO Slim → `gene_goslim_rel`
- Keywords → `gene_keyword_rel`
- EC enzymes → `gene_enzyme_rel`
- Pathways → `gene_pathway_rel`

**GO Slim Normalization with OWLTools:**

```r
# Generate GAF (Gene Association Format)
gafDf <- generateGafDf(dbConn, goStrictDf)

# Map2slim: convert specific GO terms to generic GO Slim terms
mappedGafDf <- owltoolsMap2SlimGafDf(gafDf, owl="goslim_generic.owl", subset="goslim_generic")

# Aggregate counts per proteome
summaryDf <- aggregate(dbobjectid ~ annotkwid, data = gafDf, length)
```

---

### 3. Loading: Database Insertion

#### `saveSma3sFileSet()`

The central function orchestrating the entire ETL pipeline:

1. **Insert proteome** into `taxogeno.proteome`:
   - Generates a `jobid` UUID for identification (used in the web application)
   - Saves MD5 hashes for integrity verification

2. **Insert genes** into `taxogeno.gene`:
   - Each gene receives an auto-generated `geneid`
   - Amino acid sequences are linked from the multifasta

3. **Insert 1:N relationships**:
   - GO terms, Keywords, EC, Pathways
   - Each relationship is inserted into its corresponding table

4. **Generate GO Slim summaries**:
   - Applies `owltools` to map GO to GO Slim
   - Stores results in `generated_goslim_summary`

5. **Calculate Euclidean distances**:
   - Compares the new proteome with all existing proteomes
   - Results are stored in cache tables

---

## Euclidean Distance Calculation

### Processing with Sparse Matrix and `dist()`

```r
calculateEuclideanDistances <- function(normalizedAnnotationDf) {
    # Build sparse matrix: proteomeid × annotkwid
    sparseMatrix <- xtabs(annotkwrelval ~ proteomeid + annotkwid, 
                          normalizedAnnotationDf, sparse = TRUE)
    
    # Calculate Euclidean distance between rows (proteomes)
    distObj <- dist(sparseMatrix, method = "euclidean", diag = FALSE, upper = FALSE)
    
    # Convert dist object to data frame
    distAsDf(distObj)
}
```

**Annotation Normalization:**

```sql
-- GO Slim: relative frequency per proteome
annotkwcount::float / genecount AS annotkwrelval

-- Keywords: relative frequency per proteome
COUNT(tgk.keyword) / proteome.genecount AS annotkwrelval
```

### Parallelization with `parallel`

The package uses R's `parallel` package to accelerate distance calculations:

```r
library(parallel)

# The system is prepared for parallel processing in:
# 1. Distance calculation in chunks
# 2. Processing multiple proteomes simultaneously
```

The code structure allows parallelization of calculation loops using `mclapply()` or `parLapply()`:

```r
# Each chunk of proteomes is processed independently
for(chunkVec in chunkVectorInEqualSizeFragments(otherProteomeIdVecFiltered, n = 1)) {
    # Calculate distances for a subset
    # Results accumulate in the database
}
```

Parallelization is not active by default in the sequential load script, but the code structure enables it through `mclapply()` from the `parallel` package to distribute distance calculations across multiple CPU cores.

### Memory Management: Chunked Processing

```r
chunkVectorInEqualSizeFragments <- function(x, n) {
    split(x, ceiling(seq_along(x)/n))
}

# Process in chunks of 20 proteomes (default)
for(chunkVec in chunkVectorInEqualSizeFragments(otherProteomeIdVecFiltered, n = 20)) {
    # Calculate for n=20 proteomes at a time
}
```

**Technical Rationale:** The sparse matrix generated by `xtabs(sparse = TRUE)` can exceed several gigabytes when comparing all proteomes simultaneously. Chunked processing limits memory usage to a manageable subset.

---

## Data Model

### Schema `taxogeno`

```
┌─────────────────────┐     ┌─────────────────────┐     ┌─────────────────────┐
│     taxonomy        │────▶│     proteome        │────▶│       gene          │
│─────────────────────│     │─────────────────────│     │─────────────────────│
│ ncbitaxid (PK)      │     │ proteomeid (PK)     │     │ geneid (PK)         │
│ parent_ncbitaxid    │     │ basename            │     │ proteomeid (FK)     │
│ node_rank           │     │ is_userproteome     │     │ fastaheader         │
│ scientific_name     │     │ gcaid               │     │ fastashortheader    │
│ is_ancestor         │     │ ncbitaxid (FK)      │     │ genename            │
└─────────────────────┘     │ genecount           │     │ genedescription     │
                            │ sma3sannotationmd5sum│     │ aasequence          │
                            │ multifastamd5sum    │     └─────────────────────┘
                            │ creationtimestamp   │              │
                            │ jobid (UUID)        │              ▼
                            └─────────────────────┘    ┌──────────────────────┐
                                     │                 │ gene_gop_rel         │
                                     │                 │ gene_gof_rel         │
                                     │                 │ gene_goc_rel         │
                                     │                 │ gene_goslim_rel      │
                                     │                 │ gene_keyword_rel     │
                                     │                 │ gene_enzyme_rel      │
                                     │                 │ gene_pathway_rel     │
                                     │                 └──────────────────────┘
                                     │
                                     ▼
                          ┌─────────────────────────────────────┐
                          │ generated_goslim_summary            │
                          │ euclidean_distances_generated_goslim│
                          │ euclidean_distances_keyword         │
                          └─────────────────────────────────────┘
```

### Distance Cache Tables

```sql
-- Composite primary key ensures uniqueness (proteomeid_greatest > proteomeid_least)
CREATE TABLE taxogeno.euclidean_distances_generated_goslim (
    proteomeid_greatest integer REFERENCES taxogeno.proteome(proteomeid),
    proteomeid_least integer REFERENCES taxogeno.proteome(proteomeid),
    distance numeric NOT NULL,
    CHECK (proteomeid_greatest > proteomeid_least),
    PRIMARY KEY (proteomeid_greatest, proteomeid_least)
);
```

**Cache Table Advantage:** Distances are calculated only once. The table is queried directly by the Shiny application without needing real-time recalculation.

---

## Recursive Queries (CTEs)

### `taxogeno.taxon_ancestors()`

Gets all ancestors of a given taxon using a recursive CTE:

```sql
WITH RECURSIVE rec_a (ncbitaxid, parent_ncbitaxid, node_rank, scientific_name, is_ancestor) AS (
    -- Base case: the target taxon
    SELECT ncbitaxid, parent_ncbitaxid, node_rank, scientific_name, is_ancestor
    FROM taxogeno.taxonomy
    WHERE ncbitaxid = $1
    
    UNION ALL
    
    -- Recursive case: climb to parent
    SELECT tax.ncbitaxid, tax.parent_ncbitaxid, tax.node_rank, tax.scientific_name, tax.is_ancestor
    FROM rec_a, taxogeno.taxonomy tax
    WHERE tax.ncbitaxid = rec_a.parent_ncbitaxid
      AND rec_a.ncbitaxid <> rec_a.parent_ncbitaxid
)
SELECT * FROM rec_a;
```

### `taxogeno.taxon_descendants()`

Gets all descendants of a taxon using a descending recursive CTE.

### `taxogeno.taxonomy_jsonb_children()`

Generates a JSON representation of the taxonomic tree for the Shiny application:

```sql
SELECT jsonb_build_object(
    'ncbitaxid', tax1.ncbitaxid,
    'parent_ncbitaxid', tax1.parent_ncbitaxid,
    'node_rank', tax1.node_rank,
    'scientific_name', tax1.scientific_name,
    'children', (SELECT jsonb_agg(taxogeno.taxonomy_jsonb_children(ncbitaxid))
                 FROM taxogeno.taxonomy tax2
                 WHERE tax2.parent_ncbitaxid = tax1.ncbitaxid
                   AND tax2.parent_ncbitaxid <> tax2.ncbitaxid)
)
INTO v_jsonbobj
FROM taxogeno.taxonomy tax1
WHERE ncbitaxid = p_ncbitaxid;
```

---

## Stored Procedures

### `taxogeno.insert_taxon()`

Inserts a taxon and all its ancestors using `biosql.taxon_ancestors()`:

```r
FOR v_pretaxonomy_rec IN (SELECT * FROM biosql.taxon_ancestors(p_ncbitaxid))
LOOP
    -- Determine if ancestor or target taxon
    IF v_pretaxonomy_rec.ncbi_taxon_id <> p_ncbitaxid THEN
        v_is_ancestor := TRUE;
    ELSE
        v_is_ancestor := FALSE;
    END IF;
    
    -- Insert or update with ON CONFLICT
    INSERT INTO taxogeno.taxonomy(ncbitaxid, parent_ncbitaxid, node_rank, scientific_name, is_ancestor)
    VALUES (v_ncbitaxid, v_parent_ncbitaxid, v_node_rank, v_scientific_name, v_is_ancestor)
    ON CONFLICT (ncbitaxid) DO UPDATE SET ...;
END LOOP;
```

### `taxogeno.delete_taxon()`

Deletes a taxon only if not referenced by any proteome. Uses `taxogeno.taxon_ancestors()` to identify the complete set of taxa to delete.

### `taxogeno.delete_proteome()`

Cascading deletion of a complete proteome:

```sql
DELETE FROM taxogeno.generated_goslim_summary WHERE proteomeid = p_proteomeid;
DELETE FROM taxogeno.euclidean_distances_keyword WHERE proteomeid_greatest = p_proteomeid OR proteomeid_least = p_proteomeid;
DELETE FROM taxogeno.euclidean_distances_generated_goslim WHERE proteomeid_greatest = p_proteomeid OR proteomeid_least = p_proteomeid;

FOR v_geneid IN (SELECT geneid FROM taxogeno.gene WHERE proteomeid = p_proteomeid) LOOP
    DELETE FROM taxogeno.gene_goc_rel WHERE geneid = v_geneid;
    DELETE FROM taxogeno.gene_gof_rel WHERE geneid = v_geneid;
    -- ... other relationship tables
    DELETE FROM taxogeno.gene WHERE geneid = v_geneid;
END LOOP;

DELETE FROM taxogeno.proteome WHERE proteomeid = p_proteomeid;
```

---

## External Data Integration

### Reference Database Updates

```r
updateNcbiAssemblySummaryGenbank(dbConn)  # Assembly_summary_genbank.txt
updateBiosqlNcbiTaxonomy(dbConn)          # NCBI Taxdump (via BioSQL)
updateUniprotKeyword(dbConn)              # UniProt Keywords
updateGeneOntology(dbConn)                # GO and GO Slim (OWL)
```

### `updateGeneOntology()`: OWL Parsing with XMLTABLE

Downloads Gene Ontology OWL files and parses them directly in PostgreSQL:

```sql
INSERT INTO gene_ontology.gene_ontology (goid, golabel, goaspect)
SELECT
    gene_ontology_tablified.goid,
    gene_ontology_tablified.golabel,
    gene_ontology_tablified.goaspect
FROM gene_ontology.gene_ontology_xml,
XMLTABLE(
    XMLNAMESPACES(
        'http://www.w3.org/1999/02/22-rdf-syntax-ns#' AS "rdf",
        'http://www.w3.org/2002/07/owl#' AS "owl",
        'http://www.geneontology.org/formats/oboInOwl#' AS "oboInOwl",
        'http://www.w3.org/2000/01/rdf-schema#' AS "rdfs"
    ),
    '/rdf:RDF/owl:Class'
    PASSING gene_ontology.gene_ontology_xml.xmldata
    COLUMNS
        goid     TEXT PATH 'oboInOwl:id/text()',
        golabel  TEXT PATH 'rdfs:label/text()',
        goaspect TEXT PATH 'oboInOwl:hasOBONamespace/text()'
) gene_ontology_tablified
WHERE gene_ontology.gene_ontology_xml.filename = $1
  AND gene_ontology_tablified.goid IS NOT NULL;
```

---

## Shiny Web Application

The interactive web application (`inst/webApp/webApp.R`) presents the following modules:

### Step 0: Upload Proteome
- Upload TSV and FASTA files
- Assign tags
- Generate `jobid` (UUID) for subsequent access

### Step 1: Select a Proteome
- List of available proteomes (public + user)
- Display general information: `genecount`, `ncbitaxid`, `scientific_name`
- Explore genes with their GO, Keywords, EC, Pathway annotations
- View amino acid sequences (formatted with line breaks every 60 characters)

### Step 2: Select Reference Proteomes
- Interactive taxonomic tree (`shinyTree`)
- Checkbox selection in the NCBI Taxonomy hierarchy

### Step 3: Check Annotation Info
- Interactive boxplots (`ggiraph`) of GO Slim by aspect:
  - Molecular Function
  - Biological Process
  - Cellular Component
- Tables with quartiles and interquartile range

### Step 4: Check Similarity Info
- Euclidean distance table with CSV download
- Dendrogram/heatmap of distances between selected proteomes

---

## Installation and Usage

### Installation

```r
# Install from GitHub
install.packages("devtools")
devtools::install_github("jmengom/taxogeno")

# Load the package
library(taxogeno)
```

### Database Setup

```r
# Connect to PostgreSQL
conn <- dbConnect(RPostgres::Postgres(),
                  host = "localhost",
                  dbname = "taxogeno",
                  user = "postgres",
                  password = "password")

# Initialize the database with reference data (one-time setup)
# This takes 10-30 minutes and requires internet access
initializeTaxogeno(conn)
```

### Load a Proteome

```r
# Load a single proteome from Sma3s output
jobid <- saveSma3sFileSet(
    dbConn = conn,
    sma3sAnnotationFilePath = "path/to/proteome_uniref90_go_goslim.tsv",
    multifastaFilePath = "path/to/proteome.faa",
    tagVec = c("bacteria", "example"),
    isUserProteome = FALSE,
    doUpdateNcbiAssemblyInfo = TRUE
)

print(paste("Upload job ID:", jobid))
```

---

## Database Schemas

| Schema | Contents |
|:-------|:---------|
| `taxogeno` | Main data: proteomes, genes, annotations, distances |
| `biosql` | NCBI taxonomy (BioSQL schema) with custom functions for recursive navigation |
| `ncbi` | Assembly Summary GenBank |
| `uniprot` | Keyword catalogs and UniProt proteomes |
| `gene_ontology` | GO and GO Slim (parsed OWL) |

---

## Technical Requirements

- **PostgreSQL 12+** with extensions: pl/R, XML
- **R 4.0+** with packages: DBI, RPostgres, shiny, shinyTree, ggplot2, DT, jsonlite, uuid, parallel
- **Perl** (for taxonomy loading script)
- **OWLTools** (for GO Slim mapping)
- **BioSQL** (taxonomy schema)
- **wget** (for downloading external databases)

---

## Authorship

- **R and PostgreSQL Code Development**: Javier Méndez Gómez
- **Project Supervision and Concept**: Dr. Antonio J. Pérez Pulido
- **Department of Bioinformatics**, Universidad Pablo de Olavide, Seville

The Taxogeno project was developed in 2020 as part of the training period (FCT) of the Higher Degree in Web Application Development (Royal Decree 686/2010, of May 20) at I.E.S. Hermanos Machado (Dos Hermanas, Seville). The internship was conducted at the Department of Bioinformatics, Universidad Pablo de Olavide, under the supervision of Dr. Antonio J. Pérez Pulido.

---

## References

- Sma3s: [UPOBioinfo/sma3s](https://github.com/UPOBioinfo/sma3s)
- Gene Ontology: [geneontology.org](http://geneontology.org)
- NCBI Taxonomy: [ncbi.nlm.nih.gov/taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)
- BioSQL: [biosql.org](http://biosql.org)
- OWLTools: [github.com/owlcollab/owltools](https://github.com/owlcollab/owltools)

---