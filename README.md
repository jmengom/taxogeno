# Taxogeno

**Sistema de análisis filogenético comparativo basado en anotaciones funcionales de proteomas completos**

---

## Descripción

Taxogeno es un sistema ETL (Extract-Transform-Load) para el procesamiento, almacenamiento y análisis comparativo de proteomas completos anotados funcionalmente. El sistema procesa la salida del programa [Sma3s](https://github.com/UPOBioinfo/sma3s) —desarrollado por el Dr. Antonio J. Pérez Pulido del Departamento de Bioinformática de la UPO— y construye una base de datos normalizada con:

- Anotaciones funcionales (GO, GO Slim, Keywords, EC, Pathways)
- Secuencias de proteínas
- Relaciones taxonómicas
- Matrices de distancias euclídeas entre proteomas

El sistema se implementó durante las prácticas del **Ciclo Formativo de Grado Superior en Desarrollo de Aplicaciones Web** en el Departamento de Bioinformática de la Universidad Pablo de Olavide.

---

## Arquitectura del Sistema

```
┌─────────────────────────────────────────────────────────────────────────┐
│                          TAXOGENO PIPELINE                             │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  ┌───────────────┐     ┌───────────────┐     ┌────────────────────┐    │
│  │  Sma3s TSV    │────▶│  readSma3s    │────▶│  PostgreSQL        │    │
│  │  (anotación)  │     │  Annotation   │     │  (taxogeno db)     │    │
│  └───────────────┘     └───────────────┘     └────────────────────┘    │
│                                                                         │
│  ┌───────────────┐     ┌───────────────┐     ┌────────────────────┐    │
│  │  FASTA        │────▶│  readMultif-  │────▶│  Tablas:           │    │
│  │  (proteoma)   │     │  asta         │     │  - proteome        │    │
│  └───────────────┘     └───────────────┘     │  - gene            │    │
│                                              │  - gene_goc_rel    │    │
│  ┌───────────────┐     ┌───────────────┐     │  - gene_gof_rel    │    │
│  │  OWLTools     │────▶│  map2slim     │────▶│  - gene_gop_rel    │    │
│  │  (GO Slim)    │     │               │     │  - gene_keyword_rel│    │
│  └───────────────┘     └───────────────┘     │  - gene_enzyme_rel │    │
│                                              │  - gene_pathway_rel│    │
│  ┌───────────────┐     ┌───────────────┐     │  - generated_goslim│    │
│  │  Cálculo de   │────▶│  Matriz de    │────▶│    _summary        │    │
│  │  distancias   │     │  distancias   │     │  - euclidean_      │    │
│  │  euclídeas    │     │  (caché)      │     │    distances_*     │    │
│  └───────────────┘     └───────────────┘     └────────────────────┘    │
│                                              │                         │
│                                              ▼                         │
│                              ┌───────────────────────────────┐         │
│                              │  Shiny Application (R/Shiny)  │         │
│                              │  - Consulta interactiva       │         │
│                              │  - Visualización de árboles   │         │
│                              │  - Comparación de proteomas   │         │
│                              └───────────────────────────────┘         │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Tecnologías Utilizadas

| Tecnología | Componente | Uso |
|:-----------|:-----------|:----|
| **R** | `taxogeno_library.R`, `taxogeno_shinyApp.R` | Lenguaje principal: procesamiento de datos, análisis estadístico, aplicación web |
| **PostgreSQL 12+** | `sql/taxogeno_full.sql` | Motor de base de datos con soporte JSONB, XML y funciones recursivas |
| **pl/R** | Procedimientos almacenados | Integración R-PostgreSQL para funciones estadísticas dentro de la BD |
| **Shiny** | `taxogeno_shinyApp.R` | Aplicación web interactiva para exploración y visualización |
| **OWLTools** | `bin/owltools` | Mapeo de términos GO a GO Slim |
| **BioSQL** | `sql/biosql_pgsql.sql` | Esquema para almacenamiento de taxonomía NCBI |
| **parallel** | `library(parallel)` | Procesamiento paralelo en cálculos de distancia |

---

## Estructura del Repositorio

```
taxogeno/
├── taxogeno_library.R          # Biblioteca R: todas las funciones ETL
├── taxogeno_loadProteomes.R    # Script de carga masiva (secuencial)
├── taxogeno_init.R             # Inicialización de bases de datos externas
├── taxogeno_shinyApp.R         # Aplicación Shiny (interfaz web)
├── install.sh                  # Script de instalación automatizada
├── THIRD_PARTY.md              # Licencias de terceros
├── bin/
│   ├── owltools                # Herramienta de mapeo GO
│   └── load_ncbi_taxonomy.pl   # Script Perl para carga de taxonomía NCBI
└── sql/
    ├── taxogeno_full.sql       # Esquema completo: tablas, funciones, procedimientos
    └── biosql_pgsql.sql        # Esquema BioSQL (taxonomía)
```

---

## Procesamiento de Datos: ETL

### 1. Extracción: Lectura de Datos No Estructurados

#### `readSma3sAnnotation(sma3sAnnotationPath)`

Convierte el TSV generado por Sma3s en un dataframe estructurado. El TSV contiene campos separados por tabulador que incluyen:

```r
# Columnas leídas:
fastaheader | genename | genedescription | enzyme | goc | gocname | gof | gofname | gop | gopname | keyword | pathway | goslim
```

La función extrae la cabecera corta del FASTA (`fastashortheader`) para su uso como identificador único dentro del proteoma.

#### `readMultifasta(multifastaPath)`

Procesa archivos FASTA de gran tamaño **línea a línea** para evitar cargar el contenido completo en RAM:

```r
# Lectura incremental con file() y readLines()
fileConn = file(multifastaPath, "r")
while (TRUE) {
    oneLine = readLines(fileConn, n = 1)
    # Procesamiento de cabeceras (>id) y secuencias
    # ...
}
close(fileConn)
```

**Estrategia de control de memoria:** Cada proteoma puede contener miles de secuencias. El procesamiento línea a línea permite manejar archivos de varios GB sin superar los límites de memoria del sistema.

---

### 2. Transformación: De No Estructurado a Estructurado

#### `extractColumnDfFromInsertedGeneAnnotationDf()`

Convierte campos con valores múltiples separados por `;` en relaciones 1:N:

```r
# Entrada: geneid | enzyme
#         1      | EC1;EC2;EC3
# Salida:  geneid | ec
#          1      | EC1
#          1      | EC2
#          1      | EC3
colVecList <- strsplit(insertedGeneAnnotationDf[,colName], ";")
geneidVec <- rep(insertedGeneAnnotationDf[,"geneid"], lengths(colVecList))
colVec <- unlist(colVecList)
outputDf <- data.frame(geneidVec, colVec)
```

Este patrón se aplica a:
- GO terms (goc, gof, gop) → `gene_goc_rel`, `gene_gof_rel`, `gene_gop_rel`
- GO Slim → `gene_goslim_rel`
- Keywords → `gene_keyword_rel`
- Enzimas EC → `gene_enzyme_rel`
- Pathways → `gene_pathway_rel`

**Normalización de GO Slim mediante OWLTools:**

```r
# Genera GAF (Gene Association Format)
gafDf <- generateGafDf(dbConn, goStrictDf)

# Map2slim: convierte GO terms específicos a GO Slim genéricos
mappedGafDf <- owltoolsMap2SlimGafDf(gafDf, owl="goslim_generic.owl", subset="goslim_generic")

# Agrega conteos por proteoma
summaryDf <- aggregate(dbobjectid~annotkwid, data=gafDf, length)
```

---

### 3. Carga: Inserción en Base de Datos

#### `saveSma3sFileSet()`

Función central que orquesta todo el pipeline de carga:

1. **Inserción del proteoma** en `taxogeno.proteome`:
   - Se genera un `jobid` UUID para identificación (uso en la aplicación web)
   - Se guardan hashes MD5 para verificación de integridad

2. **Inserción de genes** en `taxogeno.gene`:
   - Cada gen recibe un `geneid` serial autogenerado
   - Se vinculan las secuencias de aminoácidos desde el multifasta

3. **Inserción de relaciones 1:N**:
   - GO terms, Keywords, EC, Pathways
   - Cada relación se inserta en su tabla correspondiente

4. **Generación de resúmenes GO Slim**:
   - Se aplica `owltools` para mapear GO a GO Slim
   - Se almacena en `generated_goslim_summary`

5. **Cálculo de distancias euclídeas**:
   - Se compara el nuevo proteoma con todos los existentes
   - Los resultados se almacenan en tablas de caché

---

## Cálculo de Distancias Euclídeas

### Procesamiento con Sparse Matrix y `dist()`

```r
calculateEuclideanDistances <- function(normalizedAnnotationDf) {
    # Construye matriz dispersa: proteomeid × annotkwid
    sparseMatrix <- xtabs(annotkwrelval ~ proteomeid + annotkwid, 
                          normalizedAnnotationDf, sparse=TRUE)
    
    # Calcula distancia euclídea entre filas (proteomas)
    distObj <- dist(sparseMatrix, method="euclidean", diag=FALSE, upper=FALSE)
    
    # Convierte objeto dist a dataframe
    distAsDf(distObj)
}
```

**Normalización de anotaciones:**

```sql
-- GO Slim: frecuencia relativa por proteoma
annotkwcount::float / genecount AS annotkwrelval

-- Keywords: frecuencia relativa por proteoma
COUNT(tgk.keyword) / proteome.genecount AS annotkwrelval
```

### Paralelización con `parallel`

Taxogeno utiliza el paquete `parallel` de R para acelerar el cálculo de distancias:

```r
library(parallel)

# El sistema está preparado para procesamiento paralelo en:
# 1. Cálculo de distancias en fragmentos (chunks)
# 2. Procesamiento de múltiples proteomas simultáneamente
```

La estructura del código permite paralelizar los bucles de cálculo mediante `mclapply()` o `parLapply()`:

```r
# Cada fragmento (chunk) de proteomas se procesa independientemente
for(chunkVec in chunkVectorInEqualSizeFragments(otherProteomeIdVecFiltered, n=1)) {
    # Cálculo de distancias para un subconjunto
    # Los resultados se acumulan en la base de datos
}
```

La paralelización no está activa por defecto en el script de carga secuencial (`taxogeno_loadProteomes.R`), pero la estructura del código lo permite mediante la función `mclapply()` del paquete `parallel` para distribuir el cálculo de distancias entre múltiples núcleos de CPU.

### Gestión de Memoria: Procesamiento por Fragmentos

```r
chunkVectorInEqualSizeFragments <- function(x, n) {
    split(x, ceiling(seq_along(x)/n))
}

# Procesamiento en fragmentos de 20 proteomas (por defecto)
for(chunkVec in chunkVectorInEqualSizeFragments(otherProteomeIdVecFiltered, n=20)) {
    # Cálculo para n=20 proteomas a la vez
}
```

**Razón técnica:** La matriz dispersa generada por `xtabs(sparse=TRUE)` puede superar varios gigabytes cuando se comparan todos los proteomas simultáneamente. El procesamiento por fragmentos (chunking) limita el uso de memoria a un subconjunto manejable.

---

## Modelo de Datos

### Esquema `taxogeno`

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

### Tablas de Caché de Distancias

```sql
-- Índice primario compuesto garantiza unicidad (proteomeid_greatest > proteomeid_least)
CREATE TABLE taxogeno.euclidean_distances_generated_goslim (
    proteomeid_greatest integer REFERENCES taxogeno.proteome(proteomeid),
    proteomeid_least integer REFERENCES taxogeno.proteome(proteomeid),
    distance numeric NOT NULL,
    CHECK (proteomeid_greatest > proteomeid_least),
    PRIMARY KEY (proteomeid_greatest, proteomeid_least)
);
```

**Ventaja de la tabla de caché:** Las distancias solo se calculan una vez. La tabla es consultada directamente por la aplicación Shiny sin necesidad de recalcular en tiempo real.

---

## Funciones con Queries Recursivas (CTE)

### `taxogeno.taxon_ancestors()`

Obtiene todos los ancestros de un taxón dado mediante CTE recursiva:

```sql
WITH RECURSIVE rec_a (ncbitaxid, parent_ncbitaxid, node_rank, scientific_name, is_ancestor) AS (
    -- Caso base: el taxón objetivo
    SELECT ncbitaxid, parent_ncbitaxid, node_rank, scientific_name, is_ancestor
    FROM taxogeno.taxonomy
    WHERE ncbitaxid = $1
    
    UNION ALL
    
    -- Caso recursivo: subir al padre
    SELECT tax.ncbitaxid, tax.parent_ncbitaxid, tax.node_rank, tax.scientific_name, tax.is_ancestor
    FROM rec_a, taxogeno.taxonomy tax
    WHERE tax.ncbitaxid = rec_a.parent_ncbitaxid
      AND rec_a.ncbitaxid <> rec_a.parent_ncbitaxid
)
SELECT * FROM rec_a;
```

### `taxogeno.taxon_descendants()`

Obtiene todos los descendientes de un taxón mediante CTE recursiva descendente.

### `taxogeno.taxonomy_jsonb_children()`

Genera una representación JSON del árbol taxonómico para la aplicación Shiny:

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

## Procedimientos Almacenados

### `taxogeno.insert_taxon()`

Inserta un taxón y todos sus ancestros, utilizando `biosql.taxon_ancestors()`:

```r
FOR v_pretaxonomy_rec IN (SELECT * FROM biosql.taxon_ancestors(p_ncbitaxid))
LOOP
    -- Determina si es ancestro o el taxón objetivo
    IF v_pretaxonomy_rec.ncbi_taxon_id <> p_ncbitaxid THEN
        v_is_ancestor := TRUE;
    ELSE
        v_is_ancestor := FALSE;
    END IF;
    
    -- Inserta o actualiza con ON CONFLICT
    INSERT INTO taxogeno.taxonomy(ncbitaxid, parent_ncbitaxid, node_rank, scientific_name, is_ancestor)
    VALUES (v_ncbitaxid, v_parent_ncbitaxid, v_node_rank, v_scientific_name, v_is_ancestor)
    ON CONFLICT (ncbitaxid) DO UPDATE SET ...;
END LOOP;
```

### `taxogeno.delete_taxon()`

Elimina un taxón solo si no está referenciado por ningún proteoma. Utiliza `taxogeno.taxon_ancestors()` para identificar el conjunto completo de taxones a eliminar.

### `taxogeno.delete_proteome()`

Elimina en cascada un proteoma completo:

```sql
DELETE FROM taxogeno.generated_goslim_summary WHERE proteomeid = p_proteomeid;
DELETE FROM taxogeno.euclidean_distances_keyword WHERE proteomeid_greatest = p_proteomeid OR proteomeid_least = p_proteomeid;
DELETE FROM taxogeno.euclidean_distances_generated_goslim WHERE proteomeid_greatest = p_proteomeid OR proteomeid_least = p_proteomeid;

FOR v_geneid IN (SELECT geneid FROM taxogeno.gene WHERE proteomeid = p_proteomeid) LOOP
    DELETE FROM taxogeno.gene_goc_rel WHERE geneid = v_geneid;
    DELETE FROM taxogeno.gene_gof_rel WHERE geneid = v_geneid;
    -- ... otras tablas de relación
    DELETE FROM taxogeno.gene WHERE geneid = v_geneid;
END LOOP;

DELETE FROM taxogeno.proteome WHERE proteomeid = p_proteomeid;
```

---

## Integración de Datos Externos

### Actualización de Bases de Datos de Referencia (`taxogeno_init.R`)

```r
updateNcbiAssemblySummaryGenbank(dbConn)  # Assembly_summary_genbank.txt
updateBiosqlNcbiTaxonomy(dbConn)          # Taxdump NCBI (vía BioSQL)
updateUniprotKeyword(dbConn)              # Keywords de UniProt
updateGeneOntology(dbConn)                # GO y GO Slim (OWL)
```

### `updateGeneOntology()`: Parseo de OWL con XMLTABLE

Descarga archivos OWL de Gene Ontology y los parsea directamente en PostgreSQL:

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

## Aplicación Shiny

La aplicación web interactiva (`taxogeno_shinyApp.R`) presenta los siguientes módulos:

### Step 0: Upload Proteome
- Subida de archivos TSV y FASTA
- Asignación de tags
- Generación de `jobid` (UUID) para acceso posterior

### Step 1: Select a Proteome
- Lista de proteomas disponibles (públicos + del usuario)
- Visualización de información general: `genecount`, `ncbitaxid`, `scientific_name`
- Exploración de genes con sus anotaciones GO, Keywords, EC, Pathways
- Visualización de secuencia de aminoácidos (formato con saltos de línea cada 60 caracteres)

### Step 2: Select Reference Proteomes
- Árbol taxonómico interactivo (`shinyTree`)
- Selección mediante checkbox en la jerarquía de NCBI Taxonomy

### Step 3: Check Annotation Info
- Boxplots interactivos (`ggiraph`) de GO Slim por aspecto:
  - Molecular Function
  - Biological Process
  - Cellular Component
- Tablas con cuartiles y rango intercuartílico

### Step 4: Check Similarity Info
- Tabla de distancias euclídeas con descarga CSV
- Dendrograma (heatmap) de distancias entre proteomas seleccionados

---

## Carga Masiva (`taxogeno_loadProteomes.R`)

Procesa directorios completos con conjuntos de proteomas:

```r
dirNameVec <- c("archaeas", "bacteria", "bacteria2", "mammals", "fungi", "invert", "protozoa", "vert_other")

for(dirName in dirNameVec) {
    # Auto-detecta archivos TSV y FASTA
    baseNameVec <- gsub("^(.*)_uniref90_go_goslim[.]tsv$", "\\1", 
                        list.files(path=dirNamePath, pattern=".*_uniref90_go_goslim[.]tsv"))
    
    # Procesa cada proteoma
    for(rowNum in seq_len(nrow(filesDf))) {
        start.time <- Sys.time()
        saveSma3sFileSet(dbConn, tsvFilePath, faaFilePath, isUserProteome=FALSE)
        end.time <- Sys.time()
        timeTakenVec <- c(timeTakenVec, end.time - start.time)
        print(mean(timeTakenVec))  # Muestra media acumulada
    }
}
```

---
### Esquemas en la Base de Datos

| Esquema | Contenido |
|:--------|:----------|
| `taxogeno` | Datos principales: proteomas, genes, anotaciones, distancias |
| `biosql` | Taxonomía NCBI (esquema BioSQL) + funciones propias para la navegación recursiva del árbol taxonómico.|
| `ncbi` | Assembly Summary GenBank |
| `uniprot` | Catalogs de Keywords y Proteomas UniProt |
| `gene_ontology` | GO y GO Slim (OWL parseado) |

---

## Requisitos Técnicos

- **PostgreSQL 12+** con extensiones: pl/R, XML
- **R 4.0+** con paquetes: DBI, RPostgres, shiny, shinyTree, ggplot2, DT, jsonlite, uuid, parallel
- **Perl** (para script de carga taxonómica)
- **OWLTools** 
- **BioSQL**
- **wget** (para descarga de bases de datos externas)

---

## Autoría

- **Desarrollo del código en R y PostgreSQL**: Javier Méndez Gómez
- **Supervisión y concepto del proyecto**: Dr. Antonio J. Pérez Pulido
- **Departamento de Bioinformática**, Universidad Pablo de Olavide, Sevilla

El proyecto Taxogeno fue desarrollado en 2020 como parte del módulo de Formación en Centros de Trabajo (FCT) del Ciclo Formativo de Grado Superior en Desarrollo de Aplicaciones Web (Real Decreto 686/2010, de 20 de mayo) en el I.E.S. Hermanos Machado (Dos Hermanas, Sevilla). Las prácticas se realizaron en el Departamento de Bioinformática de la Universidad Pablo de Olavide bajo la supervisión del Dr. Antonio J. Pérez Pulido.

---

## Referencias

- Sma3s: [UPOBioinfo/sma3s](https://github.com/UPOBioinfo/sma3s)
- Gene Ontology: [geneontology.org](http://geneontology.org)
- NCBI Taxonomy: [ncbi.nlm.nih.gov/taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)
- BioSQL: [biosql.org](http://biosql.org)
- OWLTools: [github.com/owlcollab/owltools](https://github.com/owlcollab/owltools)
