# Tutorial RNA-seq de expresión diferencial con DESeq2 #
# Proyecto: Ampliseq transcriptomic array of REF seq. genes in Snai1 knock out triple negative breast cancer cells 

# ------------------------------------------------------------------------------
# 1. Configuración inicial
# ------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) # Verificamos si el paquete está instalado sin cargarlo,
  install.packages("BiocManager") # Si no está, se instala.
BiocManager::install("DESeq2") # Instalamos el programa que vamos a utilizar para realizar el análisis diferencial de genes
library(DESeq2)

# ------------------------------------------------------------------------------
# 2. Descargamos los datos
# ------------------------------------------------------------------------------

# Descargamos el conteo crudo, que va a contener los datos de expresión génica, con columnas como identificadores de genes y muestras
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
head(counts)

# Descargamos también el metadato que incluye información experimental, como genotipos o condiciones de las muestras
metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/ExperimentDesignFile.RnaSeq/experiment-design")
head(metadata)


# ------------------------------------------------------------------------------
# 3. Preparación de datos para DESeq2
# ------------------------------------------------------------------------------

# Una característica de DESeq es que espera que los conteos tengan IDs de genes como nombres de fila
head(counts)
rownames(counts) = counts$Gene.ID
head(counts)

# Y quitamos las columnas no necesarias (gene ID y gene name)
genes = counts[, c("Gene.ID", "Gene.Name")]
counts = counts[, -c(1, 2)]
head(counts)

# Por otro lado, DESeq también espera que la matriz de metadatos tenga el ID de muestra en los rownames
head(metadata)
rownames(metadata) = metadata$Run
head(metadata)

# Y simplificamos las columnas que van a ser de interés para el análisis, en este caso las condiciones del genotipo
metadata = metadata[, c("Sample.Characteristic.genotype."), drop=FALSE]
metadata

# Y renombramos la columna para que sea más manejable
colnames(metadata) = c("genotype")
metadata

# Quitamos los espacios en el nombre de las condiciones para evitar cualquier warning/aviso posterior de DESeq2
metadata$genotype[metadata$genotype == 'wild type genotype'] = 'wildtype'
metadata$genotype[metadata$genotype == 'Snai1 knockout'] = 'knockout'
metadata

# Formateamos los valores de las condiciones para convertirlos en factor
metadata$genotype = factor(metadata$genotype, levels=c("wildtype", "knockout"))
metadata$genotype


# ------------------------------------------------------------------------------
# 4. Verificación de la expresión de SNAI1
# ------------------------------------------------------------------------------

gene_id = genes$Gene.ID[ genes$Gene.Name == 'SNAI1' ] # Buscamos el ID del gen SNAI1 en los datos de los genes
gene_counts = counts[gene_id, ] # Extraemos el conteo asociado a este gen
gene_counts

gene_data = cbind(metadata, counts=as.numeric(gene_counts)) # Combinamos la información del conteo del gen, con la información del metadato
gene_data

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
ggplot(gene_data, aes(x = genotype, y = counts, fill = genotype)) + geom_boxplot() # Visualizamos con ggplot las diferencias de expresión entre genotipos


# ------------------------------------------------------------------------------
# 5. Análisis diferencial con DESeq2
# ------------------------------------------------------------------------------

# Creamos un objeto DESeqDataSet con el conteo y los metadatos, especificando el genotipo como diseño experimental
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~genotype)

# Filtramos genes con bajo conteo para mejorar la calidad del análisis
dds <- dds[rowSums(counts(dds)) > 10, ]

# Ejecutamos DESeq2
dds <- DESeq(dds)

# Extraemos la información de los resultados específicos comparando las condiciones del diseño
res = results(dds, contrast=c("genotype", "knockout", "wildtype"), alpha=1e-5) # El valor alpha es el umbral de significancia para los valores ajustados de padj. Representa la tasa máxima aceptable de falsos positivos (Error tipo I). En este caso establecemos un umbral de 0.001% para considerar un gen como significativamente diferencial, útil para datos de alta calidad o alta dimensionalidad (decenas de miles de genes a analizar, por lo que un umbral más bajo reduce el riesgo de detectar genes "falsos positivos" debido a la cantidad de pruebas realizadas)
res

### Explicación de los resultados ###

# baseMean: Promedio de conteo normalizado para un gen en todas las muestras, sin considerar el diseño experimental. Sirve como una medida básica de expresión general
# log2FoldChange: Cambio relativo en la expresión entre las dos condiciones comparadas. Un valor positivo significa que hay mayor expresión en la condición de interés, mientras un valor negativo refleja menor expresión
# lfcSE: Error estándar del log2FoldChange. Muestra variabilidad esperada del cambio relativo estimado
# stat: Estadístico de prueba usado para evaluar la hipótesis nula (que no hay diferencia en la expresión). Es el cociente entre el log2FoldChange y su error estándar.
# pvalue: Representa la probabilidad de observar un efecto tan extremo como el estimado, suponiendo que la hipótesis nula sea verdadera
# padj: Valor p ajustado para múltiples comparaciones. Este valor es crucial para determinar qué genes son diferencialmetne expresados, ya que controla la tasa de falsos descubrimientos (FDR)

## Cómo se puede filtrar más los resultados

### a. Filtrar genes con baja expresión
### b. Ajustar el valor de alpha
### c. Cambiar el umbral de log2FoldChange


# Otros ejemplos de fórmula de diseño:
# - https://www.ebi.ac.uk/gxa/experiments/E-MTAB-8031/Results
#     * RNA-seq sobre tumores pulmonares
#     * Variables:
#         1. Location: células tumorales o tejido normal cercano
#         2. Sex
#     * Formula = ~sex + location
#
# - https://www.ebi.ac.uk/gxa/experiments/E-MTAB-10041/Results
#     * RNA-seq en células de colon de ratón
#     * Variables:
#         1. Compound: medicamento contra el cancer o no
#         2. Genotype: ratón wildtype y otros 2 genotipos de mutaciones TP53
#     * Formula = ~genotype + compound  ## Estaríamos evaluando los efectos de cada variable de forma independiente, no capturaríamos cómo la interacción entre genotype y compound afecta en la expresión génica.
#     * Formula = ~group 
#         metadata$group = factor(paste(metadata$genotype, metadata$compound)) # Haciendo esto combinamos todas las variables en una sola, es útil cuando queremos analizar todas las combinaciones posibles y compararlas entre sí

# ------------------------------------------------------------------------------
# Controles adicionales
# ------------------------------------------------------------------------------

# Podemos comparar estos controles para ver si todo ha ido bien con el proyecto en el Expression Atlas
# https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5244/Results?specific=true&geneQuery=%255B%255D&filterFactors=%257B%257D&cutoff=%257B%2522foldChange%2522%253A11%252C%2522pValue%2522%253A0.01%257D&regulation=%2522UP_DOWN%2522

# Fusionar el nombre del gen con el dataframe para poder comparar los datos usando los nombres de los genes
res_df = as.data.frame(res)
head(res_df)
head(genes)
res_df = merge(res_df, genes, by='row.names')
head(res_df)

genes_to_check = c("THY1", "SFMBT2", "PASD1", "SNAI1")
res_df[res_df$Gene.Name %in% genes_to_check, ]


# ------------------------------------------------------------------------------
# 6. Visualización de los resultados
# ------------------------------------------------------------------------------

# MA plot
plotMA(res) # Los valores azules son los genes que pasan el umbral que hemos especificado del alpha

# Volcano plot
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='pvalue')

# HeatMap

if(!require(pheatmap)){
  install.packages("pheatmap")
  library(pheatmap)
}

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[(rownames(res)), ]
mat_scaled <- t(scale(t(mat)))

pheatmap(mat_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(dds)),
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50))

# Extracción de lista top genes diferencialmente expresados

# Ordenar por significancia ajustada
res_ordered <- res_df[order(res_df$padj), ]
# Seleccionar los 10 genes más significativos
top10_genes <- head(res_ordered, 10)
top10_genes

mat_subset <- mat_scaled[(head(order(res$padj), 10)), ]

pheatmap(mat_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(dds)),
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50))

# ------------------------------------------------------------------------------
# 7. Análisis adicional con biomaRt
# ------------------------------------------------------------------------------

BiocManager::install("biomaRt")
library(biomaRt)

# Buscamos los nombres en el conjunto de datos en Ensembl
ensembl <- useEnsembl(biomart="genes")
datasets = listDatasets(ensembl)
head(datasets)

dataset_nb = grep("human", datasets$description, ignore.case=TRUE)
dataset_nb

dataset = datasets$dataset[dataset_nb]
dataset

ensembl <- useDataset(dataset=dataset, mart=ensembl)

# Obtenemos las coordenadas de todos los genes
attributes <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position")
values <- list(ensembl_gene_id=c())
all.genes <- getBM(attributes=attributes, values=values, mart=ensembl)
head(all.genes)

# Renombramos la columna para que coincida con la columna de nuestro dataframe inicial
head(res_df)
colnames(all.genes)[1] = "Gene.ID"
head(all.genes)

# Fusionamos el resultado de DESeq con las coordenadas genéticas de Ensembl
merged_data <- merge(all.genes, res_df, by="Gene.ID")
head(merged_data)
