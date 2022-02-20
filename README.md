# RNA-Seq_analysis
This is a collection of scripts developed for RNA-Seq analyses.

## Pellegrinelli et al. - PEPD RNA-Seq analyses

The following scripts inside this resporoty were used in Pellegrinelli et al. 2022:

reads2counts.py
get_mapping_stats.py
combine_count_tables.py
construct_gene_expression_plot.py
construct_DESeq2_input.py
DEG_analysis20181107.R
DEG_plot.py

Pellegrinelli V., Rodriguez-Cuenca S., Rouault C., Schilbert H., Virtue S., Moreno-Navarrete J. M., Bidault G., del Carmen Vazquez Borrego M., Dias A. R., Pucker B., Dale M., Campbell M., Carobbio S., Aron-Wisnewsky J., Mora S., Masiero M., Emmanouilidou A., Mukhopadhyay S., Dougan G., den Hoed M., Loos R., Fernandez-Real J. M., Chiarugi D., Clement K. and Vidal-Puig T. Dysregulation of macrophage PEPD in obesity determines adipose tissue fibro-inflammation and insulin resistance. doi:[10.21203/rs.3.rs-57182.](https://doi.org/10.21203/rs.3.rs-57182/v1)


# construct_DESeq2_input.py
This script construct the input tables for DESeq2 analysis via DEG_analysis20181107.R

Requirements:
1) Python 2.7.x (other Python 2 versions should work as well)

```
Usage:
python construct_DESeq2_input.py
--counts <FULL_PATH_TO_OUTPUT_COUNT_TABLE>
--info <FULL_PATH_TO_SAMPLE_INFO_FILE>
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
```

# DEG_analysis20181107.R
This script performs differential gene expression analysis for RNA-seq samples.

```
Usage:
run DEG_analysis20181107.R in R
A clean data matrix and sample table (use construct_DESeq2_input.py for generation) are needed and should be provided here:
 
# --- loading sampleTable (generated by Python script) --- #
csvfile <- "clean_sample_table.txt"

# --- loading the data matrix --- #
count_data_file <- "clean_data_matrix.txt"
```

# combine_count_tables.py
This script combines all count tables given in the provided folders. It is assumed that these count tables were generated by featureCounts on the gene level. Multiple folder names can be provided separated by comma. The gene names in all count tables need to match. Raw counts and TPMs (=tags per million assigned tags) can be calculated independent of any annotation. The calculation of FPKMs is based on information about total exon length of a gene. GFF3 annotation is used for this step. The current versions assumes that the GFF3 file was downloaded from the NCBI and contains gene, transcript, and exon features. Overlapping exon features of different transcripts are merged to get the final exon length within a gene.

Requirements:
1) Python 2.7.x (other Python 2 versions should work as well)

```
Usage:
python combine_count_tables.py \
--in <FULL_PATH_INPUT_FILE(S)> (multiple paths can be provided comma-separated) \
--gff <FULL_PATH_TO_GFF3_FILE> \
--out <FULL_PATH_TO_OUTPUT_FILE>
```

Suggested citation:

this repository


# find_housekeeping_genes.py
This script identifies genes with very small variation in expression across multiple samples. Suche genes could be used as reference genes for qRT-PCR experiments. The gene expression file used as input should be in the output format of combine_count_tables.py: header line with different samples names, one row per gene starting with the gene name followed by expression values of the different samples. An annotation file can be provided to add a functional description to each gene in the output file.

Note: qPCR_gene_finder.py is an alterantive script for this function.

Requirements:
1) Python 2.7.x (other Python 2 versions should work as well)

```
Usage:
python find_housekeeping_genes.py \
--in <FULL_PATH_TO_EXPRESSION_FILE> \
--out <FULL_PATH__TO_OUTPUT_FILE>

optional:
--anno <FULL_PATH__TO_ANNOTATION_FILE> \
--cutoff <MINIMAL_EXPRESSION_PER_SAMPLE(INTEGER)>
```

Suggested citation:

this repository


# construct_gene_expression_plot.py
This script produces figures for the expression of selected genes across multiple samples. Plots can be generated for raw counts, TPMs, and FPKMs. Resulting plots can be used to analyze the expression of reference genes after normalization.

Requirements:
1) Python 2.7.x (other Python 2 versions should work as well) including matplotlib and NumPy

```
Usage:
python construct_gene_expression_plots.py \
--candidates <FULL_PATH_TO_CANDIDATE_FILE> \
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
	
at least one expression data file is required:

--counts <FULL_PATH_TO_RAW_COUNT_TABLE> \
--tpms <FULL_PATH_TO_TPM_FILE> \
--fpkms <FULL_PATH_TO_FPKM_FILE>
		
optional:

--samples <FULL_PATH_TO_SAMPLE_ORDER_FILE>
```
Suggested citation:

this repository


# DEG_plot.py
This script produces figures for differentially expressed genes sorted by the adjusted p-value to illustrate the log2FC. DESeq2 output can be used as input for this script. It is possible to add an additional column to customize gene names. Otherwise, gene IDs will be used as labels.

Requirements:
1) Python 2.7.x (other Python 2 versions should work as well) including matplotlib

```
Usage:
python DEG_plot.py \
--in <FULL_PATH_TO_INPUT_DIRECTORY> \
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
```
Suggested citation:

this repository


# get_mapping_stats.py
All reports of a STAR mapping in one folder are processed. Relevant read numbers are collected and summarized in a single table.

Requirements:
1) python 2.7

```
Usage:
python get_mapping_stats.py \
--in <FULL_PATH_TO_INPUT_DIRECTORY> \
--out <FULL_PATH_TO_OUTPUT_FILE>


Suggested citation:

this repository


# reads2counts.py
All FASTQ files in a given folder are subjected to STAR for read mapping against a given reference sequence. Based on a provided GFF3 file the expression of genes is quantified via featureCounts. 

Requirements:
1) python 2.7
2) STAR (Dobin, 2013)
3) featureCounts (Liao, 2014)

```
Usage:
python reads2counts2.py \
--fastq_file_dir <FULL_PATH_TO_DIRECTORY> \
--tmp_cluster_dir <FULL_PATH_TO_TEMPORARY_DIRECTORY_ON_CLUSTER_VOLUME> \
--result_dir <FULL_PATH_TO_RESULT_DIRECTORY> \
--ref_gff_file <FULL_PATH_TO_GFF_FILE_MATCHING_THE_PROVIDED_GENOME> \
--ref_genome_file <FULL_PATH_TO_GENOME_FILE_MATCHING_THE_PROVIDED_GFF_FILE> \
				
optional:
--dissimilarity <FLOAT, 1-identity, value between 0.0 and 1.0>[0.05] \
--length_fraction <FLOAT, value between 0.0 and 1.0>[0.9] \
--para_jobs <INTEGER, number of jobs to be processed at the compute cluster at the same time>[50]
```

Suggested citation:

this repository


# References
Pucker B., Pandey A., Weisshaar B. and Stracke R. (2020). The R2R3-MYB gene family in banana (_Musa acuminata_): genome-wide identification, classification and expression patterns.  PLoS ONE 15(10): e0239275. doi:[10.1371/journal.pone.0239275](https://doi.org/10.1371/journal.pone.0239275).

Haak, M., Vinke, S., Keller, W., Droste, J., Rückert, C., Kalinowski, J., & Pucker, B. (2018). High Quality _de novo_ Transcriptome Assembly of _Croton tiglium_. Frontiers in Molecular Biosciences, 5. doi:[10.3389/fmolb.2018.00062](https://doi.org/10.3389/fmolb.2018.00062).

Pellegrinelli V., Rodriguez-Cuenca S., Rouault C., Schilbert H., Virtue S., Moreno-Navarrete J. M., Bidault G., del Carmen Vazquez Borrego M., Dias A. R., Pucker B., Dale M., Campbell M., Carobbio S., Aron-Wisnewsky J., Mora S., Masiero M., Emmanouilidou A., Mukhopadhyay S., Dougan G., den Hoed M., Loos R., Fernandez-Real J. M., Chiarugi D., Clement K. and Vidal-Puig T. Dysregulation of macrophage PEPD in obesity determines adipose tissue fibro-inflammation and insulin resistance. doi:[10.21203/rs.3.rs-57182.](https://doi.org/10.21203/rs.3.rs-57182/v1)
