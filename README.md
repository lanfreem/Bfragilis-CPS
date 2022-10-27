Scripts used in analysis of DoTA-seq data of B fragilis CPS operons

1. Demultiplexing raw sequence files
DoTA-seq uses the Illumina I7 index read to read out the cell barcodes of the library. Thus, one needs to obtain the raw I7 read file from the sequencer in addition to the raw read 1 file. Use R4-parser.ipynb to quality filter the barcodes and associate barcodes with their amplicon reads.

2. Mapping reads to reference
Using the output of R4-parser, we use Bowtie2 to map the resulting barcode-associated reads to a reference database of expected amplicons (invertons-reference.fa), resulting in a SAM file.

3. Filtering of mapped reads
Using the outputs of Bowtie2 mapping to the database of expected amplicons, we group the reads by barcode, perform quality control and filtering on the grouped reads, and generate a summary of the output. The steps are detailed in the SAM-parser.ipynb Jupyter notebook.

4. Analysis of filtered grouped read data
The output of SAMparser is analyzed using an Rscript (CPS-analysis.R) in combination with a list of CPS combos (allcpscombos.xlsx), which filters out barcodes that do not have enough reads for each of the loci and cells with conflicting reads for any loci (both ON and OFF). It outputs an excel file, as well as a list of figures including the heatmaps used to represent the populations in figures 1 and 2. It also outputs a file called X_summary.xlsx which contains the details of the analysis, including how many barcodes passed the quality check. The results of this file contains a table of each CPS promoter variant and the number of cells sequenced, which represents the population composition. This data is extracted as a separate table with a .nodes extension and further analyzed using CPS-network-analysis.ipynb Jupyter notebook to plot the network graphs. CPS-summary-analysis.ipynb Jupyter notebook contains all other analyses excluding the model fitting analysis, which is performed in matlab and described in MCMC-scripts.zip.

5. Fitting stochastic model to timeseries single cell data
# Bfragilis-CPS
