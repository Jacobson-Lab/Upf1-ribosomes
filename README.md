# Upf1-ribosomes
 
Analysis scripts and processed data for Ganesan et al.<br/>
Preprint: [pending]

### Data availability
* Raw sequencing data generated in this study have been deposited and are available at NCBI's Gene Expression Omnibus (GEO) under accession number GSE186795.
* The mass spectrometry proteomics data have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifier PXD029577 and 10.6019/PXD029577.
* Numerical data underlying the plots and the R scripts used to generate them are in the folder _**Figures**_

### Mass spectrometry
* Processed data (exported from Scaffold) include normalized iBAQs (abundance quantification), Protein Identification Probability, and Total Spectrum Count for each sample, available in _**Processed data**/**Mass spec**_
* Criteria for valid identification of protein in a sample:
	* Protein Identification Probability >= 99%
	* Total Spectrum Count >= 2
	* normalized iBAQ > 0 for at least 2 biological replicates (if applicable)
* Scripts used to consolidate data, identify differentially recovered proteins, and plot Figure 2 are in _**Analysis scripts/Mass spec analysis**_

### RNA-Seq / Ribo-Seq sequence alignment and transcript abundance quantification
* Transcriptome used for sequence alignment is available at https://github.com/Jacobson-Lab/yeast_transcriptome_v5
* Transcript abundance results are combined and available in the folder _**Processed data**/**RSEM output**_
	* "expected_count" column from isoforms.results files
	* "FPKM" column from isoforms.results files
	* "TPM" column from isoforms.results files

### Analyses of changes in transcript abundance using DESeq2
* Differential expression (RNA-Seq) between yeast strains
  * _**Analysis scripts/Sequencing data analyses**/RNAseq_cormatrix_PCA.Rmd_ shows reprodicibility between replicates and produces Figure S10.
  * _**Analysis scripts/Sequencing data analyses**/RNAseq_analysis_DESeq2.Rmd_ produces Figure S3.
* IP vs Total ribosomes (Ribo-Seq)
  * _**Analysis scripts/Sequencing data analyses**/Riboseq_cormatrix_PCA.Rmd_ shows reprodicibility between replicates and produces Figure S9.
  * _**Analysis scripts/Sequencing data analyses**/Riboseq_analysis_DESeq2.Rmd_ produces Figure 5.
* Ribosome occupancy (Total Ribo-Seq / RNA-Seq)
	 * _**Analysis scripts/Sequencing data analyses**/RNA-vs-RiboTotal.Rmd_ produces data for Figure 5B.

### Metagene analyses
* Initial processing of bam files was done using riboWaltz package (https://github.com/LabTranslationalArchitectomics/riboWaltz). Either _reads_list_ or _reads_psite_list_ data tables were used in further analysis steps. 
* Footprint length distribution
	* Calculated by _**Analysis scripts/Metagene analyses**/rl_dist.R_. 
	* Figure 2
* Mapping of the 5’ and 3’ ends of footprints
	* Calculated by riboWaltz's `rends_heat` function. 
	* Figure 3
* Distribution of footprint abundance across the coding region
	* Metagene: 
		* Calculated by _**Analysis scripts/Metagene analyses**/binning.R_. 
		* Figure 4
	* Individual gene:
		* Calculated by _**Analysis scripts/Metagene analyses**/binning2.R_. 
		* Figure S6
* Ribo-Seq diagnostic
  * 3-nucleotide periodicity:
    * Calculated by riboWaltz's `metaprofile_psite` function.
    * Figure S5
  * Fraction of footprint's p-sites in mRNA regions:
    * Figure S8 A-B, Figure 
  * Fraction of footprint's reading frames in an mRNA region:
    * Figure S8 C-D
* Unless otherwise indicated/shown, the output table of replicate libraries were averaged. 
* The results and scripts for plotting them are available in _**Figures**/data_ and _**Figures**/scripts_, respectively.

### Codon optimality

### A-site codon occupancy

