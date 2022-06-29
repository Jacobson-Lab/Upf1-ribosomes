# Upf1-ribosomes
 
Analysis scripts and processed data for Ganesan et al.<br/>
Preprint: [pending]

### Data availability
* Raw sequencing data generated in this study have been deposited and are available at NCBI's Gene Expression Omnibus (GEO) under accession number GSE186795.
* The mass spectrometry proteomics data have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifier PXD029577 and 10.6019/PXD029577.
* Numerical data underlying the plots and the R scripts used to generate them are in the folder _**Figures**_ or within the analysis scripts (below)

### Mass spectrometry
* Processed data (exported from Scaffold) include normalized iBAQs (abundance quantification), Protein Identification Probability, and Total Spectrum Count for each sample, available in _**Processed data**/**Mass spec**_
* Criteria for valid identification of protein in a sample:
	* Protein Identification Probability >= 99%
	* Total Spectrum Count >= 2
	* normalized iBAQ > 0 for at least 2 biological replicates (if applicable)
* _**Analysis scripts/Mass spec analysis**/Mass_spec_analysis_limma.Rmd_ consolidates data, identifies differentially recovered proteins, and plots Figure 2.

### RNA-Seq and Ribo-Seq sequence alignment and transcript abundance quantification
* Transcriptome used for sequence alignment is available at https://github.com/Jacobson-Lab/yeast_transcriptome_v5
* Transcript abundance results are combined and available in the folder _**Processed data**/**RSEM output**_
	* "expected_count" column from isoforms.results files
	* "FPKM" column from isoforms.results files
	* "TPM" column from isoforms.results files

### Analyses of changes in transcript abundance using DESeq2
* Analysis scripts are in _**Analysis scripts/Sequencing data analyses**_
* Differential expression (RNA-Seq) between yeast strains
  * _RNAseq_cormatrix_PCA.Rmd_ shows reprodicibility between replicates and produces Figure S10.
  * _RNAseq_analysis_DESeq2.Rmd_ performs differential expression analysis and produces Figure S3.
* IP vs Total ribosomes (Ribo-Seq)
  * _Riboseq_cormatrix_PCA.Rmd_ shows reprodicibility between replicates and produces Figure S9.
  * _Riboseq_analysis_DESeq2.Rmd_ performs differential expression analysis, comparative analysis, and produces Figure 5.
* Ribosome occupancy (Total Ribo-Seq / RNA-Seq)
	 * _RNA-vs-RiboTotal.Rmd_ produces ribosome occupancy data used for Figure 5B.

### Metagene analyses
* Initial processing of bam files was done using riboWaltz package (https://github.com/LabTranslationalArchitectomics/riboWaltz). Either _reads_list_ or _reads_psite_list_ data tables were used in further analysis steps. 
* Analysis scripts are in _**Analysis scripts/Metagene analyses**_
* Footprint length distribution
	* _rl_dist.R_
	* Figure 2
* Mapping of the 5’ and 3’ ends of footprints
	* Calculated by riboWaltz's `rends_heat` function. 
	* Figure 3
* Distribution of footprint abundance across the coding region
	* Metagene: 
		* _binning.R_
		* Figure 4
	* Individual gene:
		* _binning2.R_
		* Figure S6
* Ribo-Seq diagnostic
  * 3-nucleotide periodicity:
    * Calculated by riboWaltz's `metaprofile_psite` function.
    * Figure S5
  * Footprint's P-sites 
  	* _psite_region_frame_fraction.R_
  	* Fraction in mRNA regions: Figure S8 A-B, Figure 7A
  	* Fraction of footprint's reading frames in an mRNA region: Figure S8 C-D, Figure 7B
* Unless otherwise indicated/shown, the output table of replicate libraries were averaged. 
* The results and scripts for plotting them are available in _**Figures**/data_ and _**Figures**/scripts_, respectively.

### Codon optimality
* Codon optimality value for each coding sequence is calculated based on https://github.com/mariodosreis/tai
* _**Analysis scripts/Codon optimality/**codon_optimality_tAI_CDS.R_ calculates codon optimality scores used for Figure 5D

### A-site codon occupancy
* Analysis scripts and required files for processing are in _**Analysis scripts/A-site codon occupancy/**_
* Calculation of mean relative occupancy:
	1. _write_rpl_by_size.R_ prepares riboWaltz's _reads_psite_list_ data tables for the analysis and export them as txt files.
	2. _codon_window_count_by_codonpos_from_rpl.py_ counts number of footprints whose A-site (or P- or E-site, as specified) are within a specified nucleotide window of a specified codon of interest.
	3. _calc_REV_window_codonpos_v2.R_ calculates mean relative occupancy.
	* _automate_codon_window_count_by_codonpos_from_rpl.sh_ automates steps 2 and 3 for multiple codons in a given codon list (such as codons_list.txt, which contains all 64 codons)
	* _automate_automate_codon_window_count_by_codonpos_from_rpl.sh_ submits the above script as multiple jobs to process multiple samples in parallel.
	* Results of all samples are combined and provided as A-site_codon_occupancy_window30_bycodon_allsamples.txt.
* _codon_occupancy_analysis.Rmd_ further analyzes the results, plots Figure 6 and Figure S7.
