# Kronos.scRT

Based on Kronosc_scRT <a href="https://www.nature.com/articles/s41467-022-30043-x" title="(Gnan et al. 2022)">(Gnan et al. 2022) </a>, Kronos.scRT is an R package containing 4 shiny apps that will lead the user through the single-cell analysis of replication timing data. 

#### Table of content
  * [Installation and preparation](#Installation-and-preparation)
  * [Shiny Apps](#Shiny-Apps)
  * [ Tutorial 1](#Tutorial-1---R-package)
    * [1.1 R-package](#1.1-R-package)
    * [1.2 Shell](#1.2-Shell)
  * [Tutorial 2](#Tutorial-2)
    * [2.1 R-package](#2.1-R-package)
    * [2.2 Shell](#2.2-Shell)
  * [Kronos scRT pipeline from shell](#Kronos-scRT-pipeline-from-shell)
  * [Citation](#Citation)
  * [Authors](#Authors)


#### Installation and preparation

To install the package, first install devtools 

    install.packages("devtools")

and then Kronos.scRT

    devtools::install_github("CL-CHEN-Lab/User_interface_for_Kronos_scRT", type = "source")

it is as well possible to use Kronos scRT from command line. To do so, it is possible to run the command

    Kronos.scRT::shell_interface()

that will return a line to add to the user's shell profile file to create the alias Kronos.

#### Shiny Apps
This package contains 4 shiny apps:

- Kronos.scRT::Pre_processing()
 <p align="center">
<img src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/Pre.png">
</p>
    -   Trims fastq files and maps demultiplexed data 
    -   Creates bins files to be used later on in the analysis
    -   Maps data and calls genomic regions Copy number
    
- Kronos.scRT::Processing()
 <p align="center">
<img src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/Pro.png">
</p>
    - Calculates scReplication Timing profiles and provides an initial overview of the data

- Kronos.scRT::scRT_tools()
 <p align="center">
<img  src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/scRTtools.png">
</p>
    - Provides multiple investigation tools for scRT data

- Kronos.scRT::RT_tools()
 <p align="center">
<img  src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/RTtools.png">
</p>
    - Provides multiple investigation tools bulk and pseudo bulk RT data
    
### R package analysis workflow

#### Tutorial 1

You can download data published in <a href="https://doi.org/10.1038/s41588-019-0474-z" title="Miura et al. 2019">Miura et al. 2019</a> from SRA using SRP130912 identifier. A selection of G1- and Mid-S- phase cells are available <a href="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/data/SRP130912_G1_MidS_cells_to_download.txt">here</a>.
#### 1.1 R-package

The first step consists in trimming our single-cell data and mapping them against a reference genome. To do so, we use the function FastqToBam.

    Kronos.scRT::FastqToBam(
      bowtie2_index = '~/location_bowtie2_Index/mm10',
      File1 = list.files(path ='~/location_SRP130912_fastq_files/',full.names = T),
      outputdir = '~/outputDirectory',
      cores = 6)

Once data have been mapped we bin our reference genome and calculate mappability and GC content per bin.

    bins_ms = Kronos.scRT::binning(
      RefGenome = '~/Path_to_Fasta_file/mm10.fa',
      bowtie2_index = '~/location_bowtie2_Index/mm10',
      directory_to_bamfiles = '~/outputDirectory/BAM/',
      cores = 6
    )
    
The next step consists in calling the CN. In this specific case we are choosing to block the mean ploidy of these cells around 2, but a range of allowed ploidy can be provided.

    Chromsize = readr::read_tsv(url('http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes'), col_names = c('chr', 'size'))

    SingleCell = Kronos.scRT::CallCNV(
      directory = '~/outputDirectory/BAM/',
      chrom_size = Chromsize,
      bins = bins_ms,
      basename = 'MsESC',
      ploidy = 2,
      cores = 6
    )

SingleCell is a list of two data frames: PerCell and CNV

The data you downloaded come with metadata attached, therefore we can provide to the software a data frame containing two columns: Cell and S_Phase. The Cell column contains all the single-cell file names while the S_Phase column contains logical values indicating whether a particular cell is in S Phase (TRUE) or not (FALSE). If you do not have this info for your data pass to the next step. For this tutorial a list has been provided.

    WhoIsWho = readr::read_tsv(url('https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/data/WhoIsWho_SRP130912.tsv'))
    
    SingleCell$PerCell = Kronos.scRT::WhoIsWho(PerCell = SingleCell$PerCell, WhoIsWho = WhoIsWho)
    
The diagnostic function can be in interactive mode (a shiny app) or not. In our specific case, after running the shiny app we want to select the manual option and proceed. For unsorted single-cell samples, it is possible to proceed with the automatic option that will be used to estimate S-phase cells. If Staging information is not available the user will have to use the Variability vs Ploidy plot to select a variability threshold to identify the S-phase.  <a href="https://github.com/CL-CHEN-Lab/Kronos_scRT" title="Kronos_scRT">As reference: figure 1a</a>

    Dagnostic_output = Kronos.scRT::diagnostic(SingleCell$PerCell)
 <p align="center">
<img src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/Diagnostic1.png">
</p>
The second stage of Diagnostic consist in adjusting the S-phase cells ploidy in order to reconstitute a continus S phase if possible.<a href="https://github.com/CL-CHEN-Lab/Kronos_scRT" title="Kronos_scRT">As reference: figure 1b</a>
 <p align="center">
<img src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/Diagnostic2.png">
</p>
The output of the diagnostic function is a list with a Settings data frame and 5 diagnostic plots. We can now proceed with the adjustment of our data using the diagnostic info.


    SingleCell$PerCell = Kronos.scRT::AdjustPerCell(PerCell = SingleCell$PerCell,
                                                    Settings = Dagnostic_output$Settings)
                                                    
                                                    
    SingleCell$CNV = Kronos.scRT::AdjustCN(PerCell = SingleCell$PerCell, scCN = SingleCell$CNV)
    
We now have to decide what will be our final resolution and create a second binning list using the function GenomeBinning

    Bins = Kronos.scRT::GenomeBinning(
      Chr_size = Chromsize,
      size = 200000,
      Chr_filter =paste0('chr', 1:19),
      Cores = 6
    )

Singe cell CN data will therefore be rebinned based on their Phase

    SingleCell$SPhase = Kronos.scRT::Rebin(PerCell = SingleCell$PerCell,scCN = SingleCell$CNV,Bins = Bins,Sphase = T)

    SingleCell$G1G2 = Kronos.scRT::Rebin(PerCell = SingleCell$PerCell,scCN = SingleCell$CNV,Bins = Bins,Sphase = F)

To call replicated and unreplicated regions we need to identify a G1G2 median profile with the function BackGround

    SingleCell$MedianG1G2=Kronos.scRT::BackGround(G1_scCN = SingleCell$G1G2)

We are now ready to call the replication state of each bin

    SingleCell$SPhase = Kronos.scRT::Replication_state(
      Samples = SingleCell$SPhase,
      background = SingleCell$MedianG1G2,
      Chr_filter = Chrom_to_keep,
      cores = 6
    )

The FilterCells function allows removing cells that are too different from the rest of the population (if present)
    
    FilterS=Kronos.scRT::FilterCells(scCN = SingleCell$SPhase,ApplyFilter = T)

this function returns a vector containing a filtered dataset FilteredData and a cell-matrix plots before and after filtering
 <p align="center">
<img src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/BeforeF.png" alt="Correlation Matrix before filter" width="400"> <img src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/AfterF.png" alt="Correlation Matrix before filter" width="400">
</p>
if needed we can overite the Sphase data with the filtered ones.

    SingleCell$SPhase = FilterS$FilteredData

We can finally calculate a pseudobulk RT from our data 

    SingleCell$pseudobulk=Kronos.scRT::pseudoBulkRT(S_scCN = SingleCell$SPhase)

as well as reformat an ESC mm10 bulk RT bedgraph to use as a control. In our data folder a liftover bulk RT profile from <a href="https://doi.org/10.1038/s41588-019-0474-z" title="Miura et al. 2019">Miura et al. 2019</a> has been provided.

    Reference=readr::read_tsv(url('https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/data/mm10_GSE108556_P239_01_02_rpm_w200ks80k_BrdUIP_Percent_q0.05.bedGraph'),col_names = c('chr','start','end','RT'))
    Reference = Kronos.scRT::RebinRT(
      RT = Reference,
      Bins = Bins,
      Chr_filter = Chrom_to_keep,
      Basename = 'ESCReference',
      Group = 'MsESC'
    )

we can now explore the data using various investigation plots such as:

Correlation plots
    
    Kronos.scRT::KCorr_plot(
      SingleCell$pseudobulk,
      Reference,
      method = 'spearman')
      
 <p align="center">    
<img src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/Corr.png" width="400">
</p>
scRT genomic regions plots

    Kronos.scRT::scRTplot(
      pseudoBulkRT = rbind(SingleCell$pseudobulk,
                           Reference),
      S_scCN = SingleCell$SPhase,
      Coordinates = list(chr = 'chr1', start = 2800000, end = 12800000) ,
      rasterized_heatmap = T
    )
    
 <p align="center">    
<img  src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/scRTplot.png">
</p>
Explore Variability based on RT or regions of interest

    Var=Kronos.scRT::Variability(S_scCN=SingleCell$SPhase,scRT=SingleCell$pseudobulk)
    
    Var=Kronos.scRT::TW_RTAnnotation(Variability=Var,RT_Groups=2)
    
    Fit_Data=Kronos.scRT::Twidth_fit_data(
                              df = Var,
                              ncores = 6
                            )
                            
     Twidth = Kronos.scRT::Twidth(Fit_Data)
     
     Kronos.scRT::Twidth_barplot(Variability = Var,Twidth = Twidth)
     
 <p align="center">
<img  src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/TW.png" width="400">
</p>

#### 1.2 Shell

The same analysis can be run from command line as follows.

To trim and map our data we used the command Kronos fastqtoBAM

    Kronos fastqtoBAM -O /dir_to_fastq/SRR6491868_1.fastq,/dir_to_fastq/SRR6491869_1.fastq,..,/dir_to_fastq/SRR6491908_1.fastq \
    -i /dir_to_bowtie2_index/mm10 \
    -o fastqtoBAM_output_directory \
    -c 6
    
inside the output folder we will find the following sub-folders:
  - FastqToBamMetrics, containing some quality control metrics
  - trimmed, containing our trimmed data
  - BAM, containing our mapped data
   
To create the bins needed for the copy number calling we use Kronos binning 

    Kronos binning -R /dir_to_reference_genome/mm10.fa \
    -o binning_output_directory \
    -i /directory_to_bowtie2_index/mm10 \
    -d /fastqtoBAM_output_directory/BAM/ \
    -c 6
    
To calculate the CNV tracks and the PerCell file we use Kronos CNV

    Kronos CNV -D /fastqtoBAM_output_directory/BAM/ \
    -B binning_output_directory/bins \
    -o CNV_output_directory \
    -C Chromosome_size_file \
    -c 6 -e MsESC -g MsESC -p 2 
    
The outputs of this command are two files: MsESC_scCNV.tsv and MsESC_PerCell.csv

As before we can impose its cycling stage to each cell using the metadata here deposeted.

    Kronos WhoIsWho -F CNV_output_directory/MsESC_PerCell.csv \
    -W Metadata.tsv \
    -o WhoIsWho_output_directory

that will create a phased_MsESC_PerCell.csv file that will be used as an input of Kronos diagnostic. Since in this case there is no interactive session, the factors to correct the S-phase progression have to be provided using the -f and -s options.

    Kronos diagnostic -f WhoIsWho_output_directory/phased_MsESC_PerCell.csv \
    -o diagnostic_output_directory -f 0.889 -s 1

We are now ready to calculate scRT tracks, pseudo-bulk RT and variability files.

    Kronos RT -F WhoIsWho_output_directory/phased_MsESC_PerCell.csv \
    -T CNV_output_directory/MsESC_scCNV.tsv \
    -C Chromosome_size_file \
    -B 200Kb \
    -o RT_output_directory \
    -f Mouse_data \
    -S diagnostic_output_directory/settings.txt

And optionally resize a reference RT for comparison

    Kronos RefRT -R ReferenceRT.tsv \
    -N Referece # the name that will be used in the plots
    -G MsESC # the same group name of your data
    -B RT_output_directory/Mouse_data_bins.rds
    -o RefRT_output_directory
    -f Mouse
    
We have now everything to obtain our plots 

- Correlation plots
    
      Kronos Corr -F RT_output_directory/MsESC_calculated_replication_timing_200Kb.tsv,RefRT_output_directory/Mouse_reference_replication_timing_200kb.tsv \
      -f Ms_bulk_vs_pseudobulk
      -o Corr_output_directory

- scPlots 

      Kronos scPlots -R RT_output_directory/MsESC_calculated_replication_timing_200Kb.tsv \
      -C RT_output_directory/MsESC_single_cells_CNV_200Kb.tsv \
      -E RefRT_output_directory/Mousere_ference_replication_timing_200kb.tsv \
      -r chr1:2800000-12800000 \
      -o scPlot_output_directory
    
- TW plots
  
      Kronos compare TW -F /Users/sgnan/Desktop/TestRT/MsESC_variability.tsv \
      -o compare_TW_output_directory

#### Tutorial 2

To give a more comprhesive idea of Kronos scRT potential a bigger dataset form <a href="https://www.nature.com/articles/s41467-022-30043-x" title="(Gnan et al. 2022)">(Gnan et al. 2022) </a> is provided with this package.

    SingleCell=Kronos.scRT::MCF7_subpop1

#### 2.1 R-package

Using the diagnostic function we can dived S-phase cells from G1/G2 cells selecting a buffer reagin between the low variability G1/G2-phase cells and the two arms made of high variability S-phase cells.

    Dagnostic_output = Kronos.scRT::diagnostic(SingleCell$PerCell)
 <p align="center">
<img  src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/MCF7_Diagnostic1.png">
</p>
Since the cells of this dataset span through most of the S phase, we can use the automatic option of Kronos scRT to fix the S-phase progression.
 <p align="center">
<img  src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/MCF7_Diagnostic2.png">
</p>
We use the parameters estimated with the diagnostic function to correct our data as before


    SingleCell$PerCell = Kronos.scRT::AdjustPerCell(PerCell = SingleCell$PerCell,Settings = Dagnostic_output$Settings)
                                                
    SingleCell$CNV = Kronos.scRT::AdjustCN(PerCell = SingleCell$PerCell, scCN = SingleCell$CNV)
    
We now have to decide what will be our final resolution and create a second binning list using the function GenomeBinning

    Chromsize = readr::read_tsv(url('http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'), col_names = c('chr', 'size'))

    Chrom_to_keep= paste0('chr', 1:22)

    Bins = Kronos.scRT::GenomeBinning(
      Chr_size = Chromsize,
      size = 200000,
      Chr_filter =Chrom_to_keep,
      Cores = 6
    )

Singe cell CN data will therefore be rebinned based on their Phase
    
    SingleCell$SPhase = Kronos.scRT::Rebin(PerCell = SingleCell$PerCell,scCN = SingleCell$CNV,Bins = Bins,Sphase = T)
    
    SingleCell$G1G2 = Kronos.scRT::Rebin(PerCell = SingleCell$PerCell,scCN = SingleCell$CNV,Bins = Bins,Sphase = F)

As before we calculate the G1G2 median profile of this data with the function BackGround

    SingleCell$MedianG1G2=Kronos.scRT::BackGround(G1_scCN = SingleCell$G1G2)

and we call the replication state of each bin is S-phase 

    SingleCell$SPhase = Kronos.scRT::Replication_state(
      Samples = SingleCell$SPhase,
      background = SingleCell$MedianG1G2,
      Chr_filter = Chrom_to_keep,
      cores = 6
    )

and in G1/G2 as well

    SingleCell$G1G2 = Kronos.scRT::Replication_state(
      Samples = SingleCell$G1G2,
      background = SingleCell$MedianG1G2,
      Chr_filter = Chrom_to_keep,
      cores = 6
    )

Remember always to check if there are cells that have to be excluded ( this can be done as well for G1/G2 Cell)
    
    FilterS=Kronos.scRT::FilterCells(scCN = SingleCell$SPhase,ApplyFilter = T,min_cor = 0.3)
    FilterG=Kronos.scRT::FilterCells(scCN = SingleCell$G1G2,ApplyFilter = T,min_cor = 0.3)

 <p align="center">
<img src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/FilterMCF7.png" alt="Correlation Matrix">
</p>
Filters can be applied as followed

    SingleCell$SPhase = FilterS$FilteredData
    SingleCell$G1G2 = FilterG$FilteredData
    
Calculate pseudobulk and load reference RT 

    SingleCell$pseudobulk=Kronos.scRT::pseudoBulkRT(S_scCN = SingleCell$SPhase)
    
    Reference = Kronos.scRT::RebinRT(
      RT = Kronos.scRT::MCF7_Reference,
      Bins = Bins,
      Chr_filter = Chrom_to_keep,
      Basename = 'MCF7 Reference',
      Group = 'MCF7 subpopulation1'
    )

Calculate correlation between Bulk and Pseudobulk
    
    Kronos.scRT::KCorr_plot(
      SingleCell$pseudobulk,
      Reference,
      method = 'spearman')
      
 <p align="center">   
<img  src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/MCF7_Corr.png" alt="Correlation Matrix" width="400">
</p>
and explore genomic regions using the following command

    Kronos.scRT::scRTplot(
      pseudoBulkRT = rbind(SingleCell$pseudobulk,
                           Reference),
      S_scCN = SingleCell$SPhase,
      Coordinates = list(chr = 'chr1', start = 40000000, end = 100000000) ,
      rasterized_heatmap = T)
      
 <p align="center">
<img  src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/MCF7_scPlot.png" alt="ScPlots">
</p>
This dataset allows to explore variability at a much higher resolution compared to the one before we can therefore up to 5 RT categories to calculate the TW.

    Var=Kronos.scRT::Variability(S_scCN=SingleCell$SPhase,scRT=SingleCell$pseudobulk)
    Var=Kronos.scRT::TW_RTAnnotation(Variability=Var,RT_Groups=5)
    Fit_Data=Kronos.scRT::Twidth_fit_data(
      df = Var,
      ncores = 6
    )
    Twidth = Kronos.scRT::Twidth(Fit_Data)
    
    Kronos.scRT::Twidth_extended_plot(Variability = Var,Fitted_data = Fit_Data,Twidth = Twidth)

    
 <p align="center">
<img  src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/TW_ext.png" alt="Correlation Matrix">
</p>

    Kronos.scRT::Twidth_barplot(Variability = Var,Twidth = Twidth)

 <p align="center">
<img  src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/TW_MCF7.png" >
</p>
We can as well test whether the differences between RT categories are statistically significant and plot the pvalues together with the barplots providing the output of Twidth_pval as an input to Twidth_barplot. N.B. this operation is highly time consuming. 

    pval=Kronos.scRT::Twidth_pval(
    variability = Var,
    twidth = Twidth,
    alternative = 'two.sided' ,
    nIterations = 10^4,
    ncores = 6, 
    adjust.methods = 'none')
    
    Kronos.scRT::Twidth_barplot(Variability = Var,Twidth = Twidth, pval = pval)

Another way to visualize RT variability is the bin probability of being replicated in function of its average replication timing in different portions of the S phase.

    BinProbS=Kronos.scRT::Prepare_S_phase_cells_forBinRepProb(S = SingleCell$SPhase,RT = SingleCell$pseudobulk)
    BinProbG=Kronos.scRT::Prepare_G1G2_phase_cells_forBinRepProb(G1.G2 = SingleCell$G1G2,RT = SingleCell$pseudobulk)
    
    Kronos.scRT::BinRepProbPlot(Variability = rbind(BinProbS,BinProbG))

<p align="center">
<img  src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/BinProb.png">
</p>

#### 2.2 Shell

To perform this analysis using your command line you need to save the MFC7 data into a folder as follow: 

    readr::write_tsv(Kronos.scRT::MCF7_subpop1$CNV,file = 'Input_dir/MCF7_subpop1_CNV.tsv')
    readr::write_csv(Kronos.scRT::MCF7_subpop1$PerCell,file = 'Input_dir/MCF7_subpop1_PerCell.csv')

Using Kronos diagnostic we can identify the S and G1/G2 populatio as well as reconstruct the S-phase progression.

    Kronos diagnostic -f Input_dir/MCF7_subpop1_PerCell.csv \
    -o diagnostic_output_directory 
    -S 0.9 -G 0.85 \
    -m 117 -C

We are now ready to calculate scRT tracks, pseudo-bulk RT and variability files.

    Kronos RT -F Input_dir/MCF7_subpop1_PerCell.csv \
    -T Input_dir/MCF7_subpop1_CNV.tsv \
    -C Chromosome_size_file \
    -B 200Kb \
    -o RT_output_directory \
    -f MCF7 \
    -S diagnostic_output_directory/settings.txt \
    --extract_G1_G2_cells

And optionally resize a reference RT for comparison

    Kronos RefRT -R ReferenceRT.tsv \
    -N Referece # the name that will be used in the plots
    -G 'MCF7 subpopulation1' # the same group name of your data
    -B RT_output_directory/MCF7_data_bins.rds
    -o RefRT_output_directory
    -f MCF7
    
We have now evetrything to obtain our plots 

- Correlation plots
    
      Kronos Corr -F RT_output_directory/MCF7_calculated_replication_timing_200Kb.tsv,RefRT_output_directory/MCF7_reference_replication_timing_200kb.tsv \
      -f MCF7_bulk_vs_pseudobulk
      -o Corr_output_directory

- scPlots 

      Kronos scPlots -R RT_output_directory/MCF7_calculated_replication_timing_200Kb.tsv \
      -C RT_output_directory/MCF7_single_cells_CNV_200Kb.tsv \
      -E RefRT_output_directory/MCF7_reference_replication_timing_200kb.tsv \
      -r chr1:40000000-100000000 \
      -o scPlot_output_directory
    
- TW plots

      Kronos compare TW -F /Users/sgnan/Desktop/TestRT/MCF7_variability.tsv \
      -o compare_TW_output_directory
    
using the option -p it is possible as well to calculate a pvalue

- BinRepProb plots

      Kronos BinRepProb -S RT_output_directory/MCF7_single_cells_CNV_200Kb.tsv \
      - G RT_output_directory/MCF7_G1_G2_single_cells_CNV_200Kb.tsv \
      -o BinRepProb_output_directory -f MCF7Sub1
    
#### Kronos scRT pipeline from shell

<p align="center">
<img  src="https://github.com/CL-CHEN-Lab/User_interface_for_Kronos_scRT/blob/main/img/schema.png">
</p>

#### Citation
If you use Kronos scRT, please cite the following paper:

Gnan, S., Josephides, J.M., Wu, X. et al. Kronos scRT: a uniform framework for single-cell replication timing analysis. Nat Commun 13, 2329 (2022). https://doi.org/10.1038/s41467-022-30043-x
#### Authors

Please contact the authors for any further questions:

[Stefano Gnan](mailto:stefano.gnan@curie.fr) and [Chunlong Chen](mailto:chunlong.chen@curie.fr) (Institut Curie)
