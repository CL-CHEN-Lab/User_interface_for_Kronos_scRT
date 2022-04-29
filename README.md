# Kronos.scRT

Based on Kronosc_scRT <a href="https://www.nature.com/articles/s41467-022-30043-x" title="(Gnan et al. 2022)">(Gnan et al. 2022) </a>, Kronos.scRT is an R package containing 3 shiny apps that will lead the user through the single-cell analysis of replication timing data. 

To install the package, first install devtools 

    install.packages("devtools")

and then Kronos.scRT

    devtools::install_github("Derfen3001/Kronos.scRT")

This package contains 3 shiny apps:

- Kronos.scRT::Pre_processing()
![](https://github.com/Derfen3001/Kronos.scRT/blob/master/img/Pre.png)
    -   Trims fastq files and maps demultiplexed data 
    -   Creates bins files to be used later on in the analysis
    -   Maps data and calls genomic regions Copy number
    
- Kronos.scRT::Processing()
![](https://github.com/Derfen3001/Kronos.scRT/blob/master/img/Pro.png)
    - Calculates scReplication Timing profiles and provides an initial overview of the data

- Kronos.scRT::Post_processing()
![](https://github.com/Derfen3001/Kronos.scRT/blob/master/img/Post.png)
    - Provides multiple investigation tools
    
### R package analysis workflow

You can download SRP130912 data from SRA

The first step consists in trimming our single-cell data and mapping them against a reference genome. To do so, we use the function FastqToBam.

    Kronos.scRT::FastqToBam(
      bowtie2_index = '~/location_bowtie2_Index/mm10',
      File1 = '~/location_SRP130912_fastq_files/',
      outputdir = '~/outputDirectory',
      cores = 6)

Once data have been mapped we bin our reference genome and calculate mappability and GC content per bin.

    bins_ms = Kronos.scRT::binning(
      RefGenome = '~/Path_to_Fasta_file/mm10.fa',
      bowtie2_index = '~/location_bowtie2_Index/mm10',
      directory_to_bamfiles = '~/outputDirectory/BAM/',
      cores = 6
    )
    
The next step consists in calling the CN
Chromsize = readr::read_tsv(url('http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes'), col_names = c('chr', 'size'))


    SingleCell = Kronos.scRT::CallCNV(
      directory = '~/outputDirectory/BAM/',
      chrom_size = Chromsize,
      bins = bins_ms,
      basename = 'MsESC',
      cores = 6
    )

SingleCell is a list of two data frames: PerCell and CNV

The data you downloaded come with metadata attached, therefore we can provide to the software a data frame containing two columns: Cell and S_Phase. The Cell column contains all the single-cell file names while the S_Phase column contains logical values indicating whether a particular cell is in S Phase (TRUE) or not (FALSE). If you do not have this info for your data pass to the next step

    SingleCell$PerCell = Kronos.scRT::WhoIsWho(PerCell = SingleCell$PerCell, WhoIsWho = WhoIsWho)
    
The diagnostic function can be in interactive mode (a shiny app) or not. In our specific case, after running the shiny app we want to select the manual option and proceed. For unsorted single-cell samples, it is possible to proceed with the automatic option that will be used to estimate S-phase cells. If Staging information is not available the user will have to use the Variability vs Ploidy plot to select a variability threshold to identify the S-phase.  <a href="https://github.com/CL-CHEN-Lab/Kronos_scRT" title="Kronos_scRT">As reference: figure 1a</a>

Dagnostic_output = Kronos.scRT::diagnostic(SingleCell$PerCell)

![](https://github.com/Derfen3001/Kronos.scRT/blob/master/img/Diagnostic1.png)

The second stage of Diagnostic consist in adjusting the S-phase cells ploidy in order to reconstitute a continus S phase if possible.<a href="https://github.com/CL-CHEN-Lab/Kronos_scRT" title="Kronos_scRT">As reference: figure 1b</a>

![](https://github.com/Derfen3001/Kronos.scRT/blob/master/img/Diagnostic2.png)

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

The FilterCells function allows removing cells that are too different from the rest of the population 
    
    FilterS=Kronos.scRT::FilterCells(scCN = SingleCell$SPhase,ApplyFilter = T)

this function returns a vector containing a filtered dataset FilteredData cell-matrix plots before and after filtering

![](https://github.com/Derfen3001/Kronos.scRT/blob/master/img/BeforeF.png) ![](https://github.com/Derfen3001/Kronos.scRT/blob/master/img/AfterF.png)



    SingleCell$SPhase = FilterS$FilteredData

We can finally calculate a pseudobulk RT from our data 

    SingleCell$pseudobulk=Kronos.scRT::pseudoBulkRT(S_scCN = SingleCell$SPhase)

as well as reformat an ESC mm10 bulk RT bedgraph to use as a control

    Reference=readr::read_tsv('~/referenceRT.tsv')
    Reference = Kronos.scRT::RebinRT(
      RT = Reference,
      Bins = Bins,
      Chr_filter = Chrom_to_keep,
      Basename = 'ESCReference',
      Group = 'MsESC'
    )

we can now explore the data using various investigation plots such as:
- correlation plots
    
    Kronos.scRT::KCorr_plot(df = rbind(
      SingleCell$pseudobulk,
      Reference
    ), method = 'spearman')


- scRT genomic regions plots

    Kronos.scRT::scRTplot(
      pseudoBulkRT = rbind(SingleCell$pseudobulk,
                           Reference),
      S_scCN = SingleCell$SPhase,
      Coordinates = list(chr = 'chr1', start = 2800000, end = 12800000) ,
      rasterized_heatmap = T
    )

- Explore Variability based on RT or regions of interest

    Var=Kronos.scRT::Variability(S_scCN=SingleCell$SPhase,scRT=SingleCell$pseudobulk)
    
    Var=Kronos.scRT::TW_RTAnnotation(Variability=Var,RT_Groups=2)
    
    Fit_Data=Kronos.scRT::Twidth_fit_data(
                              df = Var,
                              ncores = 6
                            )
                            
     Twidth = Kronos.scRT::Twidth(Fit_Data)
     
    Kronos.scRT::Twidth_extended_plot(Variability = Var,Fitted_data = Fit_Data,Twidth = Twidth)
    
![](https://github.com/Derfen3001/Kronos.scRT/blob/master/img/TW_extended.png)
    
or

    Kronos.scRT::Twidth_barplot(Variability = Var,Twidth = Twidth)

![](https://github.com/Derfen3001/Kronos.scRT/blob/master/img/TW.png)

- Visualize variability as Bin probability of Replication

    BinProb=Kronos.scRT::Prepare_S_phase_cells_forBinRepProb(S = SingleCell$SPhase,RT = SingleCell$pseudobulk)

    Kronos.scRT::BinRepProbPlot(Variability = BinProb)

![](https://github.com/Derfen3001/Kronos.scRT/blob/master/img/BinProb.png)


### Authors

Please contact the authors for any further questions:

[Stefano Gnan](mailto:stefano.gnan@curie.fr)and [Chunlong Chen](mailto:chunlong.chen@curie.fr) (Institut Curie)
