


# Admixture and association mapping pipeline

This pipeline implements admixture mapping on the output of the [AncInf](https://github.com/pmonnahan/AncInf) pipeline (i.e. RFMix output) and/or association mapping on the output of the [pre-imputation](https://github.com/pmonnahan/DataPrep) or [post-imputation](https://github.com/pmonnahan/DataPrep/postImpute) QC pipelines.  It is not completely necessary to generate the input via the mentioned pipelines, although doing so would likely reduce the chances of encountering an issue.  Prior to mapping, input data is parsed and QC'ed to identify (and optionally filter) related samples, calculate principal components, and perform genomic control matching. A more detailed description of each of these steps is provided at the bottom. 

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

 - [Requirements](#requirements)
   - [Snakemake](#snakemake)
   - [Singularity](#singularity)
 - [Running the workflow](#running-the-workflow)
 - [Pipeline Overview](#pipeline-overview)
   - [Input Data](#input-data)
   - [Data Preparation](#data-preparation)
     - [Missingness filters](#missingness-filters)
     - [LD pruning](#ld-pruning)
     - [IBD filters](#ibd-filters)
     - [Hardy-Weinberg equilibrium filters](#hardy-weinberg-equilibrium-filters)
     - [Control Matching  STOPPED HERE](#control-matching--stopped-here)
   - [Population structure and Kinship](#population-structure-and-kinship)
   - [Genome Wide Association Mapping (GWAS)](#genome-wide-association-mapping-gwas)
   - [Admixture Mapping](#admixture-mapping)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Requirements
 
### Snakemake
The pipeline is coordinated and run on an HPC (or locally) using _Snakemake_.  On UMN's HPC, snakemake can be installed by:

    module load python3/3.6.3_anaconda5.0.1
    conda install -c conda-forge -c bioconda snakemake python=3.6

The 'module load' command will likely need to be run each time prior to use of Snakemake.

Alternatively, you can try installing snakemake via _pip_:

    pip3 install --user snakemake pyaml

### Singularity

The installation of the individual programs  used throughout this pipeline can be completely avoid by utilizing a Singularity image.  This image is too large to be hosted on Github, although you can find the definitions file used to create the image [here](https://github.com/pmonnahan/AncInf/blob/master/singularity/Singularity_defs.def).  Building of images is still not currently supported at MSI, so I used a Vagrant virtual machine, which comes with Singularity pre-configured/installed (https://app.vagrantup.com/singularityware/boxes/singularity-2.4/versions/2.4).  I can also share the img file directly upon request.

However, in order to utilize the singularity image, _Singularity_ must be installed on the HPC.  Currently, the pipeline assumes that _Singularity_ will be available as a module and can be loaded into the environment via the command specified in the config.yml file, where it says 'singularity_module'.  The default setting will work for MSI at UMN.

Singularity settings in config.yml

    singularity:
      use_singularity: 'true'
      image: '/home/pmonnaha/pmonnaha/singularity/AncestryInference.sif


## Running the workflow

Clone this repository to the location where you want to store the output of the pipeline.

    git clone https://github.com/pmonnahan/AncInf.git rfmix_test
    cd rfmix_test
    
The critical files responsible for executing the pipeline are contained in the *./workflow* subdirectory contained within the cloned repo.  They are: 

* Snakefile
* config.yml
* cluster.yaml  

The **Snakefile** is the primary workhouse of _snakemake_, which specifies the dependencies of various parts of the pipeline and coordinates their submission as jobs to the MSI cluster.  No modifications to the **Snakefile** are necessary.  

However, in order for the **Snakefile** to locate all of the necessary input and correctly submit jobs to the cluster, **both** the **config.yaml** and **cluster.yaml** need to be modified. Open these files and change the required entries that are indicated with 'MODIFY'.  Other fields do not require modification, although this may be desired given the particulars of the run you wish to implement.  Details on each entry in the config file (e.g. what the program expects in each entry as well as the purpose of the entry) are provided in the _Pipeline Overview_ at the bottom.

Once these files have been modified, the entire pipeline can be run from within the cloned folder via:

    snakemake --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 32

where -j specifies the number of jobs that can be submitted at once.  

The pipeline is currently set up to run on the _small_ queue on the _mesabi_ cluster, which has a per-user submission limit of 500 jobs.  This is more than enough for the entire pipeline, so running with -n 500 will submit all necessary jobs as soon as possible.  If -j is small (e.g. 32), snakemake will submit the first 32 jobs and then submit subsequent jobs as these first jobs finish.

The attractive feature of _snakemake_ is its ability to keep track of the progress and dependencies of the different stages of the pipeline.  Specifically, if an error is encountered or the pipeline otherwise stops before the final step, _snakemake_ can resume the pipeline where it left off, avoiding redundant computation for previously completed tasks.  

To run a specific part of the pipeline, do:

    snakemake -R <rule_name> --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 20 --rerun-incomplete

where _rule\_name_ indicates the 'rule' (i.e. job) in the Snakefile that you wish to run.

All output will be contained in the parent directory labelled with the prefix provided in the config file at:

    outname: "Example" 

## Pipeline Overview

![Pipeline DAG](https://github.com/pmonnahan/admixMap/blob/master/Pipeline_DAG.png)

### Input Data

Currently, the pipeline is set up only for binary case/control phenotypes, but this could be modified in the future.

    query: "PATH_TO_PLINK_PREFIX" 
    samples: "all" 
    
The query field in the config file should be set to the PLINK prefix for the single set of PLINK files containing all data.  It is assumed that sex will be specified in the .fam file.

The samples field can be set to the path of a file containing individuals to be kept from merged query file. One sample per line and must match exactly the sample names in the RFMix output and the query .fam file.  Beware that '#' characters in sample names will likely lead to errors.  

### Data Preparation
    
After optionally subsetting the data to the samples specified above, the data is run through a series of QC steps to prepare it for mapping as well as to create a separate dataset for calculating necessary supplementary information such as population structure, relatedness, etc.

#### Missingness filters

The first step is to filter out rare variants and variants with high missingness as well as for variants in which missingness is associated with case/control status.  Thresholds for these filters can be found under the 'gwas' section of the config file:

    gwas:
      mbc_pval: '0.0001' # pvalue threshold for missingness-by-case
      missingness: '0.2' # missingness threshold
      maf: '0.01' # Minor allele frequency threshold

#### LD pruning
We seek to create a separate dataset for calculating population structure and sample relatedness.  To do so, we create a smaller dataset of well-behaving, independent markers.  To identify such markers, we first subset the data to included only controls and remove variants exhibiting high linkage disequilibrium.

    LD_pruning:
      window_size: "50"
      step_size: "5"
      r2_threshold: "0.3"
      
#### IBD filters
We then remove highly related control individuals, so that we can identify sites with significant departures from Hardy-Weinberg equilibrium (HWE).  To crudely identify related individuals, we use PLINK's Identity by Descent ([IBD](https://www.cog-genomics.org/plink/1.9/ibd)) feature.  Given the large number of pairwise comparisons of individuals, we break this calculation into multiple parallel jobs, with the number of comparisons per job specified by the 'cpj' parameter in the config file (see below).  

    IBD:
      maf: '0.05' # Only include sites with minor allele frequency above this.
      cpj: '100000' # number of comparisons per job for IBS calculation
      pi_hat: '0.1875' # IBD threshold
      filter_relateds: 'false' # This should be left as false.  

Individuals are filtered based on the 'PI_HAT' value generated by PLINK, which approximately corresponds to IBD.  [Anderson et al 2010](https://www.nature.com/articles/nprot.2010.116) suggest 0.1875 as a threshold, but higher values (~0.3) are observed in the literature.  PLINK has a tendency to overestimate these values, so using a higher value, particularly in structured populations, is justifiable.  The final parameter 'filter_relateds' pertains to case samples, and should be generally left as false.  Both of the GWAS programs can accommodate related samples, negating the need to remove them.  However, this is a necessary step in admixture mapping, so it is enforced regardless of value set here.  Related control samples are also automatically filtered, regardless of this value, and a list of filtered control samples can be found here:

    accessory/excluded_control_samples.txt

Increase the 'pi_hat' entry in the config file if too many control samples are being excluded. 

#### Hardy-Weinberg equilibrium filters
For HWE calculations, we use the structured Hardy-Weinberg ([sHWE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6827367/)) formulation, which facilitates testing of HWE when a dataset contains multiple, genetically-structured populations.  However, there are two substantial limitations to this approach.  First, it requires the user to specify the number of populations in the dataset.  Since this is generally an unknown parameter, the authors suggest calculating p-values under a range of population numbers (determined by 'pop_nums' entry in the config file; see below) and then using the maximum entropy of the p-value distribution as criteria for finding the optimal population number.  The entropy is calculated over all bins of the p-value distribution excluding the left-most bin which contains the true SNPs violating HWE (the number of bins for the p-value distribution can be adjusted by the 'bins' entry in the config file, although the default setting should generally work well).  The remainder of this distribution should be uniform if the population parameter is accurate.  Entropy increases as the distribution becomes more uniform, so max entropy is used as the criteria (see linked paper for more thorough description).

The other main limitation to the program is that it does not allow any missing data at a marker.  Since markers with no missing data could represent a biased subset of the genotype data, the pipeline performs a large series of downsampling iterations (with user-specified parameters) with the hope that, by random chance, a subset of the markers with missing data in the full dataset will have no missing data in the downsampled dataset(s).  The parameters governing this process are specified by the last 3 lines of the sHWE section in the config file listed below.  Using the example below, the pipeline will generate 2000 downsampled datasets each consisting of 1000 randomly selected individuals.  Each of these datasets will be filtered to include only markers with complete genotype data.  The first downsampled dataset is used to determine the optimum population number as described above.  Then, the sHWE test is performed on all downsampled datasets based on this optimum population number.  

After testing all downsampled datasets, the pipeline goes through a checkpoint to determine if the number of distinct markers that have been tested exceed the value specified in the 'test_threshold' entry.  If so, the pipeline continues.  If not, the pipeline halts and instructs the user to rerun the pipeline to generate an additional set of downsampled datasets (with the number of downsampled datasets again specified by 'downsample_rounds').  This process can be repeated until 'test_threshold' is reached.

sHWE parameters in config file:

    sHWE:
      maf: '0.05' # Markers with minor allele frequency below this will be removed
      pop_nums: '1,2,3,4,5,6' # The range of population numbers over which we search for optimum value
      threshold: '0.000001' # p-value threshold for filtering based on departure from sHWE
      bins: '30' # Number of bins for the p-value distribution
      downsample_rounds: '2000' # Number of downsampled (DS) datasets to create
      downsample_number: '1000' # Number of individuals in each DS dataset
      test_threshold: '500000' # Number of markers to be tested.

Only markers that are tested for sHWE will be used for calculating PCs, IBD/IBS, relatedness, etc.  While a full input dataset might contain several million SNPs, only several hundred thousand SNPs are sufficient for making these calculations.  Therefore, the default 'test_threshold' is likely sufficient.

#### Control Matching  STOPPED HERE
Control matching is optionally performed for BLINK and admixture mapping.  Logistic models implemented in GENESIS are able to handle highly imbalanced case/control ratios by way of a Saddle-Point Approximation for p-values, so control matching is not necessary.  

Control matching behavior is set by the following lines in the config file:

    match_controls: 'false'
    control_number: '4'

If 'true', then PLINK's \href{http://www.cog-genomics.org/plink/1.9/strat}{{\color{blue}{\underline{control matching program}}}} was implemented while requested the following number of controls to be retained per case sample (via -mcc flag):
```{r, ctrl_num, warning=FALSE, message=FALSE}
yaml$control_number
```

The retained samples can be found here:
```{r, matched_samples, warning=FALSE, message=FALSE}
paste(getwd(), str_replace(params$fam_file, "_sub", "_CtlMat"), sep="/")
```


### Population structure and Kinship
    
### Genome Wide Association Mapping (GWAS)
GWAS is performed via two softwares: [BLINK](http://zzlab.net/blink) and [GENESIS](https://github.com/UW-GAC/GENESIS).  BLINK implements a novel multi-locus methodology that is fast, capable of handling large datasets, and claims greater statistical power.  However, the program is incapable of generating coefficient estimates necessary for calculating odds ratios, polygenic risk scores, etc.  GENESIS scales poorly with sample size, but is well-documented, written in R (and thus more transparent), and provides much more detailed output.

Association Mapping settings in config.yml:

    gwas:
      mbc_pval: '0.0001' # pvalue threshold for missingness-by-case
      missingness: '0.2' # missingness threshold
      maf: '0.01' # Minor allele frequency threshold
      pc_num: '4' # Number of PCs to use in GWAS
      cores: '24' # Number of cores to use in GWAS.  24 is recommended
      sig_threshold: '0.000005' # pvalue threshold for annotating SNPs
      other_predictors: 'sex' # Covariates to be used in GWAS.
      blink_dir: "/usr/local/bin/BLINK" # Do not change

See the text following '#' for a description of each entry.  For the 'other_predictors' entry, the 'sex' info is taken from the .fam file specified in the 'query' entry of the config file.  Additional covariates can be listed here, but they must be provided in a separate file whose path is provided at the following entry.

    covariate_file: 'none' 
    
 This file should be tab-delimited file with a header line specifying names of covariates (excluding sex) that are listed in the 'covars' entry.  First column MUST be labeled with 'taxa' and should contain sample IDs EXACTLY as they are specified in input plink files.  Note that this feature has not been tested.

### Admixture Mapping
For the admixture mapping, this pipeline takes, as input, the output of the [RFMix](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3738819/) ancestry inference program. The model used for admixture mapping is taken from [Grinde et al 2019](https://www.cell.com/ajhg/pdfExtended/S0002-9297(19)30008-4).  In short, for every genomic window a general linear model (generalized, in the case of case/control phenotypes) is fit with the local ancestry state as the predictor variable with global ancestry (or admixture proportions) as a covariate.  For every ancestry of interest, a separate model is fit for each ancestry listed in the config.yml file at .  Note that the acronyms 

Admixture Mapping settings in config.yml:

    admixMapping:
        skip: 'false' # Skip admixture mapping entirely
        rfmix: "OUTPUT_DIRECTORY_OF RFMIX" # MODIFY; can be left alone if 'skip' is set to 'false'
        gen_since_admixture: '8'
        cores: '12'
        
The 'gen_since_admixture' entry specifies the generations since admixture began parameter that was used in the ancestry inference procedure.  This is a required parameter in RFMix, but is also used to determine the appropriate significance threshold for admixture mapping (see below).  '8' is the value frequently used in the literature for inferring local ancestry in African-American populations.  It is recommended that one either consults the literature to find a suitable value for the particular admixture scenario represented in the data.  Or, the STEAM package (see below) can be used to empirically estimate this value from the data (although I have not tested this).  
        
The Grinde et al 2019 paper also provides a framework and software (called [STEAM](https://github.com/kegrinde/STEAM)) for determining an appropriate significance threshold for multiple test correction. It is based on the genetic positions of the tested windows, the average ancestry of each individual, and the number of generations since admixture began.  

