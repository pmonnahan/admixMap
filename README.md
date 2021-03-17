# Admixture and association mapping pipeline

This pipeline implements admixture mapping on the output of the [AncInf](https://github.com/pmonnahan/AncInf) pipeline (i.e. RFMix output) and/or association mapping on the output of the [pre-imputation](https://github.com/pmonnahan/DataPrep) or [post-imputation](https://github.com/pmonnahan/DataPrep/postImpute) QC pipelines.  It is not completely necessary to generate the input via the mentioned pipelines, although doing so would likely reduce the chances of encountering an issue.  Prior to mapping, input data is parsed and QC'ed to identify (and optionally filter) related samples, calculate principal components, and perform genomic control matching. A more detailed description of each of these steps is provided at the bottom. 

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

 - [Requirements](#requirements)
   - [Snakemake](#snakemake)
   - [Singularity](#singularity)
 - [Running the workflow](#running-the-workflow)
   - [Other Notes](#other-notes)
    - [Debugging and error reports](#debugging-and-error-reports)
 - [Pipeline Overview](#pipeline-overview)
   - [Input Data](#input-data)
   - [Data Preparation](#data-preparation)
     - [Missingness filters](#missingness-filters)
     - [LD pruning](#ld-pruning)
     - [IBD filters](#ibd-filters)
     - [Hardy-Weinberg equilibrium filters](#hardy-weinberg-equilibrium-filters)
     - [Control Matching](#control-matching)
   - [Population structure and Kinship](#population-structure-and-kinship)
   - [Genome Wide Association Mapping (GWAS)](#genome-wide-association-mapping-gwas)
   - [Admixture Mapping](#admixture-mapping)
   - [Output](#output)
   - [Conditional Analysis (Optional)](#conditional-analysis-optional)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Requirements
 
### Snakemake
The pipeline is coordinated and run on an HPC (or locally) using _Snakemake_.  To install snakemake, first create a virtual environment via:
  
    module load python3/3.6.3_anaconda5.0.1
    conda install -c conda-forge mamba
    mamba create -c conda-forge -c bioconda -n <your_environment_name> snakemake
  
This will create a new virtual environment and install `snakemake`.  Then, activate this environment and perform following installations:

    conda activate <your_environment_name>
    conda install numpy yaml pandas

Anytime you need to run the pipeline, activate this environment beforehand via:

    conda activate <environment_name>

If you choose not to create an environment, you must ensure that these packages are installed and available for your python installation.

### Singularity

The installation of the individual programs used throughout this pipeline can be completely avoided by utilizing a Singularity image.  This image is too large to be hosted on Github, although you can find the definitions file used to create the image [here](https://github.com/pmonnahan/AncInf/blob/master/singularity/Singularity_defs.def).  Building of images is still not currently supported at MSI, so I used a Vagrant virtual machine, which comes with Singularity pre-configured/installed (https://app.vagrantup.com/singularityware/boxes/singularity-2.4/versions/2.4).  I can also share the img file directly upon request.

However, in order to utilize the singularity image, _Singularity_ must be installed on the HPC.  Currently, the pipeline assumes that _Singularity_ will be available as a module and can be loaded into the environment via the command specified in the config.yml file, where it says 'singularity_module'.  The default setting will work for MSI at UMN.

Singularity settings in config.yml

    singularity:
      use_singularity: 'true'
      image: '/home/pmonnaha/pmonnaha/singularity/AncestryInference.sif
      

## Running the workflow

Clone this repository to the location where you want to store the output of the pipeline.

    git clone https://github.com/pmonnahan/AncInf.git rfmix_test
    cd rfmix_test
    
The critical files responsible for executing the pipeline are contained in the _./workflow_ subdirectory contained within the cloned repo.  They are: 

* Snakefile
* config.yml
* cluster.yaml  

The _Snakefile_ is the primary workhouse of snakemake, which specifies the dependencies of various parts of the pipeline and coordinates execution.  No modifications to the _Snakefile_ are necessary.  

In order for the _Snakefile_ to locate all of the necessary input and correctly submit jobs to the cluster, **both** the _config.yaml_ and _cluster.yaml_ need to be modified. Open these files and change the required entries that are indicated with 'MODIFY'.  Other fields do not require modification, although this may be desired given the particulars of the run you wish to implement.  Details on each entry in the config file (e.g. what the program expects in each entry as well as the purpose of the entry) are provided in the _Pipeline Overview_ at the bottom.

The entire pipeline can be executed on a local machine (not recommended) or on an HPC, and the _cluster.yaml_ file is required only for the latter.  For a local run, change the `local_run` entry to `true` under the `run_settings` section of the config file, and launch snakemake from within the parent directory by the simple command:

    snakemake

However, multiple steps in the pipeline have high resource demands, and so are unlikely to be able to be run locally.  This option exists primarily for testing and troubleshooting, so the remainder of the  documentation assumes that the pipeline will be executed on an HPC.  In order to coordinate the use of the HPC, the following modifications to the snakemake command are required:

    snakemake --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 32

where -j specifies the number of jobs that can be submitted at once.  Note that the 'qsub' command is specific to the commonly-used **PBS** scheduler.  To run on a different HPC scheduler, the command would need to be modified accordingly.  For example, to coordinate submission to a **slurm** scheduler, the following command would be used:

    snakemake --cluster "sbatch --no-requeue --time={cluster.time} --mem-per-cpu={cluster.mem-per-cpu} --ntasks={cluster.ntasks} --nodes={cluster.nodes} --mail-user={cluster.mail-user} --mail-type={cluster.mail-type} -o {cluster.o} -e {cluster.e} -A {cluster.A}"" --cluster-config workflow/cluster_yale.yaml -j 32

Note also that a different _cluster.yaml_ file is required for the different scheduler.  If you open and inspect the _cluster.yaml_ file vs the _cluster_yale.yaml_ file, you will see syntax that is specific to PBS and slurm schedulers, respectively.  

One additional change in the _config.yml_ is needed in order to correctly submit jobs to the HPC.  The relevant entries are under the `run_settings` section of the config file:

    run_settings:
      local_run: 'false'
      cluster_config: 'workflow/cluster.yaml'
      scheduler: 'pbs'
      
Here, it is necessary that the `cluster_config` entry is set to the path of the cluster.yaml file that will be used in the snakemake command.  Also, the scheduler must correspond to the syntax used in the snakemake command and cluster.yaml file.  I should point out that these additional changes are needed for responsibly using PLINK within a snakemake framework, and are not directly needed for snakemake.  PLINK will attempt to auto-detect available resources upon running regardless of the resources that were requested when the job was submitted.  Therefore, we have to read and parse the requested resources in the cluster config file in order for them to be communicated to PLINK from within the Snakefile.  

### Other notes

It is recommended that _snakemake_ is run as an interactive session on an HPC.  _Snakemake_ will launch the specified number (via the -j flag) of jobs, and then will hang and wait for them to finish.  As jobs finish (and assuming no errors), _snakemake_ will launch additional jobs keeping the total running jobs at whatever -j is set for.  Although _snakemake_ should not use a lot of memory, it could have long run times, which is generally not advisable on login nodes.  

One attractive feature of _snakemake_ is its ability to keep track of the progress and dependencies of the different stages of the pipeline.  Specifically, if an error is encountered or the pipeline otherwise stops before the final step, _snakemake_ can resume the pipeline where it left off, avoiding redundant computation for previously completed tasks.  To do so, simply resubmit the original _snakemake_ command.

To run a specific part of the pipeline, do:

    snakemake -R <rule_name> --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 20 --rerun-incomplete

where _rule\_name_ indicates the 'rule' (i.e. job) in the Snakefile that you wish to run.

Also, it is often very helpful to do a 'dry-run' of the pipeline in which the different steps and dependencies are printed to screen, but no actual jobs are executed.  This can be helpful to ensure that config entries are correct, etc.  To perform a dry-run, do:

    snakemake -nrp
    
#### Debugging and error reports

Should an error be encountered in a job, snakemake will halt the pipeline and indicate in the terminal that an error has occurred.  The offending job will also be printed in red in the terminal window.  More information on why the job failed can be found in the 'stdout' and 'stderr' files that are output to the _'OandE'_ directory and will be labelled with the jobname.

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
      
See the text following '#' for a description of each entry.  

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

#### Control Matching  
Control matching is optionally performed for admixture mapping.  Logistic models implemented in GENESIS are able to handle highly imbalanced case/control ratios by way of a Saddle-Point Approximation for p-values, so control matching is not necessary.  

Control matching behavior is set by the following lines in the config file:

    match_controls: 'false'
    control_number: '4'

If 'match_controls' is set to 'true', then PLINK's [control matching](http://www.cog-genomics.org/plink/1.9/strat) feature, will be implemented such that a number of control samples will be retained for each case sample (via -mcc flag).  The number of retained controls (per case sample) is specified with the 'control_number' entry in the config file.

The retained samples can be found in the 'accessory' directory that is created at:

    accessory/samples.txt

### Population structure and Kinship

Using the sHWE-passing SNPs explained above, a first-pass kinship matrix is calculated via [KING](http://people.virginia.edu/~wc9c/KING/manual.html), which is designed to be robust in the face of population structure.  Then, these kinship estimates along with the SNP genotypes are used to calculate ancestry-informative PCs via [PC-AiR](https://rdrr.io/bioc/GENESIS/man/pcair.html). The basic idea behind PC-AiR is to calculate PCs using only unrelated individuals and then use the resulting eigenvectors/values to project PC values for the related individuals that had been excluded.  The (un)related individuals identified by PC-AiR are output to files ending in '.unrels' and '.rels' in the main directory.
    
A genetic-relatedness matrix (GRM) is then calculated via [PCRelate](https://www.rdocumentation.org/packages/GENESIS/versions/2.2.2/topics/pcrelate), which 'accounts for population structure (ancestry) among sample individuals through the use of ancestry representative principal components (PCs)'. 

Additional information on the PCs and GRM can be found in the PDF report that is automatically generated at the end of the pipeline.
    
### Genome Wide Association Mapping (GWAS)
GWAS is performed via [GENESIS](https://github.com/UW-GAC/GENESIS).  GENESIS scales poorly with sample size, but is well-documented, written in R (and thus more transparent), and provides much more detailed output.

Association Mapping settings in config.yml:

    gwas:
      mbc_pval: '0.0001' # pvalue threshold for missingness-by-case
      missingness: '0.2' # missingness threshold
      maf: '0.01' # Minor allele frequency threshold
      pc_num: '4' # Number of PCs to use in GWAS
      cores: '4' # Number of cores to use in GWAS.
      sig_threshold: '0.000005' # pvalue threshold for annotating SNPs
      other_predictors: 'sex' # Covariates to be used in GWAS.

For the 'other_predictors' entry, the 'sex' info is taken from the .fam file that is associated with the 'query' entry of the config file.  Additional covariates can be listed here, but they must be provided in a separate file whose path is provided at the following entry.

    covariate_file: 'none' 
    
 This file should be tab-delimited file with a header line specifying names of covariates (excluding sex) that are listed in the 'covars' entry.  First column MUST be labeled with 'taxa' and should contain sample IDs EXACTLY as they are specified in input plink files.  Note that this feature has not been tested.
 
 NOTE: while increasing the number of cores (via the `cores` entry) may speed up computation, doing so WITHOUT requesting additional memory (in the cluster.yaml file) may result in stalled jobs.  In GENESIS, computation is parallelized by chromosome.  Tests with a dataset containing several thousands of individuals found that GENESIS stalls with 12 cores and 150Gb of memory. 

### Admixture Mapping
For the admixture mapping, this pipeline takes, as input, the output of the [RFMix](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3738819/) ancestry inference program. The model used for admixture mapping is taken from [Grinde et al 2019](https://www.cell.com/ajhg/pdfExtended/S0002-9297(19)30008-4).  In short, for every genomic window a general linear model (generalized, in the case of case/control phenotypes) is fit with the local ancestry state as the predictor variable with global ancestry (or admixture proportions) as a covariate.  For every ancestry of interest, a separate model is fit for each ancestry listed in the config.yml file at 'ancestry_predictors'.  Note that the acronyms used here must exactly match those specified in the RFMix output.  

Admixture Mapping settings in config.yml:

    admixMapping:
        skip: 'false' # Skip admixture mapping entirely
        rfmix: "OUTPUT_DIRECTORY_OF RFMIX" # MODIFY; can be left alone if 'skip' is set to 'false'
        ancestry_predictors: 'AFR,EUR,AMR' # No spaces!
        other_predictors: 'sex'
        gen_since_admixture: '8'
        cores: '12'
        
The 'gen_since_admixture' entry specifies the generations since admixture began parameter that was used in the ancestry inference procedure.  This is a required parameter in RFMix, but is also used to determine the appropriate significance threshold for admixture mapping (see below).  '8' is the value frequently used in the literature for inferring local ancestry in African-American populations.  It is recommended that one either consults the literature to find a suitable value for the particular admixture scenario represented in the data.  Or, the STEAM package (see below) can be used to empirically estimate this value from the data (although I have not tested this).  
        
The Grinde et al 2019 paper also provides a framework and software (called [STEAM](https://github.com/kegrinde/STEAM)) for determining an appropriate significance threshold for multiple test correction. It is based on the genetic positions of the tested windows, the average ancestry of each individual, and the number of generations since admixture began.  

### Output

All output will be contained in the parent directory labelled with the prefix provided in the config file at:

    outname: "Example" 
 
 The table below provides a list of results and the corresponding file suffix.
    
| Result  | Suffix |
| ------------- | ------------- |
| GENESIS GWAS  | .genesis.txt  |
| Admixture Mapping  | .admixmap.txt  |
| PC-AiR Variance Proportions  | .pcair.varprop  |
| Phenotype Data  | .genesis.pdat  |
| PC-AiR Related Samples  | .pcair.rels  |
| PC-AiR Unrelated Samples  | .pcair.unrels  |

Additionally, a report summarizing much of the intermediate and final results is produced at the end of the pipeline and will end with the suffix 'Mapping-report.pdf'.  Figures in this report are also saved in the 'figures' directory that is created.

Lastly, SNPs whose p-value was below that specified at the 'sig_threshold' entry under 'gwas' are annotated via [SnpEff](http://snpeff.sourceforge.net/SnpEff_manual.html).  The annotated files end in 'sig.annt.txt'.  The SNPs are also annotated with the rsID using the key that is provided in the config file at:

    annotation:
        rsIDs: "/home/spectorl/pmonnaha/misc/hg19_rsID_key.txt"

This file should be a PLINK-formatted bim file containing the rsIDs.  I have already generated this file for hg19 rsIDs and can share this file upon request.

### Conditional Analysis (Optional)

An optional final step is to perform conditional analyses to locally search for additional significant variants after controlling for the effects of a particular marker.  

    conditional_analysis:
        snp_list: 'none'
        distance: '100000'

The 'snp_list' entry should specify a file path containing one marker ID (chr:position) per line.  One possibility is to use the 'significant' markers that were annotated at the final step.  For each SNP in the file, the program will perform a GWAS for all SNPs within the specified 'distance' on either side.  Importantly, the genotypes at the focal SNP are included as an additional covariate.

In order to run this step you must explicitly request the output via:

    snakemake conditional-analysis/<outname>.conditional-analysis.genesis.txt --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 1

 , where 'outname' should be replaced by the value in the 'outname' entry in the config file.  
