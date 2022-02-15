
# phASE-Extender

**Extender** for the readbackphased haplotype blocks.
***A python program to extend the ReadBackPhased haplotype blocks using markov first order transition probabilities and likelihood test.***

Developed by [Bishwa K. Giri](mailto:kirannbishwa01@gmail.com) in the [Remington Lab](https://biology.uncg.edu/people/david-remington/) at the University of North Carolina at Greensboro, Biology department.

## Citation

Giri, B. K., Remington D. L. PhaseIT - A haplotype phasing too for heterogenous and hybrid genomes and for emerging genomic models using phase-Extender and phase-Stitcher. biorxiv (2018) [not uploaded yet].

## AUTHOR/SUPPORT

Bishwa K. Giri (bkgiri@uncg.edu; kirannbishwa01@gmail.com) \
Support @ <https://groups.google.com/d/forum/phase-extender>

- [phASE-Extender](#phase-extender)
  - [Citation](#citation)
  - [AUTHOR/SUPPORT](#authorsupport)
  - [Intro to ReadBackPhasing](#intro-to-readbackphasing)
- [BACKGROUND](#background)
  - [Summary](#summary)
  - [Algorithm](#algorithm)
  - [Benefits of using phaseExtender](#benefits-of-using-phaseextender)
  - [!#f03c15 Data Requirements](#f03c15-data-requirements)
- [Tutorial](#tutorial)
  - [Prerequisites](#prerequisites)
  - [Installation and setup](#installation-and-setup)
  - [Usage:](#usage)

## Intro to ReadBackPhasing

Two heterozygote genotypes in a diploid organism genome are called to be readback phased if they are supported by aligned reads sequence. Depending upon the size and type (single end vs. paried end) of read sequence library, type of read sequence (DNAseq vs. RNAseq) readbackphased haplotype can range from the size of 2 genotypes to multiple genotypes.

**Check these links for more details on readbackphasing**

- <https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_phasing_ReadBackedPhasing.php>
- <https://github.com/secastel/phaser/tree/master/phaser>



# BACKGROUND

Haplotype phasing is a second "go to" problem in bioinformatics after read alignment. The importance of haplotype phasing applies directly to the analyses of ASE (allele specific expression), preparation of extended haplotype for EHH (extended haplotype homozygosity) test, and preparation of dipolid genome which will soon be a new standard in bioinformatics in coming years, etc. The necessity for haplotype phasing (and eventually diploid genome) increases with the increase in heterozygosity in the genome because higher hetetogeneity leads to larger alignment bias and complicates the reliability of the variants that are called using that alignment data (SAM, BAM files).

The classical approach to haplotype phasing involves application of LD (linkage disequilibrium) test between two heterozygous markers, that began with the preparation of genetic map by Alfred Sturtevant. Any existing population based haplotype phasing tool therefore uses the LD test with varying degree of complexity based on available sample size, markers, types of the marker and relationship between samples. For haplotype phasing tools like [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html), [ShapeIT](http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html), [impute2](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html), etc are used dominantly. These tools make use of the variants (SNPs, InDels) along the length of the genome by treating chain of variants along genomic positions singly and independently. So, for a diploid organism that contains chain of variants with "n" heterozygote sites there exists **"2<sup>n</sup>"** possible haplotypes. To certain extent this **"2<sup>n</sup>"** problem is approached and eased by sampling genotype data in short genomic intervals from several samples and applying identity by descent (IBD), most common haplotype method, etc. on the sampled genotypes to infer the possible haplotype in that region. However, application of these methods may not be optimal in solving phase states in organisms that have highly heterogenous genome, are hybrids and/or have very few or no reference panels and abundant genotype data.

**The main issues of the tools that deal with haplotype phasing in**"2<sup>n</sup>"**way can be summarized as:**

- Increased computation burden due to **"2<sup>n</sup>"** problem.
- Problem in phasing rare variants.
- Is mostly applicable to human and organisms with reference haplotype panel and abundant genomic data.
- Not optimal for organisms that have heterogenous genome, or are outbreds, or are hybrids, or belong to the group of organisms that have small sample of genomic resources and reference genome prepared.  
  
Therefore, using LD between two adjacent SNPs using small population of samples isn't able to provide enough resolution to solve GW (genome wide) haplotype preparation, which could result in excessive switch errors. Additionally, in heterogenous and hybrid genome, the problems arising due to switch errors in downstream analyses could be manifold.

ReadBackphasing is emerging as a new and more reliable method for preparation of short range haplotypes by joining heterozygous variants that are covered by sequence reads. These short haplotypes can be further elongated by adding sequence reads of longer length (PacBio reads) or by adding more genomic and RNAseq reads from the same individual; see tools [WhatsHap](https://whatshap.readthedocs.io/en/latest/), [hapCut](https://github.com/vibansal/hapcut), [phaser](https://github.com/secastel/phaser/blob/master/phaser/README.md) etc.  

**But, issues with existing RBphase method and tools still remain and can be summmarized as:**

   1. are mainly aimed at preparing long range haplotypes but not necessarily genome wide.
   2. existing RBphase methods are only targeted at individual or family level, i.e they require multiple inputs of "BAM" and "VCF" files for the same individual and/or trios. The increase in the size of phased haplotype blocks only depends upon multiple BAM files, or multiple sets of longer reads from the same individual, which still means additional sequencing and cost is involved.
   3. Referring to point 2 -> integration of RBPhased data with population based phasing is still missing i.e they are not able to solve phase state of two haplotype blocks in same sample by using information of the phase state of the haplotype from other samples.

With increase in the size of PE (paired end) reads from illumina and availability of longer sequences from [PacBio](https://www.pacb.com/smrt-science/smrt-sequencing/read-lengths/), it is now possible to increase the size of RBphased haplotypes considerably. Inspite of the the longer reads in the future there will always be coverage issues due ro random coverage gaps. This gaps and will break genome wide haploltype into multiple haplotype segments, thus complete RBphasing is also not an optimal solution.

## Summary

**Combination of RBphasing with population based phase extenstion reduces problem of**

- "2<sup>n</sup>" computation burden.
- requirement of large genotype samples and reference panels.

by :

  1. First preparing RBphased haplotype blocks within a sample using aligned sequence reads (BAM, SAM files).
  2. Next, the two consecutive haplotypes blocks at a break point in a sample are joined by computing the likelihood estimates of the LD (linkage disequilibrium) observed in other samples at that break point.

Since, the size of RBphased blocks increases with the increase in heterozygosity in the genome, **phaseExtender** is a tools highly suitable for the organisms with heterogenous genome and/or which have limited amount of genotype data sequenced. RBphased haplotype always has more information compared to a single SNP or InDel thereby overcoming the issues with switcherrors when preparing the long range haplotype in heterozygous population.
So, readbackphasing combined with population based phasing is able to yield higher variants per haplotype block, making **phaseExtender** a better method and tool when working with organisms with higher heterozygosity (out crossing population and hybrids).

## Algorithm

- **phASE-Extender** uses RBphased haplotype data of several individuals that belong to same family, population or species. The readbackphased VCF for a single sample can be prepared using applications like [phaser](https://github.com/secastel/phaser/tree/master/phaser), [hapcut2](https://github.com/vibansal/HapCUT2), [GATK readbackphasing](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_phasing_ReadBackedPhasing.php). Several single sample RBphased VCFs are then merged to create multisample VCF.
- The RBphased haplotype data from multisample VCF is then converted to simple tabular format (HAPLOTYPE file). The genotype in this tabular format is represented as IUPAC base and haplotype blocks using unique block index.
- Next, the two consecutive RBphased blocks in a single sample can be joined either in parallel or alternate configuration, see figure **??**.
- To join the two consecutive RBphased blocks in a single sample we use haplotype state information of other samples in the pool. The likelihood of a posssible configuration (parallel vs. alternate) is estimated as LD observed between the consecutive blocks in the other samples.
- Likelihood estimates are computed by establishing a first order markov chain between the nucleotides of two consecutive blocks. Markov chains are represented as first order transition matrix from all the nucleotides of former haplotype block to all the nucelotides in the later haplotype block and then vice versa (for the reverse chain). The observed nucleotide emission probablity and transition probablity are cumulated into maxed-sum or maxed-product value to form a meaningful likelyhood estimates for both the configuration (parallel vs. alternate).
- The phase state of two consecutive haplotype blocks is then extended if the `computed log2 (likelihood)` of either configuration is above the `threshold log2(likelihood) cutoff`.

For the **mcve** regarding the algorithm see this issue on [**stackoverflow**](https://codereview.stackexchange.com/questions/186396/solve-the-phase-state-between-two-haplotype-blocks-using-markov-transition-proba) and/or [**my blog**](https://everestialblog.wordpress.com/2018/03/27/extension-of-the-readbackphased-haplotypes-using-phase-extender/).

## Benefits of using phaseExtender

- Combination of RBphased data with populaiton based phasing therefore allows us to use small sample size to accurately predict the proper haplotype state.
- PhaseExtender provides flexibility to adjust the number of phased genotype used for building markov chains between consecutive blocks.
- PhaseExtender provide flexibility of limiting phase extension to certain bed regions.
- Ability to adjust LOD cutoff along with above discussed customization provides a means to recursively improve haplotype phasing.

## ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) Data Requirements

**phASE-Extender** can be used with the multi-sample vcf files produced by GATK pipeline or other tools that generate readbackphased haplotype blocks in the output VCF. A haplotype file is created using the RBphased VCF and then piped into **phase-Extender**. See, this example for data structure of input haplotype file [sample input haplotype file01](https://github.com/everestial/phase-Extender/blob/master/example01/input_haplotype_small.txt) - a tab separated text file with `PI` and `PG_al` value for each samples.

In several tutorial examples (test files below) I have used hapotype file prepared from RBphased VCF generated by `phaser` (<https://github.com/secastel/phaser> , <https://github.com/secastel/phaser/tree/master/phaser>). However, **phase-Extender** can be used with input haplotype file prepared from any RBphased VCF given it meets the appropriate data structure. Readbackphased VCF can be converted into `haplotype file` using an add-on tool `vcf_to_table-v3.py`.\
VCF from `phaser` consists of `Phased Genotype i.e PG` and `Phase Block Index i.e PI` values in `FORMAT` field; **PI** represents the `unique haplotype block index` and **PG** represents the `phased genotypes within that PI block`.  After converting the `RBphased VCF` to `haplotype file` **phase-Extender** uses the `Phased Genotype i.e PG` and `Phase Block Index i.e PI` values in the `haplotype file` to prepare the transition matrix probabilities and proceed with phase extension.



# Tutorial

### Prerequisites

**phase-Extender** is written in python3, so you need to have python3 installed on your system to run this code locally. If you don't have python installed then, you can install from [here](https://www.python.org/downloads/). For linux; you can get latest python3 by:

`sudo apt-get install python3`

### Installation  and setup

1. Clone this repo.

```
git clone https://github.com/everestial/phase-Extender
cd phase-Extender
```

2. Make virtual env for python and install requirements.

```
python3 -m venv myenv
source myenv/bin/activate   # for linux
myenv\Scripts\activate      # for windows
pip install -r requirements.txt
```

OR, you can install latest versions individually by:

```
pip install pandas numpy matplotlib

```

3. To run tests locally:

  ``` pip install pytest
      pytest . 
   ```

4. Run help on phase extender

```
   python3 phase_extender.py -h  
```

<pre>
Checking and importing required modules:

#######################################################################
        Welcome to phase-extender version 1
  Author: kiranNbishwa (bkgiri@uncg.edu, kirannbishwa01@gmail.com)
#######################################################################

Loading the argument variables ....
usage: phase_extender.py [-h] [--nt NT] --input INPUT --SOI SOI
                         [--output OUTPUT] [--refHap REFHAP]
                         [--useSample USESAMPLE] [--bed BED] [--snpTh SNPTH]
                         [--numHets NUMHETS] [--lods LODS] [--culLH CULLH]
                         [--writeLOD WRITELOD] [--hapStats HAPSTATS]
                         [--addMissingSites ADDMISSINGSITES]

optional arguments:
  -h, --help            show this help message and exit
  --nt NT               number of process to run -> the maximum number of
                        processes that can be run at once is the number of
                        different chromosomes (contigs) in the input haplotype
                        file.
  --input INPUT         name of the input haplotype file -> This haplotype
                        file should contain unique index represented by 'PI'
                        and phased genotype represented by 'PG_al' for several
                        samples.
  --SOI SOI             sample of interest -> Name of the sample intended for
                        haplotype extension.
  --output OUTPUT       Name of the output directory. default: 'SOI_extended'
  --refHap REFHAP       reference haplotype panel -> This file should also
                        contain 'PI' and 'PG_al' values for each sample in
                        that haplotype reference panel.default: empty
  --useSample USESAMPLE
                        list of samples -> Use phase state information only
                        from these samples while running. This is intended to
                        provide control over phase-extension by allowing
                        select samples from the pool of samples (refHap and/or
                        input). This is useful for comparing the results when
                        different individuals, populations are used in phase
                        extension process.Options:
                        'all','refHap','input','comma separated name of
                        samples'. default: 'all' - i.e all the samples in
                        (refHap + input) will be used, if refHap is given else
                        all samples only from input file is used.
  --bed BED             bed file -> Process the haplotype extension only in
                        this bed regions. This is useful to limit haplotype
                        extension only within certain regions, like - within
                        genes, exons, introns, QTL, etc. default: empty
  --snpTh SNPTH         snp threshold -> Minimum number of SNPs required in
                        both haplotype blocks before starting phase extension.
                        Option: an integer value. Default: snpTh = 3
  --numHets NUMHETS     number of heterozygote SNPs -> Maximum number of
                        heterozygote SNPs used in consecutive haplotype blocks
                        for computing likelyhood of the phase states. Option:
                        an integer value. Default: numHets = 40
  --lods LODS           log2 of Odds cut off -> Cutoff threshold to assign the
                        phase states between consecutive haplotype blocks.
                        Option: an integer value Default: 5 for parallel
                        configuration (i.e 2^5 = 32 times more likely).
  --culLH CULLH         Cumulative likelhood estimates -> The likelhoods for
                        two possible configuration can either be max-sum vs.
                        max-product. Default: maxPd i.e max-product. Options:
                        'maxPd' or 'maxSum'.
  --writeLOD WRITELOD   write log2 of Odds to the output file-> writes the
                        calculated LODs between two consecutive haplotype
                        blocks in the output file. Option: 'yes', 'no'.
                        Default: no. **Note: the 'lods-score' are printed
                        regardless if the consecutive blocks are joined or
                        not.
  --hapStats HAPSTATS   Computes the descriptive statistics and plots
                        histogram of the haplotype for input and extended
                        haplotype. Default: 'no'.Option: 'yes', 'no' .
  --addMissingSites ADDMISSINGSITES
                        write the lines that have data missing for SOI on the
                        output file. Option: yes, no
   </pre>



## Usage

1. First example:

`
python3 phase_extender.py --input example01/input_haplotype_file.txt --SOI ms02g --lods 10`
<pre>
Checking and importing required modules:

#######################################################################
        Welcome to phase-extender version 1
  Author: kiranNbishwa (bkgiri@uncg.edu, kirannbishwa01@gmail.com)
#######################################################################

Loading the argument variables ....
Assigning values to the global variables ....
  - sample of interest: "ms02g"
  - using "1" processes
  - using haplotype file "example01/input_haplotype_file.txt"
  - using log2 odds cut off of "10.0"
  - each consecutive haplotype block should have minimum of "3" SNPs
  - using maximum of "40" heterozygote sites in each consecutive blocks to compute transition probabilities
  - using "max product" to estimate the cumulative maximum likelyhood of each haplotype configuration between two consecutive blocks
  - no bed file is given.
  - no reference haplotype panel is provided
  - statistics of the haplotype before and after extension will not be prepared for the sample of interest i.e "ms02g".     Only extendent haplotype block will be prepared.
  - LOD (log 2 of odds) for consecutive block will not be written to the output file

# Reading the input haplotype file "example01/input_haplotype_file.txt"
  - Lines that have data missing for sample "ms02g" is written in the file "ms02g_extended/missingdata_ms02g.txt"

# Genomic bed file is not provided ...
  - So, phase extension will run throughout the genome.

# Haplotype reference panel is not provided ...
  So, phase extension will run using the samples available in the input haplotype file.

# Filtered the lines that have data missing for sample "ms02g"; check the file "ms02g_extended/missingdata_ms02g.txt"
  - Loaded read-backphased variants onto the memory

# Haplotype reference panel is not provided....
  - Only using the samples in the input ("example01/input_haplotype_file.txt") data.

# No bed file is given ...
  - So, grouping the haplotype file only by chromosome (contig)

# Writing initial haplotype for sample "ms02g" in the file "initial_haplotype_ms02g.txt"
  - Proceeding to phase-extension without preparing descriptive statistics of initial haplotype state.

# Starting multiprocessing using "1" processes

## Extending haplotype blocks in chromosome (contig) 2
  - Grouping the dataframe using unique "PI - phased index" values.
  - Starting MarkovChains for contig 2
  - Phase-extension completed for contig "2" in 0.017184734344482422 seconds
  - Worker maximum memory usage: 58.57 (mb)

Completed haplotype extension for all the chromosomes.time elapsed: '0.11295127868652344'
Global maximum memory usage: 79.88 (mb)
Merging dataframes together .....

Extended haplotype data for sample "ms02g" is written in the file "extended_haplotype_ms02g.txt".
Skipping the preparation of descriptive statistics of extended haplotype.

Run is complete for all the chromosomes (contigs)

writing singletons and missing sites to extended haplotype
End :)
</pre>

## Usage and Inputs

- Requires a multisample readbackphased `haplotype file` as input and returns a single sample extended haplotype file. Other results files containing statistics on the initial vs. extended haplotype are also produced.
- Optionally, haplotype reference panel (with same data structure as input haplotype) and bed file can be included to limit or improve the process of phase extension.

?? needs improvement.  
Check this detailed [step by step tutorial](https://github.com/everestial/phase-Extender/wiki/phase-Extender-Tutorial) for preparation of `input files` and know-how about running `phase-Extender`.
??



## Inputs

***haplotype file (required):*** Input haplotype file. Should contain `PI` and `PG_al` values for each sample.\
To convert the vcf file to haplotype file (from VCF to tabular format) use ??**Step 01 (a)** in the tutorial. \
The sample name should not contain "`_`" character.

***haplotype reference panel (optional):*** Unlike "haplotype reference panel" used in other phasing tools, phaseExtender requires reference panel in the same structure as a HAPLOTYPE file.
To convert the haplotype reference panel (from VCF to proper text format) use **Step 01 (b)** in the tutorial. ??
  
***bed file (optional):*** If you goal is to limit phase extension to certain genomic regions (for eg. gene, exon or QTL boundries), we suggest that you provide appropriate bed file. **phase-Extender** then exclusively limits phasing within internal boundries of the input bed regions.

    # structure of the bed file
    contig    start    end       # this header is not included    
    2         1258     199897
    2         397765   412569

    #To convert the GTF,GFF to bed file use:
    python3 gffToBed.py  --input myGTF.gtf  --output myBed.bed

**the required python file for bed file preparation is at:**
<https://github.com/everestial/SmallTools/tree/master/GtfToTable>
<https://github.com/melissacline/TCGA-GAF-source/blob/master/scripts/gffToBed.py>

Continue .............

## Arguments

### Required

* **--input** - input haplotype file. `PI` and `PG_al` should be present in header for each sample.
- **--SOI** - sample of interest. It should refer to a single sample in the haplotype the file. The sample name should not contain "`_`" character.

### Performance Related

* **--nt** *(1)* - maximum number of processes to run at once. The maximum number of processes is limited to number of chromosomes (contigs) in the input haplotype file.

### Optional

* **--python_string** *(python3)* - Calls `python 3` interpreter to run the program.
- **--output** *(SOI_extended)* - Output directory.
- **--snpTh** *(3)* - snp threshold. Minimum number of SNPs required in each consecutive haplotype block to run phase extension between two blocks.
- **--numHets** *(40)* - num of heterozygotes. Maximum number of heterozygote SNPs used from each consecutive block to compute maximum likelihood estimate of each configuration between two blocks.
- **--culLH** *(maxPd)* - cumulation of the likelihood estimates. The likelhoods for two possible configuration can either be "maxed as sum" or "maxed as product". ***Default*** is "max-product". ***Options:*** 'maxPd' or 'maxSum'.
- **--lods** *(5)* - log2 of Odds cut off threshold. The cutoff threshold used to extend consecutive haplotype blocks. **`**Note: Default value is set at (2^5 = 32 times likely). So, two consecutive blocks will be joined in parallel configuration if computed log2(likelihood) > lods threshold **
- **--useSample** *(all)* - Samples to use in the given input haplotype file (plus reference haplotype) to compute transition matrix. Options: 'all','refHap','input','comma separated name of samples'. Default: all the samples in (refHap + input) will be used.
- **--bed** - Process the haplotype extension only within this bed regions. ***This is useful if you want to limit haplotype extension only within certain regions, like - within genes, exons, introns, QTL boundries, etc.***
- **--writeLOD** *(no)* - writes the calculate LODs between two consecutive haplotype blocks when processing phase extension to the output file. **Options:** 'yes', 'no'. **`**Note: the 'lods-score' are printed regardless if the "
"consecutive blocks are joined or not.**
- **--hapStats** *(no)* - Prepare descriptive statistics, and histogram of the haplotype size distribution of the input haplotype file vs. extended haplotype for the sample of interest. **Options:** 'yes', 'no'
- **--addMissingSites** *(no)* - include the non-phased and missing genotype data from input haplotype file to the final phase extended output file. **Option:** 'yes', 'no'.




## Output Files

### initial_haplotype_*SOI*.txt & extended_haplotype_*SOI*.txt

Contains all RBphased haplotype data for the sample of interest before and after phase extension.
- 1 - **contig** - Contig name (or number).
- 2 - **pos** - Start position of haplotype (1 based).
- 3 - **ref** - Reference allele at that site.
- 4 - **all-alleles** - All the alleles represented by all the samples at that site.
- 5 - **SOI_PI** - Unique `PI` index of the haplotype blocks for sample of interest.
- 6 - **SOI_PG_al** - Phased GT (genotype) alleles at the genomic position that belong to unique `PI` indexes.
- 7 - **log2Odds** (only in **extended_haplotype_SOI.txt**) - log2Odds computed between the former and later block.



### initial_haplotype_stats_*SOI*.txt & final_haplotype_stats_*SOI*.txt

Descriptive haplotype statistics of the input haplotype file for the sample of interest.
- 1 - **contig** - Contig name (or number).
- 2 - **SOI_PI** - Comma separated list of unique `PI` index of the haplotype blocks for the sample of interest. The total number of `PI` index represents the total number of haplotype fragments that are present in the given contig in that sample.
- 3 - **num_Vars_by_PI** - Number of variants sites within each `PI` block for the sample of interest.
- 4 - **range_of_PI** - Genomic range of the each `PI` block for the sample of interest.
- 5 - **total_haplotypes** - Total number of haplotype (i.e `PI`) in the given coting for the sample of interest.
- 6 - **total_Vars** - Total number of variant sites in the given contig for the sample of interest. **Note:** The sum of (num_Vars_by_PI) = total_Vars.

**Note:** - The `SOI_PI`, and it's associated statistics are in order.



### missingdata_*SOI*.txt

Contains data from the sites that have unphased or missing GT (genotype) for the sample of interest in the input haplotype file.
**Note:** This data is merged with `extended_haplotype_SOI.txt` if `--addMissingSites` is set to "yes".



### extended_haplotype_"SOI_"allsites.txt

This file contains ReadBackPhased haplotype after phase extension concated with the missing data. This file contain equal number of row as input haplotype file and data only for sample of interest.




## Plots

**Note:** - These plots are based on the descriptive statistics generated for haplotypes before and after phase extension. It is possible to take these statistics (initial_haplotype_stats_*SOI*.txt & final_haplotype_stats_*SOI*.txt) and make custom plots in **R** or by using other methods.



### total_haps_"SOI_"initial.png  & total_haps_"SOI_"final.png

Number of haplotypes for given contig before and after phase extension.



### total_vars_"SOI_"initial.png  & total_vars_"SOI_"final.png

Number of variants for given contig before and after phase extension.



### hap_size_byVar_"SOI_"initial.png & hap_size_byVar_"SOI_"final.png

Histogram of the distribution of the haplotype size (by number of variants in the haplotype) in given contig before and after phase extension.



### hap_size_byGenomicRange_"SOI_"initial.png & hap_size_byGenomicRange_"SOI_"final.png

Histogram of the distribution of the haplotype size (by genomic range of the haplotype) in given contig before and after phase extension.




## ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)Some Q/A on phase-extender

  ***1) What kind of algorithm does phase-extender use ?***  
  phase-extender uses first-order-transition probabilities from each level of genotypes from former haplotype block to each level of genotypes to later haplotype block. This version (v1) uses **forward-1stOrder-markov chains** and **backward-1stOrder-markov chains** transition probabilities. Future versions will follow improvements by adding markov-chains of higher order.

  ***2) What is the advantage of using phase-extender ?***  
  We generally need accurate phase state with in a gene/transcript level while doing ASE, dissecting maternal-paternal effects. Long haplotypes are mostly important while preparing diploid genome, testing selective sweeps within a QTL regions etc. For emerging organism systems where genotype data are sparse CW (chromosome wide), GW (genome wide) haplotypes are more difficult to solve. Also, haplotype phasing may be more complicated in out crossing individuals and hybrids due to heterogenity.\ RBphasing actually provides an advantage with heterogenous genome because frequency of RBphased blocks increases with heterogenity in the genome. These short haplotype fragments have multiple heterozygous variants on a short haplotype block. With increase in the size of `PE` (paired end) reads the size of RBphase blocks are also increases. phase-extender comes handy at this stage; that it tries to solve phase state of two consecutive blocks from one sample at a time by using data from haplotype blocks of other samples that bridge that breakpoint. So, we can solve haplotype configuration for SOI (sample of interest) with more confidence because:
     - we have more variants within each blocks contributing to more information.
     - we only need to solve two possible phase state at one time compared to 2^n haplotype when reading one SNP at a time.
    phase-Extender also provides a more flexible and manipulative control over how to proceed with phase extension. It is also possible to control several parameters like `lods`, `snpTh`, `numHets`, `culLH`, `bed`, `useSample` to observe and compare how phase extension changes.

  ***3) Does phase-extender phase InDels ?***  
    Yes, but it is conditional. The InDels should already be reliably readbackphased to a haplotype block. That way when the haplotype is being extended for those SNPs, InDels hitchhike with it and get extended too.

  ***4) What is the minimal size of the haplotype block that is required?***  
    The bigger the two haplotypes are, the better is the likelyhood test of which haplotype is phased with which. By, default I have kept this number to 3 variants (SNPs exclusive) per haplotype block that needs extension.

  ***5) Does phase-extender do GW (genome wide) or CW (chromosome wide) haplotype phasing***?  
  There are certain situation when phase-extender is able to do GW or CW haplotype phasing.
    A) If you have lots of samples where the haplotypes breakpoint in one sample is bridged by other samples, such that breakpoint gets solved with each recursive application of `phase-Extender` then it is possible to obtain CW and GW haplotype.
    In this case we can run phase-extender for each sample there by extending the haplotype to certain extent. After this phase-extender can be applied recursively on the updated data each time, there by extending the haplotypes for each sample to full chromosome length and possibly to to full genome wide length. There is greater likelyhood of obtaining GW phase if samples are sequenced at higher coverage, increase in paired-end sequence length, availability of large sequence reads like pac-bio reads.
    B) Another situation when GW, CW phase extension might be possible is when you have at least few samples which have haplotype resolved at GW/CW level. These can include fully phased data like genome matrix file, fully phased VCF data, fully phased haplotype reference panel. For this the fully phased sample should be provided as one single blocks in the group of sample that is piped to `phase-extender`.

  ***6) Does phase-extender phase non readbackphase SNPs***?  
    No, it does not. It is a possible future update.

  ***7) Does phase-extender impute missing genotypes***?  
    No, it does not. It is a possible future update.

  ***8) Does phase-extender use haplotype reference panel***?  
    Yes, it does. Thought, the VCF (haplotype reference panel) should be convert to appropriate haplotype file.

  ***9) Does phase-extender use recombination into account***?  
    No and possibly these feature will be of least importance in phase-extender. Main objective of `Phase-Extender` is to join already phased short consecutive haplotype blocks with in a sample by using the relationship of the variants at those sites in several other samples. These haplotypes which are phased in other samples but has breakpoint in SOI are used to build transition probabilities. There is an assumption that recombination is less likely to occur exactly at that breakpoint or near it. So, most of the variation in haplotype among samples around the break point are not the result of recent recombination but only mutation.

  ***10) Does phase-extender phase rare genotypes***?  
    Yes, it does. But, the rare genotype should be the readbackphased to the short haplotype blocks. This is one of the advantage of `phase-extender` compared to other tools when it come to phasing rare genotype. When a single SNPs is used singly to phase into a haplotype, rare genotypes are really hard to phase accurately - the reason being the statistical significance of the rare genotype belonging to either two phase state is highly ambigous. But, if the rare genotype is attached to a haplotye block supported by several read-back phased genotypes, this makes phasing of rare genotypes most accurate, since likelyhoods are provided by other SNPs that are not rare.

  ***11) How fast is phase-extender***?  
    phase-extender is written in python-3, so it is comparatively slower than other tools that are built on the top of C, C++ or java.
   Coming from a pure biology background, learning python was one of the most enduring task I have taken and then building this tool was a big part of my PhD. I have optimized the part of calling VCF file using cyvcf2 (which is on average 4 times faster than old pyVCF module). phase-extender is also optimized for being able to run on multiple threads/process. But, if you are running phase-extender on big genome data and have very large number of samples, and running on laptop I suggest running on one thread, which may be time consuming but will reduce memory burden.

  ***12) Does phase-extender do trio based phase extension***?  
    No, it does not. It is a possible future update.

  ***13) What should be the relatedness of my samples***?  
    Within population, or within species level data are good.

  ***14) What is differece between phase extender and phase stitcher***?  
    `phase-Extender` is a general puprpose haplotype phasing tools. `phase-Stitcher` is specifically for F1 hybrids.

  **_15) Should I prepare my haplotype block file only using `phaser`_**?  
    `phase-Extender`, `phase-Stitcher` can be use with data generated by any RBphasing tool.




## Acknowledgement

   I have not been very fortunate to surround myself or at least get face to face help from savvy computer programmers. But, my heart is very thankful to people behind the web who have made me capable of working this problem out. **Thanks to many people on biostars, stackoverflow, seqanswer and google web searches who provided feedback on small question that were the part of `phase-Extender` project.**

   Should anyone be interested in futher improving this project via improvments on alrorithm and programming, I would be more than happy to.  




## Expected capabilities in the future (coming soon)

- Phase SNPs that are not assigned to ReadBackPhased blocks
- Genotype imputation
- Trio based phasing, Family based phasing
- Higher order markov chain capabilities
- Multiprocessing within chromosome
