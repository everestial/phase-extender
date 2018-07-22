
# phASE-Extender

**Extender** for the readbackphased haplotype blocks.\
***A python program to extend the ReadBackPhased haplotype blocks using markov first order transition probabilities.***

Developed by [Bishwa K. Giri](mailto:kirannbishwa01@gmail.com) in the [Remington Lab](website?) at the University of North Carolina at Greensboro, Biology department.

## Citation
Giri, B. K., Remington D. L. Haplotype phasing in heterozygous genome and hybrids using phase-Extender and phase-Stitcher. biorxiv (2018) [not uploaded yet].

## AUTHOR/SUPPORT
Bishwa K. Giri (bkgiri@uncg.edu; kirannbishwa01@gmail.com) \
Support @ https://groups.google.com/d/forum/phase-extender

## Intro to ReadBackPhasing
**Check these links for details on readbackphasing**
- https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_phasing_ReadBackedPhasing.php
- https://github.com/secastel/phaser/tree/master/phaser

<br>
<br>

# BACKGROUND
Haplotype phasing is the second "go to" problem in bioinformatics after read alignment. The importance of haplotype phasing applies directly to the analyses of ASE (allele specific expression), preparation of extended haplotype for EHH (extended haplotype homozygosity) test, and preparation of a diploid genome which will soon be a new standard in bioinformatics in the coming years. The necessity for haplotype phasing (and eventually diploid genome) increases with the increase in heterozygosity in the genome because higher hetetogeneity leads to larger alignment bias and complicates the reliability of the variants that are called using alignment data (SAM, BAM files).

The classical approach to haplotype phasing involves application of LD (linkage disequilibrium) test between two heterozygous markers, that began with the preparation of genetic map by Alfred Sturtevant. Any existing population based haplotype phasing tool therefore uses the LD test with varying degree of complexity based on available sample size, markers, types of the marker and relationship between samples. For haplotype phasing tools like [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html), [ShapeIT](http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html), [impute2](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html), etc are used dominantly. These tools make use of the variants (SNPs, InDels) along the length of the genome by treating chain of variants along genomic positions singularly and independently. So, for a diploid organism that contains chain of variants with "n" heterozygotic sites there exists "2^n" possible haplotypes. To a certain extent this "2^n" problem is approached and eased by sampling genotype data from short genomic regions and then applying a popular method called identity by descent (IBD) on that fragment to infer the possible haplotype. However, application of these methods may not be optimal in solving phase states in organisms with a heterozygous genome, hybrids and/or have very few or no reference panels and abundant genotype data.

**The main issues of the tools that deal with haplotype phasing in "2^n" way can be summarized as:**
  - Increased computation burden due to "2^n" problem.
  - Problem in phasing rare variants. 
  - Is mostly applicable to human and organism with reference haplotype panel and abundant genomic data. 
  - Not optimal for organisms that have heterozygous genome, or are outbreds, or are hybrids, or belong to recent group of organisms that have reference genome prepared.  
  
Therefore, using LD between two adjacent SNPs with small population of samples will not be able to provide enough resolution to solve GW haplotype preparation, which could result in excessive switch errors. Additionally, in heterozygous and hybrid genomes, problems arising due to switch errors in downstream analyses could be manifold.

ReadBackphasing is slowly emerging as a new and more reliable method for preparation of short range haplotypes by joining heterozygous variants that are covered by sequence reads. These short haplotypes can be further elongated by adding sequencing reads of longer length (PacBio reads) or by adding more genomic and RNAseq reads from the same individual; see tools [WhatsHap](https://whatshap.readthedocs.io/en/latest/), [hapCut](https://github.com/vibansal/hapcut), [phaser](https://github.com/secastel/phaser/blob/master/phaser/README.md) etc.  

**But, issues with existing RBphase method and/or tools still remain and can be summarized as:**
   - are mainly aimed at preparing long range haplotype but not necessarily genome wide. 
   - existing RBphase methods are only targeted at individual or family level, i.e they require multiple inputs of "BAM" and "VCF" files for the same individual and/or trios. The increase in the size of phased haplotype blocks only depend upon multiple BAM files, or multiple sets of longer reads from the same individual.
   - Referring to point 2 -> integration of RBPhased data with population based phasing is still missing i.e they are not able to solve phase state of two haplotype blocks in same sample by using information phase state of the haplotype from other samples.

Even with large PE reads in future there will always be coverage issues which will break genome wide haplotype in multiple segments. It is therefore imperative to establish a method that is able to integrate the RBphased data among samples to solve and prepare long range haplotypes in individual genomes. While integrating RBphase methods with population data, the concept of LD still applies but the phasing accuracy is highly increased owing to phased information contained in short haplotype blocks across several samples. Additionally, this method will put heterozygous genome at most advantage because of large RBphase blocks they are able to produce. 

For specimens that don't have multiple genomic and transcriptomic sources, reference haplotype panel the existing RBphasing tools are only able to yield short range haplotypes. Even with multiple bam inputs these tools are still not aimed at preparing genome wide haplotype, but rather are limited at preparing long streches of haplotype for certain parts of a genome. Additionally, these RBphase tools do not have method to solve phase state between two adjacent haplotype blocks that are not spanned by sequence reads (in a sample) i.e there is no method incorporated to transfer readbackphased haplotype data information across population of samples.  

With increase in the size of PE (paired end) reads from Illumina and availability of longer sequences from PacBio, it is possible to prepare longer short range haplotypes that span more than 2 heterozygous site in the genome. Thus, the issue with "2^n" problem can be overcome by 
    **1)** first preparing RBphased haplotype blocks using sequence reads, 
    **2)** Next, the two consecutive haplotypes blocks can be joined in a sample one at a time by reading phase state of the break point in other population of samples. The size of RBphased blocks also increases with increase in heterozygosity in the genome thereby making this method more suitable for highly heterozygous populations. And, there is always more information in RBphased haplotype compared with a single SNP or InDel thereby overcoming the issues with switch-errors when preparing the long range haplotype in heterozygous population. 

So, RBphasing method is able to yield higher variants per haplotype block thus making phase-Extender a better method and tool when working with genotypes (or genomes )having higher heterozygosity (out crossing population and hybrids). 

Existing haplotype phasing tools (that use RBphase and non-RBphase method) are oriented towards organism that have large genotype data and several haplotype reference panels available. For non-model organism these methods are therefore sub-optimal. RBphasing provides some phasing benefits because the phase state along short blocks that can be readily extracted. However, in emerging research models converting the short haplotype fragments to longer haplotype strings can be challenging due to a small sample size, low coverage within each sample, lack of enough genotype data, absence of reference panels. **phASE-Extender** overcomes these limitations by using the short range RBphased haplotype blocks from several samples that have several haplotype break points. **phase-Extender** is aimed at filling this technical gap by using the short fragments haplotypes from several samples, where the breakpoint of one sample is bridged by several other samples. 

**phASE-Extender** starts with already prepared RBphased haplotype blocks of several individuals at once and joins the two RBphased blocks for a single sample at once using overlapping phase states in other samples of the pool. Thus it is possible to solve the haplotype breakpoints by transferring the phase information between samples that belong to same family, population or species. *So, if we have a pool of several samples that have RBphased haplotypes, each sample will have several non-overlapping breakpoints among the samples. We can solve the phase state between two consecutive blocks (of one sample at a time) by computing the likelihood of the phase states (parallel vs. alternate configuration). Likelihood of each configuration is computed by using the observed phase state in other samples that bridge that breakpoint position.* In **phASE-Extender** we use first order markov chain to compute the relationship between two consecutive haplotype blocks by preparing first order transition matrix from all the nucleotides of former haplotype blocks to all the nucleotides in the latter haplotype block. We then compute cumulative likelihood estimates of haplotype states for both the possible configuration. The phase state of two consecutive haplotype blocks is then extended if the `computed log2 (likelihood)` of either configuration is above the `threshold log2(likelihood) cutoff`. Using RBphase method we can therefore use small sample size to accurately predict the proper haplotype state. 

## ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)Data Requirements

**phASE-Extender** can be used with the multi-sample vcf files produced by GATK pipeline or other tools that generate readbackphased haplotype blocks in the output VCF. A haplotype file is created using the RBphased VCF and then piped into **phase-Extender**. See, this example for data structure of input haplotype file [sample input haplotype file01](https://github.com/everestial/phase-Extender/blob/master/example01/input_haplotype_small.txt) - a tab separated text file with `PI` and `PG_al` value for each samples.

In several tutorial examples (test files below) I have used hapotype file prepared from RBphased VCF generated by `phaser` (https://github.com/secastel/phaser , https://github.com/secastel/phaser/tree/master/phaser). However, **phase-Extender** can be used with input haplotype file prepared from any RBphased VCF given it meets the appropriate data structure. Readbackphased VCF can be converted into `haplotype file` using an add-on tool `vcf_to_table-v3.py`.\
VCF from `phaser` consists of `Phased Genotype i.e PG` and `Phase Block Index i.e PI` values in `FORMAT` field; **PI** represents the `unique haplotype block index` and **PG** represents the `phased genotypes within that PI block`.  After converting the `RBphased VCF` to `haplotype file` **phase-Extender** uses the `Phased Genotype i.e PG` and `Phase Block Index i.e PI` values in the `haplotype file` to prepare the transition matrix probabilities and proceed with phase extension.

## Algorithm
For the **mcve** regarding the algorithm see this issue on [**stackoverflow**](https://codereview.stackexchange.com/questions/186396/solve-the-phase-state-between-two-haplotype-blocks-using-markov-transition-proba) and/or [**my blog**](https://everestialblog.wordpress.com/2018/03/27/extension-of-the-readbackphased-haplotypes-using-phase-extender/). 

<br>
<br>

# Tutorial

## Setup
**phase-Extender** is written in python3 interpreter. So, it can be directly run using the **".py"** file, given all the required modules (dependencies) are installed.

**Runs on `Python 3.x` and has the following dependencies:**

  - [pandas](http://pandas.pydata.org/)
  - [argparse](https://docs.python.org/3/library/argparse.html)
  - [collections](https://docs.python.org/3/library/collections.html?highlight=collections#module-collections)
  - [csv](https://docs.python.org/3/library/csv.html?highlight=csv)
  - [decimal](https://docs.python.org/3/library/decimal.html?highlight=decimal#module-decimal)
  - [functools](https://docs.python.org/3/library/functools.html?highlight=functools)
  - [itertools](https://docs.python.org/3/library/itertools.html?highlight=itertools)
  - [io](https://docs.python.org/3/library/io.html?highlight=io#module-io)
  - [multiprocessing](https://docs.python.org/3/library/multiprocessing.html?highlight=multiprocessing#)
  - [numpy](http://www.numpy.org/), `http://cs231n.github.io/python-numpy-tutorial/`
  - [resource](https://docs.python.org/3/library/resource.html?highlight=resource#module-resource)
  - [time](https://docs.python.org/3/library/time.html?highlight=time#module-time)
  - [matplotlib](https://matplotlib.org/2.2.2/index.html)
  
<br>
  
## Installation
```
pip3 install -r requirements.txt

#Or
sudo pip3 install -r requirements.txt
```

**If there is issues with installation while using `requirements.txt` install each dependencies individually.**\
`sudo python3 -m pip install pandas`  
  
<br>

## Usage
  - Requires a readbackphased `haplotype file` as input and returns an extended haplotype file and other results files containing statistics on the initial vs. extended haplotype. 
  - Optionally, haplotype reference panel (with same data structure as input haplotype) and bed file can be included to gain control over the outcome of phase extension.

Check this detailed [step by step tutorial](https://github.com/everestial/phase-Extender/wiki/phase-Extender-Tutorial) for the preparation of `input files` and know-how about running `phase-Extender`.
    
<br>

## Input data

***haplotype file (required):*** Input haplotype file. Should contain `PI` and `PG_al` values for each sample.\
To convert the vcf file to haplotype file (from VCF to proper text format) use **Step 01 (a)** in the tutorial. \
The sample name should not contain "`_`" character.

***haplotype reference panel (optional):*** Haplotype reference panel as text file (should have same data structure as input haplotype file). \
To convert the haplotype reference panel (from VCF to proper text format) use **Step 01 (b)** in the tutorial.
  
***bed file (optional):*** If your goal is to limit phase extension to certain genomic regions (for eg. gene, exon or QTL boundries), we suggest that you provide appropriate bed file. **phase-Extender** is exclusively limited to internal boundries of bed regions.

    # structure of the bed file
    contig    start    end       # this header is not included    
    2         1258     199897
    2         397765   412569

    #To convert the GTF,GFF to bed file use:
    python3 gffToBed.py  --input myGTF.gtf  --output myBed.bed
    
    
## Arguments
### Required
* **--input** - input haplotype file. `PI` and `PG_al` should be present in header for each sample.
* **--SOI** - sample of interest. It should refer to a single sample in the haplotype the file. The sample name should not contain "`_`" character. 

### Performance Related
* **--nt** _(1)_ - maximum number of processes to run at once. The maximum number of processes is limited to number of chromosomes (contigs) in the input haplotype file. 

### Optional
* **--python_string** _(python3)_ - Calls `python 3` interpreter to run the program.
* **--output** _(SOI_extended)_ - Output directory. 
* **--snpTh** _(3)_ - snp threshold. Minimum number of SNPs required in each consecutive haplotype block to run phase extension between two blocks.
* **--numHets** _(40)_ - num of heterozygotes. Maximum number of heterozygotic SNPs used from each consecutive block to compute maximum likelihood estimate of each configuration between two blocks.
* **--culLH** _(maxPd)_ - cumulation of the likelihood estimates. The likelihoods for two possible configuration can either be "maxed as sum" or "maxed as product". ***Default*** is "max-product". ***Options:*** 'maxPd' or 'maxSum'.
* **--lods** _(5)_ - log2 of Odds cut off threshold. The cutoff threshold used to extend consecutive haplotype blocks. **`**Note: Default value is set at (2^5 = 32 times likely). So, two consecutive blocks will be joined in parallel configuration if computed log2(likelihood) > lods threshold **
* **--useSample** _(all)_ - Samples to use in the given input haplotype file (plus reference haplotype) to compute transition matrix. Options: 'all','refHap','input','comma separated name of samples'. Default: all the samples in (refHap + input) will be used.
* **--bed** - Process the haplotype extension only within this bed region. ***This is useful if you want to limit haplotype extension only within certain regions, like - within genes, exons, introns, QTL boundries, etc.*** 
* **--writeLOD** _(no)_ - writes the calculate LODs between two consecutive haplotype blocks when processing phase extension to the output file. **Options:** 'yes', 'no'. **`**Note: the 'lods-score' are printed regardless if the "
"consecutive blocks are joined or not.**
* **--hapStats** _(no)_ - Prepare descriptive statistics, and histogram of the haplotype size distribution of the input haplotype file vs. extended haplotype for the sample of interest. **Options:** 'yes', 'no'
* **--addMissingSites** _(no)_ - include the non-phased and missing genotype data from input haplotype file to the final phase extended output file. **Option:** 'yes', 'no'.

<br>
<br>

## Output Files

### initial_haplotype_*SOI*.txt & extended_haplotype_*SOI*.txt
Contains all RBphased haplotype data for the sample of interest before and after phase extension.
* 1 - **contig** - Contig name (or number).
* 2 - **pos** - Start position of haplotype (1 based).
* 3 - **ref** - Reference allele at that site.
* 4 - **all-alleles** - All the alleles represented by all the samples at that site.
* 5 - **SOI_PI** - Unique `PI` index of the haplotype blocks for sample of interest.
* 6 - **SOI_PG_al** - Phased GT (genotype) alleles at the genomic position that belong to unique `PI` indexes.
* 7 - **log2Odds** (only in **extended_haplotype_SOI.txt**) - log2Odds computed between the former and later block.

<br>

### initial_haplotype_stats_*SOI*.txt & final_haplotype_stats_*SOI*.txt
Descriptive haplotype statistics of the input haplotype file for the sample of interest.
* 1 - **contig** - Contig name (or number).
* 2 - **SOI_PI** - Comma separated list of unique `PI` index of the haplotype blocks for the sample of interest. The total number of `PI` index represents the total number of haplotype fragments that are present in the given contig in that sample.
* 3 - **num_Vars_by_PI** - Number of variants sites within each `PI` block for the sample of interest. 
* 4 - **range_of_PI** - Genomic range of the each `PI` block for the sample of interest.
* 5 - **total_haplotypes** - Total number of haplotype (i.e `PI`) in the given coting for the sample of interest.
* 6 - **total_Vars** - Total number of variant sites in the given contig for the sample of interest. **Note:** The sum of (num_Vars_by_PI) = total_Vars.

**Note:** - The `SOI_PI`, and it's associated statistics are in order.

<br>

### missingdata_*SOI*.txt
Contains data from the sites that have unphased or missing GT (genotype) for the sample of interest in the input haplotype file.
**Note:** This data is merged with `extended_haplotype_SOI.txt` if `--addMissingSites` is set to "yes".

<br>

### extended_haplotype_"SOI_"allsites.txt
This file contains ReadBackPhased haplotype after phase extension concatenated with the missing data. This file contains equal number of row as input haplotype file and data only for sample of interest. 

<br>
<br>

## Plots
**Note:** - These plots are based on the descriptive statistics generated for haplotypes before and after phase extension. It is possible to take these statistics (initial_haplotype_stats_*SOI*.txt & final_haplotype_stats_*SOI*.txt) and make custom plots in **R** or by using other methods.

<br>

### total_haps_"SOI_"initial.png  & total_haps_"SOI_"final.png
Number of haplotypes for given contig before and after phase extension. 

<br>

### total_vars_"SOI_"initial.png  & total_vars_"SOI_"final.png
Number of variants for given contig before and after phase extension. 

<br>

### hap_size_byVar_"SOI_"initial.png & hap_size_byVar_"SOI_"final.png
Histogram of the distribution of the haplotype size (by number of variants in the haplotype) in given contig before and after phase extension. 

<br>

### hap_size_byGenomicRange_"SOI_"initial.png & hap_size_byGenomicRange_"SOI_"final.png
Histogram of the distribution of the haplotype size (by genomic range of the haplotype) in given contig before and after phase extension.

<br>
<br>

## ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)Some Q/A on phase-extender: 

  **_1) What kind of algorithm does phase-extender use ?_**  
  phase-extender uses first-order-transition probabilities from each level of genotypes from former haplotype block to each level of genotypes to later haplotype block. This version (v1) uses **forward-1stOrder-markov chains** and **backward-1stOrder-markov chains** transition probabilities. Future versions will follow improvements by adding markov-chains of higher order.    
    
  **_2) What is the advantage of using phase-extender ?_**  
  We generally need accurate phase state with in a gene/transcript level while doing ASE, dissecting maternal-paternal effects. Long haplotypes are mostly important while preparing diploid genome, testing selective sweeps within a QTL regions etc. For emerging organism systems where genotype data are sparse CW (chromosome wide), GW (genome wide) haplotypes are more difficult to solve. Also, haplotype phasing may be more complicated in out crossing individuals and hybrids due to heterogeneity.\ RBphasing actually provides an advantage with heterozygous genome because frequency of RBphased blocks increases with heterogenity in the genome. These short haplotype fragments have multiple heterozygous variants on a short haplotype block. With increase in the size of `PE` (paired end) reads the size of RBphase blocks are also increases. phase-extender comes handy at this stage; it tries to solve phase state of two consecutive blocks from one sample at a time by using data from haplotype blocks of other samples that bridge that breakpoint. So, we can solve haplotype configuration for SOI (sample of interest) with more confidence because: 
     - we have more variants within each blocks contributing to more information.
     - we only need to solve two possible phase state at one time compared to 2^n haplotype when reading one SNP at a time.
    phase-Extender also provides a more flexible and manipulative control over how to proceed with phase extension. It is also possible to control several parameters like `lods`, `snpTh`, `numHets`, `culLH`, `bed`, `useSample` to observe and compare how phase extension changes.     
     
  **_3) Does phase-extender phase InDels ?_**  
    Yes, but it is conditional. The InDels should already be reliably readbackphased to a haplotype block. That way when the haplotype is being extended for those SNPs, InDels hitchhike with it and get extended too.
    
  **_4) What is the minimal size of the haplotype block that is required?_**  
    The bigger the two haplotypes are, the better is the likelihood test of which haplotype is phased with which. By, default I have kept this number to 3 variants (SNPs exclusive) per haplotype block that needs extension.
    
  **_5) Does phase-extender do GW (genome wide) or CW (chromosome wide) haplotype phasing_**?  
  There are certain situation when phase-extender is able to do GW or CW haplotype phasing.
    A) If you have lots of samples where the haplotypes breakpoint in one sample is bridged by other samples, such that breakpoint gets solved with each recursive application of `phase-Extender` then it is possible to obtain CW and GW haplotype. 
    In this case we can run phase-extender for each sample there by extending the haplotype to certain extent. After this phase-extender can be applied recursively on the updated data each time, there by extending the haplotypes for each sample to full chromosome length and possibly to to full genome wide length. There is greater likelihood of obtaining GW phase if samples are sequenced at higher coverage, increase in paired-end sequence length, availability of large sequence reads like pac-bio reads.
    B) Another situation when GW, CW phase extension might be possible is when you have at least few samples which have haplotype resolved at GW/CW level. These can include fully phased data like genome matrix file, fully phased VCF data, fully phased haplotype reference panel. For this the fully phased sample should be provided as one single blocks in the group of sample that is piped to `phase-extender`.
       
  **_6) Does phase-extender phase non readbackphase SNPs_**?  
    No, it does not. It is a possible future update.
      
  **_7) Does phase-extender impute missing genotypes_**?  
    No, it does not. It is a possible future update. 
    
  **_8) Does phase-extender use haplotype reference panel_**?  
    Yes, it does. Though, the VCF (haplotype reference panel) should be converted to appropriate haplotype file.
   
  **_9) Does phase-extender use recombination into account_**?  
    No and possibly these feature will be of least importance in phase-extender. Main objective of `Phase-Extender` is to join already phased short consecutive haplotype blocks with in a sample by using the relationship of the variants at those sites in several other samples. These haplotypes which are phased in other samples but has breakpoint in SOI are used to build transition probabilities. There is an assumption that recombination is less likely to occur exactly at that breakpoint or near it. So, most of the variation in haplotype among samples around the break point are not the result of recent recombination but only mutation.
    
  **_10) Does phase-extender phase rare genotypes_**?  
    Yes, it does. But, the rare genotype should be the readbackphased to the short haplotype blocks. This is one of the advantage of `phase-extender` compared to other tools when it come to phasing rare genotype. When a single SNPs is used singly to phase into a haplotype, rare genotypes are really hard to phase accurately - the reason being the statistical significance of the rare genotype belonging to either two phase state is highly ambigous. But, if the rare genotype is attached to a haplotye block supported by several read-back phased genotypes, this makes phasing of rare genotypes most accurate, since likelihoods are provided by other SNPs that are not rare.
    
  **_11) How fast is phase-extender_**?  
    phase-extender is written in python-3, so it is comparatively slower than other tools that are built on the top of C, C++ or java. 
   Coming from a pure biology background, learning python was one of the most enduring task I have taken and then building this tool was a big part of my PhD. I have optimized the part of calling VCF file using cyvcf2 (which is on average 4 times faster than the old pyVCF module). phase-extender is also optimized for being able to run on multiple threads/process. But, if you are running phase-extender on big genome data and have very large number of samples, and running on laptop I suggest running on one thread, which may be time consuming but will reduce memory burden.    
    
  **_12) Does phase-extender do trio based phase extension_**?  
    No, it does not. It is a possible future update.   
   
  **_13) What should be the relatedness of my samples_**?  
    Within population, or within species level data are good.   
   
  **_14) What is differece between phase extender and phase stitcher_**?  
    `phase-Extender` is a general puprpose haplotype phasing tools. `phase-Stitcher` is specifically for F1 hybrids.   
   
  **_15) Should I prepare my haplotype block file only using `phaser`_**?  
    `phase-Extender`, `phase-Stitcher` can be use with data generated by any RBphasing tool.
   
<br>
<br>
   
   ## ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)Acknowledgement
   I have not been very fortunate to surround myself or at least get face to face help from savvy computer programmers. But, my heart is very thankful to people behind the web who have made me capable of working this problem out. **Thanks to many people on biostars, stackoverflow, seqanswer and google web searches who provided feedback on small question that were the part of `phase-Extender` project.**    
   
   Should anyone be interested in futher improving this project via improvements on algorithm and programming, I would be more than happy to.  

<br>
<br>

## ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)Expected capabilities in the future (coming soon)
### Phase SNPs that are not assigned to ReadBackPhased blocks
### Genotype imputation
### Trio based phasing, Family based phasing
### Higher order markov chain capabilities
### Multiprocessing within chromosome
