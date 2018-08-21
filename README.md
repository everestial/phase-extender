
# phASE-Extender

**Extender** for the readbackphased haplotype blocks.
***A python program to extend the ReadBackPhased haplotype blocks using markov first order transition probabilities and likelihood test.***

Developed by [Bishwa K. Giri](mailto:kirannbishwa01@gmail.com) in the [Remington Lab](https://biology.uncg.edu/people/david-remington/) at the University of North Carolina at Greensboro, Biology department.

## Citation
Giri, B. K., Remington D. L. PhaseIT - A haplotype phasing too for heterogenous and hybrid genomes and for emerging genomic models using phase-Extender and phase-Stitcher. biorxiv (2018) [not uploaded yet].

## AUTHOR/SUPPORT
Bishwa K. Giri (bkgiri@uncg.edu; kirannbishwa01@gmail.com) \
Support @ https://groups.google.com/d/forum/phase-extender

## Intro to ReadBackPhasing
Two heterozygote genotypes in a diploid organism genome are called to be readback phased if they are supported by aligned reads sequence. Depending upon the size and type (single end vs. paried end) of read sequence library, type of read sequence (DNAseq vs. RNAseq) readbackphased haplotype can range from the size of 2 genotypes to multiple genotypes.

**Check these links for more details on readbackphasing**
- https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_phasing_ReadBackedPhasing.php
- https://github.com/secastel/phaser/tree/master/phaser

<br>

# BACKGROUND
Haplotype phasing is a second "go to" problem in bioinformatics after read alignment. The importance of haplotype phasing applies directly to the analyses of ASE (allele specific expression), preparation of extended haplotype for EHH (extended haplotype homozygosity) test, and preparation of dipolid genome which will soon be a new standard in bioinformatics in coming years, etc. The necessity for haplotype phasing (and eventually diploid genome) increases with the increase in heterozygosity in the genome because higher hetetogeneity leads to larger alignment bias and complicates the reliability of the variants that are called using that alignment data (SAM, BAM files).

The classical approach to haplotype phasing involves application of LD (linkage disequilibrium) test between two heterozygous markers, that began with the preparation of genetic map by Alfred Sturtevant. Any existing population based haplotype phasing tool therefore uses the LD test with varying degree of complexity based on available sample size, markers, types of the marker and relationship between samples. For haplotype phasing tools like [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html), [ShapeIT](http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html), [impute2](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html), etc are used dominantly. These tools make use of the variants (SNPs, InDels) along the length of the genome by treating chain of variants along genomic positions singly and independently. So, for a diploid organism that contains chain of variants with "n" heterozygote sites there exists **"2<sup>n</sup>"** possible haplotypes. To certain extent this **"2<sup>n</sup>"** problem is approached and eased by sampling genotype data in short genomic intervals from several samples and applying identity by descent (IBD), most common haplotype method, etc. on the sampled genotypes to infer the possible haplotype in that region. However, application of these methods may not be optimal in solving phase states in organisms that have highly heterogenous genome, are hybrids and/or have very few or no reference panels and abundant genotype data.

**The main issues of the tools that deal with haplotype phasing in **"2<sup>n</sup>"** way can be summarized as:**
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

In several tutorial examples (test files below) I have used hapotype file prepared from RBphased VCF generated by `phaser` (https://github.com/secastel/phaser , https://github.com/secastel/phaser/tree/master/phaser). However, **phase-Extender** can be used with input haplotype file prepared from any RBphased VCF given it meets the appropriate data structure. Readbackphased VCF can be converted into `haplotype file` using an add-on tool `vcf_to_table-v3.py`.\
VCF from `phaser` consists of `Phased Genotype i.e PG` and `Phase Block Index i.e PI` values in `FORMAT` field; **PI** represents the `unique haplotype block index` and **PG** represents the `phased genotypes within that PI block`.  After converting the `RBphased VCF` to `haplotype file` **phase-Extender** uses the `Phased Genotype i.e PG` and `Phase Block Index i.e PI` values in the `haplotype file` to prepare the transition matrix probabilities and proceed with phase extension.

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
  - [numpy](http://www.numpy.org/), (http://cs231n.github.io/python-numpy-tutorial/)
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

**If there is issue with installation while using `requirements.txt` install each dependencies individually.**  
e.g: 
```
sudo python3 -m pip install pandas
``` 
  
<br>

## Usage and Inputs
  - Requires a multisample readbackphased `haplotype file` as input and returns a single sample extended haplotype file. Other results files containing statistics on the initial vs. extended haplotype are also produced. 
  - Optionally, haplotype reference panel (with same data structure as input haplotype) and bed file can be included to limit or improve the process of phase extension.

?? needs improvement.  
Check this detailed [step by step tutorial](https://github.com/everestial/phase-Extender/wiki/phase-Extender-Tutorial) for preparation of `input files` and know-how about running `phase-Extender`.
??
    
<br>

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
https://github.com/everestial/SmallTools/tree/master/GtfToTable 
https://github.com/melissacline/TCGA-GAF-source/blob/master/scripts/gffToBed.py
    
Continue .............
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
* **--numHets** _(40)_ - num of heterozygotes. Maximum number of heterozygote SNPs used from each consecutive block to compute maximum likelihood estimate of each configuration between two blocks.
* **--culLH** _(maxPd)_ - cumulation of the likelihood estimates. The likelhoods for two possible configuration can either be "maxed as sum" or "maxed as product". ***Default*** is "max-product". ***Options:*** 'maxPd' or 'maxSum'.
* **--lods** _(5)_ - log2 of Odds cut off threshold. The cutoff threshold used to extend consecutive haplotype blocks. **`**Note: Default value is set at (2^5 = 32 times likely). So, two consecutive blocks will be joined in parallel configuration if computed log2(likelihood) > lods threshold **
* **--useSample** _(all)_ - Samples to use in the given input haplotype file (plus reference haplotype) to compute transition matrix. Options: 'all','refHap','input','comma separated name of samples'. Default: all the samples in (refHap + input) will be used.
* **--bed** - Process the haplotype extension only within this bed regions. ***This is useful if you want to limit haplotype extension only within certain regions, like - within genes, exons, introns, QTL boundries, etc.*** 
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
This file contains ReadBackPhased haplotype after phase extension concated with the missing data. This file contain equal number of row as input haplotype file and data only for sample of interest. 

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
  We generally need accurate phase state with in a gene/transcript level while doing ASE, dissecting maternal-paternal effects. Long haplotypes are mostly important while preparing diploid genome, testing selective sweeps within a QTL regions etc. For emerging organism systems where genotype data are sparse CW (chromosome wide), GW (genome wide) haplotypes are more difficult to solve. Also, haplotype phasing may be more complicated in out crossing individuals and hybrids due to heterogenity.\ RBphasing actually provides an advantage with heterogenous genome because frequency of RBphased blocks increases with heterogenity in the genome. These short haplotype fragments have multiple heterozygous variants on a short haplotype block. With increase in the size of `PE` (paired end) reads the size of RBphase blocks are also increases. phase-extender comes handy at this stage; that it tries to solve phase state of two consecutive blocks from one sample at a time by using data from haplotype blocks of other samples that bridge that breakpoint. So, we can solve haplotype configuration for SOI (sample of interest) with more confidence because: 
     - we have more variants within each blocks contributing to more information.
     - we only need to solve two possible phase state at one time compared to 2^n haplotype when reading one SNP at a time.
    phase-Extender also provides a more flexible and manipulative control over how to proceed with phase extension. It is also possible to control several parameters like `lods`, `snpTh`, `numHets`, `culLH`, `bed`, `useSample` to observe and compare how phase extension changes.     
     
  **_3) Does phase-extender phase InDels ?_**  
    Yes, but it is conditional. The InDels should already be reliably readbackphased to a haplotype block. That way when the haplotype is being extended for those SNPs, InDels hitchhike with it and get extended too.
    
  **_4) What is the minimal size of the haplotype block that is required?_**  
    The bigger the two haplotypes are, the better is the likelyhood test of which haplotype is phased with which. By, default I have kept this number to 3 variants (SNPs exclusive) per haplotype block that needs extension.
    
  **_5) Does phase-extender do GW (genome wide) or CW (chromosome wide) haplotype phasing_**?  
  There are certain situation when phase-extender is able to do GW or CW haplotype phasing.
    A) If you have lots of samples where the haplotypes breakpoint in one sample is bridged by other samples, such that breakpoint gets solved with each recursive application of `phase-Extender` then it is possible to obtain CW and GW haplotype. 
    In this case we can run phase-extender for each sample there by extending the haplotype to certain extent. After this phase-extender can be applied recursively on the updated data each time, there by extending the haplotypes for each sample to full chromosome length and possibly to to full genome wide length. There is greater likelyhood of obtaining GW phase if samples are sequenced at higher coverage, increase in paired-end sequence length, availability of large sequence reads like pac-bio reads.
    B) Another situation when GW, CW phase extension might be possible is when you have at least few samples which have haplotype resolved at GW/CW level. These can include fully phased data like genome matrix file, fully phased VCF data, fully phased haplotype reference panel. For this the fully phased sample should be provided as one single blocks in the group of sample that is piped to `phase-extender`.
       
  **_6) Does phase-extender phase non readbackphase SNPs_**?  
    No, it does not. It is a possible future update.
      
  **_7) Does phase-extender impute missing genotypes_**?  
    No, it does not. It is a possible future update. 
    
  **_8) Does phase-extender use haplotype reference panel_**?  
    Yes, it does. Thought, the VCF (haplotype reference panel) should be convert to appropriate haplotype file.
   
  **_9) Does phase-extender use recombination into account_**?  
    No and possibly these feature will be of least importance in phase-extender. Main objective of `Phase-Extender` is to join already phased short consecutive haplotype blocks with in a sample by using the relationship of the variants at those sites in several other samples. These haplotypes which are phased in other samples but has breakpoint in SOI are used to build transition probabilities. There is an assumption that recombination is less likely to occur exactly at that breakpoint or near it. So, most of the variation in haplotype among samples around the break point are not the result of recent recombination but only mutation.
    
  **_10) Does phase-extender phase rare genotypes_**?  
    Yes, it does. But, the rare genotype should be the readbackphased to the short haplotype blocks. This is one of the advantage of `phase-extender` compared to other tools when it come to phasing rare genotype. When a single SNPs is used singly to phase into a haplotype, rare genotypes are really hard to phase accurately - the reason being the statistical significance of the rare genotype belonging to either two phase state is highly ambigous. But, if the rare genotype is attached to a haplotye block supported by several read-back phased genotypes, this makes phasing of rare genotypes most accurate, since likelyhoods are provided by other SNPs that are not rare.
    
  **_11) How fast is phase-extender_**?  
    phase-extender is written in python-3, so it is comparatively slower than other tools that are built on the top of C, C++ or java. 
   Coming from a pure biology background, learning python was one of the most enduring task I have taken and then building this tool was a big part of my PhD. I have optimized the part of calling VCF file using cyvcf2 (which is on average 4 times faster than old pyVCF module). phase-extender is also optimized for being able to run on multiple threads/process. But, if you are running phase-extender on big genome data and have very large number of samples, and running on laptop I suggest running on one thread, which may be time consuming but will reduce memory burden.    
    
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
   
   Should anyone be interested in futher improving this project via improvments on alrorithm and programming, I would be more than happy to.  

<br>
<br>

## ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)Expected capabilities in the future (coming soon)
### Phase SNPs that are not assigned to ReadBackPhased blocks
### Genotype imputation
### Trio based phasing, Family based phasing
### Higher order markov chain capabilities
### Multiprocessing within chromosome
