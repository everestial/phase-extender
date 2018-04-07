# phASE-Extender
**Extender** for the readbackphased haplotype blocks.
***a python program to extend the ReadBackPhased haplotype blocks using first order transition probabilities.***

Developed by [Bishwa K. Giri](mailto:kirannbishwa01@gmail.com) in the [Remington Lab](website?) at the University of North Carolina at Greensboro, Biology department.

# Citation
Giri, B. K., Remington D. L. Haplotype phase extension and preparation of diploid genome using phase-Extender. biorxiv (2018) [not uploaded yet].


# Intro to ReadBackPhasing
**Check these links for details on readbackphasing**
- https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_phasing_ReadBackedPhasing.php
- https://github.com/secastel/phaser/tree/master/phaser


**Runs on `Python 3.x` and has the following dependencies:**

[pandas](http://pandas.pydata.org/)
[argparse](https://docs.python.org/3/library/argparse.html)
[collections](https://docs.python.org/3/library/collections.html?highlight=collections#module-collections)
[csv](https://docs.python.org/3/library/csv.html?highlight=csv)
[decimal](https://docs.python.org/3/library/decimal.html?highlight=decimal#module-decimal)
[functools](https://docs.python.org/3/library/functools.html?highlight=functools)
[itertools](https://docs.python.org/3/library/itertools.html?highlight=itertools)
[io](https://docs.python.org/3/library/io.html?highlight=io#module-io)
[multiprocessing](https://docs.python.org/3/library/multiprocessing.html?highlight=multiprocessing#)
[numpy](http://www.numpy.org/, http://cs231n.github.io/python-numpy-tutorial/)
[resource](https://docs.python.org/3/library/resource.html?highlight=resource#module-resource)
[time](https://docs.python.org/3/library/time.html?highlight=time#module-time)
[matplotlib](https://matplotlib.org/2.2.2/index.html)



# BACKGROUND
Haplotype phasing is a second "go to" problem in bioinformatics after read alignment. The importance of haplotype phasing applies directly to the analyses of ASE (allele specific expression), preparation of extended haplotype for EHH (extended haplotype homozygosity) test, and preparation of dipolid genome which will soon be a new standard in bioinformatics in coming years, etc. The necessity for haplotype phasing (and diploid genome) increases with the increase in heterozygosity in the genome, because higher hetetogeneity leads to larger alignment bias and complicates the reliability of the variants that are called using that alignment data (SAM, BAM files).

For haplotype phasing tools like [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html), [ShapeIT](http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html), [impute](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html), etc are used dominantly. These tools make use of the variants (SNPs, InDels) along the length of the genome and treat each variants as single observations independent of variant in the consecutive position. So, for a chain of variants with "n" heterozygote sites there exists "2^n" possible haplotypes. Several tools treat each site separately and apply identity by descent, most common haplotype method, etc. to infer the possible haplotype for a sample. Also, the issue with "2^n" is increased computation burden.

However, with increase in the size of PE (paired end) reads from illumina and availability of longer sequences from PacBio, issue with "2^n" problem can be overcome by first preparing RBphased haplotype blocks along the genome, and next by joining the two consecutive blocks in a sample one at a time. *So, if we have a pool of several samples that have RBphased haplotypes, each sample will have several non-overlapping breakpoints among the samples. We can solve the phase state between two consecutive blocks (of one sample at a time) by computing the likely hood of the phase states (parallel vs. alternate configuration). Likelihood of each configuration is computed by using the observed phase state in other samples that bridge that breakpoint position.* In **phASE-Extender** we use first order markov chain and transition matrix to compute the relationship between two consecutive haplotype blocks by preparing first order transition matrix from all the nucleotides of former haplotype block to all the nucelotides in the later haplotype block. We then compute cumulative likelyhood estimates of haplotype states for both the possible configuration. The phase state of two consecutive haplotype blocks is then extended if the computed log2 (likelihood) of either configuration is above the threshold log2(likelihood) cutoff. 

ReadBackphasing is slowly emerging as a replacement for preparation of long range haplotypes; see tools [WhatsHap] (https://whatshap.readthedocs.io/en/latest/), [hapCut] (https://github.com/vibansal/hapcut), [phaser] (https://github.com/secastel/phaser/blob/master/phaser/README.md). The issue with existing RBphase method is that these tools prepare the RBPhased haplotype blocks based on multiple input "BAM" and "VCF" files for a single sample and/or trios. The increase in the size of phased haplotype blocks only depends upon multiple BAM files, or multiple sets of longer reads from the same individual.For specimen that don't have multiple genomic and rna sources, reference haplotype panel these method is only possible to yield short range haplotypes. Further these RBphase tool do not have method to transfer phase information across population of samples. 

Also, most of the haplotype phasing tools (that use RBphase and non-RBphase method) are oriented towards organism that have large genotype data and several haplotype reference panels available. For non-model organism these methods are therefore sub-optimal. RBphasing provides some phasing benefits because the phase state along short blocks can be readily extracted. However, converting the short haplotype fragments to longer haplotype strings can be challenging due to small sample size, lack of enough genotype data, absence of reference panels in emerging research models. 

**phASE-Extender** overcomes these limitations by using the short range RBphased haplotype blocks from several samples that have several haplotype break points. **phase-Extender** is aimed at filling this technical gap by using the short fragments haplotypes among several samples, where the breakpoint of one sample is bridged by several other samples. **phASE-Extender** starts with already prepared RBphased haplotype blocks of several individuals at once and joins the two RBphased blocks for a single sample at once using overlapping phase states in other samples of the pool. Thus it is possible to solve the haplotype breakpoints by transfering the phase information between samples that belong to same family, populaton or species. 

# Input data

**contd...**
phASE-Extender can be used with the haplotype file generate from the vcf files produced by GATK pipeline or several other tools that produces readbackphased haplotype blocks - given the input haplotype file meets the appropriate data structure. See, this example of a [sample input file](https://github.com/everestial/phase-Extender/blob/master/example01/input_haplotype_small.txt) - a tab separated text file with `PI` and `PG_al` value for each samples.

In the example file I am using vcf data produced from `phaser` that consists of `PI` and `PG` values in `FORMAT` field and is the most directly applicable input. This tool exclusively uses the `Phased Genotype i.e PG` and `Phase Block Index i.e PI`values to prepare the transition matrix. 
The 
This tool exclusively uses the `Phased Genotype i.e PG` and `Phase Block Index i.e PI` in the VCF (variant call format) generated by the program `pHASER` https://github.com/secastel/phaser , https://github.com/secastel/phaser/tree/master/phaser . PI represents the unique haplotype block index (or key) in the given VCF files and PG represents the phased genotypes within that PI block. So, it is crucial to run **phaser** on your VCF of interest before running phase extension or prepare the "haplotype file" in appropriate format.


# Algorithm
For the **mcve** regarding the algorithm see this issue on **stackoverflow** https://codereview.stackexchange.com/questions/186396/solve-the-phase-state-between-two-haplotype-blocks-using-markov-transition-proba . 
**A neat diagrammatic version will be soon included in my wordpress account **



## Tutorial

# Setup
**phase-Extender** is written in python3 interpreter. So, it can be directly run using the **".py"** file, given all the required modules are installed.

# Usage
Requires a VCF file or a haplotype file as input and returns a extended haplotype file and other results files containing statistics on the initial vs. extended haplotype. Optionally, haplotype reference panel (with same data structure as input haplotype) and bed file can be included to gain control on the phase extension.

VCF file from `phaser` can be converted to haplotype file using the python parser `vcf_to_table-v2.py`. 

    vcf_to_table-v2.py --vcf a_vcf_file.vcf --out a_haploype_file.txt
    
`--unphased yes` can be added to above command to also included the genotypes that have unphased state. This parameter will not affect the phase extension, and is only included to keep the whole data intact.
    
**Test case 01**

Using file from [example 01](https://github.com/everestial/phase-Extender/tree/master/example01)

Run command:

    python3 phase_extender_v1-final.py --input input_haplotype_small.txt --SOI ms02g 

**Test case 02**

Using file from [example 02](https://github.com/everestial/phase-Extender/tree/master/example02)

Run command: ** complete ??

    python3 phase_extender_v1-final.py --input input_haplotype_small.txt --SOI ms02g 
    
**Additional input files **

***haplotype reference panel:*** If your goal is to use reference haplotype panel, we suggest providing the haplotype reference panel with same data structure as input haplotype file.

To convert the haplotype reference panel (from VCF to proper text format) use: ** complete ?? **

    python3 hapRefVCFtoTxt.py --input hapRef.vcf --out  
  
***bed file:*** If you goal is to limit phase extension to certain genomic regions (for eg. gene, exon or QTL boundries), we suggest that you provide appropriate bed file. Remember, **phase-Extender** is exclusively limited to bed regions.   

To convert the GTF,GFF to bed file for appropriate GTF,GFF feature use:  ** complete ??

    python3 gffToBed.py  --input myGTF.gtf  --output myBed.bed
    
# Arguments
## Required
* **--bam** - Comma separated list of BAM files containing reads. Duplicates should be marked, and files should be indexed using samtools index.

# Optional
* **--python_string** _(python2.7)_ - Command to use when calling python, required for running


## Performance Related
* **--nt** _(1)_ - Maximum number of processes to run at once
snpTh
numHets


# Output Files

## *out_prefix*.haplotypes.txt

Contains all haplotypes phased along with detailed phasing statistics.

* 1 - **contig** - Contig haplotype is found on.





## **initial**.hapstats
parameters in this file ??

## **initial**.hapstats
parameters in this file ??


## **some**.plot.png





# Some Q/A on phase-extender: 

  **1) How is phase Extender different than other gw (genome wide), cw (chromosome wide) haplotype phasing program ?**
  
    The very first important thing in genomics and transcriptomics is alignment. While there has been quite some improvements 
    in alinment algorithms and which still continues, another field of genomics that has utmost importance is haplotype phasing - 
    which includes organizing the called variants accurately among each other. So, phasing represents another era of improvements
    in genomics techonology and application. 
    Most, of the haplotype phasing tools are mostly built with human genome (modern and ancient) in mind,but they do benefit 
    other ogranism which have similar level of genomic resources. They rely moslty on huge amount of genotype data and available 
    haplotype reference panel which are not readily available for other organisms. Also, these tools mainly aim for genome wide 
    haplotype, which is of least concern for emerging models. Rather solving the complete phase state within a gene, 
    transcriptome or region of interest is more crucial. 
    Most phasing tools solve the haplotype state problem by taking in SNPs data from multiple samples, then it tries to solve the 
    phase state for one sample by preparing a haplotype that can possibly give rise to all other possible haplotype 
    configurations. Recent advancement in haplotype phasing has started involving "ReadBackPhased" haplotypes. But, most of the
    tools are geared toward connecting the phase state by sequencing different tissues from sample (citation ... ?). However, 
    it is possible to solve the phase state between two haplotype blocks if the position they are broken at is covered by 
    RBPhased haplotyes in other related samples. 
    
    
  **2) What kind of algorithm does phase-extender use ?**
  
    phase-extender uses first-order-transition probabilities from each level of genotypes from former haplotype 
    block to each level of genotypes to later haplotype block. At the moment phase-extender only uses **forward-1stOrder-markov chains** transition probabilities. I soon plan to improve the algorithm by adding **backward-1stOrder-markov chains** 
    transition probabilities. These may follow improvements by adding markov-chains of higher order.
    
    
   **3) What is the advantage of using phase-extender ?**
   
    We generally need accurate phase state with in a gene/transcript level while doing ASE, dissecting maternal-paternal effects.  
    Long haplotypes are mostly important while preparing diploid genome, long range haplotype for testing selective sweeps, 
    within a QTL regions etc. 
    For emerging organism systems where CW,GW haplotypes are more difficult to solve, it is more important to solve haplotypes 
    within short gene regions. 
    But, most of RBphased Haplotype rarely cover full gene regions due to coverage issues, or due to short read sequence that 
    cannot cover whole gene, etc. In that situation, if we rather have short haplotype blocks from severals samples it can be 
    used to solve a full gene length haplotypes for each sample (one at a time). phase-extender comes handy at this stage; that 
    it takes the information from several haplotype blocks from other samples that are not broken at that
    point. So, we can solve haplotype configuration for SOI (sample of interest) pretty much easily. ***Also, with paired-end 
    read size getting increased in Illumina sequence data we can have more reliable full gene-level haplotypes for SOI.***
    
   **4) Does phase-extender phase InDels ?**
    Yes, but it is conditional. The InDels should already be phased to a reliable haplotype block. That way when the haplotype is
    being extended for those SNPs, InDels too get extended.
    
   **5) What is the minimal size of the haplotype block that is required? **  contd ....
    The bigger the two haplotypes, the better is the likelyhood test of which haplotype is phased with which. By, default
    I have kept this number to 3 variants (SNPs exclusive) per haplotype block that needs extension.
    
   Does phase-extender do GW (genome wide) or CW (chromosome wide) haplotype phasing?
    There are certain situation when phase-extender is able to do GW or CW haplotype phasing.
    A) If you have lots of samples where the haplotypes overlap between sample covers almost all chromosome length. In this
       case we can run phase-extender for each sample there by extending the haplotype. After this phase-extender can be 
       applied recursively on the updated data each time, there by extending the haplotypes for each sample to full chromosome
       length and possibly to to full genome wide length. This is going to be more and more possible due to samples being
       sequenced at higher coverage, increase in the size of paired-end sequence length, availability of large sequence reads
       like pac-bio reads. The method for doing is explained in this link/website ..... ??
       
    B) Another situation when GW/CW phase extension might be possible is when you have at least few samples which have haplotype 
       resolved at GW/CW level. These can include fully phased data like genome matrix file, fully phased VCF data, fully phased
       haplotype blocks. For this the fully phased sample should be provided as one single blocks in the group of sample that is
       fed to phase-extender. Here is an example.... ??
       
   Does phase-extender phase the SNPs that are not phased to any other SNPs?
    For, now phase-extender doesn't have algorithm built into it to do so. I am hoping to add that ASAP.
      
   Does phase-extender impute missing genotypes?
    No, it does not. I have not developed any algorithm so far to impute missing genotypes. Should anyone be interested in 
    developing an algorithm, and have data to do so; I would be more than happy to.  
    
   Does phase-extender use haplotype reference panel?
    No, until now. Though it should be possible in near future.
   
   Does phase-extender use recombination into account?
    No and possibly these feature will be of least importance in phase-extender. Phase-Extender mostly relies on already       
    phased haplotype block from other samples. These haplotype blocks which were phased in other sample but in SOI were used 
    to build transition probabilities. There is strong assumption that variation in other samples are not the result of 
    recombination but only mutation.
    
   Does phase-extender phase rare genotypes?
    Yes, it does. But, the rare genotype should be the part (phased to) of phased haplotype blocks. This is one of the good 
    thing about phase-extender compared to other tools when it come to rare genotype. When several SNPs are phased together 
    to extend the haplotype, rare genotypes are really hard to phase accurately - the reason being the statistical 
    significance of the rare genotype belonging to either two phase state is highly ambigous. But, if the rare genotype is 
    attached to a haplotye block supported by read-back phasing, this makes phasing of rare genotypes most accurate, since 
    the likely hoods are provided by other SNPs that are not rare.
    
   How fast is phase-extender?
    phase-extender is written in python-3, so it is slower than other tools that are built on the top of C, C++ or java.   
    Coming from a pure biology background, learning python was one of the most enduring task I have taking and then 
    building this tool was a big part of my PhD. I have optimized the part of calling VCF file using cyvcf2 (which is on 
    average 4 times faster than old pyVCF module). phase-extender is also optimized for being able to run on multiple 
    threads/process. But, if you are running phase-extender on big genome data and have very large number of samples, and    
    running on laptop I suggest running on one thread, which may be time consuming but will reduce memory burden.
    
   Who helped me during the preparation of this program?
    I have not been very fortunate to surround myself or at least get face to face help from savvy computer programmers.
    But, my heart is very thankful to people behind the web who have made me capable of working this problem out - Thanks 
    to many people on biostars, stackoverflow, seqanswer and google web searches.
    
   Does phase-extender do trio based phase extension?
   
   
   What should be the relatedness of my samples?
   
   
   Has phase-extender been tested?
   
   
   What is differece between phase extender and phase stitcher?
   
   
   
   Should I prepare my haplotype block file onlyusing phaser?
   
   
   
   
## Note:
- Basically this tool takes the paritally phased haplotypes
    - segregates them to maternal vs. paternal population haplotypes (using markov model)
    - stitches the haplotypes to create a genome wide haplotype
    
This program is exclusively designed to phase haplotypes in F1 hybrid for now and should work equally with data generated from genome reseq, RNAseq, WES or RAD given the individual is a hybrid derived from divergent strains and/or populations, but may be used to equally extended to phase F2 hybrids (work on progress).


# Future updates
## Phase SNPs 
## Imputation
## Trio based phasing
## higher order markov chain capabilities
## multiprocessing within chromosome


