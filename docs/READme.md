This tutorial provides a hands on for preparing required `input files` and running `phase-Extender`. **phase-Extender** has been provided with three example sets highlighting possible use cases.

## Part-A: Run phase extender.
### Step 01: Prepare input files.

  **A) Convert VCF to a haplotype file:**
```
python3 vcf_to_table-v3.py --mode VcfToHap --PI PI --PG PG --vcf RBphased_file.vcf --out haploype_file.txt
```
      
  **Note:**
  - If RBphase information is represented by FORMAT fields other than "PI" and "PG", then it will have to be placed accordingly.
  - `"--unphased yes"` switch can be added to `"vcf_to_table-v3.py"` to include the unphased genotypes. This parameter will not affect the phase extension, and is only included to keep the entire data-set in the conversion.
  - run command `python3 vcf_to_table-v3.py --help` for a detailed help section or manual to convert into the requisite tabular format from a VCF or Haplotype panel.
  
<br>

  **B) Convert haplotype reference panel (VCF) to haplotype file:**  
```
python3 vcf_to_table-v3.py --mode RefPanelToHap --PI CHROM --PG GT --vcf RefPanel.vcf --out RefPanel_haploype.txt
```

<br>

### Step 02: Run phase-Extender.
Parameters that are not called are retained at default values. 

  - **Call for help -**
```
python3 phase_extender_v1-final.py --help
```
<br>
    
  - **Example test case 01 (with minimal parameters) -**\
Use data from [example 01](https://github.com/everestial/phase-Extender/tree/master/example01)
```
python3 phase_extender_v1-final.py --input input_haplotype_file.txt --SOI ms02g --lods 10
  
# to write the computed LOD between two blocks to the output file
python3 phase_extender_v1-final.py --input input_haplotype_file.txt --SOI ms02g --writeLOD yes

#Output is stored in directory `ms02g_extended\`
```
<br>

  - **Example test case 02 (multiple cases) -**\
Use data from [example 02](https://github.com/everestial/phase-Extender/tree/master/example02)
```
# There are two steps to run case 02. 
# use 25 heterozygote SNPs (from each block) for ...
  # ... preparing transition matrix between two consecutive blocks
# set lods2 cutoff threshold at 10
# prepare descriptive statistics of the haplotype (initial vs. final) 
python3 phase_extender_v1-final.py --nt 2 --input haplotype_file_test02.txt --SOI ms02g --numHets 25 --culLH maxPd --hapStats yes --lods 10

# Output is stored in directory `ms02g_extended\`.

# to assign output to a different directory
python3 phase_extender_v1-final.py --input haplotype_file_test02.txt --SOI ms02g --output my_test
# Output is stored in directory `my_test\`.

# add "Reference Haplotype Panel" to the run
python3 phase_extender_v1-final.py --nt 2 --input haplotype_file_test02.txt --SOI ms02g --numHets 25 --culLH maxPd --hapStats yes --refHap refPanel_lyrata_test02.txt

# add "bed file" to the run
python3 phase_extender_v1-final.py --nt 1 --input haplotype_file_test02.txt --SOI ms02g --numHets 25 --culLH maxPd --hapStats yes --refHap refPanel_lyrata_test02.txt --bed bed_boundries.bed

# only use samples from "reference panel" to run phase extension
python3 phase_extender_v1-final.py --input haplotype_file_test02.txt --SOI ms02g --useSample refHap

# only use select sample to prepare transition matrix probabilities
# ** it is possible to mix sample names between reference panel and input haplotype file
python3 phase_extender_v1-final.py --input haplotype_file_test02.txt --SOI ms02g --useSample ms01e,ms02g,MA605,Sp76
```
<br>

  - **Example test case 03 (using large dataset) -**\
The data is ReadBackPhased haplotype (chr2 and chr3) from several samples of *A. lyrata* (Mayodan, Spiterstulen and F1 hybrids) used in my study. [Click here for the dataset](https://www.dropbox.com/home/public_shared/phase-Extender_example03) 
    - input file is "**inputHaplotype_chr2n3.txt**".
    - output are in directory "**ms02g_chr2n3_output**" 
```
python3 phase_extender_v1-final.py --nt 2 --input haplotype_file_RBphased_chr2n3.Vars.txt --SOI ms02g --output ms02g_chr2n3_output --numHets 25 --culLH 'maxPd' --hapStats yes --lods 10 --writeLOD yes --addMissingSites yes
```
<br>
<br>

## Part-B: Interpreting phase-Extender's output 

### **1. Interpreting phase extension between blocks (from output of test case 01) -**

**Initial phase state of sample `ms02g` in input file:**
```
ms02g_PI    ms02g_PG_al
6     C|T
6     T|C
6     A|C
6     G|T
4     T|C
4     T|C
4     T|C
4     T|A
```
<br>

**Extended phase state of sample `ms02g` (with computed LOD) in output file:**
```
contig	pos	ref	all-alleles	ms02g_PI	ms02g_PG_al	log2odds
2	15881764	.	.	6	C|T	.
2	15881767	.	.	6	T|C	.
2	15881989	.	.	6	A|C	.
2	15882091	.	.	6	G|T	.
2	15882451	.	.	6	C|T	-73.30333
2	15882454	.	.	6	C|T	-73.30333
2	15882493	.	.	6	C|T	-73.30333
2       15882505        .       .       6       A|T     -73.30333

```
<br>

**What does this tell us?**
  - The lod cutoff parametrized in the command is an absolute value. Our cutoff threshold was `10` and the computed lods between the two blocks is `-73.3033`.
  - Since, the **|computed lods|** > **lods threshold** phase-Extender proceeds with phase extension between blocks with PI (6 and 4).
  - Since, the **compute lods** is a negative value, the phase is exteded in alternate configuration.
  - So, if **|computed lods|** < **lods threshold** the phase state won't extend between two consecutive blocks.
  
<br>
<br>

### **2. Interpreting histogram plots (from output of test case 03) -**\
So, when phase-Extender runs it will join two consecutive haplotype and will increase the size of the haplotype block by number of variants within each haplotype, and also by length of the haplotype. Inversly, the total number of haplotype will decrease. This change can be compared by observing changes in the shape of the histogram (size by number of the haplotype) before vs. after phase extension.

**Histogram of the haplotype size distribution before phase extension (contig 2 & 3)**
![beforephaseextension!](https://github.com/everestial/phase-Extender/blob/master/example03/hap_size_byVar_ms02g_initial.png)

<br>
<br>

**Histogram of the haplotype size distribution after phase extension (contig 2 & 3)**
![afterphaseextension!](https://github.com/everestial/phase-Extender/blob/master/example03/hap_size_byVar_ms02g_final.png)

<br>
<br>

# Recursive application of haplotype phase extension  
  - ***phase-Extender may not be able to prepare a full length phased haplotype in one run.*** This is not the limitation but rather an intended feature in this tool. The reason is to provide flexibility and allow the user to control phase extension. A controllable haplotype extension is largely required for phase extension in emerging research models owing to few genotype samples, absence of reference panel and heterozygosity in the genome.\
  - The main idea is to first run **phase-Extender** with higher `log2Odds cut off` for several samples. Then merge the output of each sample to run another round of **phase extension** with lower `log2Odds cut off`. When applied recursively we reduce the number of haplotypes and increase the length of haplotypes in each run.
  - To achieve genome wide haplotype phasing there should be enough haplotype blocks (from other samples) bridging the haplotype breakpoint in the sample of interest.   
  - So, controllable haplotype extension is a novel feature intended in **phase-Extender**. A full length recursive phase extension is illustrated in this link [phase Extender on recursive mode](some website??) using the bash script. **(coming soon !)**
  
  
  
  
