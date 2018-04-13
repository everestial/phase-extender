
## Step 01:
Prepare required files.

  **A) Convert VCF to haplotype file:**

      python3 vcf_to_table-v3.py --mode VcfToHap --PI PI --PG PG --vcf RBphased_file.vcf --out haploype_file.txt
      
  **Note:
  - If RBphase information is represented by FORMAT fields other than "PI" and "PG", then it can be replaced accordingly.
  - `"--unphased yes"` can be added to `"vcf_to_table-v3.py"` to include the unphased genotypes. This parameter will not affect the phase extension, and is only included to keep the whole data intact.`
  - run command `python3 vcf_to_table-v3.py --help` for more details on VCF and haplotype reference panel to table conversion.
        

  **B) Convert haplotype reference panel (VCF) to haplotype file:**
  
      python3 vcf_to_table-v3.py --mode RefPanelToHap --PI CHROM --PG GT --vcf RefPanel.vcf --out RefPanel_haploype.txt
    

## Step 02:
Run phase-Extender.\
Parameters that are not called are set at default value.

**Call for help -**

    python3 phase_extender_v1-final.py --help
    
**Example test case 01 (with minimal parameters) -**\
Use data from [example 01](https://github.com/everestial/phase-Extender/tree/master/example01)

    python3 phase_extender_v1-final.py --input haplotype_file_test01.txt --SOI ms02g
    
    # write computed LOD between two blocks to the output file
    python3 phase_extender_v1-final.py --input haplotype_file_test01.txt --SOI ms02g --writeLOD yes
    
Output is stored in directory `ms02g_extended\`.


**Example test case 02 (multiple cases) -**\
Use data from [example 02](https://github.com/everestial/phase-Extender/tree/master/example02)

    # use 2 processes, 25 Het sites for transition matrix
    # lods cutoff 10 and output descriptive statistics of the haplotype    
    python3 phase_extender_v1-final.py --nt 2 --input haplotype_file_test02.txt --SOI ms02g --numHets 25 --culLH maxPd --hapStats yes --lods 10
    
    # Output is stored in directory `ms02g_extended\`.
    
    # assign output to different directory
    python3 phase_extender_v1-final.py --input haplotype_file_test02.txt --SOI ms02g --output my_test
    # Output is stored in directory `my_test\`.
    
    # add "Reference Haplotype Panel" to the run
    python3 phase_extender_v1-final.py --nt 2 --input haplotype_file_test02.txt --SOI ms02g --numHets 25 --culLH maxPd --hapStats yes --refHap refPanel_lyrata_test02.txt
    
    # add "bed file" to the run
    python3 phase_extender_v1-final.py --nt 1 --input haplotype_file_test02.txt --SOI ms02g --numHets 25 --culLH maxPd --hapStats yes --refHap refPanel_lyrata_test02.txt --bed bed_boundries.bed
    
    # only use samples from "reference panel" to run phase extension
    python3 phase_extender_v1-final.py --input haplotype_file_test02.txt --SOI ms02g --useSample refHap
        
    # only use select sample from to run phase extension
    # ** it is possible to mix sample names between reference panel and input haplotype file
    python3 phase_extender_v1-final.py --input haplotype_file_test02.txt --SOI ms02g --useSample ms01e,ms02g,MA605,Sp76
    
**Example test case 03 (using large data) -**\
Use the file in this dropbox link
The data are from chr2 and chr3 for A lyrata. 

    # some command
