#!/home/bin/python3

import argparse

# when using cyvcf on python3 activate this import
from cyvcf2 import VCF
from cyvcf2 import VCFReader
import numpy as np

# When using cyvcf on python2 activate this import
#import cyvcf as vcf

import time

print()
print('Checking required modules')
print()

''' Purpose of the program: mine the data from the VCF files. The output file consists of
ReadBackPhased genotypes as PI (block keys) and PG (genotype values). Other GT with no phased
state may be extracted as well.'''


''' function to compute allele frequencies at each position of the VCF file'''
def compute_allele_freq(alt_freq):
    # alt_freq is returned as either float or tuple of floats; so we extract data in more universal format
    if isinstance(alt_freq, tuple):
        alt_freq = [round(x, 3) for x in alt_freq]
        ref_freq = round(1 - sum(alt_freq), 3)
        all_freq = [ref_freq] + alt_freq

    elif isinstance(alt_freq, float):
        alt_freq = round(alt_freq, 3)
        ref_freq = round(1 - alt_freq, 3)
        all_freq = [ref_freq, alt_freq]

    # for the situation when ref and/or alt freq are 'NoneType' and exception arises
    elif alt_freq is None:
        alt_freq = '.'
        ref_freq = '.'
        all_freq = [ref_freq, alt_freq]

    # conver the float instances to string in "freq" data
    all_freq = [str(x) for x in all_freq]
    return all_freq



def main():

    # define argument variables
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", help="sorted vcf file", required=True)
    parser.add_argument("--out", help="name of the output file", required=True)
    parser.add_argument("--unphased", default='no',
                        help="include unphased variants in the output. "
                             "available options: yes or no")

    global args  # creating a global argument variable
    args = parser.parse_args()

    start_time01 = time.time()
    #with open("RBphased_HaploBlock02.txt", 'w') as write_block:
    with open(args.out, 'w') as write_block:
        #vcf_file = VCF('F1.phased_updated_variants.Final02.vcf')
        vcf_file = VCF(args.vcf)
        sample_ids = vcf_file.samples

        print('%i samples found' %len(sample_ids))
        print()

        output_header = ['contig', 'pos', 'ref', 'all-alleles', 'all-freq']

        for name in sample_ids:
            output_header.append(name + '_PI')
            output_header.append(name + '_PG_al')

        # write the header of the output file
        print('\t'.join(output_header), file=write_block)


        chr_on_process = ''

        ''' now, start parsing the VCF file using pyVCF '''
        print('reading the input vcf file %s' % str(args.vcf))
        print()
        for variant in vcf_file:
            contig = str(variant.CHROM)
            pos = str(variant.POS)
            ref_allele = variant.REF
            alt_alleles = variant.ALT
            all_alleles = [ref_allele] + alt_alleles

            alt_freq = variant.INFO.get('AF')

            # pass the "alt_freq" values to a function to compute "all_freq" as string
            all_freq = compute_allele_freq(alt_freq)


            # find which chr is in the process
            if chr_on_process != contig:
                print('contig %s is being processed' %str(contig))
                print()
                chr_on_process = contig


            # this is returned as array. for now it is just easy to mine data from the "all_alleles" variable
            # ** so this may be used in the future
            gt_bas = variant.gt_bases

            # start mining the phased_index (PI) and phased_genotype (PG) values for each sample
            # this outputs a list of PI and PG for all available sample in the VCF file
            PI_values = variant.format('PI')
            PG_values = variant.format('PG')

            # write data to the output file
            # ('2', 15881018, 'G', ['G', 'A', 'C'], [0.0, 1.0])
            write_block.write('\t'.join([contig, pos, ref_allele, ','.join(all_alleles),\
                                        ','.join(all_freq)]))

            # seems like "bcftools" merged variants are missing some of the PI and PG tags.
            if PI_values == None:
                PI_values = np.repeat('.', len(sample_ids), axis=0)

            if PG_values == None:
                PG_values = np.repeat('.', len(sample_ids), axis=0)

            ''' now, start extracting sample level information '''
            for n, pi_value in enumerate(PI_values):
                pg_val = PG_values[n]

                #print('n, pi and pg')
                #print(n, pi_value, pg_val)
                #print()

                if pg_val == '.':
                    pg_allele = '.'

                elif '/' in pg_val:
                    # if unphased variants are not of interest we return them as '.'
                    if args.unphased == 'no':
                        pg_allele = '.'

                    else:
                        if pg_val == './.':
                            pg_allele = '.'

                        else:
                            pg_val = pg_val.split('/')
                            pg_allele = all_alleles[int(pg_val[0])] + '/' + all_alleles[int(pg_val[1])]


                elif '|' in pg_val:
                    pg_val = pg_val.split('|')
                    pg_allele = all_alleles[int(pg_val[0])] + '|' + all_alleles[int(pg_val[1])]


                # ** for the future: in some instances there are "*" in the GT, which we currently replace as "."
                # this may change in the future
                #if '*' in pg_allele:
                    #pg_allele = '.'

                write_block.write('\t' + '\t'.join([pi_value, pg_allele]))
            write_block.write('\n')
            print()

        print('elapsed time: ', time.time() - start_time01)
        print()

if __name__ == '__main__':
    main()









