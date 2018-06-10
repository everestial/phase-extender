#!/usr/bin/python3


## Checking required modules
print("\nChecking and importing required modules: ")

import argparse
import collections
import csv
from decimal import Decimal
from functools import reduce
import itertools
from itertools import islice
from itertools import product
from io import StringIO

from multiprocessing import Pool
import numpy as np
import os
import pandas as pd
import resource
import time
import io
import shutil

# for plotting
import matplotlib.pyplot as plt


def main():

    print()
    ''' printing authorship '''

    print("#######################################################################")
    print("        Welcome to phase-extender version %i       " % 1.0)
    print("  Author: kiranNbishwa (bkgiri@uncg.edu, kirannbishwa01@gmail.com) ")
    print("#######################################################################")
    print()


    print("Loading the argument variables ....")

    ''' Define required argument for interactive mode program. '''
    parser = argparse.ArgumentParser()

    parser.add_argument("--nt",
                        help="number of process to run -> "
                             "the maximum number of processes that can be run at once is the number "
                             "of different chromosomes (contigs) in the input haplotype file.",
                        default=1, required=False)

    parser.add_argument("--input",
                        help="name of the input haplotype file -> "
                             "This haplotype file should contain unique index represented by 'PI' and "
                             "phased genotype represented by 'PG_al' for several samples. ",
                        required=True)

    parser.add_argument("--SOI",
                        help="sample of interest -> "
                             "Name of the sample intended for haplotype extension.",
                        required=True)

    parser.add_argument("--output",
                        help="Name of the output directory. default: 'SOI_extended' ",
                        default='', required=False)

    parser.add_argument("--refHap",
                        help="reference haplotype panel -> "
                             "This file should also contain 'PI' and 'PG_al' values for each "
                             "sample in that haplotype reference panel."
                             "default: empty ",
                        default='', required=False)

    parser.add_argument("--useSample",
                        help="list of samples -> "
                             "Use phase state information only from these samples while running. "
                             "This is intended to provide control over phase-extension by "
                             "allowing select samples from the pool of samples (refHap and/or input). "
                             "This is useful for comparing the results when different individuals, "
                             "populations are used in phase extension process."
                             "Options: 'all','refHap','input','comma separated name of samples'. "
                             "default: 'all' - i.e all the samples in (refHap + input) will be used, "
                             "if refHap is given else all samples only from input file is used.",
                        default='all', required=False)

    parser.add_argument("--bed",
                        help=" bed file -> "
                             "Process the haplotype extension only in this bed regions. "
                             "This is useful to limit haplotype extension only within certain "
                             "regions, like - within genes, exons, introns, QTL, etc. "
                             "default: empty ",
                        default= '', required=False)

    parser.add_argument("--snpTh",
                        help="snp threshold -> "
                             "Minimum number of SNPs required in both haplotype blocks before starting "
                             "phase extension. "
                             "Option: an integer value. "
                             "Default: snpTh = 3 ",
                        default=3, required=False)

    parser.add_argument("--numHets",
                        help="number of heterozygote SNPs -> "
                             "Maximum number of heterozygote SNPs used in consecutive haplotype blocks "
                             "for computing likelyhood of the phase states. "
                             "Option: an integer value. "
                             "Default: numHets = 40 ",
                        default=40, required=False)

    parser.add_argument("--lods",
                        help="log2 of Odds cut off -> "
                             "Cutoff threshold to assign the phase states between consecutive haplotype blocks. "
                             "Option: an integer value "
                             "Default: 5 for parallel configuration (i.e 2^5 = 32 times more likely). ",
                        default=5, required=False)
    #**Note:  positive score above the cutoff values joins two consecutive block in "
                             #"parallel configuration, negative score below the cutoff joins two consecutive block in "
                             #"alternate configuration.

    parser.add_argument("--culLH",
                        help="Cumulative likelhood estimates -> "
                             "The likelhoods for two possible configuration can either be max-sum vs. max-product. "
                             "Default: maxPd i.e max-product. "
                             "Options: 'maxPd' or 'maxSum'. ",
                        default='maxPd', required=False)

    parser.add_argument("--writeLOD",
                        help="write log2 of Odds to the output file-> "
                             "writes the calculated LODs between two consecutive haplotype blocks in the output file. "
                             "Option: 'yes', 'no'. "
                             "Default: no. "
                             "**Note: the 'lods-score' are printed regardless if the "
                             "consecutive blocks are joined or not.",
                        default='no', required=False)  # ** to do - add this feature

    parser.add_argument("--hapStats",
                        help="Computes the descriptive statistics and plots histogram of the haplotype for "
                             "input and extended haplotype. "
                             "Default: 'no'."
                             "Option: 'yes', 'no' .",
                        default='no', required=False)

    parser.add_argument("--addMissingSites",
                        help=" write the lines that have data missing for SOI on the output file. "
                             "Option: yes, no ",
                        default='no', required=False)  # ** to do - add this



    ''' create a global argument variable and declare values of the global variables.
        The values for the global variable are assigned using arguments (args.()) passed by user.'''
    # this is necessary if arguments are declared at several functions using "args.variable" parameter.
    global args;
    args = parser.parse_args()   # .. but keep it as it is.

    global output  # ** now used with "outputdir" - may change in future to control filename in output
    global outputdir
    global use_bed
    global soi
    global snp_threshold
    global num_of_hets
    global lods_cut_off
    global use_sample  # used in coherence with "sample_list"
    global maxed_as
    global writelod
    global time01       # set a start time (used to time computational run for each chromosome vs. all processes)
    global addmissingsites


    print("Assigning values to the global variables ....")
    # advantage of this method (assignment of global variable values at the beginning) ...
      # ... is that we can change the program to non-interactive mode easily and
    # point to specific values/files on this top part. This becomes helpful while debugging.

    time01 = time.time()
    soi = args.SOI
    print('  - sample of interest: "%s" ' %(soi))

    # Assign the output directory
    if args.output == '':
        outputdir = soi + '_extended'
    elif args.output != '':
        outputdir = args.output
    if os.path.exists(outputdir):
        shutil.rmtree(outputdir, ignore_errors=False, onerror=None)
    os.makedirs(outputdir, exist_ok=True)


    # assign number of process to be used
    nt = int(args.nt)  # default, nt = 1
    print('  - using "%s" processes ' % (nt))

    input_file = args.input  # the input haplotype file
    # input_file = 'allele_table_for_phase_extender.txt'
    print('  - using haplotype file "%s" ' % (input_file))

    lods_cut_off = float(args.lods)  # log_of_odds_cut_off, default = 5
    print('  - using log2 odds cut off of "%s" ' % (lods_cut_off))


    # minimum number of SNPs in a haplotype block before it can be phase-extended
    snp_threshold = args.snpTh   # default, snp_threshold = 3
    print('  - each consecutive haplotype block should have minimum of "%s" SNPs ' % (snp_threshold))

    # controls the max number of hetVars to include in likelyhood calculation
    num_of_hets = int(args.numHets)   # default, num_of_hets = 40
    print('  - using maximum of "%s" heterozygote sites in each consecutive blocks to compute '
          'transition probabilities' % (num_of_hets))

    # add argument for max sum vs. max product of likelyhood estimates before calculating the LOD-score
    maxed_as = args.culLH   # default, maxed_as = "*"
    if maxed_as == 'maxSum':
        max_is = 'max sum'
        maxed_as = '+'
    elif maxed_as == 'maxPd':
        max_is = 'max product'
        maxed_as = '*'
    print('  - using "%s" to estimate the cumulative maximum likelyhood of each haplotype configuration '
          'between two consecutive blocks ' % (max_is))


    ##set the required variables related to bed file if bed file is given.
    # use provided bed file to limit phase extension on bed regions/boundries
    if args.bed != "":      # default, use_bed = "no"
        use_bed = 'yes'
        bed_file = args.bed
        print('  - using the bed file "%s" to limit phase extension at the bed boundries '%(bed_file))
    else:
        use_bed = 'no'
        print('  - no bed file is given.')


    # if a haplotype panel is provided then the reference panel can be used as backbone or ..
      # .. meaningful data to aid in phase extension
    if args.refHap != "":  # default, use_refHap = "no"
        use_refhap = 'yes'
        refhap = args.refHap
        print('  - using the reference haplotype panel "%s" ' %(refhap))
    else:
        use_refhap = 'no'
        print('  - no reference haplotype panel is provided ')


    # which samples to use while running phase-extension, default is all (hapRef + input)
    # this variable is updated, later in the pipeline when all sample names are collected
    use_sample = args.useSample  # default, use_sample = "all"


    # print the hapstats to file and also plot histogram
    if args.hapStats == 'yes':  # default, hapstats = 'no' ** ??
        hapstats = 'yes'
        print('  - statistics of the haplotype before and after extension will '
              'be prepared for the sample of interest i.e "%s" ' %(soi))
    else:
        hapstats = 'no'
        print('  - statistics of the haplotype before and after extension will not '
          'be prepared for the sample of interest i.e "%s". '
              '    Only extendent haplotype block will be prepared.' % (soi))



    # write calculated LOD (log2 of odds) of the possible phase state between two consecutive blocks
    # default, writeLOD = "no"
    if args.writeLOD == 'yes':
        writelod = 'yes'
        print('  - LOD (log 2 of odds) for consecutive block will be written to the output file ')
    elif args.writeLOD != 'yes':
        writelod = 'no'
        print('  - LOD (log 2 of odds) for consecutive block will not be written to the output file ')

    # if missing sites are to be added to phase extended file.
    # this will output a new file instead of writing on top of phase extended file.
    if args.addMissingSites == 'yes':
        addmissingsites = 'yes'
    elif args.addMissingSites != 'yes':
        addmissingsites = 'no'



    '''Assign the number of process - this is the optimal position to start multiprocessing ! 
       **note: number of process should be declared after all the global variables are declared,
       because each pool will need to copy the variable/value of global variables. '''
    pool = Pool(processes=nt)  # number of pool to run at once; default at 1

    #### completed assignment of values to argument variables


    ''' Step 01: Read the input file and and prepare two files as output '''
    # a) One output file contains extended phase-states for the sample of interest (soi)
    # b) another output file contains the lines that have missing data for sample of interest

    ''' Step 01 - A: read the input haplotype file and prepare output files '''
    with open(input_file) as input_data, \
            open(outputdir + '/' +"missingdata_" + soi + ".txt", 'w') as missing_data:

        print()
        print()
        print('# Reading the input haplotype file "%s" '%input_data.name)
        print('  - Lines that have data missing for sample "%s" is written in the file "%s" '
              %(soi, missing_data.name))
        #print('extended haplotype data for sample "%s" will be written in the file "%s" '
              #%(soi, update_data.name))


        ''' Step 01 - B: check if "bed file" and "haplotype reference" file are given.
            - then read the "bed file" and "haplotype file" into the memory.
            - these data will be used downstream after reading the haplotype file as "good_data" '''

        # check and load bed file
        if use_bed == 'yes':
            ''' we want to extend phase state only within bed boundries.
                - so, we merge the "input haplotype-file"  with "bed-file". '''
            print()
            print('Reading the bed file "%s" ' % (bed_file))
            print('Phase extension will be limited to regions declared in the bed file')

            my_bed = pd.read_csv(bed_file, sep='\t', names=['CHROM', 'start', 'end'])
            my_bed['CHROM'] = my_bed['CHROM'].astype(str)  # setting CHROM column as string type ..
            #  this is necessary because there has been problem with groupby operations downstream

        else:
            print()
            print('# Genomic bed file is not provided ... ')
            print('  - So, phase extension will run throughout the genome.')


        # check and load "haplotype reference panel"
        if use_refhap == 'yes':
            hap_panel = pd.read_csv(refhap, sep='\t').drop(['REF', 'ALT'], axis=1)
            hap_panel['CHROM'] = hap_panel['CHROM'].astype(str)  # setting CHROM as string type data

            # also find the sample in refHap panel
            hap_panel_samples = find_samples(list(hap_panel.keys()))

        else:
            print()
            hap_panel_samples = []
            print('# Haplotype reference panel is not provided ... ')
            print('  So, phase extension will run using the samples available in the input haplotype file. ')


        ''' Step 01 - C: from the input file
        1) For the soi, remove/store the lines that have PI and allele value as empty (.) as "missing_data"
        2) Rest of the data should have meaningful data for soi, and we store it as "good_data".
           - also prepare the list of samples that will be used in phase extension.
        3) Prepare a list of samples to use while computing transition matrix.


        4) ** move this else where: write the first K1 block - after this we will only need to update K2 block values,
           when reading k1,v1 and k2,v2 consecutive blocks as pairs. '''

        ''' store the good (meaningful) portion of the input_data as variable (Good_Data).
        But, write the missing data as file (** but can also be stored in a variable) '''
        good_data = ''   # empty string

        for lines in input_data:
            if lines.startswith('CHROM'):
                head = lines.rstrip('\n').split('\t')

                # find the index positions of sample-of-interest's PI and PG_allele
                soi_PI_index = head.index(soi + ':PI')
                soi_PG_index = head.index(soi + ':PG_al')


                # update the heading for good_data and missing_data
                good_data += '\t'.join(head) + '\n'
                missing_data.write('\t'.join(head) + '\n')

                continue


            ''' Now, for the soi if PI and PG are missing (i.e, represented by '.') write it into
            "missing_file.txt", else store it as "good_data" '''
            lines = lines.rstrip('\n').split('\t')
            if len(lines) <= 1:  # to control for the last ('\n') line if present
                break


            # separate the lines with missing data (no "PI" or has "PI" but ambiguous SNP like "*")
            if lines[soi_PI_index] == '.' or '.' in lines[soi_PG_index]:
                missing_data.write('\t'.join(lines) + '\n')

            # write the good part of the RB-phased VCF
            elif lines[soi_PI_index] != '.':
                good_data += '\t'.join(lines) + '\n'

            # ** for future: this above code can be modified to include non-phased SNP variants.
            # - just remove "lines[soi:PG_index]" to put SNPs with no PI index inside "good_data"

        print()
        print('# Filtered the lines that have data missing for sample "%s"; check the file "%s" '
              %(soi, missing_data.name))
        print('  - Loaded read-backphased variants onto the memory')


        ''' Step 01 - D: Prepare the samples to use the data from. '''
        ''' Prepare a list of tuples of samples (PI, PG_al) from the input data and update it as needed.
            - **Note: the sample list should always include the soi (sample of interest)
                - this is done to include observation from soi rather than introducing a pseudo count
                  when transition is missing from some observation (n to m). '''
        samples = head.copy()
        sample_list = find_samples(samples)  # returns data from "input haplotype file"

        # update the names in "sample_list" if other samples are requested by the user:
        if use_sample == "" or use_sample == 'input':
            sample_list = sample_list

        # use all the samples from hapRefPanel and input samples
        elif use_sample == 'all':
            sample_list = sample_list + hap_panel_samples

        elif use_sample == 'refHap':
            sample_list = hap_panel_samples + [(soi + ":PI", soi +":PG_al")]  # add the self sample name to account ..
            # .. for missing observations instead of using pseudo count

        # if specific select samples are of interest, split the sample names and then prepare ..
        # .. the list of tuples of sample "PI" and "PG_al"
        else:
            sample_list = use_sample.split(',')
            sample_list = [((x + ':PI'), (x + ':PG_al')) for x in sample_list] + \
                          [(soi + ":PI", soi +":PG_al")]
        # print()
        #print('sample listing to use')
        #print(sample_list)



        ''' Step 02: pipe the data into "pandas", then:
            A) group the data by "contig" which helps in multiprocessing/threading.
              A - optional: if "bed regions" are given add the bed_regions boundries as "start_end"
            B) within each group, group again by "PI keys" of soi and then sort by
               minimum "POS" value for each "PI key"
            C) then pipe the data within each "PI key" for phase-extension computation.'''

        ''' Step 02 - A : read good part of the data into "pandas" as dataframe.'''
        good_data = pd.read_table(StringIO(good_data), delimiter='\t')
        good_data['CHROM'] = good_data['CHROM'].astype(str)  # setting CHROM as string type data # this is necessary
        # to maintain proper groupby downstream

        # ** only if "good_data" is desired as text output
        #pd.DataFrame.to_csv(good_data, 'good_data_test.txt', sep='\t', header=True, index=False)

        ''' Step 02 - A (add on - i) ** merge reference haplotype if provided '''
        print()
        if use_refhap == "yes":
            # update the "good_data" (i.e, haplotype data)
            print('Merging input haplotype data with data from the hap-reference panel')

            good_data = pd.merge(good_data, hap_panel, on=['CHROM', 'POS'],
                                     how='left').fillna('.')
            good_data.sort_values(by=['CHROM', 'POS'], inplace=True)

            # if haplotype and reference panel merged lines are desired
            #pd.DataFrame.to_csv(good_data, 'hap_and_refPanel_merged.txt', sep='\t',
                                #header=True, index=False)
            del hap_panel

        else:
            print('# Haplotype reference panel is not provided....\n'
                  '  - Only using the samples in the input ("%s") data.' %(input_data.name))


        ''' Step 02 - A (add on - ii) ** merge bed-regions if provided to limit phase extension
                                         and group the data by "contig". '''
        print()
        if use_bed == 'no':
            # group data only at "contig" level, keep the sort as it is
            print('# No bed file is given ... ')
            print('  - So, grouping the haplotype file only by chromosome (contig)')

            good_data_by_group = good_data.groupby('CHROM', sort=False)

        elif use_bed == 'yes':
            print('# Merging the bed boundries from "%s" with the input haplotype file ... "%s" '
                  % (bed_file, input_data.name))

            # merge/intersect the "bed regions" and "haplotype file"
            # then groupy "contig" and "bed regions" by passing it to function "merge_hap_with_bed()"
            good_data_by_group = merge_hap_with_bed(my_bed, good_data)
            # ** for future: we can also run multiprocessing while merging "hap file" with "bed regions"
            del my_bed


        ''' Step 02 - A (**add on - iii):
            - Write the initial haplotype data.
            - Compute the statistics of the initial phased file for SOI if required '''

        print()
        print('# Writing initial haplotype for sample "%s" in the file "%s" '
              %(soi, 'initial_haplotype_' + soi + '.txt'))

        # select the colums of interest
        initial_haplotype = good_data[['CHROM', 'POS', 'REF', 'all-alleles', soi + ':PI', soi + ':PG_al']]. \
            sort_values(by=['CHROM', 'POS'])

        # write this initial haplotype to a file
        pd.DataFrame.to_csv(initial_haplotype, outputdir + '/' + 'initial_haplotype_' + soi + '.txt',
                            sep='\t', header=True, index=False)


        if hapstats == 'yes':
            print('  - Computing the descriptive statistics of the haplotype data before phase extension')

            # pipe the data to a function to compute haplotype statistics
            compute_haplotype_stats(initial_haplotype, soi, prefix='initial')
        else:
            print('  - Proceeding to phase-extension without preparing descriptive statistics of initial haplotype state.')



        ''' Step 02 - B: - Split the data (grouped by chromosome (contig) values.
                         - Store data in disk or memory.
                         - Multiprocess each chunks separately '''
        print()
        print('# Starting multiprocessing using "%i" processes ' %(nt))

        # ** new method: create a folder to store the data to disk (rather than memory)
        # ** (see old method for comparison)
        if os.path.exists('chunked_Data_' + soi):
            shutil.rmtree('chunked_Data_' + soi, ignore_errors=False, onerror=None)
        os.makedirs('chunked_Data_' + soi + '/', exist_ok=True)


        ''' Step 02 - B (i)'''

        ################### old method - ** if possible reuse this method in future.
        # take the large dataframe that is grouped by contig and ..
        # .. keep chunks of dataframes as as OrderedDict(list of (keys, Dataframe object))
        #df_list = collections.OrderedDict()
        ########################################

        # new method - storing data to disk
        for chr_, data_by_chr in good_data_by_group:
            pd.DataFrame.to_csv(data_by_chr, 'chunked_Data_' + soi + '/' + soi + ':' + str(chr_),
                                sep='\t', index=False, header=True)


        # clear memory - does it do it's job ** ??
        initial_haplotype = None; good_data = None; input_file = None
        good_data_by_group = None; samples = None; input_data = None
        data_by_chr = None
        del initial_haplotype, good_data, input_file, good_data_by_group, samples, input_data, data_by_chr

    ''' Now, pipe the procedure to next function for multiprocessing (i.e Step 02 - C) '''
    multiproc(sample_list, pool, hapstats)

    # remove the chunked data folder ** (this can be retained if need be)
    #shutil.rmtree('chunked_Data_' + soi, ignore_errors=False, onerror=None)

    print('End :)')


def multiproc(sample_list, pool, hapstats):

    print()
    ''' Step 02 - C: Start, multiprocessing/threading - process each contig separately. '''
    ''' Step 02 - C (ii) : create a pool of process and then feed the list of dataframes to another
        function "groupby_and_read()" to run phase extension.
        - After the data is passed into that function; the steps (02-C to 9) are covered there'''

    path='chunked_Data_' + soi   # create a temp folder
    file_path = [(item, sample_list) for item in list(os.path.join(path, part) for part in os.listdir(path))]

    ## ** to do: Add "sort" method in "file_path" to read data in order. This way we can save ..
      # time/memory while doing sorting within pandas dataframe.
      # This sort method is available in "phase-Stitcher"


    result = pool.imap(groupby_and_read, file_path)
    pool.close()
    pool.join()
    pool.terminate()


    #############################################
    # old method of pooling - directly from memory
    #p = Pool(processes=nt)  # number of pool to run at once; default at 1
    # #result = pool.imap(groupby_and_read, list(df_list.items()))
    #result = pool .imap(groupby_and_read, [it + (sample_list, ) for it in list(df_list.items())])
    # ** for future: use iterator, yield or generators to reduce memory burden while multiprocessing.
    ###########################


    print()
    print("Completed haplotype extension for all the chromosomes."
          "time elapsed: '%s' " %(time.time() - time01))
    print('Global maximum memory usage: %.2f (mb)' % current_mem_usage())
    print("Merging dataframes together ....." )
    print()


    ''' Step 10: merge the returned result (extended haplotype block by each contigs)
        together and write it to an output file. '''
    result_merged = pd.concat(result).sort_values(by=['CHROM', 'POS'], ascending=[True,True])
    #result_merged = pd.DataFrame(result_merged).sort_values(by=['contig', 'pos'], ascending=[1,1])

    # write data to the final output file
    pd.DataFrame.to_csv(result_merged, outputdir + '/' + 'extended_haplotype_'+ soi +'.txt',
                        sep='\t', index=False)

    print('Extended haplotype data for sample "%s" is written in the file "%s". '
          %(soi, 'extended_haplotype_'+ soi + '.txt'))


    ''' Step 11: now, compute the descriptive stats of haplotype after phase extension. '''
    if hapstats == 'yes':
        print()
        print('Computing the descriptive statistics of the extended haplotype file.')
        compute_haplotype_stats(result_merged, soi, prefix='final')

    else:
        print('Skipping the preparation of descriptive statistics of extended haplotype.')

    print()
    print('Run is complete for all the chromosomes (contigs)')


    ''' Step 12: if phase extended data is to be merged with non phased SNPs and missing sites.
        - we can use this data to run another round of phase extension of singletons (SNP and Indels).
        - it also helps to control recursive phase extension. '''
    print()
    print('writing singletons and missing sites to extended haplotype')
    if addmissingsites == 'yes':
        missed_data_toadd = pd.read_csv(outputdir + '/' + "missingdata_" + soi + ".txt", sep='\t',
                                        usecols=['CHROM', 'POS', 'REF', 'all-alleles',
                                                 soi + ':PI', soi + ':PG_al'])
        missed_data_toadd['CHROM'] = missed_data_toadd['CHROM'].astype(str)

        # "result_merged" only contains RBphased data for SOI
        new_df = pd.concat([result_merged, missed_data_toadd]).sort_values(by=['CHROM', 'POS'], ascending=[True, True])
        new_df = new_df[['CHROM', 'POS', 'REF', 'all-alleles',	soi + ':PI', soi + ':PG_al']]

        # write data to the final output file
        pd.DataFrame.to_csv(new_df, outputdir + '/' + 'extended_haplotype_' + soi + '_allsites.txt',
                            sep='\t', index=False, na_rep='.')

        del new_df

    del result_merged





def groupby_and_read(file_path):

    print()
    good_data_by_contig = open(file_path[0], 'r')
    chr_ = good_data_by_contig.name.split(':')[-1]
    sample_list = file_path[1]
    contigs_group = pd.read_csv(StringIO(good_data_by_contig.read()), sep='\t')


    ''' After doing groupby (by chromosome) we pipe in data, for each chromosome '''

    time_chr = time.time()  # to time the process in each chromosome
    print('## Extending haplotype blocks in chromosome (contig) %s' %(chr_))

    del good_data_by_contig


    # if phase extension is to be limited to provided bed-regions
    if use_bed == 'yes':
        # start another level of grouping of the data by bed interval key
        # (start-end of the bed is used as unique key for grouping)
        print('Splitting the contig by bed regions')
        good_data_by_bed = contigs_group.groupby('start_end', sort=False)


        # clear memory
        contigs_group = None
        del contigs_group


        ''' store the dataframes split by "start_end" key in "df_list_by_bed".
            Then pass this to phase-extension process. '''

        # store dataframe/s outside bed region
        # ** - for future (convert this into collection.OrderedDict of ..
        # .. dataframes if data needs to be stored by multiple void blocks)
        data_by_bed_void = pd.DataFrame()


        # to store dataframes that are split by bed-regions but are within boundry of interest
        df_list_by_bed = []

        for bed_id, data_by_bed in good_data_by_bed:
            if 'void' in bed_id:
                # write data from these columns/row as a variable
                data_by_bed_void = data_by_bed[
                    ['CHROM', 'POS', 'REF', 'all-alleles', soi + ':PI', soi + ':PG_al']]

                # ** write void lines separately (if needed) - optimize in the future.
                #pd.DataFrame.to_csv(data_by_bed_void, 'haplotype_file_excluded_regions.txt',
                                    #sep='\t', header=None, index=False)

            else:
                # now pass this dataframe to a function "process_consecutive_blocks()"
                  # i.e, Step 02-D - that reads two consecutive keys at once
                # see **1, for future - there is possibility to multiprocess (phase-extension) this method
                phase_extended_by_bed = process_consecutive_blocks(
                    data_by_bed, soi, chr_, snp_threshold,
                    sample_list, num_of_hets, lods_cut_off)

                # append the dataframes that are returned as "phase_extended" haplotype block
                # .. thus creating a list of dataframes
                df_list_by_bed.append(phase_extended_by_bed)

        # append the regions of the dataframe that was restricted by bed-file
        df_list_by_bed.append(data_by_bed_void)

        # concat all the dataframes in the list into one big dataframe
        # this merges all the datraframe by bed-file within each chromosome/contig together
        phase_extended_by_chr = pd.concat(df_list_by_bed)


        # clear memory
        data_by_bed_void = None; df_list_by_bed = None; phase_extended_by_bed=None


        print('  - Phase-extension completed for contig "%s" in %s seconds' % (chr_, time.time() - time_chr))
        print('  - Worker maximum memory usage: %.2f (mb)' % (current_mem_usage()))
        print()

        return phase_extended_by_chr.sort_values(by=['POS'])


        ####################################   **********************
        # **1 for future - multiprocess phase extension for each bed-regions separately.
        # use pool.map() or starmap
        # create another pools of process for each chunks of dataframes
        #p2 = Pool(1)   # processing one dataframe at a time; but multiple possible
        # pass it to phase extension
        #result_by_bed = p2.map(process_consecutive_blocks, list(df_list_by_bed.values()))
        #p2.close()
        # "result_by_bed" has stored extended haplotype block for each bedregions as list.
        # so, merge those list of dataframe
        # result_by_bed_for_each_chr = pd.concat(result_by_bed)
        # merge this again to the dataframe that were excluded due to bed-regions restriction
        #result_by_bed_for_each_chr
        #return result_by_bed_for_each_chr
        #####################################  ******************************



    # if haplotype are extended at chromosome/contig level ..
    # ..just pass whole contig/chromosome at once
    else:
        phase_extended_by_chr = process_consecutive_blocks(
            contigs_group, soi, chr_, snp_threshold,
            sample_list, num_of_hets, lods_cut_off)

        # clear memory
        contigs_group = None;
        del contigs_group


    print('  - Phase-extension completed for contig "%s" in %s seconds' %(chr_, time.time()-time_chr))
    print('  - Worker maximum memory usage: %.2f (mb)' % (current_mem_usage()))
    print()

    # return the phase-extended haplotype back to the pool-process ..
    # .. which will be pd.concatenated after all pools are complete
    return phase_extended_by_chr.sort_values(by=['POS'])


''' now, read the two consecutive blocks to do phase extension. '''
def process_consecutive_blocks(contigs_group, soi, chr_, snp_threshold,
                               sample_list, num_of_hets, lods_cut_off):

    #print()
    print('  - Grouping the dataframe using unique "PI - phased index" values. ')


    ''' Step 02 - D: group dataframe again by "PI keys" of soi and then
       sort by minimum "POS" value for each "PI key".
       - This sorting is necessary because sometimes "haplotype blocks" are like 3-3-3-3  5-5-5  3-3-3-3
          - i.e there are small RBphased blocks within the boundry of larger RBphased block.
          - Not, sure what is causing this (prolly sampling difference of large vs. small chunks in PE reads)
          - This problem should go away in first round of haplotype-extension'''

    contigs_group = contigs_group. \
        assign(New=contigs_group.groupby([soi + ':PI']).
               POS.transform('min')).sort_values(['New', 'POS'])


    ''' Step 03: Now, start reading the "contigs_group" for haplotype-extension.
    A) Store the data as dictionary with 'header' values as keys. Some keys are: CHROM, POS, sample (PI, PG within sample),
       etc ... Then group the dictionary using unique "PI" values as 'keys' for grouping.
        Note: This dict-data should contain information about two adjacent haplotype blocks that needs extending.
        In this example I want to extend the haplotypes for "sample ms02g" which has two blocks 6 and 4.
        So, I read the PI and PG value for this sample. Also, data should store with some unique keys.
    B) Iterate over two consecutive Haplotype-Blocks at once.
        Note: While iterating over two blocks, initially we write the very first block of the "contig". With this
        method, now when we iterate over two consecutive blocks we can only update and write the second block.
        '''

    # covert pandas dataframe back to text like file before converting it into dictionary.
    contigs_group = pd.DataFrame.to_csv(contigs_group, sep='\t', index=False, header=True)

    ''' Step 03 - A : read the data with header as keys and groupby using that "keys" '''
    phased_dict = csv.DictReader(StringIO(contigs_group), delimiter='\t')
    phased_grouped = itertools.groupby(phased_dict, key=lambda x: x[soi + ':PI'])

    ''' Since the dictionary isn't ordered, we return the order using OrderedDictionary '''
    # ** for future: there is room for improvement in here (memory and speed)
    grouped_data = collections.OrderedDict()
    for key, grp in phased_grouped:
        grouped_data[key] = accumulate(grp)

    ''' Clear memory '''
    del phased_dict
    del phased_grouped
    del contigs_group


    #print()
    print('  - Starting MarkovChains for contig %s' % chr_)
    ''' Step 03 - B : now pipe the data for phase extension '''
    ''' Step 03 - B : And, iterate over two consecutive Haplotype-Blocks at once. This is done to obtain all
    possible Haplotype configurations between two blocks. The (keys,values) for first block is represented as
    k1,v2 and for the later block as k2,v2. '''



    ''' Step 03 - B (i): Before running consecutive blocks, we write data from the very first block to the file.
    Reason : Before we start computing and solving the haplotype phase state, we plan to write the
    data for very first block (k1, v1). So, after that, we can solve the relation between two consecutive..
    .. blocks but only write data from 2nd block each time - based on what relation comes out. '''
    very_first_block = [list(grouped_data.items())[0]]

    if len(list(grouped_data.items())) == 1:
        print('there is only one block, so skipping phase extension')

    # write header of the extended phase-block
    extended_haplotype = '\t'.join(['CHROM', 'POS', 'REF', 'all-alleles', soi +':PI', soi +':PG_al']) + '\n'

    if writelod == 'yes':  # add extra field if desired by user
        extended_haplotype = extended_haplotype.rstrip('\n') + '\tlog2odds\n'
        log2odds = ''

    # write data/values from very first block.
    for k1, v1 in very_first_block:
        for r1, vals in enumerate(v1[soi + ':PI']):
            new_line = '\t'.join([v1['CHROM'][r1], v1['POS'][r1], v1['REF'][r1], v1['all-alleles'][r1],
                                  v1[soi + ':PI'][r1], v1[soi + ':PG_al'][r1]]) + '\n'
            if writelod == 'yes':
                new_line = new_line.rstrip('\n') + '\t.\n'

            extended_haplotype += new_line


        #print('very first block end\n\n')  # marker for debugging


    ''' Step 03 - B (ii):  Starting MarkovChains.
            Now, read data from two consecutive blocks at a time.
            Note: At the end of computation write the data only from each k2 block. No need to write the data
            from k1 block of each iteration because it was written in earlier loop.'''

    ''' Step 03 - B (ii - 1) Create empty "checker variables".
        Note: This checker variables (actually multi-level boolean logic) help to carryover information from
        ealier iteration of a for-loop - i.e identify if the values from later block i.e k2, v2 were phased to
        to earlier block (k1, v1) in "parallel" vs. "alternate configuration".
        - If two consecutive blocks are phased, k2_new is now assigned k1 from earlier block; else (if not phased)
          k2_new stays empty ('').
        - So, the role of flipped variable is to keep information if k2,v2 were phased straight vs. alternate
          compared to k1, v1 in the earlier run. These checker-variables are crucial to keep the proper phase-state
          in the output file.'''

    # start checker variables
    k2_new = ''   # updates the index of k2 for each k1,v1 ; k2,v2 run
    flipped = ''  # boolean logic to check and store if the phase state flipped during extension


    ''' Step 03 - B (ii - 2): Now, read two consecutive blocks at a time'''
    for (k1, v1), (k2, v2) in zip(grouped_data.items(), islice(grouped_data.items(), 1, None)):

        ''' Step 03 - B (ii - 2-A): iterate over the first Haplotype Block, i.e the k1 block.
        The nucleotides in the left of the phased SNPs are called Block01-haplotype-A,
        and similarly on the right as Block01-haplotype-B. '''

        # iterate over the first Haplotype Block, i.e the k1 block and v1 values
        hap_block1a = [x.split('|')[0] for x in v1[soi + ':PG_al']]  # the left haplotype of block01
        hap_block1b = [x.split('|')[1] for x in v1[soi + ':PG_al']]

        # iterate over the second Haplotype Block, i.e the k2 block and v2 values
        hap_block2a = [x.split('|')[0] for x in v2[soi + ':PG_al']]
        hap_block2b = [x.split('|')[1] for x in v2[soi + ':PG_al']]


        ''' Step 03 - B (ii - 2-B) : Create possible haplotype configurations for "forward markov chain".
        Possible haplotype Configurations will be, Either :

        1) Block01-haplotype-A phased with Block02-haplotype-A,
            creating -> hapb1a-hapb2a, hapb1b-hapb2b '''
        ''' First possible configuration '''
        hapb1a_hapb2a = [hap_block1a, hap_block2a]
        hapb1b_hapb2b = [hap_block1b, hap_block2b]

        ''' Or, Second Possible Configuration
        2) block01-haplotype-A phased with Block02-haplotype-B
            creating -> hapb1a-hapb2b, hapb1b-hapb2a '''
        hapb1a_hapb2b = [hap_block1a, hap_block2b]
        hapb1b_hapb2a = [hap_block1b, hap_block2a]

        ''' Step 03 - B (ii - 2-C) :
        Create possible haplotype configurations for "reverse markov chain"
        - reverse markov chain are added to increase the confidence in likelyhood estimation. '''

        # switch the keys values for reverse markov chain
        k1_r = k2
        k2_r = k1
        v1_r = v2
        v2_r = v1

        # switch the haplotype positions for preparing the reverse markov chains
        hapb1a_hapb2a_r = [hapb1a_hapb2a[1], hapb1a_hapb2a[0]]
        hapb1b_hapb2b_r = [hapb1b_hapb2b[1], hapb1b_hapb2b[0]]

        hapb1a_hapb2b_r = [hapb1a_hapb2b[1], hapb1a_hapb2b[0]]
        hapb1b_hapb2a_r = [hapb1b_hapb2a[1], hapb1b_hapb2a[0]]


        ################################# - inactive for now- can be used for adding SNP phasing later on.
        ''' skip if one of the keys has no values - this is redundant ?? - keep it for just in case situation
        ** can also be used in the future if we want to phase the SNPs that have no assigned 'PI' values,
        i.e the "PI" will be "." '''
        if k1 == '.' or k2 == '.':
            for xi in range(len(v2[soi + ':PI'])):
                new_line = '\t'.join([v2['CHROM'][xi], v2['POS'][xi], v2['REF'][xi], v2['all-alleles'][xi],
                                      k2, hapb1a_hapb2a[1][xi] + '|' + hapb1b_hapb2b[1][xi]]) + '\n'
                if writelod == 'yes':
                    new_line = new_line.rstrip('\n') + '\t.\n'

                extended_haplotype += new_line

            # update the values of checker variables
            k2_new = ''
            flipped = ''

            continue   # to next consecutive blocks
        ######################################################


        ''' Step 03 - C : Set the threshold for the minimum number of SNPs required in haplotype block
        before continuing phase extension. '''
        ''' If all the data in soi, in either v1 or v2 are SNPs below a certain threshold we just write
        the data and continue. i.e say if a Haplotype-Block is composed of only 2 SNPs it will be less
        reliable to extend the phase-state.
        - So, this step can also be used to control the minimum number/size of the haplotypes that is required
        before it can be phase-extended.
        - by default the minimum number of SNPs (exclusive) in the soi haplotype is set to 3.
        - If minimum requirement isn't met just skip extending the phase and write it to a file and continue. '''
        number_of_snp_in_soi_v1 = len([x for x in v1[soi + ':PG_al'] if len(x) == 3])
        number_of_snp_in_soi_v2 = len([x for x in v2[soi + ':PG_al'] if len(x) == 3])

        # print('number of SNPs: ', NumSNPsInsoi_v1, NumSNPsInsoi_v2)
        if number_of_snp_in_soi_v1 < snp_threshold \
                or number_of_snp_in_soi_v2 < snp_threshold:
            for xth, vals in enumerate(v2[soi + ':PI']):
                new_line = '\t'.join([v2['CHROM'][xth], v2['POS'][xth], v2['REF'][xth], v2['all-alleles'][xth],
                                      k2, hapb1a_hapb2a[1][xth] + '|' + hapb1b_hapb2b[1][xth]]) + '\n'
                if writelod == 'yes':
                    new_line = new_line.rstrip('\n') + '\t.\n'

                extended_haplotype += new_line

            # update values of the checker variables
            # this is important, so previous k2 and flip state doesn't get carried over without purpose
            k2_new = ''
            flipped = ''

            continue   # to next consecutive blocks


        ''' Step 04: For the consecutive blocks that pass the thresholds (SNP number, have PI != '.', etc.,
        pipe the data (k1, v1 ; k2, v2) to a defined function for computation of forward and reverse
        markov chain transition probabilities for these two consecutive blocks (k1, v1; k2, v2) '''

        #### for forward chain   ########
        # ** set "orientation=reversed" to compute transition ..
          # .. from the lower tip of former block with upper tip of later block
          # .. this helps in using the closest genomic position between consecutive blocks thus ..
          # .. downsizing the effects created by recombination.
        lhfc_f, lhsc_f = \
            compute_maxLh_score(soi, sample_list, k1, k2, v1, v2, num_of_hets,
                                      hapb1a_hapb2a, hapb1b_hapb2b,
                                      hapb1a_hapb2b, hapb1b_hapb2a, orientation=reversed)

        #### for reverse chain   ########
        # set "orientation=lambda..." just passes a null value keeping orientation as it is.
        lhfc_r, lhsc_r = compute_maxLh_score \
            (soi, sample_list, k1_r, k2_r, v1_r, v2_r, num_of_hets,
             hapb1a_hapb2a_r, hapb1b_hapb2b_r,
             hapb1a_hapb2b_r, hapb1b_hapb2a_r, orientation=lambda x: x)



        ''' Step 05-06 are inside the function "compute_maxLh_score()". The values
        (lhfc_f, lhsc_f, lhfc_r, lhsc_r) returned from this function is then used in Step 07. '''



        ''' Step 07 :  previous (Step 06) returns the likelyhoods and/or LODs score for both "parallel"
        and alternate configurations (for both forward and reverse algorithm).
        - We now extend the phase states by comparing LODs score against  cutoff-values.'''

        ''' Step 07 - A(i): calculate the average of the likelyhoods, odds and then log2 of odds. '''
        # average of the likelyhooods for first vs. second configuration
        # (from both forward and reverse algorithm)
        # ** note: "maxed_as" variable doesn't apply here, because maxLH using forward vs. reverse ..
          # .. are just re-estimates. So, we simply take and average on both "maxSum" and "maxPd"
        avg_lhfc = Decimal(lhfc_f + lhfc_r) / 2
        avg_lhsc = Decimal(lhsc_f + lhsc_r) / 2

        # therefore, odds of first_vs_second_configuration is
        odds_fc_vs_sc = avg_lhfc / avg_lhsc

        ''' Step 07 - A(ii) : convert the likelyhoods to odds-ratio and then logOf 2 odds'''
        lods2_score_1st_config = Decimal(odds_fc_vs_sc).ln()/(Decimal('2').ln())
        lods2_score_2nd_config = (-lods2_score_1st_config)

        #print('logOdds')  # marker for debugging
        #print(lods2_score_1st_config)


        ''' Step 07 - B : pipe the LOD scores and write the phase state between two consecutive blocks.
                - use "lods cutoff" to decide on phase extension
                - and then store, write it to files.
         ** We can also use accumulation of this stage to run histogram building at later stage.
            - that acculated "extended_haplotype" can be all written at once - this is important while multiprocessing. '''
        k2_new, flipped, extended_haplotype = extend_phase_state(soi, k1, k2, v1, v2, k2_new, flipped,
                                             lods2_score_1st_config, lods_cut_off, extended_haplotype,
                                             hapb1a_hapb2a, hapb1b_hapb2b)

        ''' Now, go to Step 08, function "extend_phase_state" '''
        # this process udates the data in "extended_haplotype" recursively on the for-loop

    # finally return the extended haplotype as pandas dataframe
    phase_extend = extended_haplotype

    del extended_haplotype

    return pd.read_table(io.StringIO(phase_extend), sep='\t')



def merge_hap_with_bed(my_bed, good_data):

    print('Extracting bed regions and position of the haplotype file ... ')

    c1 = my_bed.CHROM.values
    s1 = my_bed.start.values
    e1 = my_bed.end.values
    c2 = good_data['CHROM'].values
    pos2 = good_data.POS.values

    # now, find the intersecting positions (between haplotype ("pos" values) and bed regions).
    overlap = (
        (c2[:, None] == c1) & ((pos2[:, None] >= s1) & (pos2[:, None] <= e1)))

    i, j = np.where(overlap)


    '''intersect the haplotype file with bed file '''

    print('Intersecting bed regions with haplotype file, and '
          'creating new haplotype boundries for each contig ... ')

    df_bed__hap_interesect = pd.DataFrame(
        np.column_stack([good_data.values[i], my_bed.values[j]]),
        columns=good_data.keys().append(my_bed.keys()))

    # drop any duplicate columns
    df_bed__hap_interesect = df_bed__hap_interesect.T.drop_duplicates().T

    df_bed__hap_interesect['start_end'] = pd.DataFrame(
        df_bed__hap_interesect.apply(lambda x: str(x.start) + '-' + str(x.end), axis=1))

    # only keep "contig, pos and intersection ("start-end")"
    df_bed__hap_interesect = df_bed__hap_interesect[['CHROM', 'POS', 'start_end']]
    # if regions of intersection are of interest
    #pd.DataFrame.to_csv(df_bed__hap_interesect, 'df_bed and hap intersect.txt', sep='\t', index=None, header=True)


    # now, merge the df(good data) with the df (intersected bed and hap) to create updated
      # haplotype file.
    # Any part that is intersected is assigned unique-keys of "start-end" values. The non-intersected
      # part are filled with "NA"
    data_frames = [good_data, df_bed__hap_interesect]

    print('Merging the bed file with the input haplotype file, '
          'and marking the bed boundries.')

    # fill the non-merging lines with "void"
    good_data_update = reduce(lambda left, right: pd.merge(
        left, right, on=['CHROM', 'POS'], how='outer'), data_frames).fillna('na')
    # ** for future: assigning an unique identifier like "NA01, NA02" for non-intersecting lines as blocks.


    # if updated "haplotype file" with intersected "haplotype file" and "bed file" is of interest
    #pd.DataFrame.to_csv(good_data_update, 'df_bed and hap merged.txt', sep='\t', index=None,
                        #header=True)


    print('Grouping the haplotype file by chromosome (contigs) .... ')
    good_data_by_group = good_data_update.groupby('CHROM', sort=False)

    # clear memory
    del c1, c2, s1, e1, pos2, i, j, good_data, good_data_update, df_bed__hap_interesect, my_bed

    return good_data_by_group



''' function to compute transition probs between two consecutive blocks. '''
def compute_maxLh_score(soi, sample_list, k1, k2, v1, v2, num_of_hets,
                              hapb1a_hapb2a, hapb1b_hapb2b,
                              hapb1a_hapb2b, hapb1b_hapb2a,orientation):



    #print()
    ''' Step 05 : Start preping the first order markov transition matrix and probabilities for the
    consecutive haplotypes between two blocks.
    - To store the sum-values of the product of the transition probabilities.
    These sum are added as the product-of-transition is returned by nested for-loop;
    from the loop "for m in range(....)" '''

    # empty variable for storing "sum of the product of transition probabilities" between two blocks.
    # we can either "max sum" or "max product" the likelyhoods, default is "max-product".
    # cul -> cumulative likelyhood estimates (either summed (+) or product(*))
    if maxed_as == '+':
        cul_of_pt_hapb1a_b2a = 0
        cul_of_pt_hapb1b_b2b = 0

        cul_of_pt_hapb1a_b2b = 0
        cul_of_pt_hapb1b_b2a = 0

    elif maxed_as == '*':
        cul_of_pt_hapb1a_b2a = 1
        cul_of_pt_hapb1b_b2b = 1

        cul_of_pt_hapb1a_b2b = 1
        cul_of_pt_hapb1b_b2a = 1


    '''Step 05 - A : We start calculating the nucelotide probabilities (aka "emission probabilities"),
    and then transition probabilities from each level of v1[soi] to each level of v2[soi].
    - we are using reversed() to actually make haplotype configurations starting from the lower tip (i.e last n-th)
    element of former block with the upper tip (first m-th) element of the later block.
    - using reverse() also helps us control if we just want trans-prob calculation between first 200 Het_SNPs.
        - the reason being that when haplotype becomes larger and large it introduces computation burden.
        - also as the haplotype becomes larger the chances of LD at the tips of the haplotypes increase,
            thereby changing the likelyhoods in another direction.
        - orientation : this parameter helps us to control which tip (top vs. bottom) of the former haplotype is
            used with which tip of later haplotype when preparing transition matrix.
            # make this explanation figurative in thesis/dissertation. '''


    # to control the number of Heterozygote site we using for phase-extenstion
    num_of_het_site = 0

    # n-ranges from block01
    for n in orientation(range(len(v1[soi + ':PI']))):

        ''' Skip computation if n'th item is InDel or has "*" allele in the genotype.
        - These indels are only removed from computation part but are added in the final output file.
        - The Indels are phased based on which SNPs they hitchhike with. '''
        if len(v1[soi + ':PG_al'][n]) > 3 \
            or '*' in v1[soi + ':PG_al'][n]:
            continue

        ''' only use certain number of HetVars to compute likely hood. This saves computation burden.
            - set at default at 50 heterozygous variants. '''
        num_of_het_site += 1
        if num_of_het_site >= num_of_hets:
            continue

        #print('num_of_het_site: ', num_of_het_site, n)  # marker for debugging


        ''' Step 05 - A (i) : compute "nucleotide probabilities" aka "emission probabilities"
        - create dictionary to store "nucleotide counts".
        - Set the counts of each nucleotides at zero (0).
        - the counts are updated at "n-th" level after counting the existing nucleotides (either A,T,G,C)
          at "n-th" position for all the samples. '''

        # setting dictionary for "emission/nucleotide counts"
        nucleotide_count_dict = {'A': 0, 'T': 0, 'G': 0, 'C': 0}

        # dictionary for "emission" probabilities
        nucleotide_prob_dict = nucleotide_count_dict.copy()   # makes a shallow copy



        ''' Creating variables to store the product-values of the transition probabilities.
        These are updated for each level of "n" paired with each level of "m". '''
        # potp -> product of transition probability

        ##** need to change here ??

        ############# deprecated ??
        #potp_hapb1a_b2a = 1
        #potp_hapb1b_b2b = 1

        #potp_hapb1a_b2b = 1
        #potp_hapb1b_b2a = 1
        #########################

        if maxed_as == '+':
            tp_hapb1a_b2a = 0
            tp_hapb1b_b2b = 0

            tp_hapb1a_b2b = 0
            tp_hapb1b_b2a = 0

        elif maxed_as == '*':
            tp_hapb1a_b2a = 1
            tp_hapb1b_b2b = 1

            tp_hapb1a_b2b = 1
            tp_hapb1b_b2a = 1




        ''' Step 05 - A (ii) :
        - Now, calculate the initial nucleotide counts at "n-th" level,
        before computing the transition counts.
        - only calculated from "v1" at n-th level and only once for each parse/iteration '''
        for (x, y) in sample_list:
            for nucleotide in 'ATGC':
                nucleotide_count_dict[nucleotide] += v1[y][n].split('|').count(nucleotide)


        # now, compute emission probabilities
        total_nucl = sum(list(nucleotide_count_dict.values()))

        for ky, vy in nucleotide_count_dict.items():
            #print(ky, vy)
            nucleotide_prob_dict[ky] = vy/total_nucl


        ''' Step 05 - B : Count number of transition from each nucleotide (n-th) to each nucleotide (m-th).
        - Now, we read nucleotides (ATGC) at each level of "m"
        to compute the transition from each level of "n". '''

        # to control certain number of hetVars to include in computation
        num_of_het_site_at_m = 0
        for m in range(len(v2[soi + ':PI'])):    # m-ranges from block02

            ''' Like at "n-th" level, skip if InDel present at this "m-th" level.
            But InDel will be phased along with hitchhiking SNPs. '''
            if len(v2[soi + ':PG_al'][m]) > 3\
                or '*' in v2[soi + ':PG_al'][m]:
                continue

            ''' only use certain number of HetVars to compute likely hood at "later block".
                - This saves computation burden.
                - set at default at 200 heterozygous variants. '''
            num_of_het_site_at_m += 1
            if num_of_het_site_at_m >= num_of_hets:
                continue


            ''' Creating an empty dictionary to store transition counts from "n-th" level of V1 to "m-th" level of V2.
            - ('A', 'T') represents transition from 'A' to 'T'.
            - ** for future: probably upgrade using numpy '''
            transition_count_dict = {
                ('A', 'A'): 0, ('A', 'T'): 0, ('A', 'G'): 0, ('A', 'C'): 0,
                ('T', 'A'): 0, ('T', 'T'): 0, ('T', 'G'): 0, ('T', 'C'): 0,
                ('G', 'A'): 0, ('G', 'T'): 0, ('G', 'G'): 0, ('G', 'C'): 0,
                ('C', 'A'): 0, ('C', 'T'): 0, ('C', 'G'): 0, ('C', 'C'): 0, }

            ''' Dictionary to store above counts as probabilities. '''
            # shallow copy should be fine  # **for future - do deep copy if need be
            transition_prob_dict = transition_count_dict.copy()

            ''' Step 05 - B(i) : Now, loop through each sample to compute transition counts (from_ , to)
            for each nucleotide combination'''
            for x, y in sample_list:

                ''' Skip data mining if values are missing for the sample.
                Continues to next sample, since there is no point in computing transition.
                - We also skip any variant that has "*" allele (i.e deletion spanning variants) in it. '''
                if v1[y][n] == '.' or v2[y][m] == '.' \
                        or '*' in v1[y][n] or '*' in v2[y][m]:
                    continue

                #############################################
                ## ** for future: add min number of samples called in a line, before proceeding
                #############################################

                ''' Mine nucleotides at n'th and m'th position in v1 and v2 for each samples.
                - soi (sample of interest) is also included in the counting, rather than adding pseudo counts.'''
                nucl_B1 = (v1[y][n]).split('|')  # nucleotides at n'th of Block01
                nucl_B2 = (v2[y][m]).split('|')  # nucleotides at m'th of Block02


                ''' Step 05 - B(ii) : Create zipped vs. product configuration of the
                haploptyes from two positions .
                - i.e create possible haplotype configurations between ("n-th" and "m-th") nucleotides.
                 If the index (PI value) are same we create zip, if index (PI value) are
                 different we create product (see "possible haplotypes" below).
                 E.g        same PI-values     vs.   different PI-values
                            POS  PI   PG_al          POS  PI   PG_al
                    n'th    21   4     A|T           21   4    A|T
                    m'th    35   4     C|G           35   7    C|G

                 possible haplotypes:
                            A-C and T-G,       vs.   A-C, A-G, T-C, T-G
                '''

                ''' Create zipped configuration'''
                # if PI values are same/equal prepare zip between phased nucleotides
                if v1[x][n] == v2[x][m]:
                    ziped_nuclb1b2 = list(zip(nucl_B1, nucl_B2))

                    ''' then count the number of transitions '''
                    for (from_, to) in transition_count_dict:
                        transition_count_dict[(from_, to)] += ziped_nuclb1b2.count((from_, to))

                # if PI values are not same then we create product between phased nucleotides
                elif v1[x][n] != v2[x][m]:
                    prod_nuclb1b2 = list(product(nucl_B1, nucl_B2))

                    ''' then count the number of transitions '''
                    for (from_, to) in transition_count_dict:
                        transition_count_dict[(from_, to)] += prod_nuclb1b2.count((from_, to)) / 2


            '''Step 06 : Compute likelyhood of "parallel" vs. "alternate" configuration of two consecutive haplotype blocks.
            A) calculate transition probabilities from "n-th" to "m-th" level.
            B) then, calculate the product of the transition from one level of "n" to several levels of "m".
            C) then, calculate joint probability (summed or product) for transitions from several levels of "n"
               to several levels of "m".
            D) then, calculate overall likelyhood and log2Odds of "parallel" vs. "alternate" configuration
               using both forward and backward algorith. '''

            ''' Step 06 - A : using nucleotide counts and transition counts (at "n-th" and "m-th" of each sample),
                              compute transition probabilities.
            ** Note (for future): At the end this all possible transition probabilities should sum to 1.
            But, at some places transitions may be missing due to situation like ('A', '.')
            i.e from 'A' to '.' empty.
            ** future - This problem may be adjusted in the future'''

            ''' Step 06 - A (i) : compute all observed transition probabilities '''
            for (from_, to) in transition_prob_dict:
                transition_prob_dict[(from_, to)] = \
                    compute_transition_probs(transition_count_dict[(from_, to)],
                                             nucleotide_count_dict[from_])


            ''' Step 06 - A (ii) : find observed configuration for soi at "n-th" and "m-th" level
            i.e (from_, to). '''
            hapb1a_hapb2a_transition = (hapb1a_hapb2a[0][n], hapb1a_hapb2a[1][m])
            hapb1b_hapb2b_transition = (hapb1b_hapb2b[0][n], hapb1b_hapb2b[1][m])

            hapb1a_hapb2b_transition = (hapb1a_hapb2b[0][n], hapb1a_hapb2b[1][m])
            hapb1b_hapb2a_transition = (hapb1b_hapb2a[0][n], hapb1b_hapb2a[1][m])

            ''' Step 06 - B : computes the "product of transition probabilities" from "n-th" to
             several levels of "m" using for-loop.
             - ** because the below code is indented one level below "for m in range(len(v2['ms02g:PI']))".
             - ** no need to add pseudo-counts, because if no haplotypes are observed in any samples except soi,
               the prob(from_, to) for each configuration will be 1/4 there by nullifying the likelyhoods to "1".
            '''

            ######### deprecated ??
            #potp_hapb1a_b2a *= transition_prob_dict[hapb1a_hapb2a_transition]
            #potp_hapb1b_b2b *= transition_prob_dict[hapb1b_hapb2b_transition]
            #potp_hapb1a_b2b *= transition_prob_dict[hapb1a_hapb2b_transition]
            #potp_hapb1b_b2a *= transition_prob_dict[hapb1b_hapb2a_transition]


            # new additions ****
            # adding emission probabilities to the method to get more fair estimate of the likelihoods
            # p(X given Y) = p(X) * p(XtY - transitions)
            # ** for future: this if/else can be rather fixed by establishing a new function
            ## we are multiplying p(transition)*p(emission)  and doing maxSum or maxProd depending upon what is required
            ## the codes below can be dramatically improved.

            if maxed_as == '+':
                tp_hapb1a_b2a += transition_prob_dict[hapb1a_hapb2a_transition] \
                                 * Decimal(nucleotide_prob_dict[hapb1a_hapb2a_transition[0]])
                tp_hapb1b_b2b += transition_prob_dict[hapb1b_hapb2b_transition] \
                                 * Decimal(nucleotide_prob_dict[hapb1b_hapb2b_transition[0]])
                tp_hapb1a_b2b += transition_prob_dict[hapb1a_hapb2b_transition] \
                                 * Decimal(nucleotide_prob_dict[hapb1a_hapb2b_transition[0]])
                tp_hapb1b_b2a += transition_prob_dict[hapb1b_hapb2a_transition] \
                                 * Decimal(nucleotide_prob_dict[hapb1b_hapb2a_transition[0]])

            elif maxed_as == '*':
                tp_hapb1a_b2a *= transition_prob_dict[hapb1a_hapb2a_transition] \
                                 * Decimal(nucleotide_prob_dict[hapb1a_hapb2a_transition[0]])
                tp_hapb1b_b2b *= transition_prob_dict[hapb1b_hapb2b_transition] \
                                 * Decimal(nucleotide_prob_dict[hapb1b_hapb2b_transition[0]])
                tp_hapb1a_b2b *= transition_prob_dict[hapb1a_hapb2b_transition] \
                                 * Decimal(nucleotide_prob_dict[hapb1a_hapb2b_transition[0]])
                tp_hapb1b_b2a *= transition_prob_dict[hapb1b_hapb2a_transition] \
                                 * Decimal(nucleotide_prob_dict[hapb1b_hapb2a_transition[0]])
            ##



        ''' Step 06 - C : compute the max sum or max product of the transition probabilities
        across several levels of "n" to several levels of "m".
        ** Note: we can either do "max sum" or "max product".
        So, "cul_of_pt_hapb1a_b2a" is the sum of the likelyhoods of hapBlock1A being phased with hapBlock2A'''

        if maxed_as == '+':
            cul_of_pt_hapb1a_b2a += tp_hapb1a_b2a
            cul_of_pt_hapb1b_b2b += tp_hapb1b_b2b

            cul_of_pt_hapb1a_b2b += tp_hapb1a_b2b
            cul_of_pt_hapb1b_b2a += tp_hapb1b_b2a

        elif maxed_as == '*':
            cul_of_pt_hapb1a_b2a *= (tp_hapb1a_b2a)
            cul_of_pt_hapb1b_b2b *= (tp_hapb1b_b2b)

            cul_of_pt_hapb1a_b2b *= (tp_hapb1a_b2b)
            cul_of_pt_hapb1b_b2a *= (tp_hapb1b_b2a)


    ### Test markers
    #print()
    #print(orientation)
    #print('cumulative scores')
    #print(cul_of_pt_hapb1a_b2a)
    #print(cul_of_pt_hapb1b_b2b)
    #print(cul_of_pt_hapb1a_b2b)
    #print(cul_of_pt_hapb1b_b2a)


    ''' Step 06 - D : Now, compute the "Odds ratio" and "log2 of the Odds"
        of each possible haplotype configuration, for both 1st vs. 2nd configurations '''

    ''' Step 06 - D(i) : compute the likely hood of first configuration (lhfc) vs. second configuration (lhsc)
    First Configuration = hapb1a with hapb2a, and hapb1b with hapb2b.
    Second Configuration = hapb1a with hapb2b, and hapb1b with hapb2a. '''
    if maxed_as == '+':
        lhfc = Decimal(cul_of_pt_hapb1a_b2a + cul_of_pt_hapb1b_b2b)
        lhsc = Decimal(cul_of_pt_hapb1a_b2b + cul_of_pt_hapb1b_b2a)

    elif maxed_as == '*':
        lhfc = Decimal(cul_of_pt_hapb1a_b2a * cul_of_pt_hapb1b_b2b)
        lhsc = Decimal(cul_of_pt_hapb1a_b2b * cul_of_pt_hapb1b_b2a)

    #lhfc = Decimal(cul_of_pt_hapb1a_b2a * cul_of_pt_hapb1b_b2b)
    #lhsc = Decimal(cul_of_pt_hapb1a_b2b * cul_of_pt_hapb1b_b2a)

    return lhfc, lhsc



    ##################################################
    ''' ** Note - for the future: The likelyhood can also be tested with each configuration split into two parts.
    And, we can check if Hapb1a has higher affinity to Hapb2a compared to Hapb2b vs. if Hapb1b has higher
    affinity to Hapb2b compared to Hapb2a. This two way testing can provide important insights if recombination
    may have happened with in gene, or if there is problem with paralogous alignment (paralog detection is
    possible if homVar haplotype is phased to hetVar haplotype; and hap_1a shows equal affinity to hap_2a and 2b).
    In such situation hap1a should have been homVar for all phased SNPs but due to paralog sequence aligment it's
    varaints were called as hetVar.
    - but keeping this problem for the future.

    likelyhood_hapb1a_with_b2a_Vs_b2b = LH_hapb1awb2avsb2b = \
        likelyhood(cul_of_pt_hapb1a_b2a, cul_of_pt_hapb1a_b2b)

    .. and similarly others
    '''
    #############################################################


def extend_phase_state(soi, k1, k2, v1, v2, k2_new, flipped,
                       lods2_score_1st_config, lods_cut_off,
                       extended_haplotype,
                       hapb1a_hapb2a, hapb1b_hapb2b):

    #print()
    hmmm = 1  # ** - just added to highlight the docstring in pyCharm; delete in future

    ''' Step 08 : extend the phase states between two consecutive haplotypes based on the evidence.
     - default cut-off is set at log2 score of 4 (i.e 2^4 times likely), which can be changed by user
     - this function also returns the updated values of checker-variable i.e "k2_new" and "flipped",
       these checker-variables are very important to extend phase in a proper state.
       - "k2_new" keeps the log of "if the haplotype was extended"
       - "flipped" keeps the log of "if it was "parallel" vs. "alternate" configuration".
     - we write the data to a file, and also store in another variable to compute "haplotype stats" if need be.
    '''

    # when "k2_new" is not updated
      # this happens for the very first consecutive block of a contig
      # this happens when the previous consecutive blocks failed to extend
    if k2_new == '':

        # if the |lods score| is above cut-off the phase needs extension in parallel configuration
        if lods2_score_1st_config > lods_cut_off:  # i.e 2^4 times more likely
            k2_new = k1     # now, "k1" value gets carried over for next consecutive run
            flipped = 'no'  # since it was phased in parallel configuration, there for "no flip"

            for xth in range(len(v2[soi + ':PI'])):
                new_line = '\t'.join([v2['CHROM'][xth], v2['POS'][xth],  v2['REF'][xth], v2['all-alleles'][xth],
                                      k2_new, hapb1a_hapb2a[1][xth] + '|' + hapb1b_hapb2b[1][xth]]) + '\n'
                if writelod == 'yes':
                    new_line = new_line.rstrip('\n') + '\t'+ str(round(lods2_score_1st_config, 5)) +'\n'

                extended_haplotype += new_line


        # if the |-lods score| is above cut-off, phase needs extension in "alternate" configuration
        elif lods2_score_1st_config < -lods_cut_off:  # i.e 2^4 times less likely
            k2_new = k1
            flipped = 'yes'

            for xth in range(len(v2[soi + ':PI'])):
                new_line = '\t'.join([v2['CHROM'][xth], v2['POS'][xth],  v2['REF'][xth], v2['all-alleles'][xth],
                                      k2_new, hapb1b_hapb2b[1][xth] + '|' + hapb1a_hapb2a[1][xth]]) + '\n'
                if writelod == 'yes':
                    new_line = new_line.rstrip('\n') + '\t'+ str(round(lods2_score_1st_config, 5)) +'\n'

                extended_haplotype += new_line

        # if the |lods score| does not pass cut-off threshold in either direction, there is no phase extension
        # so, reset the "checker-variables" and write/store data in original configuration as input data
        else:  # i.e 4 times less likely
            k2_new = ''
            flipped = ''

            # print('no phase extension, no flip')  # marker for debugging
            for xth in range(len(v2[soi + ':PI'])):
                new_line = '\t'.join([v2['CHROM'][xth], v2['POS'][xth], v2['REF'][xth], v2['all-alleles'][xth],
                                      k2, hapb1a_hapb2a[1][xth] + '|' + hapb1b_hapb2b[1][xth]]) + '\n'
                if writelod == 'yes':
                    new_line = new_line.rstrip('\n') + '\t'+ str(round(lods2_score_1st_config, 5)) +'\n'

                extended_haplotype += new_line


    # when previous consecutive blocks were extended,
      # "k2_new" now carries a value from previous (k2) which is actually (k1) in this run
      # "flipped" is also carrying an updated value
    elif k2_new != '':
        if flipped == 'no':

            # previous config was extended but not flipped, so it must stay in "parallel" config this time.
            if lods2_score_1st_config > lods_cut_off:  # i.e 4 times more likely
                k2_new = k2_new
                flipped = 'no'

                #print (1_A is phased with 2_A, no flip, k2:%s \n' %k2_new)  # marker for debugging
                for xth in range(len(v2[soi + ':PI'])):
                    new_line = '\t'.join([v2['CHROM'][xth], v2['POS'][xth],  v2['REF'][xth], v2['all-alleles'][xth],
                                          k2_new, hapb1a_hapb2a[1][xth] + '|' + hapb1b_hapb2b[1][xth]]) + '\n'
                    if writelod == 'yes':
                        new_line = new_line.rstrip('\n') + '\t' + str(round(lods2_score_1st_config, 5)) + '\n'

                    extended_haplotype += new_line

            # previous config was extended but not flipped. So, now when the phase extends
              # the configuration is kept "alternate" because lods2 score is negative.
            elif lods2_score_1st_config < -lods_cut_off:  # i.e 4 times less likely
                k2_new = k2_new
                flipped = 'yes'

                #print ('1_A is phased with 2_B, yes flipped, k2:%s \n' % k2_new) # marker for debugging
                for xth in range(len(v2[soi + ':PI'])):
                    new_line = '\t'.join([v2['CHROM'][xth], v2['POS'][xth],  v2['REF'][xth], v2['all-alleles'][xth],
                                          k2_new, hapb1b_hapb2b[1][xth] + '|' + hapb1a_hapb2a[1][xth]]) + '\n'
                    if writelod == 'yes':
                        new_line = new_line.rstrip('\n') + '\t' + str(round(lods2_score_1st_config, 5)) + '\n'

                    extended_haplotype += new_line

            # if phase state isn't extending, write the data as is and reset the "checker-variables"
            else:
                k2_new = ''
                flipped = ''

                #print ('no phase extension, no flip, k2 %s\n' % k2)  # marker for debugging
                for xth in range(len(v2[soi + ':PI'])):
                    new_line = '\t'.join([v2['CHROM'][xth], v2['POS'][xth], v2['REF'][xth], v2['all-alleles'][xth],
                                          k2, hapb1a_hapb2a[1][xth] + '|' + hapb1b_hapb2b[1][xth]]) + '\n'
                    if writelod == 'yes':
                        new_line = new_line.rstrip('\n') + '\t' + str(round(lods2_score_1st_config, 5)) + '\n'

                    extended_haplotype += new_line

        # if previous block were extended but flipped
        elif flipped == 'yes':

            # previous config was extended but flipped, so to maintain it's phase state
              # it must stay in "parallel" config this time by flipping with previous extension.
            if lods2_score_1st_config > lods_cut_off:  # i.e 4 times more likely
                k2_new = k2_new
                flipped = 'yes'

                #print('\n1_A is phased with 2_A, but yes flipped, k2:%s \n' % k2_new)  # marker for debug
                for xth in range(len(v2[soi + ':PI'])):
                    new_line = '\t'.join([v2['CHROM'][xth], v2['POS'][xth], v2['REF'][xth], v2['all-alleles'][xth],
                                          k2_new, hapb1b_hapb2b[1][xth] + '|' + hapb1a_hapb2a[1][xth]]) + '\n'
                    if writelod == 'yes':
                        new_line = new_line.rstrip('\n') + '\t' + str(round(lods2_score_1st_config, 5)) + '\n'

                    extended_haplotype += new_line


            # previous config was extended but flipped, so to maintain it's phase state
              # it must stay in "alternate" config this time by not flipping with previous extension.
            elif lods2_score_1st_config < -lods_cut_off:  # i.e 4 times less likely
                k2_new = k2_new
                flipped = 'no'

                #print('\n1_A is phased with 2_B, but not flipped, k2:%s \n' % k2_new)  # debug
                for xth in range(len(v2[soi + ':PI'])):
                    new_line = '\t'.join([v2['CHROM'][xth], v2['POS'][xth], v2['REF'][xth], v2['all-alleles'][xth],
                                          k2_new, hapb1a_hapb2a[1][xth] + '|' + hapb1b_hapb2b[1][xth]]) + '\n'
                    if writelod == 'yes':
                        new_line = new_line.rstrip('\n') + '\t' + str(round(lods2_score_1st_config, 5)) + '\n'

                    extended_haplotype += new_line

            # previous extension was phased, but this one didn't
            # so, write the data as it is and reset the "checker-variables"
            else:
                k2_new = ''
                flipped = ''

                #print('\nno phase extension, no flip (reset), k2 %s\n' % k2)  # marker for debugging
                for xth in range(len(v2[soi + ':PI'])):
                    new_line = '\t'.join([v2['CHROM'][xth], v2['POS'][xth], v2['REF'][xth], v2['all-alleles'][xth],
                                                 k2, hapb1a_hapb2a[1][xth] + '|' + hapb1b_hapb2b[1][xth]]) + '\n'
                    if writelod == 'yes':
                        new_line = new_line.rstrip('\n') + '\t' + str(round(lods2_score_1st_config, 5)) + '\n'

                    extended_haplotype += new_line


    return k2_new, flipped, extended_haplotype


''' function that computes several stats from haplotype block before and after phase extension.
These data can be used to compare the improvements in phase extension, and plotting histogram.
This can include the size based on genome co-ordinate of the "POS" in each block.
This can also include the average number of variants per haplotype before/after phase extension.
and make a plot out of it - using pyPlot, MatlibPlot or R. '''
def compute_haplotype_stats(hap_data, soi, prefix) :

    print()
    ''' this function computes the stats of the haplotype file both before and after phase extension.
        - we pass in "prefix" as "initial" vs. "final" for the haplotype block before and after
          extension to compute statistics for respective phased file. '''

    ''' Step 09 - A : write this haplotype information in a more concise format.
        - in this step we are restructuring the data so proper statistics can be handled properly. '''
    hap_stats = hap_data.groupby(['CHROM', soi +':PI'], sort=False)\
        .size().rename('num_Vars_by_PI')\
        .reset_index(level=1).astype(str)\
        .groupby(level=0).agg(','.join).reset_index()

    ''' add other required columns '''
    hap_stats['range_of_PI'] = hap_data.groupby(['CHROM', soi +':PI'], sort=False)['POS']\
        .apply(lambda g: g.max() - g.min()).rename('range').reset_index(level=1)\
        .astype(str).groupby(level=0).agg(','.join).reset_index()['range']

    hap_stats['total_haplotypes'] = hap_stats\
        .apply(lambda row: len(row[soi +':PI'].split(',')), axis=1)

    hap_stats['total_Vars'] = hap_stats \
        .apply(lambda row: sum([int (i) for i in row.num_Vars_by_PI.split(',')]), axis=1)

    ''' write this concise statistics to a file. '''
    pd.DataFrame.to_csv(hap_stats, outputdir + '/' + prefix+'_haplotype_stats_' + soi + '.txt',
                        sep='\t', header=True, index=False)


    ## ** for future - plot haplotype STATs - may need improvement.
    ''' Step 09 - B : Now, use the above dataframe to plot several statistical plots from here on.
        - plot the haplotype statistics using matplotlib, pylab or ggplot (python).
        - start with "with open() .... " to open a file, write the plot and then close. '''
    # some good links for hist-plots: http://pandas.pydata.org/pandas-docs/version/0.18/visualization.html

    # plots total number of variants per-chromosome
    with open(outputdir + '/' + 'total_vars_'+ soi +'_'+ prefix+'.png', 'wb') as fig_initial:
        hap_stats.plot(x='CHROM', y='total_Vars', kind='bar')
        plt.xlabel('chromosomes')
        plt.ylabel('number of variants')
        plt.suptitle('number of variants for each chromosome')
        plt.savefig(fig_initial)
        # plt.show()  # optional: to show the plot to the console

    # plots total number of haplotypes per-chromosome
    with open(outputdir + '/' + 'total_haps_' + soi + '_' + prefix + '.png', 'wb') as fig_initial:
        hap_stats.plot(x='CHROM', y='total_haplotypes', kind='bar')
        plt.xlabel('chromosomes')
        plt.ylabel('number of haplotypes')
        plt.suptitle('number of haplotypes for each chromosome')
        plt.savefig(fig_initial)


    # plot histogram of haplotype size (by number of variants) per-chromosome
    # make different plots for each chromosome, but X-axis is shared.
    with open(outputdir + '/' + 'hap_size_byVar_' + soi + '_' + prefix + '.png', 'wb') as fig_initial:

        ## making histogram plot.
        # raising if-else when number of chromosome is 1 vs. more than 1. ** This can be optimize in the future.
        fig, ax = plt.subplots(nrows=len(hap_stats), sharex=True)

        # when we have only one chromosome
        if len(hap_stats) == 1:
            # fig, ax = plt.subplots(nrows=len(hap_stats), sharex=True, squeeze=False)
            data_i = hap_stats['num_Vars_by_PI'].str.split(',')
            data_i = pd.Series([int(x) for x in data_i[0]])
            data_i.plot(kind='hist', label=str(hap_stats['CHROM']), alpha=0.5)
            plt.ylabel('frequency of the haplotypes')

        elif len(hap_stats) > 1:
            for i, data in hap_stats.iterrows():
                # first convert data to list of integers
                data_i = [int(x) for x in data['num_Vars_by_PI'].split(',')]
                # print(data)
                # print(data_i)
                ax[i].hist(data_i, label=str(data['CHROM']), alpha=0.5)
                ax[i].legend()

            fig.text(.05, .5, 'frequency of the haplotypes', ha='center', va='center', rotation='vertical')

        # add a common x-label
        plt.xlabel('size of the haplotype (by number of variants)')
        # plt.ylabel('frequency of the haplotypes')
        plt.suptitle(prefix + ' histogram of size of the haplotype (by number of variants) \n'
                     'for each chromosome')
        plt.savefig(fig_initial)



    # plot histogram of haplotype size (by genomic ranges of haplotype) per-chromosome
    # make different plots for each chromosome, but X-axis is shared.
    with open(outputdir + '/' + 'hap_size_byGenomicRange_' + soi + '_' + prefix + '.png', 'wb') as fig_initial:

        ## making histogram plot.
        # raising if-else when number of chromosome is 1 vs. more than 1. ** This can be optimize in the future.
        fig, ax = plt.subplots(nrows=len(hap_stats), sharex=True)

        # when we have only one chromosome
        if len(hap_stats) == 1:
            # fig, ax = plt.subplots(nrows=len(hap_stats), sharex=True, squeeze=False)
            data_i = hap_stats['range_of_PI'].str.split(',')
            data_i = pd.Series([int(x) for x in data_i[0]])
            data_i.plot(kind='hist', label=str(hap_stats['CHROM']), alpha=0.5)
            plt.ylabel('frequency of the haplotypes')

        # when the dataframe has more than one chromosome
        elif len(hap_stats) > 1:
            for i, data in hap_stats.iterrows():
                # first convert data to list of integers
                data_i = [int(x) for x in data['range_of_PI'].split(',')]
                # print(data)
                # print(data_i)
                ax[i].hist(data_i, label=str(data['CHROM']), alpha=0.5)
                ax[i].legend()

            # adding y-label separately if there are multiple chromosomes
            fig.text(.05, .5, 'frequency of the haplotypes', ha='center', va='center', rotation='vertical')

        # writing a common x-label
        plt.xlabel('size of the haplotype (by Genomic Distance)')
        #plt.ylabel('frequency of the haplotypes')
        plt.suptitle(prefix + ' histogram of size of the haplotype (by Genomic Distance)\n'
                     'for each chromosome')
        plt.savefig(fig_initial)


    #print('Global maximum memory usage: %.2f (mb)' % current_mem_usage())

    # clear memory
    del hap_data, hap_stats



''' function to return the transitions probabilities from transition counts'''
def compute_transition_probs(pX_Y, pX):
    # returns transition probs from "X" to "Y"
    try:
        return Decimal(pX_Y / pX)

    except ZeroDivisionError:
        return 0

    # ** Note: the exception is raised because all alleles (A,T,G,C) may not be represented at "n-th" level
    # so, "pX" for allele "A" might be zero, hence raising the error.
    # so, transition from "A" to any (A,T,G,C) is technically zero,
        # but theoretically it should be at least 1/16.
    # ** future: this theoretical estimation may be optimized in the future.



''' Merge ordered dictionaries, mapping each key to lists of all values mapped to in order of occurence.'''
def accumulate(data):
    acc = collections.OrderedDict()
    for dataum in data:
        for kg, vg in dataum.items():
            acc.setdefault(kg, []).append(vg)
    return acc


''' prepare list of all available samples (index, value) in the input file
    ** note: we also control output of specific sample if user desires. '''
def find_samples(samples):
    # the below method preserves the order of samples.
    seen = set()
    sample_list = [x.split(':')[0] for x in samples if ':' in x]
    sample_list = [x for x in sample_list if not (x in seen or seen.add(x))]
    sample_list = [((x + ':PI'), (x + ':PG_al')) for x in sample_list]


    return sample_list

    # ** Note: this returns data as:
    # sample_list = [('ms01e:PI', 'ms01e:PG_al'), ('ms02g:PI', 'ms02g:PG_al'),
               #('ms03g:PI', 'ms03g:PG_al'), ('ms04h:PI', 'ms04h:PG_al')]


''' to monitor memory '''
def current_mem_usage():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.



if __name__ == '__main__':
    main()