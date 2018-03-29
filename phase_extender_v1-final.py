#!/home/bin/python3


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
import pandas as pd
import resource
import time
import io

# for plotting
import matplotlib.pyplot as plt



def main():

    print()
    print("loading the argument variables and assigning values to them....")
    ''' Define required argument for interactive mode program. '''
    parser = argparse.ArgumentParser()

    parser.add_argument("--nt",
                        help="number of threads -> the maximum number of threads that can be run is the maximum number "
                             "of different chromosomes/chromosomes in the input haplotype file.",
                        default=1, required=False)

    parser.add_argument("--input",
                        help="name of the input haplotype file -> This haplotype file should contain phased genotype "
                             "and haplotype blocks for several samples. The phased haplotype block index should be "
                             "represented by 'PI' and should be an unique value in that sample. The phased genotype "
                             "for this unique block should be represented by 'PG_al' in that sample.",
                        required=True)

    parser.add_argument("--SOI",
                        help="Sample of interest -> type in the name of sample for which you want to "
                             "extend the haplotype.",
                        required=True)

    parser.add_argument("--output",
                        help="Prefix for the name of the output file. ** default: '' ",
                        default='', required=False)

    parser.add_argument("--refHap",
                        help="reference haplotype panel -> This file should also contain 'PI' and 'PG_al' values "
                             "for each sample in that haplotype reference panel."
                             "** default: no reference haplotype file. ",
                        default='', required=False)

    parser.add_argument("--useSample",
                        help="list of samples -> use phase state information only from these samples while running "
                             "phase-extenstion for the SOI (sample of interest). "
                             "This helps you control the phase-extension on SOI by using select samples from "
                             "the pool of samples (of refHap and/or input). "
                             "this becomes helpful when you want to see if including samples from particular .."
                             " population increases the phase extension."
                             "Available options: 'all','refHap','input','comma separated name of samples'. "
                             "** default: all the samples in (refHap + input) will be used.",
                        default='all', required=False)

    parser.add_argument("--bed",
                        help=" bed file -> process the haplotype extension only in this bed regions. "
                             "This is useful if you want to limit haplotype extension only within certain "
                             "regions, like - within genes, exons, introns, QTL, etc. ** default: no bed file ",
                        default= '', required=False)

    parser.add_argument("--snpTh",
                        help="snp threshold -> minimum number of SNPs required in both haplotype blocks before starting "
                             "phase extension. Option: an integer value. ** default snpTh = 3 ",
                        default=3, required=False)

    parser.add_argument("--numHets",
                        help="number of heterozygote SNPs -> maximum number of heterozygote SNPs used in consecutive "
                             "haplotype blocks for computing likelyhood of the phase states. "
                             "Option: type an integer value. ** default numHets = 40 ",
                        default=40, required=False)

    parser.add_argument("--lods",
                        help="log2 of Odds cut off -> cutoff threshold to assign the phase state between consecutive "
                             "haplotype blocks. Option: type an integer value. eg: 3. "
                             "**Note: Default value is set at (2^5 = 32). So, two consecutive blocks will be joined "
                             "in parallel configuration if default lods-score is used.",
                        default=5, required=False)
    #**Note:  positive score above the cutoff values joins two consecutive block in "
                             #"parallel configuration, negative score below the cutoff joins two consecutive block in "
                             #"alternate configuration.

    parser.add_argument("--culLH",
                        help="cumulation of the likelhood estimates. The likelhoods for two possible configuration can "
                             "either be max-sum vs. max-product. Default is max-product. Options: 'maxPd' or 'maxSum'. ",
                        default='maxPd', required=False)

    parser.add_argument("--writeLOD",
                        help="print log2 of Odds to the output file-> writes the calculate LODs between two "
                             "consecutive haplotype blocks when processing phase extension. "
                             "Options: 'yes', 'no'. **Note: the 'lods-score' are printed regardless if the "
                             "consecutive blocks are joined or not.",
                        default='no', required=False)  # ** to do - add this feature

    parser.add_argument("--hapStats",
                        help="computes the statistics of the of the input haplotype file vs. extended haplotype "
                             "for the sample of interest. Also, produces the histgram of the haplotype.",
                        default='no', required=False)

    parser.add_argument("--onlyPhasedSites",  help=" write the lines that have data missing for SOI on extended "
                                                   "haplotype block."
                        ,default='yes', required=False)  # ** to do - add this



    ''' create a global argument variable and declare values of the global variables.
        The values for the global variable are assigned using arguments (args.()) passed by user.'''
    # this is necessary if arguments are declared at several functions using "args.variable" parameter.
    global args;
    args = parser.parse_args()   # .. but keep it as it is.

    global output
    global use_bed
    global soi
    global snp_threshold
    global num_of_hets
    global lods_cut_off
    global use_sample
    global sample_list  # used in coherence with "use_sample"
    global maxed_as


    global update_data  # local ??
    global extended_haplotype  # probably local and deprecated ?? ??
    global time01       # set time for each chromosome vs. all chromosome; deprecated or local ??


    ''' assign the values for the global variables '''
    # advantage of this method is that we can change the program to non-interactive mode easily and
    # point to specific values/files on this top part. This becomes helpful during debugging.

    time01 = time.time()

    soi = args.SOI
    print('using sample "%s" ' %(soi))

    output = args.output

    # assign number of process to be used
    nt = int(args.nt)  # default, nt = 1
    print('using "%s" multiple processes ' % (nt))

    input_file = args.input  # the input haplotype file
    # input_file = 'allele_table_for_phase_extender.txt'
    print('using haplotype file "%s" ' % (input_file))

    lods_cut_off = args.lods  # log_of_odds_cut_off, default = 5
    print('using log2 odds cut off of "%s" ' % (lods_cut_off))

    # minimum number of SNPs in a haplotype block before it can be phase-extended
    snp_threshold = args.snpTh   # default, snp_threshold = 3
    print('consecutive haplotype block each should have minimum of "%s" SNPs ' % (snp_threshold))

    # controls the max number of hetVars to include in likelyhood calculation
    num_of_hets = int(args.numHets)   # default, num_of_hets = 40
    print('using maximum of "%s" heterozygote sites in each consecutive blocks to compute '
          'transition probabilities' % (num_of_hets))

    # add argument for max sum vs. max product of likelyhood estimates before calculating the LOD-score
    maxed_as = args.culLH   # default, maxed_as = "*"
    if maxed_as == 'maxSum':
        max_is = 'max sum'
        maxed_as = '+'
    elif maxed_as == 'maxPd':
        max_is = 'max product'
        maxed_as = '*'
    print('using "%s" to estimate the cumulative maximum likelyhood of the transition '
          'probabilities between two consecutive blocks ' % (max_is))


    ##set the required variables related to bed file if bed file is given.
    # use provided bed file to limit phase extension on bed regions/boundries
    if args.bed != "":      # default, use_bed = "no"
        use_bed = 'yes'
        bed_file = args.bed
        print('using the bed file "%s" to limit phase extension at the bed boundries '%(bed_file))
    else:
        use_bed = 'no'
        print('no bed file is given.'
              '\nextending haplotype block through out the contig boundries')


    # if a haplotype panel is provided then the reference panel can be used as backbone or ..
      # .. meaningful data to aid in phase extension
    if args.refHap != "":  # default, use_refHap = "no"
        use_refhap = 'yes'
        refhap = args.refHap
        print('using the reference haplotype panel "%s" ' %(refhap))
    else:
        use_refhap = 'no'
        print('no reference haplotype panel is provided ')


    # which samples to use while running phase-extension, default is all (hapRef + input)
    # this variable is updated, later in the pipeline when all sample names are collected
    use_sample = args.useSample  # default, use_sample = "all"


    # print the hapstats to file and also plot histogram
    if args.hapStats == 'yes':  # default, hapstats = 'yes'
        hapstats = 'yes'
        print('statistics of the haplotype before and after extension will '
              'be prepared for the sample of interest i.e "%s" ' %(soi))
    else:
        hapstats = 'no'
        print('statistics of the haplotype before and after extension will not '
          'be prepared for the sample of interest i.e "%s".'
              '\nOnly extendent haplotype block will be prepared.' % (soi))



    # write calculated LOD (log2 of odds) of the possible phase state between two consecutive blocks
    writeLOD = args.writeLOD    # default, writeLOD = "no"
    if writeLOD == 'yes':
        writelod = 'yes'
        print('LOD (log 2 of odds) for consecutive block will be written to the file ')
    elif writeLOD != 'yes':
        writelod = 'no'


    #### completed assignment of values to argument variables




    ''' printing authorship '''
    print()
    print("##################################################")
    print("         Welcome to phase-extender version %i       "%1.0)
    print("     Author: kiranNbishwa (bkgiri@uncg.edu, kirannbishwa01@gmail.com) ")
    print("##################################################")
    print()



    ''' Step 01: Read the input file and and prepare two files as output '''
    # a) One output file contains extended phase-states for the sample of interest (soi)
    # b) another output file contains the lines that have missing data for sample of interest

    ''' Step 01 - A: read the input haplotype file and prepare output files '''
    with open(input_file) as input_data, \
            open("missingdata_" + soi + ".txt", 'w') as missing_data, \
        open('extended_haplotype_'+ soi +'.txt', 'w') as update_data:

        print('reading the input haplotype file "%s" '%input_data.name)
        print('lines that have data missing for sample "%s" is written in the file "%s" '
              %(soi, missing_data.name))
        print('extended haplotype data for sample "%s" will be written in the file "%s" '
              %(soi, update_data.name))


        ''' Step 01 - B: check if "bed file" and "haplotype reference" file are given.
            - then read the "bed file" and "haplotype file" into the memory.
            - these data will be used downstream after reading the haplotype file as "good_data" '''

        # check and load bed file
        if use_bed == 'yes':
            ''' we want to extend phase state only within bed boundries.
                - so, we merge the "input haplotype-file"  with "bed-file". '''
            print()
            print('reading the bed file "%s" ' % (bed_file))
            print('phase extension will be limited to regions declared in the bed file')

            my_bed = pd.read_csv(bed_file, sep='\t', names=['contig', 'start', 'end'])

        else:
            print()
            print('bed file is not provided')
            print('so the phase extension will run throughout the genome')


        # check and load "haplotype reference panel"
        if use_refhap == 'yes':
            hap_panel = pd.read_csv(refhap, sep='\t').drop(['ref', 'alt'], axis=1)

            # also find the sample in refHap panel
            hap_panel_samples = find_samples(list(hap_panel.keys()))

        else:
            print()
            hap_panel_samples = []
            print('haplotype reference panel is not provided')
            print('so the phase extension will run using only the samples in the input file ')


        ''' Step 01 - C: from the input file
        1) For the soi, remove/store the lines that have PI and allele value as empty (.) as "missing_data"
        2) Rest of the data should have meaningful data for soi, and we store it as "good_data".
           - also prepare the list of samples that will be used in phase extension
        3) write the first K1 block - after this we will only need to update K2 block values,
           when reading k1,v1 and k2,v2 consecutive blocks as pairs. '''

        ''' store the good (meaningful) portion of the input_data as variable (Good_Data).
        But, write the missing data as file (** but can also be stored in a variable) '''
        good_data = ''   # empty string

        for lines in input_data:
            if lines.startswith('contig'):
                header = lines.strip('\n').split('\t')

                # find the index positions of sample-of-interest's PI and PG_allele
                soi_PI_index = header.index(soi + '_PI')
                soi_PG_index = header.index(soi + '_PG_al')

                # update the header for good_data and missing_data
                good_data += '\t'.join(header) + '\n'
                missing_data.write('\t'.join(header) + '\n')

                continue


            ''' Now, for the soi if PI and PG are missing (i.e, represented by '.') write it into
            "missing_file.txt", else store it as "good_data" '''
            lines = lines.strip('\n').split('\t')

            # separate the lines with missing data (no "PI" or has "PI" but ambiguous SNP like "*")
            if lines[soi_PI_index] == '.' or lines[soi_PG_index] == '.':
                missing_data.write('\t'.join(lines) + '\n')

            # write the good part of the RB-phased VCF
            elif lines[soi_PI_index] != '.':
                good_data += '\t'.join(lines) + '\n'

            # ** for future: this code can be modified to include non-phased SNP variants.
            # - just remove "lines[soi_PG_index]" to put SNPs with no PI index inside "good_data"

        print()
        print('filtered the lines that have data missing for sample "%s"; check the file "%s" '
              %(soi, missing_data.name))
        print('loaded read-backphased variants onto the memory')


        ''' Prepare a list of tuples of samples from the input data and update it as needed.
            - **Note: the sample list should always include the soi (sample of interest)
                - this is done to include soi rather than introducing a pseudo count. '''
        # ** to do: fix this to include samples from input, vs. hapRef panel vs. custom samples ??
        samples = header.copy()
        sample_list = find_samples(samples)  # these are from input haplotype file

        # update the names in "sample_list" if requested by the user:
        if use_sample == "" or use_sample == 'input':
            sample_list = sample_list

        # use all the samples from hapRefPanel and input samples
        elif use_sample == 'all':
            sample_list = sample_list + hap_panel_samples

        elif use_sample == 'refHap':
            sample_list = hap_panel_samples + [(soi + "_PI", soi +"_PG_al")]

        # use the input sample names, split it and then prepare ..
        # .. the list of tuples of sample "PI" and "PG_al"
        else:
            sample_list = use_sample.split(',')
            sample_list = [((x + '_PI'), (x + '_PG_al')) for x in sample_list] + \
                          [(soi + "_PI", soi +"_PG_al")]

        #print('sample listing to use')
        #print(sample_list)
        #print(len(sample_list))



        ''' Step 02: pipe the data into "pandas", then:
            A) group the data by "contig" which helps in multiprocessing/threading.
              A - optional: if "bed regions" are given add the bed_regions boundries as "start_end"
            B) within each group, group again by "PI keys" of soi and then sort by
               minimum "pos" value for each "PI key"
            C) then pipe the data within each "PI key" for phase-extension computation.'''

        ''' Step 02 - A : read data into "pandas" as dataframe.'''
        good_data = pd.read_table(StringIO(good_data), delimiter='\t')

        # ** only if "good_data" is desired as text output
        #pd.DataFrame.to_csv(good_data, 'good_data_test.txt', sep='\t', header=True, index=False)

        ''' Step 02 - A (** only if reference haplotype is provided) '''
        if use_refhap == "yes":
            # update the "good_data" (i.e, haplotype data)
            print()
            print('merging input data with data from the hap-Reference panel')

            good_data = pd.merge(good_data, hap_panel, on=['contig', 'pos'],
                                     how='left').fillna('.')

            # if haplotype and reference panel merged lines are desired
            #pd.DataFrame.to_csv(good_data, 'hap_and_refPanel_merged.txt', sep='\t',
                                #header=True, index=False)
            del hap_panel

        else:
            print()
            print('Haplotype reference panel is not provided.'
                  '\nOnly using the samples in input (haplotype file) data.')


        '''** only, if bed file is given add the bed reiongs to limit phase extension.
            and group the data by "contig"'''
        if use_bed == 'no':
            # group data only at "contig" level, keep the sort as it is
            print()
            print('no bed file given')
            print('grouping the haplotype file by chromosome/contig')

            good_data_by_group = good_data.groupby('contig', sort=False)

        elif use_bed == 'yes':
            print()
            print('merging the bed boundries from "%s" with the haplotype file "%s" '
                  % (bed_file, input_data.name))

            # merge/intersect the "bed regions" and "haplotype file"
            # then groupy "contig" and "bed regions" by passing it to function "merge_hap_with_bed()"
            good_data_by_group = merge_hap_with_bed(my_bed, good_data)
            # ** for future: we can also run multiprocessing while merging "hap file" with "bed regions"

            del my_bed


        ''' **Optional: compute the statistics of the initial phased file for soi if required '''
        if hapstats == 'yes':
            print()
            print('computing the haplotype statistics of the RBphased data before phase extension')
            # select the colums of interest
            initial_haplotype = good_data[['contig', 'pos', 'ref',
                                           'all-alleles', soi + '_PI', soi + '_PG_al']]

            ''' write this initial haplotype to a file '''
            pd.DataFrame.to_csv(initial_haplotype, 'initial_haplotype_' + soi + '.txt',
                                sep='\t', header=True, index=False)

            # pipe the data to a function to compute haplotype statistics
            compute_haplotype_stats(initial_haplotype, soi, prefix='initial')


        else:
            print()
            print('proceeding to phase-extension without preparing haplotype statistics. ')
            # select the columns of interest
            initial_haplotype = good_data[['contig', 'pos', 'ref',
                                           'all-alleles', soi + '_PI', soi + '_PG_al']]

            ''' write this initial haplotype to a file '''
            pd.DataFrame.to_csv(initial_haplotype, 'initial_haplotype_' + soi + '.txt',
                                sep='\t', header=True, index=False)


        # clear memory
        del good_data, initial_haplotype



        ''' Step 02 - B: Start, multiprocessing/threading - process each contig separately. '''
        print()
        print('Starting multiprocessing using "%i" number of processes ' %(nt))


        ''' Step 02 - B (i)'''
        # take the large dataframe that is grouped by contig and ..
        # .. keep chunks of dataframes as as OrderedDict(list of (keys, Dataframe object))
        df_list = collections.OrderedDict()
        for chr_, data_by_chr in good_data_by_group:
            df_list[chr_] = data_by_chr

        del good_data_by_group  #clearing memory

        print('testing memory ')
        print(type(list(df_list.items())))
        #print(list(df_list.items()))
        print(type(itertools.product(list(df_list.items()))))



        ''' Step 02 - B (ii) : create a pool of process and then feed the list of dataframes to another
            function "groupby_and_read()" to run phase extension.
            - After the data is passed into that function; the steps (02-C to 9) are covered there'''

        p = Pool(nt)  # number of pool to run at once; default at 1
        result = p.map(groupby_and_read, list(df_list.items()))
        # ** for future: use iterator, yield or generators to reduce memory burden while multiprocessing.


        p.close()
        p.join()

        print()
        print("completed haplotype extension for all the chromosomes."
              "time elapsed: '%s' " %(time.time() - time01))
        print('Global maximum memory usage: %.2f (mb)' % current_mem_usage())
        print("merging dataframes together ....." )


        ''' Step 10: merge the returned result (extended haplotype block by each contigs)
            together and write it to an output file. '''
        result_merged = pd.concat(result)
        # write data to the final output file
        pd.DataFrame.to_csv(result_merged, 'extended_haplotype_'+ soi + output + '.txt',
                            sep='\t', index=False)


        ''' Step 11: now, compute the haplotype stats of the phase-extended file '''
        if hapstats == 'yes':
            print()
            print('computing the haplotype statistics of the extended RBphased file.')
            compute_haplotype_stats(result_merged, soi, prefix='final')


    print('end :)')





def groupby_and_read(good_data_by_group):

    print()
    ''' After doing groupby (by chromosome) we pipe in data, for each chromosome '''
    # from earlier pool.map() process the data from each contig/chromosome passes ..
    # .. here as tuple(chromosome, dataframe of that chromosome)
    # we reassign to new variable and clear memory of old variable
    chr_ = good_data_by_group[0]
    contigs_group = good_data_by_group[1]

    time_chr = time.time()
    print('extending haplotype blocks in chromosome/contig %s' %(chr_))

    del good_data_by_group


    # if phase extension is to be limited to provided bed-regions
    if use_bed == 'yes':
        # start another level of grouping of the data by bed interval key
        # (start-end of the bed is used as unique key for grouping)
        good_data_by_bed = contigs_group.groupby('start_end', sort=False)

        print('splitting the contig by bed regions')


        del contigs_group  # clear memory

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
                    ['contig', 'pos', 'ref', 'all-alleles', soi + '_PI', soi + '_PG_al']]

                # ** write void lines separately (if needed)
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

        print('phase-extension completed for contig "%s" in %s seconds' % (chr_, time.time() - time_chr))
        print('\tWorker maximum memory usage: %.2f (mb)' % (current_mem_usage()))
        print()

        return phase_extended_by_chr.sort_values(by=['pos'])


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

    print('phase-extension completed for contig "%s" in %s seconds' %(chr_, time.time()-time_chr))
    print('\tWorker maximum memory usage: %.2f (mb)' % (current_mem_usage()))
    print()

    # return the phase-extended haplotype back to the pool-process ..
    # .. which will be pd.concatenated after all pools are complete
    return phase_extended_by_chr.sort_values(by=['pos'])


''' now, read the two consecutive blocks to do phase extension. '''
def process_consecutive_blocks(contigs_group, soi,
                               chr_, snp_threshold, sample_list, num_of_hets,
                               lods_cut_off):

    print()
    print('now grouping the dataframe using "phased index" values. ')


    ''' Step 02 - D: group dataframe again by "PI keys" of soi and then
       sort by minimum "pos" value for each "PI key".
       - This sorting is necessary because sometimes "haplotype blocks" are like 3-3-3-3  5-5-5  3-3-3-3
          - i.e there are small RBphased blocks within the boundry of larger RBphased block.
          - Not, sure what is causing this (prolly sampling difference of large vs. small chunks in PE reads)
          - This problem should go away in first round of haplotype-extension'''

    contigs_group = contigs_group. \
        assign(New=contigs_group.groupby([soi + '_PI']).
               pos.transform('min')).sort_values(['New', 'pos'])


    ''' Step 03: Now, start reading the "contigs_group" for haplotype-extension.
    A) Store the data as dictionary with 'header' values as keys. Some keys are: chr, pos, sample (PI, PG within sample),
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
    phased_grouped = itertools.groupby(phased_dict, key=lambda x: x[soi + '_PI'])

    ''' Since the dictionary isn't ordered, we return the order using OrderedDictionary '''
    # ** for future: there is room for improvement in here (memory and speed)
    grouped_data = collections.OrderedDict()
    for key, grp in phased_grouped:
        grouped_data[key] = accumulate(grp)

    ''' Clear memory '''
    del phased_dict
    del phased_grouped
    del contigs_group


    print()
    print('Starting MarkovChains for contig %s' % chr_)
    ''' Step 03 - B : now pipe the data for phase extension '''
    ''' Step 03 - B : And, iterate over two consecutive Haplotype-Blocks at once. This is done to obtain all
    possible Haplotype configurations between two blocks. The (keys,values) for first block is represented as
    k1,v2 and for the later block as k2,v2. '''

    ##########################################
    # ** to do -
    initial_haplotype_stats = collections.OrderedDict()  #
    final_haplotype_stats = collections.OrderedDict()  # update using pandas
    #########################################

    ''' Step 03 - B (i): Before running consecutive blocks, we write data from the very first block to the file.
    Reason : Before we start computing and solving the haplotype phase state, we plan to write the
    data for very first block (k1, v1). So, after that, we can solve the relation between two consecutive..
    .. blocks but only write data from 2nd block each time - based on what relation comes out. '''
    very_first_block = [list(grouped_data.items())[0]]


    if len(list(grouped_data.items())) == 1:
        print('there is only one block, so skipping phase extension')

    # write header of the extended phase-block
    extended_haplotype = '\t'.join(['contig', 'pos', 'ref', 'all-alleles', soi +'_PI', soi +'_PG_al']) + '\n'

    # write data/values from very first block.
    for k1, v1 in very_first_block:
        for r1, vals in enumerate(v1[soi + '_PI']):
            new_line = '\t'.join([v1['contig'][r1], v1['pos'][r1], v1['ref'][r1], v1['all-alleles'][r1],
                                  v1[soi + '_PI'][r1], v1[soi + '_PG_al'][r1]]) + '\n'
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
        hap_block1a = [x.split('|')[0] for x in v1[soi + '_PG_al']]  # the left haplotype of block01
        hap_block1b = [x.split('|')[1] for x in v1[soi + '_PG_al']]

        # iterate over the second Haplotype Block, i.e the k2 block and v2 values
        hap_block2a = [x.split('|')[0] for x in v2[soi + '_PG_al']]
        hap_block2b = [x.split('|')[1] for x in v2[soi + '_PG_al']]


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
            for xi in range(len(v2[soi + '_PI'])):
                new_line = '\t'.join([v2['contig'][xi], v2['pos'][xi], v2['ref'][xi], v2['all-alleles'][xi],
                                      k2, hapb1a_hapb2a[1][xi] + '|' + hapb1b_hapb2b[1][xi]]) + '\n'
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
        number_of_snp_in_soi_v1 = len([x for x in v1[soi + '_PG_al'] if len(x) == 3])
        number_of_snp_in_soi_v2 = len([x for x in v2[soi + '_PG_al'] if len(x) == 3])

        # print('number of SNPs: ', NumSNPsInsoi_v1, NumSNPsInsoi_v2)
        if number_of_snp_in_soi_v1 < snp_threshold \
                or number_of_snp_in_soi_v2 < snp_threshold:
            for xth, vals in enumerate(v2[soi + '_PI']):
                new_line = '\t'.join([v2['contig'][xth], v2['pos'][xth], v2['ref'][xth], v2['all-alleles'][xth],
                                      k2, hapb1a_hapb2a[1][xth] + '|' + hapb1b_hapb2b[1][xth]]) + '\n'
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
          # .. downsizing the effect of recombination.
        lhfc_f, lhsc_f = \
            compute_transitions_probs(soi, sample_list, k1, k2, v1, v2, num_of_hets,
                                      hapb1a_hapb2a, hapb1b_hapb2b,
                                      hapb1a_hapb2b, hapb1b_hapb2a, orientation=reversed)

        #### for reverse chain   ########
        # set "orientation=lambda..." just passes a null value keeping orientation as it is.
        lhfc_r, lhsc_r = compute_transitions_probs \
            (soi, sample_list, k1_r, k2_r, v1_r, v2_r, num_of_hets,
             hapb1a_hapb2a_r, hapb1b_hapb2b_r,
             hapb1a_hapb2b_r, hapb1b_hapb2a_r, orientation=lambda x: x)



        ''' Step 05-06 are inside the function "compute_transitions_probs()". The values
        (lhfc_f, lhsc_f, lhfc_r, lhsc_r) returned from this function is then used in Step 07. '''


        ''' Step 07 :  previous (Step 06) returns the likelyhoods and/or LODs score for both "parallel"
        and alternate configurations (for both forward and reverse algorithm).
        - We now extend the phase states by comparing LODs score against  cutoff-values.'''

        ''' Step 07 - A(i): calculate the average of the likelyhoods, odds and then log2 of odds. '''
        # average of the likelyhooods for first vs. second configuration
        # (from both forward and reverse algorithm)
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

    print()
    print('extracting bed regions and position of the haplotype file')

    c1 = my_bed.contig.values
    s1 = my_bed.start.values
    e1 = my_bed.end.values
    c2 = good_data['contig'].values
    pos2 = good_data.pos.values

    # now, find the intersecting positions (between haplotype ("pos" values) and bed regions).
    overlap = (
        (c2[:, None] == c1) & ((pos2[:, None] >= s1) & (pos2[:, None] <= e1)))

    i, j = np.where(overlap)


    '''intersect the haplotype file with bed file '''

    print('intersecting bed regions with haplotype file, and '
          'creating new haplotype boundries for each contig ... ')

    df_bed__hap_interesect = pd.DataFrame(
        np.column_stack([good_data.values[i], my_bed.values[j]]),
        columns=good_data.keys().append(my_bed.keys()))

    # drop any duplicate columns
    df_bed__hap_interesect = df_bed__hap_interesect.T.drop_duplicates().T

    df_bed__hap_interesect['start_end'] = pd.DataFrame(
        df_bed__hap_interesect.apply(lambda x: str(x.start) + '-' + str(x.end), axis=1))

    # only keep "contig, pos and intersection ("start-end")"
    df_bed__hap_interesect = df_bed__hap_interesect[['contig', 'pos', 'start_end']]
    # if regions of intersection are of interest
    #pd.DataFrame.to_csv(df_bed__hap_interesect, 'df_bed and hap intersect.txt', sep='\t', index=None, header=True)


    # now, merge the df(good data) with the df (intersected bed and hap) to create updated
      # haplotype file.
    # Any part that is intersected is assigned unique-keys of "start-end" values. The non-intersected
      # part are filled with "NA"
    data_frames = [good_data, df_bed__hap_interesect]

    print('merging bed file with haplotype file')

    # fill the non-merging lines with "void"
    good_data_update = reduce(lambda left, right: pd.merge(
        left, right, on=['contig', 'pos'], how='outer'), data_frames).fillna('na')
    # ** for future: assigning an unique identifier like "NA01, NA02" for non-intersecting lines as blocks.


    # if updated "haplotype file" with intersected "haplotype file" and "bed file" is of interest
    #pd.DataFrame.to_csv(good_data_update, 'df_bed and hap merged.txt', sep='\t', index=None,
                        #header=True)


    print('grouping the haplotype file by chromosome/contigs .... ')
    good_data_by_group = good_data_update.groupby('contig', sort=False)

    # clear memory
    del c1, c2, s1, e1, pos2, i, j, good_data, good_data_update, df_bed__hap_interesect, my_bed

    return good_data_by_group



''' function to compute transition probs between two consecutive blocks. '''
def compute_transitions_probs(soi, sample_list, k1, k2, v1, v2, num_of_hets,
                              hapb1a_hapb2a, hapb1b_hapb2b,
                              hapb1a_hapb2b, hapb1b_hapb2a,orientation):



    print()
    ''' Step 05 : Start preping the first order markov transition matrix and probabilities for the
    consecutive haplotypes between two blocks.
    - To store the sum-values of the product of the transition probabilities.
    These sum are added as the product-of-transition is returned by nested for-loop;
    from the loop "for m in range(....)" '''

    # empty variable for storing "sum of the product of transition probabilities" between two blocks.
    # ** to do - we can either "max sum" or "max product" the likelyhoods, default is "max-product".
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
    for n in orientation(range(len(v1[soi + '_PI']))):

        ''' Skip computation if n'th item is InDel or has "*" allele in the genotype.
        - These indels are only removed from computation part but are added in the final output file.
        - The Indels are phased based on which SNPs they hitchhike with. '''
        if len(v1[soi + '_PG_al'][n]) > 3 \
            or '*' in v1[soi + '_PG_al'][n]:
            continue

        ''' only use certain number of HetVars to compute likely hood. This saves computation burden.
            - set at default at 200 heterozygous variants. '''
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

        ''' Creating variables to store the product-values of the transition probabilities.
        These are updated for each level of "n" paired with each level of "m". '''
        # potp -> product of transition probability
        potp_hapb1a_b2a = 1
        potp_hapb1b_b2b = 1

        potp_hapb1a_b2b = 1
        potp_hapb1b_b2a = 1


        ''' Step 05 - A (ii) :
        - Now, calculate the initial nucleotide counts at "n-th" level,
        before computing the transition counts.
        - only calculated from "v1" at n-th level and only once for each parse/iteration '''
        for (x, y) in sample_list:
            for nucleotide in 'ATGC':
                nucleotide_count_dict[nucleotide] += v1[y][n].count(nucleotide)


        ''' Step 05 - B : Count number of transition from each nucleotide (n-th) to each nucleotide (m-th).
        - Now, we read nucleotides (ATGC) at each level of "m"
        to compute the transition from each level of "n". '''

        # to control certain number of hetVars to include in computation
        num_of_het_site_at_m = 0
        for m in range(len(v2['ms02g_PI'])):    # m-ranges from block02

            ''' Like at "n-th" level, skip if InDel present at this "m-th" level.
            But InDel will be phased along with hitchhiking SNPs. '''
            if len(v2[soi + '_PG_al'][m]) > 3\
                or '*' in v2[soi + '_PG_al'][m]:
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
                            pos  PI   PG_al          pos  PI   PG_al
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
                    compute_transition_probs(transition_count_dict[(from_, to)], nucleotide_count_dict[from_])


            ''' Step 06 - A (ii) : find observed configuration for soi at "n-th" and "m-th" level
            i.e (from_, to). '''
            hapb1a_hapb2a_transition = (hapb1a_hapb2a[0][n], hapb1a_hapb2a[1][m])
            hapb1b_hapb2b_transition = (hapb1b_hapb2b[0][n], hapb1b_hapb2b[1][m])

            hapb1a_hapb2b_transition = (hapb1a_hapb2b[0][n], hapb1a_hapb2b[1][m])
            hapb1b_hapb2a_transition = (hapb1b_hapb2a[0][n], hapb1b_hapb2a[1][m])

            ''' Step 06 - B : computes the "product of transition probabilities" from "n-th" to
             several levels of "m" using for-loop.
             - ** because the below code is indented one level below "for m in range(len(v2['ms02g_PI']))".
             - ** no need to add pseudo-counts, because if no haplotypes are observed in any samples except soi,
               the prob(from_, to) for each configuration will be 1/4 there by nullifying the likelyhoods to "1".
            '''
            potp_hapb1a_b2a *= transition_prob_dict[hapb1a_hapb2a_transition]
            potp_hapb1b_b2b *= transition_prob_dict[hapb1b_hapb2b_transition]
            potp_hapb1a_b2b *= transition_prob_dict[hapb1a_hapb2b_transition]
            potp_hapb1b_b2a *= transition_prob_dict[hapb1b_hapb2a_transition]


        ''' Step 06 - C : compute the max sum or max product of the transition probabilities
        across several levels of "n" to several levels of "m".
        ** Note: we can either do "max sum" or "max product".
        So, "cul_of_pt_hapb1a_b2a" is the sum of the likelyhoods of hapBlock1A being phased with hapBlock2A'''

        if maxed_as == '+':
            cul_of_pt_hapb1a_b2a += potp_hapb1a_b2a
            cul_of_pt_hapb1b_b2b += potp_hapb1b_b2b

            cul_of_pt_hapb1a_b2b += potp_hapb1a_b2b
            cul_of_pt_hapb1b_b2a += potp_hapb1b_b2a

        elif maxed_as == '*':
            cul_of_pt_hapb1a_b2a *= (potp_hapb1a_b2a)
            cul_of_pt_hapb1b_b2b *= (potp_hapb1b_b2b)

            cul_of_pt_hapb1a_b2b *= (potp_hapb1a_b2b)
            cul_of_pt_hapb1b_b2a *= (potp_hapb1b_b2a)


    ''' Step 06 - D : Now, compute the "Odds ratio" and "log2 of the Odds"
        of each possible haplotype configuration, for both 1st vs. 2nd configurations '''

    ''' Step 06 - D(i) : compute the likely hood of first configuration (lhfc) vs. second configuration (lhsc)
    First Configuration = hapb1a with hapb2a, and hapb1b with hapb2b.
    Second Configuration = hapb1a with hapb2b, and hapb1b with hapb2a. '''
    lhfc = Decimal(cul_of_pt_hapb1a_b2a * cul_of_pt_hapb1b_b2b)
    lhsc = Decimal(cul_of_pt_hapb1a_b2b * cul_of_pt_hapb1b_b2a)

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

    print()
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

            for xth in range(len(v2[soi + '_PI'])):
                new_line = '\t'.join([v2['contig'][xth], v2['pos'][xth],  v2['ref'][xth], v2['all-alleles'][xth],
                                      k2_new, hapb1a_hapb2a[1][xth] + '|' + hapb1b_hapb2b[1][xth]]) + '\n'
                extended_haplotype += new_line


        # if the |-lods score| is above cut-off, phase needs extension in "alternate" configuration
        elif lods2_score_1st_config < -lods_cut_off:  # i.e 2^4 times less likely
            k2_new = k1
            flipped = 'yes'

            for xth in range(len(v2[soi + '_PI'])):
                new_line = '\t'.join([v2['contig'][xth], v2['pos'][xth],  v2['ref'][xth], v2['all-alleles'][xth],
                                      k2_new, hapb1b_hapb2b[1][xth] + '|' + hapb1a_hapb2a[1][xth]]) + '\n'
                extended_haplotype += new_line

        # if the |lods score| does not pass cut-off threshold in either direction, there is no phase extension
        # so, reset the "checker-variables" and write/store data in original configuration as input data
        else:  # i.e 4 times less likely
            k2_new = ''
            flipped = ''

            # print('no phase extension, no flip')  # marker for debugging
            for xth in range(len(v2[soi + '_PI'])):
                new_line = '\t'.join([v2['contig'][xth], v2['pos'][xth], v2['ref'][xth], v2['all-alleles'][xth],
                                      k2, hapb1a_hapb2a[1][xth] + '|' + hapb1b_hapb2b[1][xth]]) + '\n'
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
                for xth in range(len(v2[soi + '_PI'])):
                    new_line = '\t'.join([v2['contig'][xth], v2['pos'][xth],  v2['ref'][xth], v2['all-alleles'][xth],
                                          k2_new, hapb1a_hapb2a[1][xth] + '|' + hapb1b_hapb2b[1][xth]]) + '\n'
                    extended_haplotype += new_line

            # previous config was extended but not flipped. So, now when the phase extends
              # the configuration is kept "alternate" because lods2 score is negative.
            elif lods2_score_1st_config < -lods_cut_off:  # i.e 4 times less likely
                k2_new = k2_new
                flipped = 'yes'

                #print ('1_A is phased with 2_B, yes flipped, k2:%s \n' % k2_new) # marker for debugging
                for xth in range(len(v2[soi + '_PI'])):
                    new_line = '\t'.join([v2['contig'][xth], v2['pos'][xth],  v2['ref'][xth], v2['all-alleles'][xth],
                                          k2_new, hapb1b_hapb2b[1][xth] + '|' + hapb1a_hapb2a[1][xth]]) + '\n'
                    extended_haplotype += new_line

            # if phase state isn't extending, write the data as is and reset the "checker-variables"
            else:
                k2_new = ''
                flipped = ''

                #print ('no phase extension, no flip, k2 %s\n' % k2)  # marker for debugging
                for xth in range(len(v2[soi + '_PI'])):
                    new_line = '\t'.join([v2['contig'][xth], v2['pos'][xth], v2['ref'][xth], v2['all-alleles'][xth],
                                          k2, hapb1a_hapb2a[1][xth] + '|' + hapb1b_hapb2b[1][xth]]) + '\n'
                    extended_haplotype += new_line

        # if previous block were extended but flipped
        elif flipped == 'yes':

            # previous config was extended but flipped, so to maintain it's phase state
              # it must stay in "parallel" config this time by flipping with previous extension.
            if lods2_score_1st_config > lods_cut_off:  # i.e 4 times more likely
                k2_new = k2_new
                flipped = 'yes'

                #print('\n1_A is phased with 2_A, but yes flipped, k2:%s \n' % k2_new)  # marker for debug
                for xth in range(len(v2[soi + '_PI'])):
                    new_line = '\t'.join([v2['contig'][xth], v2['pos'][xth], v2['ref'][xth], v2['all-alleles'][xth],
                                          k2_new, hapb1b_hapb2b[1][xth] + '|' + hapb1a_hapb2a[1][xth]]) + '\n'
                    extended_haplotype += new_line


            # previous config was extended but flipped, so to maintain it's phase state
              # it must stay in "alternate" config this time by not flipping with previous extension.
            elif lods2_score_1st_config < -lods_cut_off:  # i.e 4 times less likely
                k2_new = k2_new
                flipped = 'no'

                #print('\n1_A is phased with 2_B, but not flipped, k2:%s \n' % k2_new)  # debug
                for xth in range(len(v2[soi + '_PI'])):
                    new_line = '\t'.join([v2['contig'][xth], v2['pos'][xth], v2['ref'][xth], v2['all-alleles'][xth],
                                          k2_new, hapb1a_hapb2a[1][xth] + '|' + hapb1b_hapb2b[1][xth]]) + '\n'
                    extended_haplotype += new_line

            # previous extension was phased, but this one didn't
            # so, write the data as it is and reset the "checker-variables"
            else:
                k2_new = ''
                flipped = ''

                #print('\nno phase extension, no flip (reset), k2 %s\n' % k2)  # marker for debugging
                for xth in range(len(v2[soi + '_PI'])):
                    new_line = '\t'.join([v2['contig'][xth], v2['pos'][xth], v2['ref'][xth], v2['all-alleles'][xth],
                                                 k2, hapb1a_hapb2a[1][xth] + '|' + hapb1b_hapb2b[1][xth]]) + '\n'
                    extended_haplotype += new_line


    return k2_new, flipped, extended_haplotype


''' function that computes several stats from haplotype block before extension.
These data can be used to compare the improvements in phase extension, and plotting histogram.
**Note - to do ?? : Write a code (here) to compute the size of the haplotype before phasing.
This can include the size based on genome co-ordinate of the pos in each block.
This can also include the average number of variants per haplotype before/after phase extension.
and make aplot out of it - using pyPlot, MatlibPlot or R. '''
def compute_haplotype_stats(hap_data, soi, prefix) :
    print()
    ''' this function computes the stats of the haplotype file both before and after phase extension.
        - we pass in "prefix" as "initial" vs. "final" for the haplotype block before and after
          extension to compute statistics for respective phased file. '''

    ''' Step 09 - A : write this haplotype information in a more concise format.
        - in this step we are restructuring the data so proper statistics can be handled properly. '''
    hap_stats = hap_data.groupby(['contig', soi +'_PI'], sort=False)\
        .size().rename('num_Vars_by_PI')\
        .reset_index(level=1).astype(str)\
        .groupby(level=0).agg(','.join).reset_index()

    ''' add other required columns '''
    hap_stats['range_of_PI'] = hap_data.groupby(['contig', soi +'_PI'], sort=False)['pos']\
        .apply(lambda g: g.max() - g.min()).rename('range').reset_index(level=1)\
        .astype(str).groupby(level=0).agg(','.join).reset_index()['range']

    hap_stats['total_haplotypes'] = hap_stats\
        .apply(lambda row: len(row[soi +'_PI'].split(',')), axis=1)

    hap_stats['total_Vars'] = hap_stats \
        .apply(lambda row: sum([int (i) for i in row.num_Vars_by_PI.split(',')]), axis=1)

    ''' write this concise statistics to a file. '''
    pd.DataFrame.to_csv(hap_stats,  prefix+'_haplotype_stats_' + soi + '.txt',
                        sep='\t', header=True, index=False)


    ## ** - plot haplotype STATs - may need improvement.
    ''' Step 09 - B : Now, use the above dataframe to plot several statistical plots from here on.
        - plot the haplotype statistics using matplotlib, pylab or ggplot (python).
        - start with "with open() .... " to open a file, write the plot and then close. '''
    # some good links for hist-plots: http://pandas.pydata.org/pandas-docs/version/0.18/visualization.html

    # plots total number of variants per-chromosome
    with open('total_vars_'+ soi +'_'+ prefix+'.png', 'wb') as fig_initial:
        hap_stats.plot(x='contig', y='total_Vars', kind='bar')
        plt.xlabel('chromosomes')
        plt.ylabel('number of variants')
        plt.suptitle('number of variants for each chromosome')
        plt.savefig(fig_initial)
        # plt.show()  # optional: to show the plot to the console


    # plots total number of haplotypes per-chromosome
    with open('total_haps_'+ soi +'_'+ prefix+'.png', 'wb') as fig_initial:
        hap_stats.plot(x='contig', y='total_haplotypes', kind='bar')
        plt.xlabel('chromosomes')
        plt.ylabel('number of haplotypes')
        plt.suptitle('number of haplotypes for each chromosome')
        plt.savefig(fig_initial)


    # plot histogram of haplotype size (by number of variants) per-chromosome
        # make different plots for each chromosome, but X-axis is shared.
    with open('hap_size_byVar_'+ soi +'_'+ prefix+'.png', 'wb') as fig_initial:
        fig, ax = plt.subplots(nrows=len(hap_stats), sharex=True)
        for i, data in hap_stats.iterrows():
            # first convert data to list of integers
            data_i = [int(x) for x in data['num_Vars_by_PI'].split(',')]
            ax[i].hist(data_i, label=str(data['contig']), alpha=0.5)
            ax[i].legend()
        plt.savefig(fig_initial)


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
    # so, transition from "A" to any (A,T,G,C) is technically zero, but theoretically it should be 1/16.
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
    sample_list = [x.split('_')[0] for x in samples if '_' in x]
    sample_list = [x for x in sample_list if not (x in seen or seen.add(x))]
    sample_list = [((x + '_PI'), (x + '_PG_al')) for x in sample_list]


    return sample_list

    # ** Note: this returns data as:
    # sample_list = [('ms01e_PI', 'ms01e_PG_al'), ('ms02g_PI', 'ms02g_PG_al'),
               #('ms03g_PI', 'ms03g_PG_al'), ('ms04h_PI', 'ms04h_PG_al')]


''' to monitor memory '''
def current_mem_usage():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.



if __name__ == '__main__':
    main()