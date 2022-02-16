import os
import collections
import csv
import time
import itertools
import shutil


from decimal import Decimal
from functools import reduce, partial
from io import StringIO
from multiprocessing import Pool

import numpy as np
import pandas as pd

from phase_extender.hapstats import compute_haplotype_stats
from phase_extender.compute_score import compute_maxLh_score, extend_phase_state
from phase_extender.utils import accumulate, current_mem_usage


def phase_converter(
    soi,
    outputdir,
    nt,
    input_file,
    lods_cut_off,
    snp_threshold,
    num_of_hets,
    maxed_as,
    bed_file,
    refhap,
    use_sample,
    hapstats,
    writelod,
    addmissingsites,
):

    """Assign the number of process - this is the optimal position to start multiprocessing !
    **note: number of process should be declared after all the global variables are declared,
    because each pool will need to copy the variable/value of global variables."""
    pool = Pool(processes=nt)  # number of pool to run at once; default at 1

    """ Step 01: Read the input file and and prepare two files as output """
    # a) One output file contains extended phase-states for the sample of interest (soi)
    # b) another output file contains the lines that have missing data for sample of interest
    data = pd.read_csv(input_file, sep="\t")
    data_header = list(data.columns)
    pg_al_set = {al for al in data_header if al.endswith(":PG_al")}
    pi_set = {pi for pi in data_header if pi.endswith(":PI")}
    soi_PI_index = soi + ":PI"
    soi_PG_index = soi + ":PG_al"

    # check if soi is in header
    if not soi_PI_index in pi_set:
        assert False, "soi pi index is not found"

    if not soi_PG_index in pg_al_set:
        assert False, "soi pg index is not found"

    os.makedirs(outputdir, exist_ok=True)
    missing_fpath = outputdir + "/" + "missingdata_" + soi + ".txt"

    missing = data[(data[soi_PI_index] == ".") | (data[soi_PG_index] == ".")]
    missing.to_csv(
        missing_fpath,
        sep="\t",
        index=False,
    )
    good_data = pd.concat([data, missing]).drop_duplicates(keep=False)

    """ Step 01 - B: check if "bed file" and "haplotype reference" file are given.
        - then read the "bed file" and "haplotype file" into the memory.
        - these data will be used downstream after reading the haplotype file as "good_data" """

    # check and load bed file
    if bed_file:
        """we want to extend phase state only within bed boundries.
        - so, we merge the "input haplotype-file"  with "bed-file"."""
        my_bed = pd.read_csv(bed_file, sep="\t", names=["CHROM", "start", "end"])
        my_bed["CHROM"] = my_bed["CHROM"].astype(
            str
        )  # setting CHROM column as string type ..
        #  this is necessary because there has been problem with groupby operations downstream

    else:
        print("# Genomic bed file is not provided ... ")

    # check and load "haplotype reference panel"
    if refhap:
        hap_panel = pd.read_csv(refhap, sep="\t").drop(["REF", "ALT"], axis=1)
        hap_panel["CHROM"] = hap_panel["CHROM"].astype(
            str
        )  # setting CHROM as string type data

        # also find the sample in refHap panel
        hap_panel_samples = find_samples(list(hap_panel.keys()))

    else:
        hap_panel_samples = []
        print("# Haplotype reference panel is not provided ... ")
        print(
            "  So, phase extension will run using the samples available in the input haplotype file. "
        )

    """ Step 01 - D: Prepare the samples to use the data from. """
    """ Prepare a list of tuples of samples (PI, PG_al) from the input data and update it as needed.
        - **Note: the sample list should always include the soi (sample of interest)
            - this is done to include observation from soi rather than introducing a pseudo count
              when transition is missing from some observation (n to m). """
    sample_list = find_samples(data_header)  # returns data from "input haplotype file"

    # update the names in "sample_list" if other samples are requested by the user:
    if use_sample == "" or use_sample == "input":
        sample_list = sample_list

    # use all the samples from hapRefPanel and input samples
    elif use_sample == "all":
        sample_list = sample_list + hap_panel_samples

    elif use_sample == "refHap":
        sample_list = hap_panel_samples + [
            (soi + ":PI", soi + ":PG_al")
        ]  # add the self sample name to account ..
        # .. for missing observations instead of using pseudo count

    # if specific select samples are of interest, split the sample names and then prepare ..
    # .. the list of tuples of sample "PI" and "PG_al"
    else:
        sample_list = use_sample.split(",")
        sample_list = [((x + ":PI"), (x + ":PG_al")) for x in sample_list] + [
            (soi + ":PI", soi + ":PG_al")
        ]

    """ Step 02: pipe the data into "pandas", then:
        A) group the data by "contig" which helps in multiprocessing/threading.
            A - optional: if "bed regions" are given add the bed_regions boundries as "start_end"
        B) within each group, group again by "PI keys" of soi and then sort by
            minimum "POS" value for each "PI key"
        C) then pipe the data within each "PI key" for phase-extension computation."""

    """ Step 02 - A : read good part of the data into "pandas" as dataframe."""
    # good_data = pd.read_table(StringIO(good_data), delimiter='\t')
    good_data["CHROM"] = good_data["CHROM"].astype(
        str
    )  # setting CHROM as string type data # this is necessary
    # to maintain proper groupby downstream

    # ** only if "good_data" is desired as text output
    # pd.DataFrame.to_csv(good_data, 'good_data_test.txt', sep='\t', header=True, index=False)

    """ Step 02 - A (add on - i) ** merge reference haplotype if provided """
    if refhap:
        # update the "good_data" (i.e, haplotype data)
        print("Merging input haplotype data with data from the hap-reference panel")

        good_data = pd.merge(
            good_data, hap_panel, on=["CHROM", "POS"], how="left"
        ).fillna(".")
        good_data.sort_values(by=["CHROM", "POS"], inplace=True)

        # if haplotype and reference panel merged lines are desired
        # pd.DataFrame.to_csv(good_data, 'hap_and_refPanel_merged.txt', sep='\t',
        # header=True, index=False)
        del hap_panel

    else:
        print(
            "# Haplotype reference panel is not provided....\n"
            '  - Only using the samples in the input ("%s") data.' % (input_file)
        )

    """ Step 02 - A (add on - ii) ** merge bed-regions if provided to limit phase extension
                                        and group the data by "contig". """
    if not bed_file:
        # group data only at "contig" level, keep the sort as it is
        print("# No bed file is given ... ")
        print("  - So, grouping the haplotype file only by chromosome (contig)")

        good_data_by_group = good_data.groupby("CHROM", sort=False)

    elif bed_file:
        print(
            '# Merging the bed boundries from "%s" with the input haplotype file ... "%s" '
            % (bed_file, input_file)
        )

        # merge/intersect the "bed regions" and "haplotype file"
        # then groupy "contig" and "bed regions" by passing it to function "merge_hap_with_bed()"
        good_data_by_group = merge_hap_with_bed(my_bed, good_data)
        # ** for future: we can also run multiprocessing while merging "hap file" with "bed regions"
        del my_bed

    """ Step 02 - A (**add on - iii):
        - Write the initial haplotype data.
        - Compute the statistics of the initial phased file for SOI if required """

    print(
        '# Writing initial haplotype for sample "%s" in the file "%s" '
        % (soi, "initial_haplotype_" + soi + ".txt")
    )

    # select the colums of interest
    initial_haplotype = good_data[
        ["CHROM", "POS", "REF", "all-alleles", soi + ":PI", soi + ":PG_al"]
    ].sort_values(by=["CHROM", "POS"])

    # write this initial haplotype to a file
    initial_haplotype.to_csv(
        outputdir + "/" + "initial_haplotype_" + soi + ".txt",
        sep="\t",
        header=True,
        index=False,
    )

    if hapstats == "yes":
        print(
            "  - Computing the descriptive statistics of the haplotype data before phase extension"
        )

        # pipe the data to a function to compute haplotype statistics
        compute_haplotype_stats(initial_haplotype, soi, "initial", outputdir)
    else:
        print(
            "  - Proceeding to phase-extension without preparing descriptive statistics of initial haplotype state."
        )

    """ Step 02 - B: - Split the data (grouped by chromosome (contig) values.
                        - Store data in disk or memory.
                        - Multiprocess each chunks separately """
    print()
    print('# Starting multiprocessing using "%i" processes ' % (nt))

    # ** new method: create a folder to store the data to disk (rather than memory)
    # ** (see old method for comparison)
    # if os.path.exists('chunked_Data_' + soi):
    #     shutil.rmtree('chunked_Data_' + soi, ignore_errors=False, onerror=None)
    # os.makedirs('chunked_Data_' + soi + '/', exist_ok=True)

    """ Step 02 - B (i)"""

    ################### old method - ** if possible reuse this method in future.
    # take the large dataframe that is grouped by contig and ..
    # .. keep chunks of dataframes as as OrderedDict(list of (keys, Dataframe object))
    # df_list = collections.OrderedDict()
    ########################################

    # # new method - storing data to disk
    # for chr_, data_by_chr in good_data_by_group:
    #     chunked_path = 'chunked_Data_' + soi + '/' + soi + ':' + str(chr_)
    #     data_by_chr.to_csv(chunked_path,sep='\t', index=False, header=True)

    # clear memory - does it do it's job ** ??
    # initial_haplotype = None; good_data = None; input_file = None
    # # good_data_by_group = None; samples = None
    # data_by_chr = None
    # del initial_haplotype, good_data, input_file, good_data_by_group, samples, data_by_chr

    """ Now, pipe the procedure to next function for multiprocessing (i.e Step 02 - C) """
    multiproc(
        sample_list,
        pool,
        hapstats,
        soi,
        outputdir,
        addmissingsites,
        bed_file,
        snp_threshold,
        num_of_hets,
        lods_cut_off,
        maxed_as,
        writelod,
        good_data_by_group,
    )

    # remove the chunked data folder ** (this can be retained if need be)
    # shutil.rmtree('chunked_Data_' + soi, ignore_errors=False, onerror=None)

    print("End :)")


def multiproc(
    sample_list,
    pool,
    hapstats,
    soi,
    outputdir,
    addmissingsites,
    bed_file,
    snp_threshold,
    num_of_hets,
    lods_cut_off,
    maxed_as,
    writelod,
    good_df_by_grp,
):

    print()
    """ Step 02 - C: Start, multiprocessing/threading - process each contig separately. """
    """ Step 02 - C (ii) : create a pool of process and then feed the list of dataframes to another
        function "groupby_and_read()" to run phase extension.
        - After the data is passed into that function; the steps (02-C to 9) are covered there"""

    # path='chunked_Data_' + soi   # create a temp folder
    # file_path = [(item, sample_list) for item in list(os.path.join(path, part) for part in os.listdir(path))]
    # print(file_path)
    # breakpoint()
    # file_path = [item for item in list(os.path.join(path, part) for part in os.listdir(path))]

    ## ** to do: Add "sort" method in "file_path" to read data in order. This way we can save ..
    # time/memory while doing sorting within pandas dataframe.
    # This sort method is available in "phase-Stitcher"
    df_list = (good_df_by_grp.get_group(x) for x in good_df_by_grp.groups)

    partial_group = partial(
        groupby_and_read,
        bed_file=bed_file,
        soi=soi,
        snp_threshold=snp_threshold,
        sample_list=sample_list,
        num_of_hets=num_of_hets,
        lods_cut_off=lods_cut_off,
        maxed_as=maxed_as,
        writelod=writelod,
    )

    result = pool.imap(partial_group, df_list)
    # result = pool.starmap( groupby_and_read, *(file_path, use_bed, soi, snp_threshold, sample_list, num_of_hets, lods_cut_off, maxed_as) )
    pool.close()
    pool.join()
    pool.terminate()

    print("Global maximum memory usage: %.2f (mb)" % current_mem_usage())
    print("Merging dataframes together .....")

    """ Step 10: merge the returned result (extended haplotype block by each contigs)
        together and write it to an output file. """
    result_merged = pd.concat(result).sort_values(
        by=["CHROM", "POS"], ascending=[True, True]
    )
    # result_merged = pd.DataFrame(result_merged).sort_values(by=['contig', 'pos'], ascending=[1,1])

    # write data to the final output file
    pd.DataFrame.to_csv(
        result_merged,
        outputdir + "/" + "extended_haplotype_" + soi + ".txt",
        sep="\t",
        index=False,
    )

    print(
        'Extended haplotype data for sample "%s" is written in the file "%s". '
        % (soi, "extended_haplotype_" + soi + ".txt")
    )

    """ Step 11: now, compute the descriptive stats of haplotype after phase extension. """
    if hapstats == "yes":
        print()
        print("Computing the descriptive statistics of the extended haplotype file.")
        compute_haplotype_stats(result_merged, soi, "final", outputdir)

    else:
        print(
            "Skipping the preparation of descriptive statistics of extended haplotype."
        )

    print()
    print("Run is complete for all the chromosomes (contigs)")

    """ Step 12: if phase extended data is to be merged with non phased SNPs and missing sites.
        - we can use this data to run another round of phase extension of singletons (SNP and Indels).
        - it also helps to control recursive phase extension. """
    print()
    print("writing singletons and missing sites to extended haplotype")
    if addmissingsites == "yes":
        missed_data_toadd = pd.read_csv(
            outputdir + "/" + "missingdata_" + soi + ".txt",
            sep="\t",
            usecols=["CHROM", "POS", "REF", "all-alleles", soi + ":PI", soi + ":PG_al"],
        )
        missed_data_toadd["CHROM"] = missed_data_toadd["CHROM"].astype(str)

        # "result_merged" only contains RBphased data for SOI
        new_df = pd.concat([result_merged, missed_data_toadd]).sort_values(
            by=["CHROM", "POS"], ascending=[True, True]
        )
        new_df = new_df[
            ["CHROM", "POS", "REF", "all-alleles", soi + ":PI", soi + ":PG_al"]
        ]

        # write data to the final output file
        pd.DataFrame.to_csv(
            new_df,
            outputdir + "/" + "extended_haplotype_" + soi + "_allsites.txt",
            sep="\t",
            index=False,
            na_rep=".",
        )

        del new_df

    del result_merged


def groupby_and_read(
    file_path,
    bed_file,
    soi,
    snp_threshold,
    sample_list,
    num_of_hets,
    lods_cut_off,
    maxed_as,
    writelod,
):
    # good_data_by_contig = open(file_path[0], 'r')
    # chr_ = good_data_by_contig.name.split(':')[-1]
    # sample_list = file_path[1]
    # contigs_group = pd.read_csv(StringIO(good_data_by_contig.read()), sep='\t')
    contigs_group = file_path
    chr_ = contigs_group["CHROM"].values[0]

    """ After doing groupby (by chromosome) we pipe in data, for each chromosome """

    time_chr = time.time()  # to time the process in each chromosome
    # print('## Extending haplotype blocks in chromosome (contig) %s' %(chr_))

    # del good_data_by_contig

    # if phase extension is to be limited to provided bed-regions
    if bed_file:
        # start another level of grouping of the data by bed interval key
        # (start-end of the bed is used as unique key for grouping)
        print("Splitting the contig by bed regions")
        good_data_by_bed = contigs_group.groupby("start_end", sort=False)

        # clear memory
        # contigs_group = None
        # del contigs_group

        """ store the dataframes split by "start_end" key in "df_list_by_bed".
            Then pass this to phase-extension process. """

        # store dataframe/s outside bed region
        # ** - for future (convert this into collection.OrderedDict of ..
        # .. dataframes if data needs to be stored by multiple void blocks)
        data_by_bed_void = pd.DataFrame()

        # to store dataframes that are split by bed-regions but are within boundry of interest
        df_list_by_bed = []

        for bed_id, data_by_bed in good_data_by_bed:
            if "void" in bed_id:
                # write data from these columns/row as a variable
                data_by_bed_void = data_by_bed[
                    ["CHROM", "POS", "REF", "all-alleles", soi + ":PI", soi + ":PG_al"]
                ]

                # ** write void lines separately (if needed) - optimize in the future.
                # pd.DataFrame.to_csv(data_by_bed_void, 'haplotype_file_excluded_regions.txt',
                # sep='\t', header=None, index=False)

            else:
                # now pass this dataframe to a function "process_consecutive_blocks()"
                # i.e, Step 02-D - that reads two consecutive keys at once
                # see **1, for future - there is possibility to multiprocess (phase-extension) this method
                # phase_extended_by_bed = process_consecutive_blocks(
                #     data_by_bed, soi, chr_, snp_threshold,
                #     sample_list, num_of_hets, lods_cut_off, maxed_as)
                phase_extended_by_bed = process_consecutive_blocks(
                    contigs_group,
                    soi,
                    chr_,
                    snp_threshold,
                    sample_list,
                    num_of_hets,
                    lods_cut_off,
                    writelod,
                    maxed_as,
                )

                # append the dataframes that are returned as "phase_extended" haplotype block
                # .. thus creating a list of dataframes
                df_list_by_bed.append(phase_extended_by_bed)

        # append the regions of the dataframe that was restricted by bed-file
        df_list_by_bed.append(data_by_bed_void)

        # concat all the dataframes in the list into one big dataframe
        # this merges all the datraframe by bed-file within each chromosome/contig together
        phase_extended_by_chr = pd.concat(df_list_by_bed)

        # clear memory
        data_by_bed_void = None
        df_list_by_bed = None
        phase_extended_by_bed = None

        print(
            '  - Phase-extension completed for contig "%s" in %s seconds'
            % (chr_, time.time() - time_chr)
        )
        print("  - Worker maximum memory usage: %.2f (mb)" % (current_mem_usage()))
        print()

        return phase_extended_by_chr.sort_values(by=["POS"])

    # if haplotype are extended at chromosome/contig level ..
    # ..just pass whole contig/chromosome at once
    else:
        phase_extended_by_chr = process_consecutive_blocks(
            contigs_group,
            soi,
            chr_,
            snp_threshold,
            sample_list,
            num_of_hets,
            lods_cut_off,
            writelod,
            maxed_as,
        )

        # clear memory
        contigs_group = None
        del contigs_group

    print(
        '  - Phase-extension completed for contig "%s" in %s seconds'
        % (chr_, time.time() - time_chr)
    )
    print("  - Worker maximum memory usage: %.2f (mb)" % (current_mem_usage()))
    print()

    # return the phase-extended haplotype back to the pool-process ..
    # .. which will be pd.concatenated after all pools are complete
    return phase_extended_by_chr.sort_values(by=["POS"])


""" now, read the two consecutive blocks to do phase extension. """


def process_consecutive_blocks(
    contigs_group,
    soi,
    chr_,
    snp_threshold,
    sample_list,
    num_of_hets,
    lods_cut_off,
    writelod,
    maxed_as,
):

    # print()
    print('  - Grouping the dataframe using unique "PI - phased index" values. ')

    """ Step 02 - D: group dataframe again by "PI keys" of soi and then
       sort by minimum "POS" value for each "PI key".
       - This sorting is necessary because sometimes "haplotype blocks" are like 3-3-3-3  5-5-5  3-3-3-3
          - i.e there are small RBphased blocks within the boundry of larger RBphased block.
          - Not, sure what is causing this (prolly sampling difference of large vs. small chunks in PE reads)
          - This problem should go away in first round of haplotype-extension"""

    contigs_group = contigs_group.assign(
        New=contigs_group.groupby([soi + ":PI"]).POS.transform("min")
    ).sort_values(["New", "POS"])

    """ Step 03: Now, start reading the "contigs_group" for haplotype-extension.
    A) Store the data as dictionary with 'header' values as keys. Some keys are: CHROM, POS, sample (PI, PG within sample),
       etc ... Then group the dictionary using unique "PI" values as 'keys' for grouping.
        Note: This dict-data should contain information about two adjacent haplotype blocks that needs extending.
        In this example I want to extend the haplotypes for "sample ms02g" which has two blocks 6 and 4.
        So, I read the PI and PG value for this sample. Also, data should store with some unique keys.
    B) Iterate over two consecutive Haplotype-Blocks at once.
        Note: While iterating over two blocks, initially we write the very first block of the "contig". With this
        method, now when we iterate over two consecutive blocks we can only update and write the second block.
        """

    # covert pandas dataframe back to text like file before converting it into dictionary.
    contigs_group = pd.DataFrame.to_csv(
        contigs_group, sep="\t", index=False, header=True
    )

    """ Step 03 - A : read the data with header as keys and groupby using that "keys" """
    phased_dict = csv.DictReader(StringIO(contigs_group), delimiter="\t")
    phased_grouped = itertools.groupby(phased_dict, key=lambda x: x[soi + ":PI"])

    """ Since the dictionary isn't ordered, we return the order using OrderedDictionary """
    # ** for future: there is room for improvement in here (memory and speed)
    grouped_data = collections.OrderedDict()
    for key, grp in phased_grouped:
        grouped_data[key] = accumulate(grp)

    """ Clear memory """
    del phased_dict
    del phased_grouped
    del contigs_group

    # print()
    print("  - Starting MarkovChains for contig %s" % chr_)
    """ Step 03 - B : now pipe the data for phase extension """
    """ Step 03 - B : And, iterate over two consecutive Haplotype-Blocks at once. This is done to obtain all
    possible Haplotype configurations between two blocks. The (keys,values) for first block is represented as
    k1,v2 and for the later block as k2,v2. """

    """ Step 03 - B (i): Before running consecutive blocks, we write data from the very first block to the file.
    Reason : Before we start computing and solving the haplotype phase state, we plan to write the
    data for very first block (k1, v1). So, after that, we can solve the relation between two consecutive..
    .. blocks but only write data from 2nd block each time - based on what relation comes out. """
    very_first_block = [list(grouped_data.items())[0]]

    if len(list(grouped_data.items())) == 1:
        print("there is only one block, so skipping phase extension")

    # write header of the extended phase-block
    extended_haplotype = (
        "\t".join(["CHROM", "POS", "REF", "all-alleles", soi + ":PI", soi + ":PG_al"])
        + "\n"
    )

    if writelod == "yes":  # add extra field if desired by user
        extended_haplotype = extended_haplotype.rstrip("\n") + "\tlog2odds\n"
        log2odds = ""

    # write data/values from very first block.
    for k1, v1 in very_first_block:
        for r1, vals in enumerate(v1[soi + ":PI"]):
            new_line = (
                "\t".join(
                    [
                        v1["CHROM"][r1],
                        v1["POS"][r1],
                        v1["REF"][r1],
                        v1["all-alleles"][r1],
                        v1[soi + ":PI"][r1],
                        v1[soi + ":PG_al"][r1],
                    ]
                )
                + "\n"
            )
            if writelod == "yes":
                new_line = new_line.rstrip("\n") + "\t.\n"

            extended_haplotype += new_line

        # print('very first block end\n\n')  # marker for debugging

    """ Step 03 - B (ii):  Starting MarkovChains.
            Now, read data from two consecutive blocks at a time.
            Note: At the end of computation write the data only from each k2 block. No need to write the data
            from k1 block of each iteration because it was written in earlier loop."""

    """ Step 03 - B (ii - 1) Create empty "checker variables".
        Note: This checker variables (actually multi-level boolean logic) help to carryover information from
        ealier iteration of a for-loop - i.e identify if the values from later block i.e k2, v2 were phased to
        to earlier block (k1, v1) in "parallel" vs. "alternate configuration".
        - If two consecutive blocks are phased, k2_new is now assigned k1 from earlier block; else (if not phased)
          k2_new stays empty ('').
        - So, the role of flipped variable is to keep information if k2,v2 were phased straight vs. alternate
          compared to k1, v1 in the earlier run. These checker-variables are crucial to keep the proper phase-state
          in the output file."""

    # start checker variables
    k2_new = ""  # updates the index of k2 for each k1,v1 ; k2,v2 run
    flipped = ""  # boolean logic to check and store if the phase state flipped during extension

    """ Step 03 - B (ii - 2): Now, read two consecutive blocks at a time"""
    for (k1, v1), (k2, v2) in zip(
        grouped_data.items(), itertools.islice(grouped_data.items(), 1, None)
    ):

        """Step 03 - B (ii - 2-A): iterate over the first Haplotype Block, i.e the k1 block.
        The nucleotides in the left of the phased SNPs are called Block01-haplotype-A,
        and similarly on the right as Block01-haplotype-B."""

        # iterate over the first Haplotype Block, i.e the k1 block and v1 values
        hap_block1a = [
            x.split("|")[0] for x in v1[soi + ":PG_al"]
        ]  # the left haplotype of block01
        hap_block1b = [x.split("|")[1] for x in v1[soi + ":PG_al"]]

        # iterate over the second Haplotype Block, i.e the k2 block and v2 values
        hap_block2a = [x.split("|")[0] for x in v2[soi + ":PG_al"]]
        hap_block2b = [x.split("|")[1] for x in v2[soi + ":PG_al"]]

        """ Step 03 - B (ii - 2-B) : Create possible haplotype configurations for "forward markov chain".
        Possible haplotype Configurations will be, Either :

        1) Block01-haplotype-A phased with Block02-haplotype-A,
            creating -> hapb1a-hapb2a, hapb1b-hapb2b """
        """ First possible configuration """
        hapb1a_hapb2a = [hap_block1a, hap_block2a]
        hapb1b_hapb2b = [hap_block1b, hap_block2b]

        """ Or, Second Possible Configuration
        2) block01-haplotype-A phased with Block02-haplotype-B
            creating -> hapb1a-hapb2b, hapb1b-hapb2a """
        hapb1a_hapb2b = [hap_block1a, hap_block2b]
        hapb1b_hapb2a = [hap_block1b, hap_block2a]

        """ Step 03 - B (ii - 2-C) :
        Create possible haplotype configurations for "reverse markov chain"
        - reverse markov chain are added to increase the confidence in likelyhood estimation. """

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
        """ skip if one of the keys has no values - this is redundant ?? - keep it for just in case situation
        ** can also be used in the future if we want to phase the SNPs that have no assigned 'PI' values,
        i.e the "PI" will be "." """
        if k1 == "." or k2 == ".":
            for xi in range(len(v2[soi + ":PI"])):
                new_line = (
                    "\t".join(
                        [
                            v2["CHROM"][xi],
                            v2["POS"][xi],
                            v2["REF"][xi],
                            v2["all-alleles"][xi],
                            k2,
                            hapb1a_hapb2a[1][xi] + "|" + hapb1b_hapb2b[1][xi],
                        ]
                    )
                    + "\n"
                )
                if writelod == "yes":
                    new_line = new_line.rstrip("\n") + "\t.\n"

                extended_haplotype += new_line

            # update the values of checker variables
            k2_new = ""
            flipped = ""

            continue  # to next consecutive blocks
        ######################################################

        """ Step 03 - C : Set the threshold for the minimum number of SNPs required in haplotype block
        before continuing phase extension. """
        """ If all the data in soi, in either v1 or v2 are SNPs below a certain threshold we just write
        the data and continue. i.e say if a Haplotype-Block is composed of only 2 SNPs it will be less
        reliable to extend the phase-state.
        - So, this step can also be used to control the minimum number/size of the haplotypes that is required
        before it can be phase-extended.
        - by default the minimum number of SNPs (exclusive) in the soi haplotype is set to 3.
        - If minimum requirement isn't met just skip extending the phase and write it to a file and continue. """
        number_of_snp_in_soi_v1 = len([x for x in v1[soi + ":PG_al"] if len(x) == 3])
        number_of_snp_in_soi_v2 = len([x for x in v2[soi + ":PG_al"] if len(x) == 3])

        # print('number of SNPs: ', NumSNPsInsoi_v1, NumSNPsInsoi_v2)
        if (
            number_of_snp_in_soi_v1 < snp_threshold
            or number_of_snp_in_soi_v2 < snp_threshold
        ):
            for xth, vals in enumerate(v2[soi + ":PI"]):
                new_line = (
                    "\t".join(
                        [
                            v2["CHROM"][xth],
                            v2["POS"][xth],
                            v2["REF"][xth],
                            v2["all-alleles"][xth],
                            k2,
                            hapb1a_hapb2a[1][xth] + "|" + hapb1b_hapb2b[1][xth],
                        ]
                    )
                    + "\n"
                )
                if writelod == "yes":
                    new_line = new_line.rstrip("\n") + "\t.\n"

                extended_haplotype += new_line

            # update values of the checker variables
            # this is important, so previous k2 and flip state doesn't get carried over without purpose
            k2_new = ""
            flipped = ""

            continue  # to next consecutive blocks

        """ Step 04: For the consecutive blocks that pass the thresholds (SNP number, have PI != '.', etc.,
        pipe the data (k1, v1 ; k2, v2) to a defined function for computation of forward and reverse
        markov chain transition probabilities for these two consecutive blocks (k1, v1; k2, v2) """

        #### for forward chain   ########
        # ** set "orientation=reversed" to compute transition ..
        # .. from the lower tip of former block with upper tip of later block
        # .. this helps in using the closest genomic position between consecutive blocks thus ..
        # .. downsizing the effects created by recombination.
        lhfc_f, lhsc_f = compute_maxLh_score(
            soi,
            sample_list,
            k1,
            k2,
            v1,
            v2,
            num_of_hets,
            hapb1a_hapb2a,
            hapb1b_hapb2b,
            hapb1a_hapb2b,
            hapb1b_hapb2a,
            maxed_as,
            orientation=reversed,
        )

        #### for reverse chain   ########
        # set "orientation=lambda..." just passes a null value keeping orientation as it is.
        lhfc_r, lhsc_r = compute_maxLh_score(
            soi,
            sample_list,
            k1_r,
            k2_r,
            v1_r,
            v2_r,
            num_of_hets,
            hapb1a_hapb2a_r,
            hapb1b_hapb2b_r,
            hapb1a_hapb2b_r,
            hapb1b_hapb2a_r,
            maxed_as,
            orientation=lambda x: x,
        )

        """ Step 05-06 are inside the function "compute_maxLh_score()". The values
        (lhfc_f, lhsc_f, lhfc_r, lhsc_r) returned from this function is then used in Step 07. """

        """ Step 07 :  previous (Step 06) returns the likelyhoods and/or LODs score for both "parallel"
        and alternate configurations (for both forward and reverse algorithm).
        - We now extend the phase states by comparing LODs score against  cutoff-values."""

        """ Step 07 - A(i): calculate the average of the likelyhoods, odds and then log2 of odds. """
        # average of the likelyhooods for first vs. second configuration
        # (from both forward and reverse algorithm)
        # ** note: "maxed_as" variable doesn't apply here, because maxLH using forward vs. reverse ..
        # .. are just re-estimates. So, we simply take and average on both "maxSum" and "maxPd"
        avg_lhfc = Decimal(lhfc_f + lhfc_r) / 2
        avg_lhsc = Decimal(lhsc_f + lhsc_r) / 2

        # therefore, odds of first_vs_second_configuration is
        odds_fc_vs_sc = avg_lhfc / avg_lhsc

        """ Step 07 - A(ii) : convert the likelyhoods to odds-ratio and then logOf 2 odds"""
        lods2_score_1st_config = Decimal(odds_fc_vs_sc).ln() / (Decimal("2").ln())
        lods2_score_2nd_config = -lods2_score_1st_config

        # print('logOdds')  # marker for debugging
        # print(lods2_score_1st_config)

        """ Step 07 - B : pipe the LOD scores and write the phase state between two consecutive blocks.
                - use "lods cutoff" to decide on phase extension
                - and then store, write it to files.
         ** We can also use accumulation of this stage to run histogram building at later stage.
            - that acculated "extended_haplotype" can be all written at once - this is important while multiprocessing. """
        k2_new, flipped, extended_haplotype = extend_phase_state(
            soi,
            k1,
            k2,
            v1,
            v2,
            k2_new,
            flipped,
            lods2_score_1st_config,
            lods_cut_off,
            extended_haplotype,
            hapb1a_hapb2a,
            hapb1b_hapb2b,
            writelod,
        )

        """ Now, go to Step 08, function "extend_phase_state" """
        # this process udates the data in "extended_haplotype" recursively on the for-loop

    # finally return the extended haplotype as pandas dataframe
    phase_extend = extended_haplotype

    del extended_haplotype

    return pd.read_table(StringIO(phase_extend), sep="\t")


def merge_hap_with_bed(my_bed, good_data):

    print("Extracting bed regions and position of the haplotype file ... ")

    c1 = my_bed.CHROM.values
    s1 = my_bed.start.values
    e1 = my_bed.end.values
    c2 = good_data["CHROM"].values
    pos2 = good_data.POS.values

    # now, find the intersecting positions (between haplotype ("pos" values) and bed regions).
    overlap = (c2[:, None] == c1) & ((pos2[:, None] >= s1) & (pos2[:, None] <= e1))

    i, j = np.where(overlap)

    """intersect the haplotype file with bed file """

    print(
        "Intersecting bed regions with haplotype file, and "
        "creating new haplotype boundries for each contig ... "
    )

    df_bed__hap_interesect = pd.DataFrame(
        np.column_stack([good_data.values[i], my_bed.values[j]]),
        columns=good_data.keys().append(my_bed.keys()),
    )

    # drop any duplicate columns
    df_bed__hap_interesect = df_bed__hap_interesect.T.drop_duplicates().T

    df_bed__hap_interesect["start_end"] = pd.DataFrame(
        df_bed__hap_interesect.apply(lambda x: str(x.start) + "-" + str(x.end), axis=1)
    )

    # only keep "contig, pos and intersection ("start-end")"
    df_bed__hap_interesect = df_bed__hap_interesect[["CHROM", "POS", "start_end"]]
    # if regions of intersection are of interest
    # pd.DataFrame.to_csv(df_bed__hap_interesect, 'df_bed and hap intersect.txt', sep='\t', index=None, header=True)

    # now, merge the df(good data) with the df (intersected bed and hap) to create updated
    # haplotype file.
    # Any part that is intersected is assigned unique-keys of "start-end" values. The non-intersected
    # part are filled with "NA"
    data_frames = [good_data, df_bed__hap_interesect]

    print(
        "Merging the bed file with the input haplotype file, "
        "and marking the bed boundries."
    )

    # fill the non-merging lines with "void"
    good_data_update = reduce(
        lambda left, right: pd.merge(left, right, on=["CHROM", "POS"], how="outer"),
        data_frames,
    ).fillna("na")
    # ** for future: assigning an unique identifier like "NA01, NA02" for non-intersecting lines as blocks.

    # if updated "haplotype file" with intersected "haplotype file" and "bed file" is of interest
    # pd.DataFrame.to_csv(good_data_update, 'df_bed and hap merged.txt', sep='\t', index=None,
    # header=True)

    print("Grouping the haplotype file by chromosome (contigs) .... ")
    good_data_by_group = good_data_update.groupby("CHROM", sort=False)

    # clear memory
    del (
        c1,
        c2,
        s1,
        e1,
        pos2,
        i,
        j,
        good_data,
        good_data_update,
        df_bed__hap_interesect,
        my_bed,
    )

    return good_data_by_group


""" prepare list of all available samples (index, value) in the input file
    ** note: we also control output of specific sample if user desires. """


def find_samples(samples):
    """the below method preserves the order of samples.

    ** Note: this returns data as:
    sample_list = [('ms01e:PI', 'ms01e:PG_al'), ('ms02g:PI', 'ms02g:PG_al'),
               ('ms03g:PI', 'ms03g:PG_al'), ('ms04h:PI', 'ms04h:PG_al')]
    """
    seen = set()
    sample_list = [x.split(":")[0] for x in samples if ":" in x]
    sample_list = [x for x in sample_list if not (x in seen or seen.add(x))]
    sample_list = [((x + ":PI"), (x + ":PG_al")) for x in sample_list]

    return sample_list
