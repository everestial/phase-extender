import os
import shutil


def args_to_val(args):
    soi = args.SOI
    print('  - sample of interest: "%s" ' % (soi))

    # Assign the output directory
    outputdir = args.output if args.output != "" else soi + "_extended"
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
    snp_threshold = args.snpTh  # default, snp_threshold = 3
    print(
        '  - each consecutive haplotype block should have minimum of "%s" SNPs '
        % (snp_threshold)
    )

    # controls the max number of hetVars to include in likelyhood calculation
    num_of_hets = int(args.numHets)  # default, num_of_hets = 40
    print(
        '  - using maximum of "%s" heterozygote sites in each consecutive blocks to compute '
        "transition probabilities" % (num_of_hets)
    )

    # add argument for max sum vs. max product of likelyhood estimates before calculating the LOD-score
    maxed_as = "+" if args.culLH in ("+", "maxSum") else "*"  # default, maxed_as = "*"
    print(
        '  - using "%s" to estimate the cumulative maximum likelyhood of each haplotype configuration '
        "between two consecutive blocks " % (maxed_as)
    )

    ##set the required variables related to bed file if bed file is given.
    # use provided bed file to limit phase extension on bed regions/boundries
    bed_file = args.bed
    if bed_file:
        print(
            '  - using the bed file "%s" to limit phase extension at the bed boundries '
            % (bed_file)
        )
    else:
        print("  - no bed file is given.")

    # if a haplotype panel is provided then the reference panel can be used as backbone or ..
    # .. meaningful data to aid in phase extension
    refhap = args.refHap
    if refhap:
        print('  - using the reference haplotype panel "%s" ' % (refhap))
    else:
        print("  - no reference haplotype panel is provided ")

    # which samples to use while running phase-extension, default is all (hapRef + input)
    # this variable is updated, later in the pipeline when all sample names are collected
    use_sample = args.useSample  # default, use_sample = "all"

    # print the hapstats to file and also plot histogram
    if args.hapStats == "yes":  # default, hapstats = 'no' ** ??
        hapstats = "yes"
        print(
            "  - statistics of the haplotype before and after extension will "
            'be prepared for the sample of interest i.e "%s" ' % (soi)
        )
    else:
        hapstats = "no"
        print(
            "  - statistics of the haplotype before and after extension will not "
            'be prepared for the sample of interest i.e "%s". '
            "    Only extendent haplotype block will be prepared." % (soi)
        )

    # write calculated LOD (log2 of odds) of the possible phase state between two consecutive blocks
    # default, writeLOD = "no"
    if args.writeLOD == "yes":
        writelod = "yes"
        print(
            "  - LOD (log 2 of odds) for consecutive block will be written to the output file "
        )
    else:
        writelod = "no"
        print(
            "  - LOD (log 2 of odds) for consecutive block will not be written to the output file "
        )



    # if missing sites are to be added to phase extended file.
    # this will output a new file instead of writing on top of phase extended file.
    if args.addMissingSites == "yes":
        addmissingsites = "yes"
    else:
        addmissingsites = "no"

    return (
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
    )
