import argparse


def get_args(parser):
    """Define required argument for interactive mode program."""

    parser.add_argument(
        "--nt",
        help="number of process to run -> "
        "the maximum number of processes that can be run at once is the number "
        "of different chromosomes (contigs) in the input haplotype file.",
        default=1,
        required=False,
    )

    parser.add_argument(
        "--input",
        help="name of the input haplotype file -> "
        "This haplotype file should contain unique index represented by 'PI' and "
        "phased genotype represented by 'PG_al' for several samples. ",
        required=True,
    )

    parser.add_argument(
        "--SOI",
        help="sample of interest -> "
        "Name of the sample intended for haplotype extension.",
        required=True,
    )

    parser.add_argument(
        "--output",
        help="Name of the output directory. default: 'SOI_extended' ",
        default="",
        required=False,
    )

    parser.add_argument(
        "--refHap",
        help="reference haplotype panel -> "
        "This file should also contain 'PI' and 'PG_al' values for each "
        "sample in that haplotype reference panel."
        "default: empty ",
        default=None,
        required=False,
    )

    parser.add_argument(
        "--useSample",
        help="list of samples -> "
        "Use phase state information only from these samples while running. "
        "This is intended to provide control over phase-extension by "
        "allowing select samples from the pool of samples (refHap and/or input). "
        "This is useful for comparing the results when different individuals, "
        "populations are used in phase extension process."
        "Options: 'all','refHap','input','comma separated name of samples'. "
        "default: 'all' - i.e all the samples in (refHap + input) will be used, "
        "if refHap is given else all samples only from input file is used.",
        default="all",
        required=False,
    )

    parser.add_argument(
        "--bed",
        help=" bed file -> "
        "Process the haplotype extension only in this bed regions. "
        "This is useful to limit haplotype extension only within certain "
        "regions, like - within genes, exons, introns, QTL, etc. "
        "default: empty ",
        default=None,
        required=False,
    )

    parser.add_argument(
        "--snpTh",
        help="snp threshold -> "
        "Minimum number of SNPs required in both haplotype blocks before starting "
        "phase extension. "
        "Option: an integer value. "
        "Default: snpTh = 3 ",
        default=3,
        required=False,
    )

    parser.add_argument(
        "--numHets",
        help="number of heterozygote SNPs -> "
        "Maximum number of heterozygote SNPs used in consecutive haplotype blocks "
        "for computing likelyhood of the phase states. "
        "Option: an integer value. "
        "Default: numHets = 40 ",
        default=40,
        required=False,
    )

    parser.add_argument(
        "--lods",
        help="log2 of Odds cut off -> "
        "Cutoff threshold to assign the phase states between consecutive haplotype blocks. "
        "Option: an integer value "
        "Default: 5 for parallel configuration (i.e 2^5 = 32 times more likely). ",
        default=5,
        required=False,
    )
    # **Note:  positive score above the cutoff values joins two consecutive block in "
    # "parallel configuration, negative score below the cutoff joins two consecutive block in "
    # "alternate configuration.

    parser.add_argument(
        "--culLH",
        help="Cumulative likelhood estimates -> "
        "The likelhoods for two possible configuration can either be max-sum vs. max-product. "
        "Default: maxPd i.e max-product. "
        "Options: 'maxPd' or 'maxSum'. ",
        default="maxPd",
        required=False,
    )

    parser.add_argument(
        "--writeLOD",
        help="write log2 of Odds to the output file-> "
        "writes the calculated LODs between two consecutive haplotype blocks in the output file. "
        "Option: 'yes', 'no'. "
        "Default: no. "
        "**Note: the 'lods-score' are printed regardless if the "
        "consecutive blocks are joined or not.",
        default="no",
        required=False,
    )  # ** to do - add this feature

    parser.add_argument(
        "--hapStats",
        help="Computes the descriptive statistics and plots histogram of the haplotype for "
        "input and extended haplotype. "
        "Default: 'no'."
        "Option: 'yes', 'no' .",
        default="no",
        required=False,
    )

    parser.add_argument(
        "--addMissingSites",
        help=" write the lines that have data missing for SOI on the output file. "
        "Option: yes, no ",
        default="no",
        required=False,
    )  # ** to do - add this

    args = parser.parse_args()
    return args


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


# parser.add_argument("--nice", type=str2bool, nargs='?',
#                         const=True, default=False,
#                         help="Activate nice mode.")
