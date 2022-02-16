from itertools import product
from decimal import Decimal


""" function to compute transition probs between two consecutive blocks. """


def compute_maxLh_score(
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
    orientation,
):

    # print()
    '''Step 05 : Start preping the first order markov transition matrix and probabilities for the
    consecutive haplotypes between two blocks.
    - To store the sum-values of the product of the transition probabilities.
    These sum are added as the product-of-transition is returned by nested for-loop;
    from the loop "for m in range(....)"'''

    # empty variable for storing "sum of the product of transition probabilities" between two blocks.
    # we can either "max sum" or "max product" the likelyhoods, default is "max-product".
    # cul -> cumulative likelyhood estimates (either summed (+) or product(*))
    if maxed_as == "+":
        cul_of_pt_hapb1a_b2a = 0
        cul_of_pt_hapb1b_b2b = 0

        cul_of_pt_hapb1a_b2b = 0
        cul_of_pt_hapb1b_b2a = 0

    elif maxed_as == "*":
        cul_of_pt_hapb1a_b2a = 1
        cul_of_pt_hapb1b_b2b = 1

        cul_of_pt_hapb1a_b2b = 1
        cul_of_pt_hapb1b_b2a = 1

    """Step 05 - A : We start calculating the nucelotide probabilities (aka "emission probabilities"),
    and then transition probabilities from each level of v1[soi] to each level of v2[soi].
    - we are using reversed() to actually make haplotype configurations starting from the lower tip (i.e last n-th)
    element of former block with the upper tip (first m-th) element of the later block.
    - using reverse() also helps us control if we just want trans-prob calculation between first 200 Het_SNPs.
        - the reason being that when haplotype becomes larger and large it introduces computation burden.
        - also as the haplotype becomes larger the chances of LD at the tips of the haplotypes increase,
            thereby changing the likelyhoods in another direction.
        - orientation : this parameter helps us to control which tip (top vs. bottom) of the former haplotype is
            used with which tip of later haplotype when preparing transition matrix.
            # make this explanation figurative in thesis/dissertation. """

    # to control the number of Heterozygote site we using for phase-extenstion
    num_of_het_site = 0

    # n-ranges from block01
    for n in orientation(range(len(v1[soi + ":PI"]))):

        """Skip computation if n'th item is InDel or has "*" allele in the genotype.
        - These indels are only removed from computation part but are added in the final output file.
        - The Indels are phased based on which SNPs they hitchhike with."""
        if len(v1[soi + ":PG_al"][n]) > 3 or "*" in v1[soi + ":PG_al"][n]:
            continue

        """ only use certain number of HetVars to compute likely hood. This saves computation burden.
            - set at default at 50 heterozygous variants. """
        num_of_het_site += 1
        if num_of_het_site >= num_of_hets:
            continue

        # print('num_of_het_site: ', num_of_het_site, n)  # marker for debugging

        """ Step 05 - A (i) : compute "nucleotide probabilities" aka "emission probabilities"
        - create dictionary to store "nucleotide counts".
        - Set the counts of each nucleotides at zero (0).
        - the counts are updated at "n-th" level after counting the existing nucleotides (either A,T,G,C)
          at "n-th" position for all the samples. """

        # setting dictionary for "emission/nucleotide counts"
        nucleotide_count_dict = {"A": 0, "T": 0, "G": 0, "C": 0}

        # dictionary for "emission" probabilities
        nucleotide_prob_dict = nucleotide_count_dict.copy()  # makes a shallow copy

        """ Creating variables to store the product-values of the transition probabilities.
        These are updated for each level of "n" paired with each level of "m". """
        # potp -> product of transition probability

        ##** need to change here ??

        ############# deprecated ??
        # potp_hapb1a_b2a = 1
        # potp_hapb1b_b2b = 1

        # potp_hapb1a_b2b = 1
        # potp_hapb1b_b2a = 1
        #########################

        if maxed_as == "+":
            tp_hapb1a_b2a = 0
            tp_hapb1b_b2b = 0

            tp_hapb1a_b2b = 0
            tp_hapb1b_b2a = 0

        elif maxed_as == "*":
            tp_hapb1a_b2a = 1
            tp_hapb1b_b2b = 1

            tp_hapb1a_b2b = 1
            tp_hapb1b_b2a = 1

        """ Step 05 - A (ii) :
        - Now, calculate the initial nucleotide counts at "n-th" level,
        before computing the transition counts.
        - only calculated from "v1" at n-th level and only once for each parse/iteration """
        for (x, y) in sample_list:
            for nucleotide in "ATGC":
                nucleotide_count_dict[nucleotide] += (
                    v1[y][n].split("|").count(nucleotide)
                )

        # now, compute emission probabilities
        total_nucl = sum(list(nucleotide_count_dict.values()))

        for ky, vy in nucleotide_count_dict.items():
            # print(ky, vy)
            nucleotide_prob_dict[ky] = vy / total_nucl

        """ Step 05 - B : Count number of transition from each nucleotide (n-th) to each nucleotide (m-th).
        - Now, we read nucleotides (ATGC) at each level of "m"
        to compute the transition from each level of "n". """

        # to control certain number of hetVars to include in computation
        num_of_het_site_at_m = 0
        for m in range(len(v2[soi + ":PI"])):  # m-ranges from block02

            """Like at "n-th" level, skip if InDel present at this "m-th" level.
            But InDel will be phased along with hitchhiking SNPs."""
            if len(v2[soi + ":PG_al"][m]) > 3 or "*" in v2[soi + ":PG_al"][m]:
                continue

            """ only use certain number of HetVars to compute likely hood at "later block".
                - This saves computation burden.
                - set at default at 200 heterozygous variants. """
            num_of_het_site_at_m += 1
            if num_of_het_site_at_m >= num_of_hets:
                continue

            """ Creating an empty dictionary to store transition counts from "n-th" level of V1 to "m-th" level of V2.
            - ('A', 'T') represents transition from 'A' to 'T'.
            - ** for future: probably upgrade using numpy """
            transition_count_dict = {
                ("A", "A"): 0,
                ("A", "T"): 0,
                ("A", "G"): 0,
                ("A", "C"): 0,
                ("T", "A"): 0,
                ("T", "T"): 0,
                ("T", "G"): 0,
                ("T", "C"): 0,
                ("G", "A"): 0,
                ("G", "T"): 0,
                ("G", "G"): 0,
                ("G", "C"): 0,
                ("C", "A"): 0,
                ("C", "T"): 0,
                ("C", "G"): 0,
                ("C", "C"): 0,
            }

            """ Dictionary to store above counts as probabilities. """
            # shallow copy should be fine  # **for future - do deep copy if need be
            transition_prob_dict = transition_count_dict.copy()

            """ Step 05 - B(i) : Now, loop through each sample to compute transition counts (from_ , to)
            for each nucleotide combination"""
            for x, y in sample_list:

                """Skip data mining if values are missing for the sample.
                Continues to next sample, since there is no point in computing transition.
                - We also skip any variant that has "*" allele (i.e deletion spanning variants) in it."""
                if (
                    v1[y][n] == "."
                    or v2[y][m] == "."
                    or "*" in v1[y][n]
                    or "*" in v2[y][m]
                ):
                    continue

                #############################################
                ## ** for future: add min number of samples called in a line, before proceeding
                #############################################

                """ Mine nucleotides at n'th and m'th position in v1 and v2 for each samples.
                - soi (sample of interest) is also included in the counting, rather than adding pseudo counts."""
                nucl_B1 = (v1[y][n]).split("|")  # nucleotides at n'th of Block01
                nucl_B2 = (v2[y][m]).split("|")  # nucleotides at m'th of Block02

                """ Step 05 - B(ii) : Create zipped vs. product configuration of the
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
                """

                """ Create zipped configuration"""
                # if PI values are same/equal prepare zip between phased nucleotides
                if v1[x][n] == v2[x][m]:
                    ziped_nuclb1b2 = list(zip(nucl_B1, nucl_B2))

                    """ then count the number of transitions """
                    for (from_, to) in transition_count_dict:
                        transition_count_dict[(from_, to)] += ziped_nuclb1b2.count(
                            (from_, to)
                        )

                # if PI values are not same then we create product between phased nucleotides
                elif v1[x][n] != v2[x][m]:
                    prod_nuclb1b2 = list(product(nucl_B1, nucl_B2))

                    """ then count the number of transitions """
                    for (from_, to) in transition_count_dict:
                        transition_count_dict[(from_, to)] += (
                            prod_nuclb1b2.count((from_, to)) / 2
                        )

            """Step 06 : Compute likelyhood of "parallel" vs. "alternate" configuration of two consecutive haplotype blocks.
            A) calculate transition probabilities from "n-th" to "m-th" level.
            B) then, calculate the product of the transition from one level of "n" to several levels of "m".
            C) then, calculate joint probability (summed or product) for transitions from several levels of "n"
               to several levels of "m".
            D) then, calculate overall likelyhood and log2Odds of "parallel" vs. "alternate" configuration
               using both forward and backward algorith. """

            """ Step 06 - A : using nucleotide counts and transition counts (at "n-th" and "m-th" of each sample),
                              compute transition probabilities.
            ** Note (for future): At the end this all possible transition probabilities should sum to 1.
            But, at some places transitions may be missing due to situation like ('A', '.')
            i.e from 'A' to '.' empty.
            ** future - This problem may be adjusted in the future"""

            """ Step 06 - A (i) : compute all observed transition probabilities """
            for (from_, to) in transition_prob_dict:
                transition_prob_dict[(from_, to)] = compute_transition_probs(
                    transition_count_dict[(from_, to)], nucleotide_count_dict[from_]
                )

            """ Step 06 - A (ii) : find observed configuration for soi at "n-th" and "m-th" level
            i.e (from_, to). """
            hapb1a_hapb2a_transition = (hapb1a_hapb2a[0][n], hapb1a_hapb2a[1][m])
            hapb1b_hapb2b_transition = (hapb1b_hapb2b[0][n], hapb1b_hapb2b[1][m])

            hapb1a_hapb2b_transition = (hapb1a_hapb2b[0][n], hapb1a_hapb2b[1][m])
            hapb1b_hapb2a_transition = (hapb1b_hapb2a[0][n], hapb1b_hapb2a[1][m])

            """ Step 06 - B : computes the "product of transition probabilities" from "n-th" to
             several levels of "m" using for-loop.
             - ** because the below code is indented one level below "for m in range(len(v2['ms02g:PI']))".
             - ** no need to add pseudo-counts, because if no haplotypes are observed in any samples except soi,
               the prob(from_, to) for each configuration will be 1/4 there by nullifying the likelyhoods to "1".
            """

            ######### deprecated ??
            # potp_hapb1a_b2a *= transition_prob_dict[hapb1a_hapb2a_transition]
            # potp_hapb1b_b2b *= transition_prob_dict[hapb1b_hapb2b_transition]
            # potp_hapb1a_b2b *= transition_prob_dict[hapb1a_hapb2b_transition]
            # potp_hapb1b_b2a *= transition_prob_dict[hapb1b_hapb2a_transition]

            # new additions ****
            # adding emission probabilities to the method to get more fair estimate of the likelihoods
            # p(X given Y) = p(X) * p(XtY - transitions)
            # ** for future: this if/else can be rather fixed by establishing a new function
            ## we are multiplying p(transition)*p(emission)  and doing maxSum or maxProd depending upon what is required
            ## the codes below can be dramatically improved.

            if maxed_as == "+":
                tp_hapb1a_b2a += transition_prob_dict[
                    hapb1a_hapb2a_transition
                ] * Decimal(nucleotide_prob_dict[hapb1a_hapb2a_transition[0]])
                tp_hapb1b_b2b += transition_prob_dict[
                    hapb1b_hapb2b_transition
                ] * Decimal(nucleotide_prob_dict[hapb1b_hapb2b_transition[0]])
                tp_hapb1a_b2b += transition_prob_dict[
                    hapb1a_hapb2b_transition
                ] * Decimal(nucleotide_prob_dict[hapb1a_hapb2b_transition[0]])
                tp_hapb1b_b2a += transition_prob_dict[
                    hapb1b_hapb2a_transition
                ] * Decimal(nucleotide_prob_dict[hapb1b_hapb2a_transition[0]])

            elif maxed_as == "*":
                tp_hapb1a_b2a *= transition_prob_dict[
                    hapb1a_hapb2a_transition
                ] * Decimal(nucleotide_prob_dict[hapb1a_hapb2a_transition[0]])
                tp_hapb1b_b2b *= transition_prob_dict[
                    hapb1b_hapb2b_transition
                ] * Decimal(nucleotide_prob_dict[hapb1b_hapb2b_transition[0]])
                tp_hapb1a_b2b *= transition_prob_dict[
                    hapb1a_hapb2b_transition
                ] * Decimal(nucleotide_prob_dict[hapb1a_hapb2b_transition[0]])
                tp_hapb1b_b2a *= transition_prob_dict[
                    hapb1b_hapb2a_transition
                ] * Decimal(nucleotide_prob_dict[hapb1b_hapb2a_transition[0]])
            ##

        """ Step 06 - C : compute the max sum or max product of the transition probabilities
        across several levels of "n" to several levels of "m".
        ** Note: we can either do "max sum" or "max product".
        So, "cul_of_pt_hapb1a_b2a" is the sum of the likelyhoods of hapBlock1A being phased with hapBlock2A"""

        if maxed_as == "+":
            cul_of_pt_hapb1a_b2a += tp_hapb1a_b2a
            cul_of_pt_hapb1b_b2b += tp_hapb1b_b2b

            cul_of_pt_hapb1a_b2b += tp_hapb1a_b2b
            cul_of_pt_hapb1b_b2a += tp_hapb1b_b2a

        elif maxed_as == "*":
            cul_of_pt_hapb1a_b2a *= tp_hapb1a_b2a
            cul_of_pt_hapb1b_b2b *= tp_hapb1b_b2b

            cul_of_pt_hapb1a_b2b *= tp_hapb1a_b2b
            cul_of_pt_hapb1b_b2a *= tp_hapb1b_b2a

    """ Step 06 - D : Now, compute the "Odds ratio" and "log2 of the Odds"
        of each possible haplotype configuration, for both 1st vs. 2nd configurations """

    """ Step 06 - D(i) : compute the likely hood of first configuration (lhfc) vs. second configuration (lhsc)
    First Configuration = hapb1a with hapb2a, and hapb1b with hapb2b.
    Second Configuration = hapb1a with hapb2b, and hapb1b with hapb2a. """
    if maxed_as == "+":
        lhfc = Decimal(cul_of_pt_hapb1a_b2a + cul_of_pt_hapb1b_b2b)
        lhsc = Decimal(cul_of_pt_hapb1a_b2b + cul_of_pt_hapb1b_b2a)

    elif maxed_as == "*":
        lhfc = Decimal(cul_of_pt_hapb1a_b2a * cul_of_pt_hapb1b_b2b)
        lhsc = Decimal(cul_of_pt_hapb1a_b2b * cul_of_pt_hapb1b_b2a)

    # lhfc = Decimal(cul_of_pt_hapb1a_b2a * cul_of_pt_hapb1b_b2b)
    # lhsc = Decimal(cul_of_pt_hapb1a_b2b * cul_of_pt_hapb1b_b2a)

    return lhfc, lhsc

    ##################################################
    """ ** Note - for the future: The likelyhood can also be tested with each configuration split into two parts.
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
    """
    #############################################################


def extend_phase_state(
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
):

    # print()
    hmmm = 1  # ** - just added to highlight the docstring in pyCharm; delete in future

    """ Step 08 : extend the phase states between two consecutive haplotypes based on the evidence.
     - default cut-off is set at log2 score of 4 (i.e 2^4 times likely), which can be changed by user
     - this function also returns the updated values of checker-variable i.e "k2_new" and "flipped",
       these checker-variables are very important to extend phase in a proper state.
       - "k2_new" keeps the log of "if the haplotype was extended"
       - "flipped" keeps the log of "if it was "parallel" vs. "alternate" configuration".
     - we write the data to a file, and also store in another variable to compute "haplotype stats" if need be.
    """
    # breakpoint()

    # when "k2_new" is not updated
    # this happens for the very first consecutive block of a contig
    # this happens when the previous consecutive blocks failed to extend
    if k2_new == "":

        # if the |lods score| is above cut-off the phase needs extension in parallel configuration

        if lods2_score_1st_config > lods_cut_off:  # i.e 2^4 times more likely
            k2_new = k1  # now, "k1" value gets carried over for next consecutive run
            flipped = "no"  # since it was phased in parallel configuration, there for "no flip"

            for xth in range(len(v2[soi + ":PI"])):
                new_line = (
                    "\t".join(
                        [
                            v2["CHROM"][xth],
                            v2["POS"][xth],
                            v2["REF"][xth],
                            v2["all-alleles"][xth],
                            k2_new,
                            hapb1a_hapb2a[1][xth] + "|" + hapb1b_hapb2b[1][xth],
                        ]
                    )
                    + "\n"
                )
                if writelod == "yes":
                    new_line = (
                        new_line.rstrip("\n")
                        + "\t"
                        + str(round(lods2_score_1st_config, 5))
                        + "\n"
                    )

                extended_haplotype += new_line

        # if the |-lods score| is above cut-off, phase needs extension in "alternate" configuration
        elif lods2_score_1st_config < -lods_cut_off:  # i.e 2^4 times less likely
            k2_new = k1
            flipped = "yes"

            for xth in range(len(v2[soi + ":PI"])):
                new_line = (
                    "\t".join(
                        [
                            v2["CHROM"][xth],
                            v2["POS"][xth],
                            v2["REF"][xth],
                            v2["all-alleles"][xth],
                            k2_new,
                            hapb1b_hapb2b[1][xth] + "|" + hapb1a_hapb2a[1][xth],
                        ]
                    )
                    + "\n"
                )
                if writelod == "yes":
                    new_line = (
                        new_line.rstrip("\n")
                        + "\t"
                        + str(round(lods2_score_1st_config, 5))
                        + "\n"
                    )

                extended_haplotype += new_line
        # if the |lods score| does not pass cut-off threshold in either direction, there is no phase extension
        # so, reset the "checker-variables" and write/store data in original configuration as input data
        else:  # i.e 4 times less likely
            k2_new = ""
            flipped = ""

            # print('no phase extension, no flip')  # marker for debugging
            for xth in range(len(v2[soi + ":PI"])):
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
                    new_line = (
                        new_line.rstrip("\n")
                        + "\t"
                        + str(round(lods2_score_1st_config, 5))
                        + "\n"
                    )

                extended_haplotype += new_line

    # when previous consecutive blocks were extended,
    # "k2_new" now carries a value from previous (k2) which is actually (k1) in this run
    # "flipped" is also carrying an updated value
    elif k2_new != "":
        if flipped == "no":

            # previous config was extended but not flipped, so it must stay in "parallel" config this time.
            if lods2_score_1st_config > lods_cut_off:  # i.e 4 times more likely
                k2_new = k2_new
                flipped = "no"

                # print (1_A is phased with 2_A, no flip, k2:%s \n' %k2_new)  # marker for debugging
                for xth in range(len(v2[soi + ":PI"])):
                    new_line = (
                        "\t".join(
                            [
                                v2["CHROM"][xth],
                                v2["POS"][xth],
                                v2["REF"][xth],
                                v2["all-alleles"][xth],
                                k2_new,
                                hapb1a_hapb2a[1][xth] + "|" + hapb1b_hapb2b[1][xth],
                            ]
                        )
                        + "\n"
                    )
                    if writelod == "yes":
                        new_line = (
                            new_line.rstrip("\n")
                            + "\t"
                            + str(round(lods2_score_1st_config, 5))
                            + "\n"
                        )

                    extended_haplotype += new_line

            # previous config was extended but not flipped. So, now when the phase extends
            # the configuration is kept "alternate" because lods2 score is negative.
            elif lods2_score_1st_config < -lods_cut_off:  # i.e 4 times less likely
                k2_new = k2_new
                flipped = "yes"

                # print ('1_A is phased with 2_B, yes flipped, k2:%s \n' % k2_new) # marker for debugging
                for xth in range(len(v2[soi + ":PI"])):
                    new_line = (
                        "\t".join(
                            [
                                v2["CHROM"][xth],
                                v2["POS"][xth],
                                v2["REF"][xth],
                                v2["all-alleles"][xth],
                                k2_new,
                                hapb1b_hapb2b[1][xth] + "|" + hapb1a_hapb2a[1][xth],
                            ]
                        )
                        + "\n"
                    )
                    if writelod == "yes":
                        new_line = (
                            new_line.rstrip("\n")
                            + "\t"
                            + str(round(lods2_score_1st_config, 5))
                            + "\n"
                        )

                    extended_haplotype += new_line

            # if phase state isn't extending, write the data as is and reset the "checker-variables"
            else:
                k2_new = ""
                flipped = ""

                # print ('no phase extension, no flip, k2 %s\n' % k2)  # marker for debugging
                for xth in range(len(v2[soi + ":PI"])):
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
                        new_line = (
                            new_line.rstrip("\n")
                            + "\t"
                            + str(round(lods2_score_1st_config, 5))
                            + "\n"
                        )

                    extended_haplotype += new_line

        # if previous block were extended but flipped
        elif flipped == "yes":

            # previous config was extended but flipped, so to maintain it's phase state
            # it must stay in "parallel" config this time by flipping with previous extension.
            if lods2_score_1st_config > lods_cut_off:  # i.e 4 times more likely
                k2_new = k2_new
                flipped = "yes"

                # print('\n1_A is phased with 2_A, but yes flipped, k2:%s \n' % k2_new)  # marker for debug
                for xth in range(len(v2[soi + ":PI"])):
                    new_line = (
                        "\t".join(
                            [
                                v2["CHROM"][xth],
                                v2["POS"][xth],
                                v2["REF"][xth],
                                v2["all-alleles"][xth],
                                k2_new,
                                hapb1b_hapb2b[1][xth] + "|" + hapb1a_hapb2a[1][xth],
                            ]
                        )
                        + "\n"
                    )
                    if writelod == "yes":
                        new_line = (
                            new_line.rstrip("\n")
                            + "\t"
                            + str(round(lods2_score_1st_config, 5))
                            + "\n"
                        )

                    extended_haplotype += new_line

            # previous config was extended but flipped, so to maintain it's phase state
            # it must stay in "alternate" config this time by not flipping with previous extension.
            elif lods2_score_1st_config < -lods_cut_off:  # i.e 4 times less likely
                k2_new = k2_new
                flipped = "no"

                # print('\n1_A is phased with 2_B, but not flipped, k2:%s \n' % k2_new)  # debug
                for xth in range(len(v2[soi + ":PI"])):
                    new_line = (
                        "\t".join(
                            [
                                v2["CHROM"][xth],
                                v2["POS"][xth],
                                v2["REF"][xth],
                                v2["all-alleles"][xth],
                                k2_new,
                                hapb1a_hapb2a[1][xth] + "|" + hapb1b_hapb2b[1][xth],
                            ]
                        )
                        + "\n"
                    )
                    if writelod == "yes":
                        new_line = (
                            new_line.rstrip("\n")
                            + "\t"
                            + str(round(lods2_score_1st_config, 5))
                            + "\n"
                        )

                    extended_haplotype += new_line

            # previous extension was phased, but this one didn't
            # so, write the data as it is and reset the "checker-variables"
            else:
                k2_new = ""
                flipped = ""

                # print('\nno phase extension, no flip (reset), k2 %s\n' % k2)  # marker for debugging
                for xth in range(len(v2[soi + ":PI"])):
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
                        new_line = (
                            new_line.rstrip("\n")
                            + "\t"
                            + str(round(lods2_score_1st_config, 5))
                            + "\n"
                        )

                    extended_haplotype += new_line

    return k2_new, flipped, extended_haplotype


""" function to return the transitions probabilities from transition counts"""


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
