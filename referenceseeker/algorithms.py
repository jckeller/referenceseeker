
import math


def calculate(args, ref_id_values, common_references, cohort_results):
    """Calculates common ANIs and conDNAs based on choosen algorithm."""

    raw_scores = {}
    print("ref_id_values: \n", ref_id_values, "\ncommon_references:\n", common_references, "\ncohort_results:\n", cohort_results)
    if args.bidirectional:
        for ref_id in common_references:
            for results in cohort_results:
                ref_id_values[ref_id][0] *= results[ref_id][0][0] * results[ref_id][1][
                    0]  # Calculating bidirectional ANI
                ref_id_values[ref_id][1] *= results[ref_id][1][0] * results[ref_id][1][
                    1]  # Calculating bidirectional conDNA
    else:
        for ref_id in common_references:
            for results in cohort_results:
                ref_id_values[ref_id][0] *= results[ref_id][0][0]  # Calculating ANI
                ref_id_values[ref_id][1] *= results[ref_id][0][1]  # Calculating conDNA

    # ANIconDNA calculation: Product
    if args.algorithm == "product":
        for ref_id in common_references:
            anicondna = ref_id_values[ref_id][0] * ref_id_values[ref_id][1]
            ref_id_values[ref_id].append(anicondna)

    # ANIconDNA calculation: algorithmic mean
    if args.algorithm == "mean":
        for ref_id in common_references:
            anicondna = (ref_id[ref_id][0] + ref_id_values[ref_id][1])/2
            ref_id_values[ref_id].append(anicondna)

    # ANIconDNA calculation: geometric mean
    if args.algorithm == "geometric":
        for ref_id in common_references:
            anicondna = math.sqrt(ref_id_values[ref_id][0] * ref_id_values[ref_id][1])
            ref_id_values[ref_id].append(anicondna)

    # ANIconDNA calculation: harmonic mean
    if args.algorithm == "harmonic":
        for ref_id in common_references:
            anicondna = 2/((1/ref_id[ref_id][0]) + (1/ref_id_values[ref_id][1]))
    return ref_id_values
