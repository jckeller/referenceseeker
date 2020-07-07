
import math


def calculate(args, ref_id_values, common_references, cohort_results):
    """Calculates common ANIs and conDNAs based on choosen algorithm."""
    # Product
    if args.algorithm == "product":
        if args.bidirectional:
            for ref_id in common_references:
                for results in cohort_results:
                    ref_id_values[ref_id][0] *= results[ref_id][0][0] * results[ref_id][1][0]  # Calculating forth and back ANI
                    ref_id_values[ref_id][1] *= results[ref_id][1][0] * results[ref_id][1][1]  # Calculating forth and back conDNA
        else:
            for ref_id in common_references:
                for results in cohort_results:
                    ref_id_values[ref_id][0] *= results[ref_id][0]  # Calculating forth and back ANI
                    ref_id_values[ref_id][1] *= results[ref_id][1]  # Calculating forth and back conDNA

    # Mean
    if args.algorithm == "mean":
        if args.bidirectional:
            for ref_id in common_references:
                sum_list_ani = []
                sum_list_condna = []
                for results in cohort_results:
                    sum_list_ani.append(results[ref_id][0][0] and results[ref_id][1][0])  # Calculating sum of all forth and back ANI
                    sum_list_condna.append(results[ref_id][1][0] and results[ref_id][1][1])  # Calculating sum of all forth and back conDNA
                ref_id_values[ref_id][0] = sum(sum_list_ani)/(len(sum_list_ani)*2)
                ref_id_values[ref_id][1] = sum(sum_list_condna)/(len(sum_list_condna)*2)

        else:
            for ref_id in common_references:
                sum_list_ani = []
                sum_list_condna = []
                for results in cohort_results:
                    sum_list_ani.append(results[ref_id][0])  # Calculating sum of all forth and back ANI
                    sum_list_condna.append(results[ref_id][1][0])  # Calculating sum of all forth and back conDNA
                ref_id_values[ref_id][0] = sum(sum_list_ani) / len(sum_list_ani)
                ref_id_values[ref_id][1] = sum(sum_list_condna) / len(sum_list_condna)

    # Root-Mean-Square
    if args.algorithm == "rms":
        if args.bidirectional:
            for ref_id in common_references:
                sum_list_ani = []
                sum_list_condna = []
                for results in cohort_results:
                    sum_list_ani.append((results[ref_id][0][0])^2 and (results[ref_id][1][0])^2)  # Calculating sum of all square forth and back ANI
                    sum_list_condna.append((results[ref_id][1][0])^2 and (results[ref_id][1][1])^2)  # Calculating sum of all square forth and back conDNA
                ref_id_values[ref_id][0] = math.sqrt(sum(sum_list_ani)/(len(sum_list_ani)*2))
                ref_id_values[ref_id][1] = math.sqrt(sum(sum_list_condna)/(len(sum_list_condna)*2))
        else:
            for ref_id in common_references:
                sum_list_ani = []
                sum_list_condna = []
                for results in cohort_results:
                    sum_list_ani.append(results[ref_id][0])  # Calculating sum of all forth and back ANI
                    sum_list_condna.append(results[ref_id][1][0])  # Calculating sum of all forth and back conDNA
                ref_id_values[ref_id][0] = math.sqrt(sum(sum_list_ani) / len(sum_list_ani))
                ref_id_values[ref_id][1] = math.sqrt(sum(sum_list_condna) / len(sum_list_condna))
    return ref_id_values
