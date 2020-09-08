
import math


def calculate(args, ref_id_values, common_references, cohort_results):
    """Calculates common ANIs and conDNAs based on choosen algorithm."""
    # Product
    if args.algorithm == "product":
        if args.bidirectional:
            for ref_id in common_references:
                for results in cohort_results:
                    ref_id_values[ref_id][0] *= results[ref_id][0][0] * results[ref_id][1][0]  # Calculating bidirectional ANI
                    ref_id_values[ref_id][1] *= results[ref_id][1][0] * results[ref_id][1][1]  # Calculating bidirectional conDNA
        else:
            for ref_id in common_references:
                for results in cohort_results:
                    ref_id_values[ref_id][0] *= results[ref_id][0][0]  # Calculating ANI
                    ref_id_values[ref_id][1] *= results[ref_id][0][1]  # Calculating conDNA
        # ANIconDNA calculation
        for ref_id in common_references:
            anicondna = ref_id_values[ref_id][0] * ref_id_values[ref_id][1]
            ref_id_values[ref_id].append(anicondna)

    # Mean
    if args.algorithm == "mean":
        if args.bidirectional:
            for ref_id in common_references:
                sum_list_ani = []
                sum_list_condna = []
                for results in cohort_results:
                    sum_list_ani.append(results[ref_id][0][0] and results[ref_id][1][0])  # Calculating sum of all bidirectional ANI
                    sum_list_condna.append(results[ref_id][1][0] and results[ref_id][1][1])  # Calculating sum of all bidirectional conDNA
                ref_id_values[ref_id][0] = sum(sum_list_ani)/(len(sum_list_ani)*2)
                ref_id_values[ref_id][1] = sum(sum_list_condna)/(len(sum_list_condna)*2)

        else:
            for ref_id in common_references:
                sum_list_ani = []
                sum_list_condna = []
                for results in cohort_results:
                    sum_list_ani.append(results[ref_id][0])  # Calculating sum ANI
                    sum_list_condna.append(results[ref_id][1][0])  # Calculating sum conDNA
                ref_id_values[ref_id][0] = sum(sum_list_ani) / len(sum_list_ani)  # Calculating mean ANI
                ref_id_values[ref_id][1] = sum(sum_list_condna) / len(sum_list_condna) #Calculating mean conDNA

    # Root-Mean-Square
    if args.algorithm == "rms":
        if args.bidirectional:
            for ref_id in common_references:
                sum_list_ani = []
                sum_list_condna = []
                for results in cohort_results:
                    sum_list_ani.append((results[ref_id][0][0])^2 and (results[ref_id][1][0])^2)  # Calculating bidirectional root_mean-square ANI
                    sum_list_condna.append((results[ref_id][1][0])^2 and (results[ref_id][1][1])^2)  # Calculating bidirectional root-mean-square conDNA
                ref_id_values[ref_id][0] = math.sqrt(sum(sum_list_ani)/(len(sum_list_ani)*2))
                ref_id_values[ref_id][1] = math.sqrt(sum(sum_list_condna)/(len(sum_list_condna)*2))
        else:
            for ref_id in common_references:
                sum_list_ani = []
                sum_list_condna = []
                for results in cohort_results:
                    sum_list_ani.append(results[ref_id][0])  # Calculating sum ANI
                    sum_list_condna.append(results[ref_id][1][0])  # Calculating sum conDNA
                ref_id_values[ref_id][0] = math.sqrt(sum(sum_list_ani) / len(sum_list_ani))  # Calculating rms ANI
                ref_id_values[ref_id][1] = math.sqrt(sum(sum_list_condna) / len(sum_list_condna))  # Calculating rms conDNA

    return ref_id_values
