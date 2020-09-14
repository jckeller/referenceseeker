import json


def calculate(args, ref_id_values, common_references, cohort_results, query_genomes):
    """Calculates common ANIs and conDNAs based on choosen algorithm."""
    # Export the cohort results
    with open("cohort_results.txt", "w+") as cohort_results_file:
        cohort_results_file.write("#cohort results of the analysis\n")
        cohort_results_file.write(json.dumps(cohort_results, sort_keys=True, indent=4))
    
    # Calculate rawscores of each reference-genome to each query genome
    raw_scores = {}
    if args.bidirectional:
        for ref_id in common_references:
            dic = {}
            for name in query_genomes:
                for results in cohort_results:
                    ani_score = results[ref_id][0][0] * results[ref_id][1][0] # Calculating bidirectional ANI
                    condna_score = results[ref_id][0][1] * results[ref_id][1][1]  # Calculating bidirectional conDNA
                    anicondna_score = ani_score*condna_score
                    dic[name] = (ani_score, condna_score, anicondna_score)
            raw_scores[ref_id] = dic

    else:
        for ref_id in common_references:
            dic = {}
            for name in query_genomes:
                for results in cohort_results:
                    ani_score = results[ref_id][0][0] # Calculating bidirectional ANI
                    condna_score = results[ref_id][0][1]  # Calculating bidirectional conDNA
                    anicondna_score = ani_score * condna_score
                    dic[name] = (ani_score, condna_score, anicondna_score)
            raw_scores[ref_id] = dic

    # Export the rawscores
    with open("rawscores.txt", "w+") as rawscores_file:
        rawscores_file.write("#rawscores of ANI, conDNA and ANIconDNA of each reference-genome to each query-genome\n")
        rawscores_file.write(json.dumps(raw_scores, sort_keys=True, indent=4))

    # ANIconDNA calculation: Product
    if args.algorithm == "product":
        for ref_id_key, ref_id_dic in raw_scores.items():
            ani_score_list = []
            condna_score_list = []
            anicondna_score_list = []
            for key2, query in ref_id_dic.items():
                ani_score_list.append(query[0])
                condna_score_list.append(query[1])
                anicondna_score_list.append(query[2])
            ani = 1
            condna = 1
            anicondna = 1
            for score in ani_score_list:
                ani *= score
            for score in condna_score_list:
                condna *= score
            for score in anicondna_score_list:
                anicondna *= score

            ref_id_values[ref_id_key][0] = ani
            ref_id_values[ref_id_key][1] = condna
            ref_id_values[ref_id_key][2] = anicondna

    # ANIconDNA calculation: algorithmic mean
    if args.algorithm == "mean":
        for ref_id_key, ref_id_dic in raw_scores.items():
            ani_score_list = []
            condna_score_list = []
            anicondna_score_list = []
            for key2, query in ref_id_dic.items():
                ani_score_list.append(query[0])
                condna_score_list.append(query[1])
                anicondna_score_list.append(query[2])
            ani = 0
            condna = 0
            anicondna = 0
            for score in ani_score_list:
                ani += score
            ani = ani / len(ani_score_list)
            for score in condna_score_list:
                condna += score
            condna = condna / len(condna_score_list)
            for score in anicondna_score_list:
                anicondna += score
            anicondna = anicondna / len(anicondna_score_list)

            ref_id_values[ref_id_key][0] = ani
            ref_id_values[ref_id_key][1] = condna
            ref_id_values[ref_id_key][2] = anicondna

    # ANIconDNA calculation: geometric mean
    if args.algorithm == "geometric":
        for ref_id_key, ref_id_dic in raw_scores.items():
            ani_score_list = []
            condna_score_list = []
            anicondna_score_list = []
            for key2, query in ref_id_dic.items():
                ani_score_list.append(query[0])
                condna_score_list.append(query[1])
                anicondna_score_list.append(query[2])
            ani = 1
            condna = 1
            anicondna = 1
            for score in ani_score_list:
                ani *= score
            ani = ani**(1/len(ani_score_list))
            for score in condna_score_list:
                condna *= score
            condna = condna**(1/len(condna_score_list))
            for score in anicondna_score_list:
                anicondna *= score
            anicondna = anicondna**(1/len(anicondna_score_list))

            ref_id_values[ref_id_key][0] = ani
            ref_id_values[ref_id_key][1] = condna
            ref_id_values[ref_id_key][2] = anicondna

    # ANIconDNA calculation: harmonic mean
    if args.algorithm == "harmonic":
        for ref_id_key, ref_id_dic in raw_scores.items():
            ani_score_list = []
            condna_score_list = []
            anicondna_score_list = []
            for key2, query in ref_id_dic.items():
                ani_score_list.append(query[0])
                condna_score_list.append(query[1])
                anicondna_score_list.append(query[2])
            ani = 1
            condna = 1
            anicondna = 1
            for score in ani_score_list:
                ani += (1/score)
            ani = len(ani_score_list)/ani
            for score in condna_score_list:
                condna += (1/score)
            condna = len(condna_score_list)/condna
            for score in anicondna_score_list:
                anicondna += (1/score)
            anicondna = len(anicondna_score_list)/anicondna

            ref_id_values[ref_id_key][0] = ani
            ref_id_values[ref_id_key][1] = condna
            ref_id_values[ref_id_key][2] = anicondna

    # Export final table
    with open("final_scores.txt", "w+") as final_scores_file:
        final_scores_file.write("#final scores of ANI, conDNA and ANIconDNA of each reference-genome\n")
        final_scores_file.write(json.dumps(raw_scores, sort_keys=True, indent=4))

    return ref_id_values
