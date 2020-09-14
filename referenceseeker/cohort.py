import shutil
import sys
import ntpath

import referenceseeker.util as util
import referenceseeker.ani as rani
import referenceseeker.mash as mash
import concurrent.futures as cf
import referenceseeker.algorithms as algo


def cohort(args, config):
    """Allows cohort genome analysis."""
    all_paths = []
    for path in args.cohort_genomes:
        try:
            all_paths.append(util.check_path(path))
        except FileNotFoundError:
            sys.exit('ERROR: genome file %s is not readable!' % path)
        except PermissionError:
            sys.exit('ERROR (permission): genome file %s is not accessible' % path)
        except OSError:
            sys.exit('ERROR: genome file %s is empty!' % path)
    if all_paths:
        config['genome_path'] = all_paths
    else:
        sys.exit("ERROR: an unexpected path error has occured. Please check your query paths.")

    # mashing
    if args.verbose:
        print('\nEstimate genome distances...')
    mash_output_path = config['tmp'].joinpath('mash.out')

    if args.verbose:
        print('\nEstimate genome distances...')
    mash.exec_mash(config, mash_output_path)

    # Parse mash results
    mash_results, filtered_ids, mash_distances_list = mash.parse_mash_cohort(config, mash_output_path)

    # get genomes from RefSeq by accessions
    ref_genomes = util.read_reference_genomes(config)  # read database reference genomes
    screened_ref_genomes = {k: v for k, v in ref_genomes.items() if k in filtered_ids}

    # build DNA fragments
    dna_fragments_path = config['tmp'].joinpath('dna-fragments.fasta')
    dna_fragments_list = []
    for path in config["genome_path"]:
        dna_fragments = util.build_dna_fragments(path, dna_fragments_path)
        dna_fragments_list.append(dna_fragments)

    # align query fragments to reference genomes and compute ANI/conserved DNA
    cohort_results = []
    query_genomes = []
    if args.verbose:
        print('\nCompute ANIs...')
    with cf.ThreadPoolExecutor(max_workers=args.threads) as tpe:
        futures = []
        results = {}
        for genome_path, dna_fragments in zip(config['genome_path'], dna_fragments_list):
            for identifier, ref_genome in screened_ref_genomes.items():
                futures.append(tpe.submit(rani.align_query_genome, config, dna_fragments_path, dna_fragments, identifier))
            for f in futures:
                ref_genome_id, ani, conserved_dna = f.result()
                results[ref_genome_id] = [(ani, conserved_dna)]

            # align reference genomes fragments to query genome and compute ANI/conserved DNA
            if args.bidirectional:
                if args.verbose:
                    print('\nCompute reverse ANIs...')
                futures = []
                for identifier, ref_genome in screened_ref_genomes.items():
                    futures.append(tpe.submit(rani.align_reference_genome, config, genome_path, identifier))
                for f in futures:
                    ref_genome_id, ani, conserved_dna = f.result()
                    results[ref_genome_id].append((ani, conserved_dna))
            query_genomes.append(ntpath.basename(genome_path).split(".", 1)[0])
            cohort_results.append(results)

    # remove tmp dir
    shutil.rmtree(str(config['tmp']))

    # filter and sort results
    filtered_reference_ids_list = []
    for results in cohort_results:
        filtered_reference_ids = []
        for ref_genome_id, result in results.items():
            if args.unfiltered:
                filtered_reference_ids.append(ref_genome_id)
            else:
                if args.bidirectional:
                    query_ref = result[0]
                    ref_query = result[1]
                    if ((query_ref[0] >= config['ani']) and (query_ref[1] >= config['conserved_dna'])
                            and (ref_query[0] >= config['ani']) and (ref_query[1] >= config['conserved_dna'])):
                        filtered_reference_ids.append(ref_genome_id)
                else:
                    ani, conserved_dna = result[0]
                    if (conserved_dna >= config['conserved_dna']) and (ani >= config['ani']):
                        filtered_reference_ids.append(ref_genome_id)
        filtered_reference_ids_list.append(filtered_reference_ids)

    # Find common Reference genomes
    duplicate_check_list = []
    common_references = []
    for filtered_reference_ids in filtered_reference_ids_list:
        for filtered_ref_id in filtered_reference_ids:
            duplicate_check_list.append(filtered_ref_id)
    for elem in duplicate_check_list:
        if duplicate_check_list.count(elem) == len(filtered_reference_ids_list):
            if elem not in common_references:
                common_references.append(elem)

    # Calculate and print results based on ANI and conDNA
    if args.bidirectional:
        ref_id_values = {r: [1, 1, 1] for r in common_references}
        ref_id_values = algo.calculate(args, ref_id_values, common_references, cohort_results, query_genomes)  # Calculating ANI and conDNA

        common_references = sorted(common_references, key=lambda k: ref_id_values[k][0], reverse=True)

        # printing results
        print('#ID\tMash Distance\tANI\tCon. DNA\tANIconDNA-coefficient\tTaxonomy ID\tAssembly Status\tOrganism')  # "Aniconda?"
        for id in common_references:  # print results to STDOUT
            ref_genome = ref_genomes[id]
            result = ref_id_values[id]

            print(
                '%s\t%1.5f\t%2.2f\t%2.2f\t%2.2f\t%s\t%s\t%s' %
                (
                    id,
                    mash_distances_list[0][id],
                    result[0] * 100,
                    result[1] * 100,
                    result[2] * 100,
                    ref_genome['tax'],
                    ref_genome['status'],
                    ref_genome['name']
                )
                )
    else:
        ref_id_values = {r: [1, 1, 1] for r in common_references}
        ref_id_values = algo.calculate(args, ref_id_values, common_references, cohort_results, query_genomes)  # Calculating ANI and conDNA

        common_references = sorted(common_references, key=lambda k: ref_id_values[k][2], reverse=True)

        # printing results
        print('#ID\tMash Distance\tANI\tCon. DNA\tANIconDNA-coefficient\tTaxonomy ID\tAssembly Status\tOrganism')  # "Aniconda?"
        for id in common_references:  # print results to STDOUT
            ref_genome = ref_genomes[id]
            result = ref_id_values[id]

            print(
                '%s\t%1.5f\t%2.2f\t%2.2f\t%2.2f\t%s\t%s\t%s' %
                (
                    id,
                    mash_distances_list[0][id],
                    result[0] * 100,
                    result[1] * 100,
                    result[0] * result[1],
                    ref_genome['tax'],
                    ref_genome['status'],
                    ref_genome['name']
                )
                )
