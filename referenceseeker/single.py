import shutil
import referenceseeker.mash as mash
import referenceseeker.ani as rani
import referenceseeker.util as util


def single(args, config):
    """allows single genome analysis"""
    # mash out best hits
    screened_ref_genome_ids, mash_distances = mash.run_mash(args, config)

    # get genomes from RefSeq by accessions
    ref_genomes = util.read_reference_genomes(config)
    screened_ref_genomes = {k: v for k, v in ref_genomes.items() if k in screened_ref_genome_ids}

    # build dna fragments
    dna_fragments_path = config['tmp'].joinpath('dna-fragments.fasta')
    dna_fragments = util.build_dna_fragments(config['genome_path'], dna_fragments_path)

    results = rani.align(args, config, screened_ref_genomes, dna_fragments_path, dna_fragments)

    # remove tmp dir
    shutil.rmtree(str(config['tmp']))

    # filter and sort results
    filtered_reference_ids = []
    for ref_genome_id, result in results.items():
        if args.unfiltered:
            filtered_reference_ids.append(ref_genome_id)
        else:
            if args.bidirectional:
                query_ref = result[0]
                ref_query = result[1]
                if((query_ref[0] >= config['ani']) and (query_ref[1] >= config['conserved_dna'])
                        and (ref_query[0] >= config['ani']) and (ref_query[1] >= config['conserved_dna'])):
                    filtered_reference_ids.append(ref_genome_id)
            else:
                (ani, conserved_dna) = result[0]
                if (conserved_dna >= config['conserved_dna']) and (ani >= config['ani']):
                    filtered_reference_ids.append(ref_genome_id)

    # sort and print results according to ANI * conserved DNA values
    if args.bidirectional:
        filtered_reference_ids = sorted(filtered_reference_ids, key=lambda k: (results[k][0][0] * results[k][0][1] * results[k][1][0] * results[k][1][1]), reverse=True)
        if args.verbose:
            print('')
        print('#ID\tMash Distance\tQR ANI\tQR Con. DNA\tRQ ANI\tRQ Con. DNA\tTaxonomy ID\tAssembly Status\tOrganism')
        for id in filtered_reference_ids:  # print results to STDOUT
            ref_genome = ref_genomes[id]
            result = results[id]
            print(
                '%s\t%1.5f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%s\t%s\t%s' %
                (
                    id,
                    mash_distances[id],
                    result[0][0] * 100,
                    result[0][1] * 100,
                    result[1][0] * 100,
                    result[1][1] * 100,
                    ref_genome['tax'],
                    ref_genome['status'],
                    ref_genome['name']
                )
            )
    else:
        filtered_reference_ids = sorted(filtered_reference_ids, key=lambda k: (results[k][0][0] * results[k][0][1]), reverse=True)
        if args.verbose:
            print('')
        print('#ID\tMash Distance\tANI\tCon. DNA\tTaxonomy ID\tAssembly Status\tOrganism')
        for id in filtered_reference_ids:  # print results to STDOUT
            ref_genome = ref_genomes[id]
            result = results[id][0]
            print(
                '%s\t%1.5f\t%2.2f\t%2.2f\t%s\t%s\t%s' %
                (
                    id,
                    mash_distances[id],
                    result[0] * 100,
                    result[1] * 100,
                    ref_genome['tax'],
                    ref_genome['status'],
                    ref_genome['name']
                )
            )