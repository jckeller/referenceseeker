
import subprocess as sp
import sys

import referenceseeker.constants as rc


def exec_mash(config, mash_output_path):
    with mash_output_path.open(mode='w') as fh:
        cmd = [
            'mash',
            'dist',
            '-d', rc.UNFILTERED_MASH_DIST if config['unfiltered'] else rc.MAX_MASH_DIST,
            '-p', str(config['threads']),
            str(config['db_path'].joinpath('db.msh')),
        ]
        for genome_path in config["genome_path"]:
            cmd.append(str(genome_path))

        proc = sp.run(
            cmd,
            cwd=str(config['tmp']),
            env=config['env'],
            stdout=fh,
            stderr=sp.PIPE,
            universal_newlines=True
        )
        if proc.returncode != 0:
            sys.exit("ERROR: failed to execute Mash!\nexit=%d\ncmd=%s" % (proc.returncode, cmd))


def parse_mash_results(config, mash_output_path):
    accession_ids = []
    mash_distances = {}
    with mash_output_path.open() as fh:
        for line in fh:
            cols = line.rstrip().split()
            accession_ids.append(cols[0])
            mash_distances[cols[0]] = float(cols[2])
    return accession_ids, mash_distances


def parse_mash_cohort(config, mash_output_path):
    """Parses mash cohort output"""
    single_fasta = []
    mash_results = []
    with mash_output_path.open("r") as fh:
        all_lines = fh.readlines()
        for i, line in enumerate(all_lines):
            line = line.split()
            if i == len(all_lines)-1:  # catch last iteration
                single_fasta.append(line)
                mash_results.append(single_fasta)
            elif not single_fasta:  # catch first iteration
                single_fasta.append(line)
            elif line[1] == single_fasta[0][1]:  # if the line is part of the same fasta
                single_fasta.append(line)
            else:  # if next fasta
                mash_results.append(single_fasta)
                single_fasta = [line]

    # filter mash results for duplicates, save intersection reference hits, delete singles and reduce size of list
    top_mash_results = []
    for mash in mash_results:  # filter for best 100 mash results
        mash = sorted(mash, key=lambda k: k[2], reverse=True)
        mash = mash[0:config['n_mash_results']]
        top_mash_results.append(mash)
    id_list = [ref_genome[0] for ref_genome in top_mash_results[0]]

    # filter ids for intersection mash results of every query genome
    filtered_ids = []
    for id in id_list:
        for mash_result in top_mash_results:
            for (next_ref_genome, _, _, _, _) in mash_result:
                if next_ref_genome == id:
                    break  # id is in mash_result
            else:
                break  # id is not in mash result
        else:
            filtered_ids.append(id)

    filtered_mash_results = []
    for mash_result in top_mash_results:  # sublist containing mash results of query genomes with duplicate ref-genome IDs
        mash = []
        for next_ref_genome in mash_result:
            if next_ref_genome[0] in filtered_ids:
                mash.append(next_ref_genome)
        filtered_mash_results.append(mash)

    # build specific output lists for further processing
    screened_ref_genomes_ids_list = []
    mash_distances_list = []
    for single_fasta in filtered_mash_results:
        mash_distances = {}
        for entry in single_fasta:
            mash_distances[entry[0]] = float(entry[2])  # 0 = ID, 2 = mash distance
        mash_distances_list.append(mash_distances)
    return filtered_mash_results, filtered_ids, mash_distances_list


def run_mash(args, config, mash_output_path):
    """calculates genome distances with mash, extracts the hits and filters for the best hits"""
    # calculate genome distances via Mash
    if args.verbose:
        print('\nEstimate genome distances...')
    exec_mash(config, mash_output_path)

    # extract hits and store dist
    screened_ref_genome_ids, mash_distances = parse_mash_results(config, mash_output_path)
    if args.verbose:
        print("\tscreened %d potential reference genome(s)" % len(screened_ref_genome_ids))

    # reduce Mash output to best hits (args.crg)
    if len(screened_ref_genome_ids) > args.crg:
        if args.verbose:
            print("\treduce to best %d hits..." % args.crg)
        tmp_screened_ref_genome_ids = sorted(screened_ref_genome_ids, key=lambda k: mash_distances[k])
        screened_ref_genome_ids = tmp_screened_ref_genome_ids[:args.crg]
    return screened_ref_genome_ids, mash_distances
