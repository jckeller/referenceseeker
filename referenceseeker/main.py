
import argparse
import multiprocessing as mp
import sys

import referenceseeker
import referenceseeker.constants as rc
import referenceseeker.util as util
import referenceseeker.single as single


def main():
    # parse options and arguments
    parser = argparse.ArgumentParser(
        prog='referenceseeker',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Rapid determination of appropriate reference genomes.',
        epilog="Citation:\n%s\n\nGitHub:\nhttps://github.com/oschwengers/referenceseeker" % rc.CITATION,
        add_help=False
    )
    parser.add_argument('db', metavar='<database>', help='ReferenceSeeker database path')
    parser.add_argument('genome', metavar='<genome>', help='target draft genome in fasta format')
    group_workflow = parser.add_argument_group('Filter options / thresholds', 'These options control the filtering and alignment workflow.')
    group_workflow.add_argument('--crg', '-r', action='store', type=int, default=100, help='Max number of candidate reference genomes to pass kmer prefilter (default = 100)')
    group_workflow.add_argument('--ani', '-a', action='store', type=float, default=0.95, help='ANI threshold (default = 0.95)')
    group_workflow.add_argument('--conserved-dna', '-c', action='store', dest='conserved_dna', type=float, default=0.69, help='Conserved DNA threshold (default = 0.69)')
    group_workflow.add_argument('--unfiltered', '-u', action='store_true', help='Set kmer prefilter to extremely conservative values and skip species level ANI cutoffs (ANI >= 0.95 and conserved DNA >= 0.69')
    group_workflow.add_argument('--bidirectional', '-b', action='store_true', help='Compute bidirectional ANI/conserved DNA values (default = False)')

    group_runtime = parser.add_argument_group('Runtime & auxiliary options')
    group_runtime.add_argument('--help', '-h', action='help', help='Show this help message and exit')
    group_runtime.add_argument('--version', '-V', action='version', version='%(prog)s ' + referenceseeker.__version__)
    group_runtime.add_argument('--verbose', '-v', action='store_true', help='Print verbose information')
    group_runtime.add_argument('--threads', '-t', action='store', type=int, default=mp.cpu_count(), help='Number of used threads (default = number of available CPU cores)')

    subparsers = parser.add_subparsers(dest='subcommand', help='sub-command help')
    # add "single" sub-command option
    parser_single = subparsers.add_parser('single', help='start reference genome search for single genome')
    # add "cohort" sub-command option
    parser_cohort = subparsers.add_parser('cohort', help='start reference genome search for genome cohort')

    args = parser.parse_args()

    # setup global configuration
    config = util.setup_configuration(args)
    util.test_binaries(config)

    # check parameters
    try:
        config['db_path'] = util.check_path(args.db)
    except FileNotFoundError:
        sys.exit('ERROR: database directory is not readable!')
    except PermissionError:
        sys.exit('ERROR (permission): database directory is not accessible')
    except OSError:
        sys.exit('ERROR: database directory is empty')

    try:
        config['genome_path'] = util.check_path(args.genome)
    except FileNotFoundError:
        sys.exit('ERROR: genome file is not readable!')
    except PermissionError:
        sys.exit('ERROR (permission): genome file is not accessible')
    except OSError:
        sys.exit('ERROR: genome file (%s) is empty!' % config['genome_path'])

    # print verbose information
    if args.verbose:
        print("ReferenceSeeker v%s" % referenceseeker.__version__)
        print('Options, parameters and arguments:')
        print("\tuse bundled binaries: %s" % str(config['bundled-binaries']))
        print("\tdb path: %s" % str(config['db_path']))
        print("\tgenome path: %s" % str(config['genome_path']))
        print("\ttmp path: %s" % str(config['tmp']))
        print("\tunfiltered: %s" % str(config['unfiltered']))
        print("\tbidirectional: %s" % str(config['bidirectional']))
        print("\tANI: %0.2f" % config['ani'])
        print("\tconserved DNA: %0.2f" % config['conserved_dna'])
        print("\t# CRG: %d" % config['crg'])
        print("\t# threads: %d" % config['threads'])

    if args.subcommand == 'single':
        single.single(args, config)

    elif args.subcommand == 'cohort':
        pass

    else:
        parser.print_help()
        sys.exit("Error: no subcommand provided!")


if __name__ == '__main__':
    main()
