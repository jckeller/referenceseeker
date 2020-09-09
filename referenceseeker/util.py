
import os
import subprocess as sp
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO

import referenceseeker.constants as rc


def read_reference_genomes(config):
    ref_genomes = {}
    with open(str(config['db_path'].joinpath('db.tsv')), 'r') as fh:
        for line in fh:
            if line[0] != '#':
                cols = line.strip().split('\t')
                accession_id = cols[0]
                ref_genomes[accession_id] = {
                    'id': accession_id,
                    'tax': cols[1],
                    'status': cols[2],
                    'name': cols[3]
                }
    return ref_genomes


def build_dna_fragments(genome_path, dna_fragments_path):
    """Build DNA fragments.

    :param genome_path: Path to source DNA Fasta file.
    :param dna_fragments_path: Path to DNA fragments output Fasta file.

    :rtype {idx, dna_fragment}: A dict containing index and DNA fragment objects.
    """

    dna_fragments = {}
    dna_fragment_idx = 0
    with dna_fragments_path.open(mode='w') as fh:
        try:
            for record in SeqIO.parse(str(genome_path), 'fasta'):
                sequence = str(record.seq)
                while len(sequence) > (rc.FRAGMENT_SIZE + rc.MIN_FRAGMENT_SIZE):  # forestall fragments shorter than MIN_FRAGMENT_SIZE
                    dna_fragment = sequence[:rc.FRAGMENT_SIZE]
                    dna_fragment_idx += 1
                    fh.write('>')
                    fh.write(str(dna_fragment_idx))
                    fh.write('\n')
                    fh.write(str(dna_fragment))
                    fh.write('\n')
                    dna_fragments[dna_fragment_idx] = {
                        'id': dna_fragment_idx,
                        'length': len(dna_fragment)
                    }
                    sequence = sequence[rc.FRAGMENT_SIZE:]
                dna_fragment = sequence
                dna_fragment_idx += 1
                fh.write('>')
                fh.write(str(dna_fragment_idx))
                fh.write('\n')
                fh.write(str(dna_fragment))
                fh.write('\n')
                dna_fragments[dna_fragment_idx] = {
                    'id': dna_fragment_idx,
                    'length': len(dna_fragment)
                }
                # sequence = sequence[rc.FRAGMENT_SIZE:]
        except ImportError:
            sys.exit("ERROR: Genome file is not in fasta format")
    return dna_fragments


def setup_configuration(args):
    """Test environment and build a runtime configuration."""
    try:
        config = {
            'tmp': Path(tempfile.mkdtemp()),
            'bundled-binaries': False,
            'threads': args.threads,
            'unfiltered': args.unfiltered,
            'bidirectional': args.bidirectional,
            'crg': args.crg,
            'ani': args.ani,
            'conserved_dna': args.conserved_dna,
            'n_mash_results': int(args.n_mash_results)
        }
    except ValueError:
        sys.exit('Error: n_mash_results must be a number ("integer")')
    set_path(config)

    base_dir = Path(__file__).parent.parent
    db_path = base_dir.joinpath('db')
    if os.access(str(db_path), os.R_OK & os.X_OK):
        config['db'] = db_path
    return config


def set_path(config):
    config['env'] = os.environ.copy()
    base_dir = Path(__file__).parent.parent
    share_dir = base_dir.joinpath('share')
    if os.access(str(share_dir), os.R_OK & os.X_OK):
        config['env']['PATH'] = str(share_dir) + ':' + config['env']['PATH']
        config['bundled-binaries'] = True


def test_binaries(config):
    """Test the proper installation of necessary 3rd party executables."""

    # test Mash
    try:
        sp.check_call(
            ['mash', 'dist', '-h'],
            env=config['env'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'Mash\' was not found!')
    except:
        sys.exit('ERROR: \'Mash\' was not exeutable!')

    # test nucmer
    try:
        sp.check_call(
            ['nucmer', '--help'],
            env=config['env'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'nucmer\' was not found!')
    except:
        sys.exit('ERROR: \'nucmer\' was not exeutable!')

    # test delta-filter
    try:
        sp.check_call(
            ['delta-filter', '-h'],
            env=config['env'],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL
        )
    except FileNotFoundError:
        sys.exit('ERROR: \'delta-filter\' was not found!')
    except:
        sys.exit('ERROR: \'delta-filter\' was not exeutable!')
    return


def check_path(path):
    """checks if a given path is a viable directory"""
    path = Path(path)
    if not os.path.exists(path):
        raise FileNotFoundError
    if not os.access(str(path), os.R_OK):
        raise PermissionError
    if path.stat().st_size == 0:
        raise OSError
    path = path.resolve()
    return path
