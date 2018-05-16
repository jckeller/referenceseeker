# ReferenceSeeker: Fast determination of reference genomes.
Author: Oliver Schwengers (oliver.schwengers@computational.bio.uni-giessen.de)


## Contents
- Description
- Input & Output
- Installation
- Usage
- Examples
- Databases
- Dependencies
- Citation


## Description
ReferenceSeeker determines closely related and finished reference genomes from
RefSeq (<https://www.ncbi.nlm.nih.gov/refseq>) following a hierarchical approach
combining a kmer based lookup and ANI calculations.

ReferenceSeeker computes kmer based genome distances between a query genome and
and a database built on finished RefSeq genomes via Mash (Ondov et al. 2016).
Currently, ReferenceSeeker offers bacterial, archeael and viral databases.
For subsequent candidates ReferenceSeeker computes ANI (average nucleotide identity)
values picking genomes meeting community standard thresholds (Goris, Konstantinos et al. 2007)
ranked by ANI and conserved DNA. Additionally, ReferenceSeeker can use MeDuSa
(Bosi, Donati et al. 2015) to build scaffolds based on the 20 closest reference genomes.


## Input & Output
Input:
draft or finished genomes in fasta format

Output:
tab separated to STDOUT comprising the following columns:
- RefSeq ID
- ANI
- conserved DNA
- NCBI Taxonomy ID
- Organism (genus species strain)


## Installation
To setup ReferenceSeeker just do the following:
1. clone the latest version of the repository
2. set REFERENCE_SEEKER_HOME environment variable pointing to the repository directory
3. download and extract the databases or create one yourself

Example:
```
git clone git@github.com:oschwengers/referenceseekr.git
export REFERENCE_SEEKER_HOME=referenceseekr
wget db...
tar -xzf db.tar.gz
rm db.tar.gz
```

Alternatively, just use the aforementioned Docker image (oschwengers/referenceseekr) in order to ease the setup process.


## Usage
Usage:
    --db          ReferenceSeeker database path
    --genome      Target draft genome
    --scaffolds   Build scaffolds via MeDuSa (Bosi, Donati et al. 2015) based on detected references
    --output      Output fasta file for built scaffolds
    --cpus        Number of CPUs to use (default = all available)
    --help        Show this help message and exit
    --verbose     Print verbose information


## Examples
Simple:
```
referenceseekr.py --db <REFERENCE_SEEKER_DB> --genome <GENOME>
```

Expert: creating scaffolds with verbose output using defined # of CPUs:
```
referenceseekr.py --db <REFERENCE_SEEKER_DB> --genome <GENOME> --scaffolds --output scaffolds.fasta --verbose --cpus 8
```

With Docker:
```
docker pull oschwengers/referenceseekr:latest
docker run --rm -v <REFERENCE_SEEKER_DB>:/db -v <DATA_DIR>:/data oschwengers/referenceseekr:latest --genome <GENOME>
```

With Docker shell script:
```
docker pull oschwengers/referenceseekr:latest
referenceseekr.sh <REFERENCE_SEEKER_DB> <GENOME>
```


## Databases
ReferenceSeeker depends on custom databases based on complete NCBI RefSeq genomes comprising 
kmer hash subsets as well as fasta files.
These databases (RefSeq release 87) can be downloaded HERE: (type, # complete genomes, size zipped, size unzipped)
- bacteria, 9443, 12 Gb, 37 Gb
- archaea, 255, 205 Mb, 642 Mb
- viral, 7532, 508 Mb, 762 Mb

The latest versions can be built using a custom nextflow pipeline.
Valid values for `DB_TYPE` are:
- 'archaea'
- 'bacteria'
- 'viral'

Download and install Nextflow:
```
curl -fsSL get.nextflow.io | bash
```

Build database:
```
export REFERENCE_SEEKER_HOME=<REFERENCE_SEEKER_DIR>
sh build-db.sh <DB_TYPE>
```

## Dependencies
ReferenceSeeker depends on the following packages:
- Python (3.5.2) and BioPython (1.66)
- Mash (2.0) <https://github.com/marbl/Mash>
- MUMmer (4.0.0-beta) <https://github.com/gmarcais/mummer>
- MeDuSa (1.6) <https://github.com/combogenomics/medusa>

ReferenceSeeker has been tested against aforementioned versions.


## Citation
Manuscript is in preparation.
