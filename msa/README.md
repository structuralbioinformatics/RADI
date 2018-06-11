This folder contains scripts `multiple sequence alignment` (MSA) construction.

## Dependencies
* [`FAMSA`](http://sun.aei.polsl.pl/REFRESH/famsa)
* [`MMseqs2`](https://github.com/soedinglab/mmseqs2)
* [`Python 2.7`](https://www.python.org/download/releases/2.7/)

## Requirements
Note that `MMseqs2` requires ~500GB of memory to search the UniRef100 DB (refer to [`MMseqs2` documentation](https://github.com/soedinglab/mmseqs2/wiki#memory-consumption) for more details on memory consumption).

## Usage
The script `buildmsa.py` provides multiple sequence alignments for use as input by `raDI`. It requires the following inputs:
* The full path to the directory containing the `FAMSA` binaries (option `-f`)
* An input sequence (*i.e.* query) in [`FASTA`](https://en.wikipedia.org/wiki/FASTA_format) format (option `-i`)
* The full path to the directory containing the `MMseqs2` binaries (option `-m`)
* The full path to where the bash script `uniref.sh` was executed (option `-u`)

The script builds a profile of the query by searching a non-redundant UniRef DB (option `-n`; by default `uniref50`) with `MMseqs2` for 4 iterations, after which it switches to searching a (more) redundant UniRef DB (option `-r`; by default `uniref100`) with that profile. `MMseqs2`-formatted [UniRef 100, 90 and 50 databases](https://www.uniprot.org/help/uniref) can be obtained by executing the shell script `uniref.sh` located in the `./uniref/` folder.

Other non-mandatory options include:
* The full path to a dummy directory for MMseqs2 (option `-dummy`)
* The full path to the directory to where to output the MSA (option `-o`)
* The max. number of sequences to include in the MSA (option `-s`; by default is set to 100,000 sequences)
* The total number of cores to be used for the MSA construction (option `-t`; by default uses 1 core)

Note that: 1) the running times of both `FAMSA` and `MMseqs2` will decrease with the number of cores used (**8 or more cores is highly recommended**); and 2) the running time of `FAMSA` will [increase **exponentially**](https://www.nature.com/articles/srep33964/figures/3) with the number of sequences.