This folder contains the script `buildmsa.py` for multiple sequence alignment (MSA) construction.

## Dependencies
It requires the following dependencies:
* [`FAMSA`](http://sun.aei.polsl.pl/REFRESH/famsa)
* [`MMseqs2`](https://github.com/soedinglab/mmseqs2)
* [`Python 2.7`](https://www.python.org/download/releases/2.7/)
* [UniRef](https://www.uniprot.org/help/uniref)100, 90 and 50 databases

## Usage
And the following inputs:
* The full path to the directory containing the `FAMSA` binaries (option `-f`)
* An input sequence (*i.e.* query) in [`FASTA`](https://en.wikipedia.org/wiki/FASTA_format) format (option `-i`)
* The full path to the directory containing the `MMseqs2` binaries (option `-m`)
* The full path to where the bash script `uniref.sh` was executed (option `-u`)

The bash script `uniref.sh` creates `MMseqs2`-formatted UniRef databases.

Moreover, users can specify different `UniRef` databases using the options `-n` (*i.e.* non-redundant) and `-r` (*i.e.* redundant). By default, the script sets the options `-n` and `-r` to `UniRef50` and `UniRef100`, respectively.

Other non-mandatory options include:
* The full path to a dummy directory for `MMseqs2` (option `-dummy`)
* The full path to the directory to where to output the MSA (option `-o`)
* The maximum number of sequences to include in the MSA (option `-s`; by default is set to 100,000 sequences)
* The total number of cores to be used during the process (option `-t`; by default the script uses 1 core)

The construction of the MSA is as follows: First, the script builds a profile of the query searching for similar sequences in the non-redundant database with `MMseqs2` for 4 iterations. Next, it uses the query profile to search for more sequence relatives in the redundant database. Then, the script builds a MSA of the query and the identified sequences (up to 100,000) with `FAMSA`. Finally, it removes any columns that would cause insertions to the query from the MSA. `MMseqs2` is executed with options `-s 7.5` and `--max-seq-id 1.0` for a more sensitive search.

**Important** usage notes:
1) When searching the `UniRef100` database, `MMseqs2` requires ~500GB of memory (refer to [`MMseqs2` documentation](https://github.com/soedinglab/mmseqs2/wiki#memory-consumption) for more details on memory requirements).
2) The execution times of `FAMSA` and `MMseqs2` decrease substantially with the number of cores used (**8 or more cores is highly recommended**).
3) The execution time of `FAMSA` will [increase with the number of sequences to be aligned](https://www.nature.com/articles/srep33964/figures/3) (*i.e.* option `-s`).