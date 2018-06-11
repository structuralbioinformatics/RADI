# RADI
RADI (Reduced Alphabet Direct Information) is a variation of the direct-coupling analysis (DCA) algorithm that groups physico-chemically equivalent amino acids prior to the computation of DCA. RADI achieves similar results than the original DCA while reducing the computation time by more than 30-fold for sequences with length of ~1,000 residues. Moreover, RADI overcomes the quality of predictions based on mutual information (MI); the accuracies and distributions of contacts predicted by RADI are higher and more useful than those predicted using MI, while reducing the time of computation more than 180-fold.

Please cite: *TO BE FILLED*

## Content
The repository is organized as follows:
* The `bin` folder contains something else
* The `msa` folder contains scripts for multiple sequence alignment (MSA) construction
* The `src` folder contains the RADI program

## Dependencies
RADI requires the following dependencies:
* [`FAMSA`](http://sun.aei.polsl.pl/REFRESH/famsa)
* [`MMseqs2`](https://github.com/soedinglab/mmseqs2)
* [`Python 2.7`](https://www.python.org/download/releases/2.7/) with the [`Biopython`](http://biopython.org), [`CoreAPI`](http://www.coreapi.org), [`numpy`](http://www.numpy.org), [`tqdm`](https://pypi.org/project/tqdm/) and [`UniProt`](https://github.com/boscoh/uniprot) libraries
* The [`UCSC binaries`](http://hgdownload.cse.ucsc.edu/admin/exe/) for standalone command-line use

## Configuration
Say that everything is easier with Miniconda and Homebrew...