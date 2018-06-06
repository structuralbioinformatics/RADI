# Compiling and running raDI

## Compiling

**raDI** is written in C code. Once in the ``src`` folder, it should be as straight forward as:

```bash
make raDI
```

This will create a executable both in the source folder as well as in the ``bin`` folder called ``raDI``.

## Running

Execution of ``raDI`` depends on the user being allowed to use an unlimited number of processors, a feature normally blocked in unix-based systems.

This property can be re-written with:

```bash
ulimit -s unlimited
```

for **bash** or:

```csh
limit  stacksize unlimited
```

for **csh**.

Without this, the executable will rise a ``segmentation fault``.

As this parameters do not affect normal functionality of the system, we recommend to add the command to the `~/.bashrc` or `~/.cshrc` configuration file.

Call `raDI -h` to see the options of the executable:

```
./raDI -h

#          Program for Mutual Information and Direct Information between positions of a Multiple Sequence Alignment(MSA)
#          It includes the possibility to check contact-pairs using a PDB file
#          Authors:
#                 Bernat Anton, Mireia Besalú, Gemma de las Cuevas,
#                 Oriol Fornes, Jaume Bonet &  Baldo Oliva
#          Structural Bioinformatics lab (GRIB)
#          Universitat Pompeu Fabra
#          Catalonia
#          Contact: baldo.oliva@upf.edu  (+34 933160509)
#          Last Update: March 2018


	 -msa 		 File with Multiple Sequence Alignments (MSA)
	 -ssa 		 File with Secondary Structure (SS) aligned with the seed in the MSA file
	 -xsa 		 File with sequence of PDB structure aligned with the seed in the MSA file
	 -pdb 		 PDB File of the SINGLE chain structure of the seed
	 -sdis		 Minimum distance in sequence between residue-pairs (default 10)
	 -lfrg		 Fragment length to avoid redundant results (selecting the best pair between two fragments, default 5)
	 -ra  		 Reduced alphabet of AA (default is no-reduction:0)
	      		    0 {-}{A}{C}{D}{E}{F}{G}{H}{I}{K}{L}{M}{N}{P}{Q}{R}{S}{T}{V}{W}{Y}
	      		    1 {-}{RKH}{DE}{STNQ}{AVLIM}{FWY}{C}{G}{P}
	      		    2 {-}{RKHDESTNQC}{AVLIMFWY}{G}{P}
	      		    3 {-}{RKHDESTNQCG}{AVLIMFWYP}
	 -o   		 Root name of the outputs with MI, DI and Contact-Map data for GNU plot.
	 -swap		 Swap the AA (or raAA) without affecting gaps
	 -v   		 Verbose
	 -h   		 This help (if MSA file is null help is also shown)
```