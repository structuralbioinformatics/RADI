# PUBLICATION DATA

This folder contains all the data necessary to reproduce the data presented in the final manuscript.

## `input_SecStr.tar.gz`

This file compiles all the sequence - secondary structure assignations for each protein of the benchmark.

## `input_msa.tar.gz`

This file is not included *"as is"* due to is size.

To generate it, run:

```bash
cat input_msa_splita* > input_msa.tar.gz
```

## `*_contact_maps.tar.gz`

Each one of this compressed files contains all the `.gif` outputs obtained according to the alphabet definition selected.


