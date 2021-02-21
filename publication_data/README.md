# PUBLICATION DATA

## S1 Contact Maps of the proteins in the benchmark and predictions at RA0.

Files containing all plots of the contact maps. Real contacts are shown in grey for regions sufficiently covered by sequences in the MSA and pink for those without enough information. Red crosses above the diagonal (upper left) show the 40 top pairs with higher MI values. Under the diagonal (bottom right) are shown the 40 top pairs with higher DI values (in blue stars). MI and DI values are calculated using the effective sequence of the MSA using the standard alphabet of 20 residues (RA0).

**filename:** `ra0_contact_maps.tar.gz`

## S2 Contact Maps of the proteins in the benchmark and predictions at RA1.

Files containing all plots of the contact maps. Real contacts are shown in grey for regions sufficiently covered by sequences in the MSA and pink for those without enough information. Red crosses above the diagonal (upper left) show the 40 top pairs with higher MI values. Under the diagonal (bottom right) are shown the 40 top pairs with higher DI values (in blue stars). MI and DI values are calculated using the effective sequence of the MSA after grouping residues into a new residue alphabet defined as RA1.

**filename:** `ra1_contact_maps.tar.gz`

## S3 Contact Maps of the proteins in the benchmark and predictions at RA2.

Files containing all plots of the contact maps. Real contacts are shown in grey for regions sufficiently covered by sequences in the MSA and pink for those without enough information. Red crosses above the diagonal (upper left) show the 40 top pairs with higher MI values. Under the diagonal (bottom right) are shown the 40 top pairs with higher DI values (in blue stars). MI and DI values are calculated using the effective sequence of the MSA after grouping residues into a new residue alphabet defined as RA2.

**filename:** `ra2_contact_maps.tar.gz`

## S4 Contact Maps of the proteins in the benchmark and predictions at RA3.

Files containing all plots of the contact maps. Real contacts are shown in grey for regions sufficiently covered by sequences in the MSA and pink for those without enough information. Red crosses above the diagonal (upper left) show the 40 top pairs with higher MI values. Under the diagonal (bottom right) are shown the 40 top pairs with higher DI values (in blue stars). MI and DI values are calculated using the effective sequence of the MSA after grouping residues into a new residue alphabet defined as RA3.

**filename:** `ra3_contact_maps.tar.gz`

## S5 Input Sequence and Secondary Structure Prediction

Files containing sequence data and secondary structure assignation in multiFASTA format.

**filename:** `input_SecStr.tar.gz`

## S6 Multiple Sequence Alignments

**filename:** `input_msa.tar.gz`

This file is not included *"as is"* due to is size.

To generate it, run:

```bash
cat input_msa_splita* > input_msa.tar.gz
```

The final compiled file contains all the **M**ultiple **S**equence **A**lignments necessary to run **raDI** for the benchmark proteins.

## S7 Models of a benchmark

**filename:** `Supplementary_files_1.zip` and `Supplementary_files_2.zip`

Supplementary files of contact maps, outputs of RADI with alphabets RA0, RA1, RA2 and RA3, alignments with sMotif templates, template structures and models of a benchmark with 50 folds, distributed in folders per fold (under PDB code)
1N9L, 1H98, 1FR3, 1LSS, 1C02, 153L, 1LJ9, 1FCA, 1KU7, 1EZM, 1ATG, 1IG1, 1ID0, 1LUC, 1LR0, 1IR6, 1OR7, 1I52, 1GDT, 1D4A, 1FP6, 1QSA, 1K38, 1HW1, 2A0B, 1DI6, 1I8O, 1K20, 1G8K, 1M6K, 1F1U, 1J5Y, 1BOO, 1G72, 1KQ3, 1PJR, 1EK9, 1D5Y, 1M7J, 1G6O, 1IXC, 1MOQ, 1DD9,1GG4, 7REQ, 1G60, 1FEP, 1B7E, 1A0P, 1QKS

## S8 CCMPred results in the benchmark

**filename:** `CCMPred_1.zip` and `CCMPred_2.zip`

Supplementary files of contact maps and outputs of CCMPred of a benchmark with 50 folds, distributed in folders per fold (under PDB code)
1N9L, 1H98, 1FR3, 1LSS, 1C02, 153L, 1LJ9, 1FCA, 1KU7, 1EZM, 1ATG, 1IG1, 1ID0, 1LUC, 1LR0, 1IR6, 1OR7, 1I52, 1GDT, 1D4A, 1FP6, 1QSA, 1K38, 1HW1, 2A0B, 1DI6, 1I8O, 1K20, 1G8K, 1M6K, 1F1U, 1J5Y, 1BOO, 1G72, 1KQ3, 1PJR, 1EK9, 1D5Y, 1M7J, 1G6O, 1IXC, 1MOQ, 1DD9,1GG4, 7REQ, 1G60, 1FEP, 1B7E, 1A0P, 1QKS

## S9 python script to run CCMPred

Given a multiple sequence alignment in FastA format, for example 1TGA_A.msa, and the address where CCMPred is installed (e.g. /usr/local/share/ccmpred ), we apply a python script to obtain a set of residue-residue predicted contacts with CCMPred, that will be stored in the file “1ATG_A.top”.

The script is:

\# Path of CCMpred executable
path=”/usr/local/shared/ccmpred”

\# Data files
target=”1TGA_A.msa”
outputali=”1TGA_A.ali”
outputmat=”1ATG_A.mat”
output=”1ATG_A.top”

\# Change the format of the alignment
align = "python “+path+”/scripts/convert_alignment.py "+target+" fasta "+outputali
os.system(align)

\# CCMPred execution and time evaluation
start_time = time.time()
execute = path+“/bin/ccmpred "+outputali+" " +outputmat
os.system(execute)
print("Done for " + target + " in --- %s seconds ---" % (time.time() - start_time))

\# Selection of top pairs
select= "python “+path+”/scripts/top_couplings.py " + outputmat + “ > ”+output
os.system(select)


## S10 Table S1

**filename:** `Table_S1.xlsx`

Legend to Supplementary Table S1. Comparison of all models with the crystallographic structures. Sub-Table “Data” compares the quality of the models and the number of effective sequences in the MSA (“size MSA”). Columns show the data as in Table 2, using restraints derived with all alphabets (RA0, RA1, RA2 and RA3), for different measures of quality (TM-scores, RMSD and Z-score of Prosa2003), the number of effective sequences in the MSA, the length of the sequence, the coverage by sMotifs and the number of sequences in the original MSA. The last columns have included the CPU time, in seconds, that each approach has taken. Sub-Table “Number of sequences in MSA” compares the total of effective sequences in the MSA of each target with the original alignment obtained with the sequences of all potential homologs found with MMSeq2. Sub-Table “Graphs” analyzes the correlations of selected datasets of the quality of the models.



