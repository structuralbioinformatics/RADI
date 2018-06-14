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
