# Example using a known structure (1ATG_A.pdb)

# Executable Path
REPO=$(git rev-parse --show-toplevel)

# Input Files
MSALIGN=1ATG_A.mfa
DSSP=1ATG_A_2str.fa
SEQ=1ATG_pdb.fa
PDB=1ATG_A.pdb

# Output Prefix
OPREFIX=1ATG_A_with_pdb

${REPO}/bin/raDI -msa ${MSALIGN} -ssa ${DSSP} -xsa ${SEQ} -pdb ${PDB} -lfrg 9 -sdis 10 -v -ra 3 -o ${OPREFIX} >& ${OPREFIX}.log
