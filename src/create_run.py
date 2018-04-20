import os, sys, re
import gzip

msa_folder="/home/boliva/PROJECTS/DIRECT_INFORMATION/bernat/0_MSA/buildali"
msa_files=[x.strip(".mfa") for x in os.listdir(msa_folder) if x.endswith("mfa")]

for p in msa_files:
    prot,chain=p.split("_")
    code = prot[1:3].lower()
    sys.stdout.write("csh run_raDI_exe.sh %s %s %s %s\n"%(prot,chain,prot.lower(),code))


