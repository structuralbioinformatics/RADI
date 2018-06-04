import os,sys,re
import argparse
import ConfigParser
import shutil
import string
# Add "." to sys.path #
src_path =  os.path.abspath(os.path.dirname(__file__))
sys.path.append(src_path)

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(src_path, "config.ini")
config.read(config_file)

# Read SBI path
sbi_path = os.path.join(src_path, config.get("Paths", "sbi_path"))
sys.path.append(sbi_path)
from SBI.structure import *

def main():

 options = parse_user_arguments()
 split(options)

def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Split the structure of a macrocomplex in chains",
        epilog      = "@oliva's lab 2015")
    parser.add_argument('-i','--pdb_file',dest='pdb',action = 'store',
                        help = 'PDB Input file')
    parser.add_argument('-o','--output_file',dest='out',action = 'store', default='output',
                        help = 'Output rootname for files (default is output)')
    parser.add_argument('-v','--verbose',dest='show',action = 'store_true',
                        help = 'Verbose execution ')

    options=parser.parse_args()

    return options


def fileExist (file):
 if file is not None:
  return os.path.exists(file) and os.path.isfile(file)
 else:
  return False



def printverbose(f,flag,message):
 """Define verbose to print in file 'f', if 'flag'=True, message given as 'message'"""
 if flag: f.write("%s"%(message))


def split(options):

 if fileExist(options.pdb):
   pdb=PDB(options.pdb)
 else:
  sys.stderr.write("Missing PDB file %s \n"%options.pdb)


 for code in pdb.chain_identifiers:
   chain=pdb.get_chain_by_id(code)
   pdb_new=PDB()
   pdb_new.add_chain(chain)
   new_name=options.out+code+".pdb"
   new_fasta=options.out+code+".fa"
   name=options.out.split("/")[-1]+code
   if chain.chaintype == "N": sequence=chain.nucleotide_sequence()
   if chain.chaintype == "P": sequence=chain.protein_sequence
   pdb_new.clean
   pdb_new.write(new_name)
   fasta=open(new_fasta,"w")
   fasta.write(">%s\n%s\n"%(name,sequence))
   fasta.close()
 
 

if  __name__ == "__main__":
    main()



