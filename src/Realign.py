import string
import sys
import os
import re
from Bio import SeqIO
from Bio import ExPASy
from Bio import AlignIO
from Bio.Align.Applications import *
from urllib import *
import argparse
import subprocess


def main():
 options = parse_user_arguments()
 try:
  alignss(options)
 except IOError as e:
  print("I/O error (%s): %s" %(e.errno, e.strerror))

def parse_user_arguments(*args, **kwds):
 parser = argparse.ArgumentParser(
  description = "Realign the MSA with the secondary structure of the seed",
  epilog      = "@oliva's lab 2015")
 parser.add_argument('-ssa','--SS_align_file',dest='ssa',action = 'store',default=None,
  help = 'Input file of the alignment between the seed sequence and its secondary structure')
 parser.add_argument('-msa','--MSA_input_file',dest='msa',action = 'store',default=None,
  help = 'Input file of the multiple sequence alignment (the seed MUST BE the first sequence)')
 parser.add_argument('-o','--fasta_output_file',dest='output',action = 'store',default=sys.stdout,
  help = 'Output file with the correct alignment of the seed and the secondary structure')
 parser.add_argument('-exe','--clustalw_exe',dest='exe',action = 'store',default='/soft/bio/sequence/clustalw-2.0.12/bin/clustalw2',
  help = 'Executable file of t_coffee (i.e.  /soft/bio/sequence/clustalw-2.0.12/bin/clustalw2)')
# parser.add_argument('-shell','--shell_environment',dest='shell',action = 'store',default='/usr/bash',
#  help = 'Shell environment to run execution of multiple sequence alignment (i.e. /usr/bash)')
 parser.add_argument('-v','--verbose',dest='show',action = 'store_true',default=False,
  help = 'Verbose execution')

 options=parser.parse_args()
 return options

def fileExist (file):
    if file is not None: return os.path.exists(file) and os.path.isfile(file)
    else: return False

def printverbose(f,flag,message):
 """Define verbose to print in file 'f', if 'flag'=True, message given as 'message'"""
 if flag: f.write("%s"%(message))

def realign(a,b,template):
 """Realign the 'b' gapped sequence aligned with 'a' upon a new 'template' of 'a' """
 c=''
 k=0
 j=0
 while (k<len(a) and j<len(template)):
  if (template[j] ==  a[k] and not template[j] == '-'):
    c+=b[k]
    k=k+1
    j=j+1
  elif  (template[j] == '-'):
    c+='-'
    j=j+1
  elif  (a[k] == '-'):
    k=k+1
  else:
    c+='-'
    k=k+1
    j=j+1
 return c

def alignss(options):

 if not fileExist(options.ssa):
   print("Missing Secondary Structure Alignment file\n")
   sys.exit(10)

 if not fileExist(options.msa):
   print("Missing MSA file\n")
   sys.exit(10)

 log=sys.stderr

 ssa=AlignIO.read(options.ssa,'fasta')
 msa=AlignIO.read(options.msa,'fasta')

 fd=open("tmp.fa","w")
 fd.write(">{0:s}\n{1:s}\n".format(ssa[0].name,ssa[0].seq))
 fd.write(">{0:s}\n{1:s}\n".format(msa[0].name,msa[0].seq))
 fd.close()
 msa_exe=options.exe
 #msa_cline=ClustalwCommandline(msa_exe,infile="tmp.fa",outfile="tmp.aln")
 #child = subprocess.Popen(str(msa_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell="/bin/bash")
 #child.communicate()
 os.system("%s tmp.fa >& tmp_clustalw.log"%msa_exe)
 cmp=AlignIO.read("tmp.aln",'clustal')
 seed_seq_str=ssa[0].seq
 seed_seq_msa=msa[0].seq
 seed_ss_str=ssa[1].seq
 seed_cmp_str=cmp[0].seq
 seed_cmp_msa=cmp[1].seq

 seed_ss_cmp=''
 k=0
 j=0
 while (k<len(seed_seq_str) and j<len(seed_cmp_str)):
  #print("Check %s[%d] vs %s[%d] \n"%(seed_cmp_str[j],j,seed_seq_str[k],k))
  if (seed_cmp_str[j] ==  seed_seq_str[k] and not seed_cmp_str[j] == '-'):
    seed_ss_cmp+=seed_ss_str[k]
    k=k+1
    j=j+1
  elif  (seed_cmp_str[j] == '-'):
    seed_ss_cmp+='-'
    j=j+1
  elif  (seed_seq_str[k] == '-'):
    if (k+1<len(seed_seq_str)):
     printverbose(log,options.show,("Skip secondary structure known gap %d (.. %s ..)\n"%(k,seed_seq_str[k-1]+seed_seq_str[k]+seed_seq_str[k+1])))
    else:
     printverbose(log,options.show,("Skip secondary structure known gap %d (.. %s ..)\n"%(k,seed_seq_str[k-1]+seed_seq_str[k])))
    #seed_ss_cmp+=seed_ss_str[k]
    k=k+1
  else:
    print("Error to align SS: sequence unmatch %s [%d] vs %s [%d]\n"%(seed_cmp_str[j],j,seed_seq_str[k],k))
    seed_ss_cmp+='-'
    k=k+1
    j=j+1


 fd=open("tmp.ss_modified.fa","w")
 fd.write(">{0:s}\n{1:s}\n".format(ssa[0].name,seed_cmp_str))
 fd.write(">{0:s}\n{1:s}\n".format(ssa[1].name,seed_ss_cmp))
 fd.close()


 seed_ss_msa=''
 k=0
 j=0
 while (j<len(seed_seq_msa) and k<len(seed_cmp_msa)  ):
  if (seed_cmp_msa[k] ==  seed_seq_msa[j]):
    if (k<len(seed_ss_cmp)):
      seed_ss_msa+=seed_ss_cmp[k]
    else:
      seed_ss_msa+="-"
    k=k+1
    j=j+1
  elif  (seed_seq_msa[j] == '-'):
    seed_ss_msa+='-'
    if (j+1<len(seed_seq_msa)):
     printverbose(log,options.show,("Skip secondary structure in MSA %d (.. %s ..)\n"%(j,seed_seq_msa[j-1]+seed_seq_msa[j]+seed_seq_msa[j+1])))
    else:
     printverbose(log,options.show,("Skip secondary structure in MSA %d (.. %s ..)\n"%(j,seed_seq_msa[j-1]+seed_seq_msa[j])))
    j=j+1
  elif  (seed_cmp_msa[k] == '-'):
    if (k+1<len(seed_cmp_msa)):
     printverbose(log,options.show,("Skip secondary structure in TMP %d (.. %s ..)\n"%(k,seed_cmp_msa[k-1]+seed_cmp_msa[k]+seed_cmp_msa[k+1])))
    else:
     printverbose(log,options.show,("Skip secondary structure in TMP %d (.. %s ..)\n"%(k,seed_cmp_msa[k-1]+seed_cmp_msa[k])))
    k=k+1
  else:
    print("Error to assign SS: sequence unmatch %s [%d] vs %s [%d]\n"%(seed_cmp_str[j],j,seed_seq_str[k],k))
    seed_ss_msa+='-'
    k=k+1
    j=j+1

 fd=options.output
 if not options.output ==  sys.stdout: fd=open(options.output,"w")

 fd.write(">{0:s}\n{1:s}\n".format(msa[0].name,seed_seq_msa))
 fd.write(">{0:s}\n{1:s}\n".format(ssa[1].name,seed_ss_msa))
 fd.close()

 if not options.show:
  delete_tmp=subprocess.Popen(["rm", "-f","tmp.ss_modified.fa", "tmp.fa", "tmp.aln", "tmp_clustalw.log"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  for line in delete_tmp.stdout:
    print(line.strip())
  for line in delete_tmp.stderr:
    print(line.strip())
 else:
  print("Check temporary files 'tmp.fa' 'tmp.ss_modified.fa' and 'tmp.aln' for the cross-alignments\n")

if __name__=="__main__":
  main()

