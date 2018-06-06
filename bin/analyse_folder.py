import os, sys, re
import gzip

def main():

  files_DI={}
  files_DI.setdefault("ra0",[x for x in os.listdir("./") if x.endswith("_ra0_DI.out")])
  files_DI.setdefault("ra1",[x for x in os.listdir("./") if x.endswith("_ra1_DI.out")])
  files_DI.setdefault("ra2",[x for x in os.listdir("./") if x.endswith("_ra2_DI.out")])
  files_DI.setdefault("ra3",[x for x in os.listdir("./") if x.endswith("_ra3_DI.out")])

  data_output={}
  for ra,list_of_files in files_DI.iteritems():
    for file_DI in list_of_files:
      d=parser_output(file_DI)
      if d is not None:
         data_output.setdefault(ra,[]).append(d)

# Plotting time

  for ra,list_of_data in data_output.iteritems():
      plot_MI = [(x["eLength"],x["time_MI"]) for x in list_of_data]
      plot_DI = [(x["eLength"],x["time_DI"]) for x in list_of_data]
 
      fM = open("plot_time_MI_"+ra+".dat","w")
      fD = open("plot_time_DI_"+ra+".dat","w")
      
      for x,y in plot_MI:
          fM.write("%f\t%f\n"%(x,y))
      for x,y in plot_DI:
          fD.write("%f\t%f\n"%(x,y))

      fM.close()
      fD.close()

  gnu_file_MI = open("plot_time_MI.gnu","w")
  gnu_file_MI.write("set term gif\n")
  gnu_file_MI.write("set xrange [0:]\n")
  gnu_file_MI.write("set yrange [0:]\n")
  gnu_file_MI.write("set output 'plot_time_MI.gif' \n")
  gnu_file_MI.write("plot 'plot_time_MI_ra0.dat' title 'ra0' lt rgb 'black' , 'plot_time_MI_ra1.dat' title 'ra1' lt rgb 'blue', 'plot_time_MI_ra2.dat' title 'ra2' lt rgb 'green' , 'plot_time_MI_ra3.dat' title 'ra3' lt rgb 'red' ")
  gnu_file_MI.close()

  gnu_file_DI = open("plot_time_DI.gnu","w")
  gnu_file_DI.write("set term gif\n")
  gnu_file_DI.write("set xrange [0:]\n")
  gnu_file_DI.write("set yrange [0:]\n")
  gnu_file_DI.write("set output 'plot_time_DI.gif' \n")
  gnu_file_DI.write("plot 'plot_time_DI_ra0.dat' title 'ra0' lt rgb 'black' , 'plot_time_DI_ra1.dat' title 'ra1' lt rgb 'blue', 'plot_time_DI_ra2.dat' title 'ra2' lt rgb 'green' , 'plot_time_DI_ra3.dat' title 'ra3' lt rgb 'red' ")
  gnu_file_DI.close()

# Boxplots to compare ra1,ra2,ra3 with ra0

  compare_alphabet("MI_similar_pairs","pairs_MI",data_output)
  compare_alphabet("DI_similar_pairs","pairs_DI",data_output)

# Boxplots of Success, coverage and Effective number of sequences

  make_boxplot_files("coverage","Coverage",data_output)
  make_boxplot_files("Effective","eNseq",data_output)
  make_boxplot_files("TP_MI_CA","tp_MI_ca",data_output)
  make_boxplot_files("TP_MI_CB","tp_MI_cb",data_output)
  make_boxplot_files("TP_MI_MIN","tp_MI_min",data_output)
  make_boxplot_files("TP_MI","tp_MI",data_output)
  make_boxplot_files("TP_DI_CA","tp_DI_ca",data_output)
  make_boxplot_files("TP_DI_CB","tp_DI_cb",data_output)
  make_boxplot_files("TP_DI_MIN","tp_DI_min",data_output)
  make_boxplot_files("TP_DI","tp_DI",data_output)


def compare_alphabet(title,key,data):
   pairs={}
   box_file=open("box_"+title+".dat","w")
   for ra,list_of_data in data.iteritems():
       pairs.setdefault(ra,{})
       assign={}
       for x in list_of_data:
           assign.setdefault(x["PDB"],x[key])
       pairs[ra]=assign
   pairs_original=pairs["ra0"]
   for ra in data.iterkeys():
      if ra != "ra0":
         name=os.path.abspath(".")+"/"+title+"_"+ra+".dat"
         fc=open(name,"w")
         pairs_test=pairs[ra]
         for p,list_of_pairs in pairs_test.iteritems():
           if pairs_original.has_key(p):
             original_set_of_pairs = set(pairs_original[p])
             n=check_similar_pairs(original_set_of_pairs,list_of_pairs,9)
             fc.write("%s\t%f\n"%(p,(float)(n)))
           else:
             print("%s has no list of RA0 pairs\n"%p)
         fc.close()
         box_file.write("%s\t%s\n"%(ra,name))
   box_file.close()
       

def make_boxplot_files(title,key,data):
   box_file=open("box_"+title+".dat","w")
   for ra,list_of_data in data.iteritems():
       name=os.path.abspath(".")+"/"+title+"_"+ra+".dat"
       y=[(x["PDB"],x[key]) for x in list_of_data]
       fc=open(name,"w")
       for p,c in y:
         fc.write("%s\t%f\n"%(p,c))
       fc.close()
       box_file.write("%s\t%s\n"%(ra,name))
   box_file.close()


def check_similar_pairs(set_of_pairs,list_of_pairs,interval):
   pairs_to_test=set()
   origina_pull=set_of_pairs
   n=0
   for (x,y) in list_of_pairs:
        pairs_to_test=set()
        pairs_to_test=set([(a,b) for a in xrange(x-interval/2,x+interval/2+1) for b in xrange(y-interval/2,y+interval/2+1)])
        pairs_to_test.union(set([(a,b) for a in xrange(y-interval/2,y+interval/2+1) for b in xrange(x-interval/2,x+interval/2+1)]))
        found=0
        for (a,b) in pairs_to_test:
          if  not found:
           for (i,j) in origina_pull:
             if (a,b) == (i,j):
                n=n+1
                found=1
                remove_pair=(i,j)
                break
          if found:
             origina_pull.remove(remove_pair)
             break
   return n        
                
       


def parser_output(file_DI):
  print("open %s"%file_DI)
  fd = open(file_DI,"r")
  pdb=os.path.basename(file_DI)[0:6]
  read_file=0
  read_rank=0
  pairs_MI=[]
  pairs_DI=[]
  for line in fd:
      if not line.startswith("#") and not read_file: continue
      if read_rank>=40: break    
      if line.startswith("#Length:"): Length = line.strip().split()[-1]
      if line.startswith("#Effective Length:"): eLength = line.strip().split()[-1]
      if line.startswith("#Coverage:"):  Coverage = line.strip().split()[-1]
      if line.startswith("#Num. Sequences:"):  Nseq = line.strip().split()[-1]
      if line.startswith("#Num. effective:"):  eNseq = line.strip().split()[-1]
      if line.startswith("#Execution time MI:"):  time_MI = line.strip().split()[-1]
      if line.startswith("#Execution time DI:"):  time_DI = line.strip().split()[-1]
      if line.startswith("#True Posititves MI:"):  
             word = line.strip().split()
             tp_MI_ca  =  word[3]
             tp_MI_cb  =  word[5]
             tp_MI_min =  word[7]
             tp_MI     =  word[9]
      if line.startswith("#True Posititves DI:"):  
             word = line.strip().split()
             tp_DI_ca  =  word[3]
             tp_DI_cb  =  word[5]
             tp_DI_min =  word[7]
             tp_DI     =  word[9]
      if line.startswith("#Rank"): 
         read_file=1
         continue
      if read_file and read_rank < 40:
         word = line.strip().split()
         x_MI = word[1]
         y_MI = word[2]
         if x_MI != "-" and y_MI != "-":
            pairs_MI.append(((int)(x_MI),(int)(y_MI)))
         x_DI = word[4]
         y_DI = word[5]
         if x_DI != "-" and y_DI != "-":
            pairs_DI.append(((int)(x_DI),(int)(y_DI)))
         read_rank = read_rank + 1
         continue
  fd.close()

  if read_file == 1:
     data= { "PDB":pdb,"Length":(float)(Length), "eLength":(float)(eLength), "Coverage":(float)(Coverage), "Nseq":(float)(Nseq) , "eNseq":(float)(eNseq) , 
          "time_MI":(float)(time_MI) , "time_DI":(float)(time_DI) , "pairs_MI":pairs_MI,  "pairs_DI":pairs_DI,
          "tp_MI_ca":(float)(tp_MI_ca) , "tp_MI_cb":(float)(tp_MI_cb) , "tp_MI_min":(float)(tp_MI_min) , "tp_MI":(float)(tp_MI) , 
          "tp_DI_ca":(float)(tp_DI_ca) , "tp_DI_cb":(float)(tp_DI_cb) , "tp_DI_min":(float)(tp_DI_min) , "tp_DI":(float)(tp_DI) }
  else:
     data=None
  return data


if __name__ == '__main__':
    main()


