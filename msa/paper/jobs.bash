cat marks_et_al_2011.fasta | perl -e '%seqs=();while(<>){chomp;if(/^>(.+)$/){$header=$1;$seqs{$header}=""}else{$seqs{$header}.=$_}}foreach $header(keys %seqs){open($f,">","$header.fa");print $f "$header\n$seqs{$header}\n";close($f);}'
echo "#!/bin/bash" > buildmsa.pbs
echo "#SBATCH --account=rrg-wyeth" >> buildmsa.pbs
echo "#SBATCH --job-name=buildmsa" >> buildmsa.pbs
echo "#SBATCH --mail-type=FAIL" >> buildmsa.pbs
echo "#SBATCH --mail-user=oriol@cmmt.ubc.ca" >> buildmsa.pbs
echo "#SBATCH --mem=500GB" >> buildmsa.pbs
echo "#SBATCH --nodes=1" >> buildmsa.pbs
echo "#SBATCH --ntasks=32" >> buildmsa.pbs
echo "#SBATCH --output=buildmsa.log" >> buildmsa.pbs
echo "#SBATCH --time=7-00:00" >> buildmsa.pbs
ls *.fa | perl -e '$job=1;while(<>){chomp;$header=substr($_,0,5);open($f,">","$job.sh");system("echo \""python /home/ofornes/RADI/msa/buildmsa.py --dummy=/home/ofornes/scratch/RADI/paper/$header/tmp/ -i /home/ofornes/RADI/msa/paper/$header.fa -o /home/ofornes/scratch/RADI/paper/$header/\n"\" >> buildmsa.pbs")$job++;}'
