cat marks_et_al_2011.fasta | perl -e '%seqs=();while(<>){chomp;if(/^>(.+)$/){$header=$1;if(length($header)==5){$seqs{$header}="";}}elsif(length($header)==5){$seqs{$header}.=$_}}foreach $header(keys %seqs){open($f,">","$header.fa");print $f "$header\n$seqs{$header}\n";close($f);}'
ls *.fa | perl -e '$job=1;while(<>){chomp;$header=substr($_,0,5);open($f,">","$job.sh");print $f "python /home/ofornes/RADI/msa/buildmsa.py --dummy=/home/ofornes/scratch/RADI/paper/$header/tmp/ -i /home/ofornes/scratch/RADI/paper/$header.fa -o /home/ofornes/scratch/RADI/paper/$header/\n";close($f);$job++;}'

