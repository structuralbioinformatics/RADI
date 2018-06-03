#!/usr/bin/perl 

#use lib '/sc2/rob/stamp.4.3/perl';
#use Palign;

$format{"c"} = "Clustal";
$format{"m"} = "MSF";
$format{"b"} = "AMPS BLOCK";
$format{"f"} = "FASTA";
$format{"p"} = "PIR";
$format{"e"} = "MEntA";
$format{"s"} = "Pfam";
$format{"h"} = "HMMalign";
$format{"y"} = "Phylip";

$from = "x"; 
$to   = "c"; # anything to clustalw by default
$nogap = 0;

for($i=0; $i<=$#ARGV; ++$i) {
   if($ARGV[$i] eq "-in") { # conversion specified
 	if(($i+1)>($#ARGV)) { exit_error(); }
	$from = substr($ARGV[$i+1],0,1);
	if($from !~ /[cmbfpsyh]/) { exit_error(); }
	$i++;
   } elsif($ARGV[$i] eq "-out") { # conversion specified
 	if(($i+1)>($#ARGV)) { exit_error(); }
	$to   = substr($ARGV[$i+1],0,1);
	if($to !~ /[cmbfpsye]/) { exit_error(); } 
	$i++;
   } elsif($ARGV[$i] eq "-nogap") { # no gaps
	$nogap = 1;
   } elsif($ARGV[$i] eq "-separate") { # no gaps
	$separate = 1;
   } else {
	exit_error();
   }
}

if(($nogap==1) && ($to !~ /[fp]/)) {
	die "Error: can only do nogap in FASTA or PIR output (-out f or -out p)\n";
}
if(($separate==1) && ($to !~ /[f]/)) {
	die "Error: can only do separate in FASTA output (-out f)\n";
}
#printf("Converting from %s to %s (stdin/stdout)\n",$format{$from},$format{$to});

# Read data

@data = read_align("-",$type);
if($from eq "x") {
    $from = $type;
}
if($from eq "x") {
   die "Couldn't figure out alignment type: please specify with \"-in <type>\"\n";
}

if($from eq "m") { 
	$align = get_msf(@data); 
} elsif($from eq "c") {
	 $align = get_clustal(@data); 
} elsif($from eq "b") {
	 $align = get_block(@data); 
} elsif($from eq "f") {
	 $align = get_afasta(@data); 
} elsif($from eq "p") {
	 $align = get_apir(@data); 
} elsif($from eq "e") {
	 $align = get_amenta(@data); 
} elsif($from eq "s") {
	 $align = get_pfam(@data); 
} elsif($from eq "h") {
	 $align = get_hmmalign(@data); 
} elsif($from eq "y") {
	 $align = get_phylip(@data); 
}



# Write data 
if($to eq "m") { 
	write_msf($align,"-"); 
} elsif($to eq "c") {
	write_clustal($align,"-");
} elsif($to eq "b") {
	write_block($align,"-");
} elsif($to eq "f") {
	write_afasta($align,"-",$nogap,$separate);
} elsif($to eq "p") {
	write_apir($align,"-",$nogap);
} elsif($to eq "e") {
	write_amenta($align,"-",$nogap);
} elsif($from eq "y") {
	 write_phylip(@data); 
}


sub exit_error {
die" format: aconvert [-in <type> -out <type>] < file1 > file2
        where <type> is one of:
          m = MSF format
          c = Clustal format
          b = AMPS BLOCK format
          f = FASTA format
          p = PIR format
	  s = PFAM (Sanger) format
        defaults are in: any out: 'b' \n";
}

# Rob's Perl alignment reading/modification/etc. routines

# Generally speaking, most routines deal with a reference
#  That stores alignments:
#
# $align->{file} - name of file that the data came from
# $align->{ids}  - hash containing data for each entry in the alignment
#  e.g  $align->{ids}{"MYG_HUMAN"}{seq}       - sequence data
#  e.g  $align->{ids}{"MYG_HUMAN"}{prebumpf}  - text data
#  e.g  $align->{ids}{"MYG_HUMAN"}{postbumpf} - text data
#		(the above is mostly to do with MSF for gmat, pre/post
#                refer to the bits before and after the accession in the
#                descriptor lines; for other formats, post should be
#                redundant, or they are the same)
#  e.g. $align->{ids}{"MYG_HUMAN"}{start}     - start field, if present
#  e.g. $align->{ids}{"MYG_HUMAN"}{end}       - end field, if present
#  e.g. $align->{ids}{"MYG_HUMAN"}{ranges}    - "start"-"end" string
#
# $align->{list} - array of strings giving ids in the order they were read
#                 in to the file.
#
# $align->{nseq} - number of sequences in the alignment
# $align->{alen} - length of the alignment
#
# Change (June 1999):
# Now "gaps" are always treated as spaces internally
# So get/write_clustal  will change "-" to " "
#    get/write_msf      will change "." to " "
#
# Change (August 1999):
# new array added to reference to deal with information stored
# in block (i.e. vertical) formatted alignments
#
# $align->{p}[0..N]{pre}   - text before the alignment
# $align->{p}[0..N]{post}  - text after the alignment 
#  [where 0..N covers the positions in the alignment]
#
# STAMP specific change
# looks for STAMP fields in get_block
# These will be stored as:
# $align->{stamp}{$t} where $t is the letter
# The data will simply be in the {post} field above (accessed by a split function)

#package Palign;
#require Exporter;
#@ISA = qw(Exporter);
#@EXPORT = qw(read_align get_msf write_msf get_clustal write_clustal get_block write_block get_afasta write_afasta 
#	    get_apir write_apir get_pfam get_hmmalign merge_align do_clustalw_align extract_align get_phylip write_phylip
#	    read_exchange_matrix  get_print_string prettify write_one_fasta);

#use Search::Dict;

sub read_align { # Just read in text and modify format if necessary

#
# Likely format is determined by looking for signs of one format over another
#  and storing them in %votes.  This won't always work, of course, but
#  generally the user can force one format over another, check to see
#  if the second variable passed to the routine is "x", and complain
#  if so. 
#

    my(%votes) = ( "m" => 0, 
                   "c" => 0, 
                   "b" => 0, 
                   "f" => 0, 
                   "p" => 0,
                   "s" => 0,
		   "h" => 0,  
		   "e" => 0,
                   "y" => 0 );

   my($file) = $_[0];
   my(@data);
   my($i,$type);
   my($winner,$highest);

   @data=();


   open(IN,$file) || die "$file wouldn't open in read_align\n";

   $block_start=0; $block_end=0; 
   while(<IN>) { push(@data,$_); }
   close(IN);

   for($i=0; $i<=$#data; ++$i) {
	$_ = $data[$i];
        if(($i==0) && ($_  =~ / *[0-9]+ +[0-9]+ */)) { $votes{"y"}+=1000; }
	elsif(($_ =~ /^ *Name:/) || ($_ =~ /pileUp/) || ($_ =~ /MSF.*Check.*\.\./)) { $votes{"m"}++; last; }
	elsif($_ =~ /^CLUSTAL/) { $votes{"c"}+=10; last; }
	elsif($_ =~ /^>P1;/) { $votes{"p"}+=10;  last; }
   	elsif($_ =~ /HMMER/) { $votes{"h"}+=10; last; }
   	elsif(($_ =~ /^#=SQ/) || ($_ =~ /^#=RF/)) { $votes{"h"}++; }
	elsif($_ =~ /^>/) { $votes{"f"}++; $votes{"b"}++;  }
	elsif($_ =~ /^ *\* iteration [0-9]*/) { $votes{"b"}++; $block_start++; }
	elsif($_ =~ /^ *\*/) { $votes{"b"}++; $block_end++; }
	elsif(($_ =~ /^ID  /) || ($_ =~ /^CC  /) || ($_ =~ /^AC  /) || ($_ =~ /^SE  /)) { $votes{"s"}++; }
        elsif(($_ =~ /^HSSP .*HOMOLOGY DERIVED SECONDARY STRUCTURE OF PROTEINS/)) { $votes{"e"}+=1000; }

   }

   # Block and FASTA are quite hard to tell apart in a quick parse,
   #  This hack tries to fix this
   if(($votes{"f"} == $votes{"b"}) && ($votes{"f"}>0)) {
	if($block_start==0) { $votes{"f"}++; }
    	if($block_end==0) { $votes{"f"}++; }
   }

   $winner = "x";
   $highest = 0;
   foreach $type (keys %votes) {
#	print $type," ", $votes{$type},"\n";
	if($votes{$type}>$highest) { $winner = $type; $highest = $votes{$type}; }
   }
#   print "File is apparently of type $winner\n";
   $_[1] = $winner;
   return @data;
}
	
   


sub get_msf {

# Reads MSF.  Not too finickity, but you do need at least something like "MSF:.* .." at the
#  top of the file  to denote the start of the names, etc.
#
# Presumes that "." characters are always gaps and treats them as spaces internally
#  anything else will be left (e.g. "-" characters)

   my(@data) = @_;
   my($align,$i,$j,$k,$n);
   my($in_names,$in_align,$label,$name_count);
   my($seqdat);
   my($id);


   $align = {};

   $align->{ids}  = ();
   $align->{list} = ();
   $align->{nseq} = 0;

   $in_names = 0; $in_align = 0;
   $name_count = 0;
   for($n=0; $n<=$#data; ++$n) {
      $_ = $data[$n];
      chop;
      if((/Type.* Check:.* \.\.$/) || (/MSF:.* \.\./)) { $in_names = 1; $in_align=0; }
      if(/\/\//) { $in_align = 1; $in_names = 0; }
      if(($in_names==1) && (/Name: /)) {
        $_ =~ s/^ *//;
        $label = (split(/ +/))[1];
	if(defined($align->{ids}{$label})) {
		$label = $label . "-" . $name_count;
	}
	$align->{list}[$name_count]=$label;
	$align->{ids}{$label}{seq}   = "";
	($align->{ids}{$label}{prebumpf},$label,$align->{ids}{$label}{postbumpf}) = split(/ +/,$_,3);
	$name_count++;
      } elsif(($in_align==1) && (/[A-Za-z\.]/)) { 
        $_ =~ s/^ *//;
	($label,$seqdat) = split(/ +/,$_,2);
	$seqdat =~ s/ //g;
	$seqdat =~ s/\./ /g;
	if(defined($align->{ids}{$label})) {
	   $align->{ids}{$label}{seq} .= $seqdat;
        } else {
            print "Error: sequence in alignment not found in list: ",$label,"\n";
	    for($i=0; $i<$name_count; ++$i) {
		print "ID ",$i+1," = ",$align->{list}[$i],"\n";
	    }
	    die;
        }
     }
   }
#   foreach $id (keys %{$align->{ids}}) {
#		print "Read $id\n";
#   }
   $align->{nseq} = $name_count;
   $align->{alen} = length($align->{ids}{$align->{list}[0]}{seq});
   $align->{id_prs} = get_print_string($align);
   return $align;
}

sub write_msf {

# Writes MSF, duh.

	my($align) = $_[0];
	my($outfile) = $_[1];
	my($i,$j,$k);
	my($ts,$id);
        my($aa);
  	my($tmp) = \*OUT;

        open(OUT,">$outfile") || die "Error opening output file $outfile\n";



	print OUT "pileUp of: *.Gcg

Symbol comparison table: GenRunData:PileUpPep.Cmp  CompCheck: 1254

                   GapWeight: 6.0 
             GapLengthWeight: 0.2 

aconvert.msf MSF:  416  Type: P today Check: 100 ..
";

	for($i=0; $i<$align->{nseq}; ++$i) {
		$id = $align->{list}[$i];
		if(defined($align->{ids}{$id}{start})) {
			$align->{ids}{$id}{newid} = $id . "/" . ($align->{ids}{$id}{start}+1) . "-" . ($align->{ids}{$id}{end}+1);
			if(defined($align->{ids}{$id}{ranges})) {
				$align->{ids}{$id}{newid} .= $align->{ids}{$id}{ranges};
			}
		} else {
			$align->{ids}{$id}{newid} = $id;
		}
		if((!defined($align->{ids}{$id}{prebumpf})) || ($align->{ids}{$id}{prebumpf} !~ /Name:/)) {
			$align->{ids}{$id}{prebumpf} =  "Name: ";
		}
		if((!defined($align->{ids}{$id}{postbumpf})) || (($align->{ids}{$id}{postbumpf} !~ /Len:/) && ($align->{ids}{$id}{postbumpf} !~ /Check:/) && ($align->{ids}{$id}{postbumpf} !~ /Weight:/))) {
			$align->{ids}{$id}{postbumpf} =  " Len:  321 Check: 100  Weight: 1.00";
		}
		$ts = "  %s %-".$align->{id_prs}."s  %s\n";
		printf(OUT $ts,$align->{ids}{$id}{prebumpf},$align->{ids}{$id}{newid},$align->{ids}{$id}{postbumpf});
	}

	print OUT "\n\/\/\n";

	$i=0;
	while($i<$align->{alen}) {	
		for($j=0; $j<$align->{id_prs}+2; ++$j) { print OUT " "; }
		printf(OUT "%-3d",$i+1);
		for($j=0; $j<59; ++$j) {
			last if(($i+$j)>=($align->{alen}-1));
			print OUT " ";
		}
		printf(OUT "%3d\n",$i+$j+1);
	 	for($k=0; $k<$align->{nseq}; ++$k) {
	 		$id = $align->{list}[$k];
			$ts = "%-".$align->{id_prs}."s  ";
			printf(OUT $ts,$align->{ids}{$id}{newid});
			for($j=0; $j<60; ++$j) {
				last if(($i+$j)>=$align->{alen});
				$aa = substr($align->{ids}{$id}{seq},($i+$j),1);
				if($aa eq " ") { print OUT "."; }
				else { print OUT $aa; }
				if((($i+$j+1)%10)==0) { print OUT " "; }
			}
			printf(OUT "\n");
		}
		$i+=60;
		printf(OUT "\n");
	}
	close(OUT);
}

sub get_clustal {

# Presumes that "-" characters are always gaps and treats them as spaces internally
#  anything else will be left as is (e.g. "." characters)

   my(@data) = @_;
   my($name_count);
   my($label,$seq);
   my($align,$i,$j,$k,$n);

   $align = {};

   $align->{ids}  = ();
   $align->{nseq} = 0;

   $name_count = 0;
   for($n=0; $n<=$#data; ++$n) {
      $_ = $data[$n];
      chop;
      if((!/^CLUSTAL/) && (/^[^ ]/) && (length($_)>1)) {
	($label,$seq) = (split(/ +/));
	$seq =~ s/-/ /g;
	if(defined($align->{ids}{$label})) {
		$align->{ids}{$label}{seq} .= $seq;
	} else {
		$align->{list}[$name_count]=$label;
		$align->{ids}{$label}{seq} = $seq;
		$name_count++;
	}
      }
   }
   $align->{nseq} = $name_count;
   $align->{alen} = length($align->{ids}{$align->{list}[0]}{seq});
#   print "Alignment length is $align->{alen}\n";
   $align->{id_prs} = get_print_string($align);

   return $align;
}


sub write_clustal {

	$| = 1;
	my($align) = $_[0];
	my($i,$j,$k,$id,$ts);

	my($outfile) = $_[1];

  	my($tmp) = \*OUT;

        open(OUT,">$outfile") || die "Error opening output file $outfile\n";



	print OUT "CLUSTAL W(1.60) multiple sequence alignment\n\n";

	foreach $id (keys %{$align->{ids}}) {
		if(defined($align->{ids}{$id}{start})) {
			$align->{ids}{$id}{newid} = $id . "/" . ($align->{ids}{$id}{start}+1) . "-" . ($align->{ids}{$id}{end}+1);
			if(defined($align->{ids}{$id}{ranges})) {
				$align->{ids}{$id}{newid} .= $align->{ids}{$id}{ranges};
			}
		} else {
			$align->{ids}{$id}{newid} = $id;
		}
	}
	$i=0;
	while($i<$align->{alen}) {	
	 	for($k=0; $k<$align->{nseq}; ++$k) {
	 		$id = $align->{list}[$k];
			$ts = "%-".$align->{id_prs}."s  ";
#			print "id is $id\n";
			printf(OUT $ts,$align->{ids}{$id}{newid});
			for($j=0; $j<60; ++$j) {
				last if(($i+$j)>=$align->{alen});
				$aa = substr($align->{ids}{$id}{seq},($i+$j),1);
				if($aa eq " ") { print OUT "-"; }
				else { print OUT $aa; }
			}
			printf(OUT "\n");
		}
		$i+=60;
		printf(OUT "\n");
	}
	close(OUT);
}

sub get_block {
  
# Only " " characters are gaps.  "_" and "." characters are left as is (block format)
#
# Now saves pre and post bumpf around each position in the alignment
# align->{p}[0..N]{pre}
# align->{p}[0..N]{post}

   my(@data) = @_;
   my($align,$i,$j,$k,$n);

   $align = {};

   $align->{ids}  = ();
   $align->{list} = ();
   $align->{nseq} = 0;

   my($in_names) = 1; 
   my($in_stamp) = 0; 
   my($in_align) = 0;
   my($name_count) = 0;
   my($start_field);
   my($stamp_letter,$nstamp);
   my($p)=0;

   my($prebumpf,$postbumpf);
   my($label);

   for($n=0; $n<=$#data; ++$n) {
      $_ = $data[$n];
      chop;
      if(($in_names==1) && (/>/)) {
	  ($prebumpf,$label) = split(/>/);
	  $postbumpf = $label;
	  $label =~ s/ .*//;
	  $postbumpf =~ s/$label//;
	  if(defined($align->{ids}{$label})) {
		$label = $label . "-" . $name_count;
	  }
	  $align->{list}[$name_count]=$label;
	  $align->{ids}{$label}{seq}   = "";
	  $align->{ids}{$label}{prebumpf} = $prebumpf;
	  $align->{ids}{$label}{postbumpf} = $postbumpf;
  	  $name_count++;
      } elsif(($in_align==1) && (/[^ ]/)) {
	last if(/\*/);
	$align->{p}[$p]{pre}  = substr($_,0,$start_field);
	$align->{p}[$p]{post} = substr($_,$start_field+$name_count);
#	if($p==0) { print "Pre and post are $align->{p}[$p]{pre} and $align->{p}[$p]{post}\n"; } 
	$p++;
        $seqdat = substr($_,$start_field);

#	No - don't do the below.  Dashes mean something
#	$seqdat =~ s/[-\.]/ /g;


	for($i=0; $i<$name_count; ++$i) {
		$label = $align->{list}[$i];
		$align->{ids}{$label}{seq} .= substr($seqdat,$i,1);
#		if(substr($seqdat,$i,1) eq "-") { print "Dash encountered in sequeunce $label\n"; }
	}
     } elsif(/^\#/) {
##T   -- '1' = equivalenced residues
	$stamp_letter = substr($_,1,1);
	$align->{stamp}{$stamp_letter} = $nstamp;
	$nstamp++;
#	print " FOUND STAMP FIELD $_ $stamp_letter\n";
     } elsif(/\*/) {
	  $in_names = 0;
	  $in_align = 1;
	  $start_field=0; 
	  while(substr($_,$start_field,1) ne "*") { $start_field++; }
     }
   }
#   foreach $id (keys %{$align->{ids}}) {
#		print "Read $id ",$align->{ids}{$id}{seq},"\n";
#   }
   $align->{nseq} = $name_count;
   $align->{alen} = length($align->{ids}{$align->{list}[0]}{seq});
#   print "Details from get_block: nseq = ",$align->{nseq}," alen = ",$align->{alen},"\n";
#   die;
   $align->{id_prs} = get_print_string($align);
#   if(defined($align->{p}[0]{pre})) {
#	print "PRE is defined\n";
#   }
   return $align;
}


sub write_block {

	# Now writes pre and post bumpf as discussed above

	my($align) = $_[0];
	my($i,$j,$k,$id);
	my($outfile) = $_[1];
  	my($tmp) = \*OUT;

        open(OUT,">$outfile") || die "Error opening output file $outfile\n";


	foreach $id (keys %{$align->{ids}}) {
		if(defined($align->{ids}{$id}{start})) {
			$align->{ids}{$id}{newid} = $id . "/" . ($align->{ids}{$id}{start}+1) . "-" . ($align->{ids}{$id}{end}+1);
			if(defined($align->{ids}{$id}{ranges})) {
				$align->{ids}{$id}{newid} .= $align->{ids}{$id}{ranges};
			}
		} else {
			$align->{ids}{$id}{newid} = $id;
		}
	}
 	for($k=0; $k<$align->{nseq}; ++$k) {
 		$id = $align->{list}[$k];
		printf( OUT ">%-30s ",$align->{ids}{$id}{newid});
		if(defined($align->{ids}{$id}{prebumpf})) { print OUT  $align->{ids}{$id}{prebumpf}," "; }
		if(defined($align->{ids}{$id}{postbumpf})) { print OUT  $align->{ids}{$id}{postbumpf}," "; }
		print OUT  "\n";
#		$align->{ids}{$id}{seq} =~ s/[-\.]/ /g;
	}
	if(defined($align->{p}[0]{pre})) {
		for($i=0; $i<length($align->{p}[0]{pre}); ++$i) { print OUT " "; }
	}
	print OUT  "* iteration 1\n";
	for($i=0; $i<$align->{alen}; ++$i) {
#		if($i==0) { print "Pre and post are $align->{p}[$i]{pre} and $align->{p}[$i]{post}\n"; }
		if(defined($align->{p}[$i]{pre})) { print OUT $align->{p}[$i]{pre}; }
 		for($k=0; $k<$align->{nseq}; ++$k) {
 			$id = $align->{list}[$k];
			print OUT  substr($align->{ids}{$id}{seq},$i,1);
		}
		if(defined($align->{p}[$i]{post})) { print OUT $align->{p}[$i]{post}; }
		print OUT  "\n";
	}
	if(defined($align->{p}[0]{pre})) {
		for($i=0; $i<length($align->{p}[0]{pre}); ++$i) { print OUT " "; }
	}
	print OUT  "*\n";
	close(OUT);
	
}

sub get_afasta {


# Get aligned FASTA data (i.e. afasta = aligned fasta)
#
# Assumes that only "-" characters are gaps (treats as spaces internally)
#  (actually, of course this means spaces will be gaps, or maybe they should be ignored
#  strange dilemma that I have got myself into.  For the minute they will just be kept as is
#  (spaces are gaps as well).  This could lead to some problems.

   my(@data) = @_;
   my($align,$i,$j,$k,$n);

   $align = {};

   $align->{ids}  = ();
   $align->{list} = ();
   $align->{nseq} = 0;

   my($name_count) = 0;
   for($n=0; $n<=$#data; ++$n) {
	$_ = $data[$n];
	chop;
	if(/^>/) {
		$label = substr($_,1);
		$label =~ s/ .*//;
		$postbumpf = $_; $postbumpf =~ s/>[^ ]* //;
	  	if(defined($align->{ids}{$label})) {
			$label = $label . "-" . $name_count;
	  	}
	  	$align->{list}[$name_count]=$label;
		$align->{ids}{$label}{seq} = "";
		$align->{ids}{$label}{postbumpf} = $postbumpf;
		$name_count ++;
	} else {
		if(defined($align->{ids}{$label})) {
			$_ =~ s/-/ /g;
			$align->{ids}{$label}{seq} .= $_;
			#print "Adding $_ to $label\n";
		}
	}
   }
   $align->{nseq} = $name_count;
   $align->{alen} = length($align->{ids}{$align->{list}[0]}{seq});
   $align->{id_prs} = get_print_string($align);
   return $align;
}

sub write_afasta {

	my($align) = $_[0];
	my($nogap);
	my($outfile) = $_[1];
	if(defined($_[2])) { $nogap = $_[2]; }
	else { $nogap = 0; }
	if(defined($_[3])) { $separate = $_[3]; }
	my($i,$j,$k,$id);
  	my($tmp) = \*OUT;

        open(OUT,">$outfile") || die "Error opening output file $outfile\n";


	foreach $id (keys %{$align->{ids}}) {
		if(defined($align->{ids}{$id}{start})) {
			$align->{ids}{$id}{newid} = $id . "/" . ($align->{ids}{$id}{start}+1) . "-" . ($align->{ids}{$id}{end}+1);
			if(defined($align->{ids}{$id}{ranges})) {
				$align->{ids}{$id}{newid} .= $align->{ids}{$id}{ranges};
			}
		} else {
			$align->{ids}{$id}{newid} = $id;
		}
	}

	for($i=0; $i<$align->{nseq}; ++$i) {
	   $id = $align->{list}[$i];
           if($separate==1) { 
                $filename = $id . ".fasta";
                $filename =~ s/\//-/g;
     	  	print "Filename is $filename\n";
 		$seq = $align->{ids}{$id}{seq};
		if($nogap==1) { $seq =~ s/ //g; }
     		write_one_fasta($filename,$id,$align->{ids}{$id}{prebumpf},$seq);
           } else {
		print OUT ">",$id," ";
		if(defined($align->{ids}{$id}{prebumpf})) { print OUT $align->{ids}{$id}{prebumpf}," "; }
		if(defined($align->{ids}{$id}{postbumpf})) { print OUT $align->{ids}{$id}{postbumpf}," "; }
		print OUT "\n";
		$k=0;
		for($j=0; $j<length($align->{ids}{$id}{seq}); ++$j) {
			if(($nogap==0) || (substr($align->{ids}{$id}{seq},$j,1) ne " ")) { 
				$aa = substr($align->{ids}{$id}{seq},($j),1);
				if($aa eq " ") { print OUT "-"; }
				else { print OUT $aa; }
			#	print OUT substr($align->{ids}{$id}{seq},$j,1);
				if((($k+1)%60)==0) { print OUT "\n"; }
				$k++;
			}
		}
		print OUT "\n";
           }
	}
	close(OUT);
}
		
sub get_apir {

# Get aligned PIR data (i.e. afasta = aligned PIR)
#
# Assumes that only "-" characters are gaps (treats as spaces internally)
#  (actually, of course this means spaces will be gaps, or maybe they should be ignored
#  strange dilemma that I have got myself into.  For the minute they will just be kept as is
#  (spaces are gaps as well).  This could lead to some problems.

   my(@data) = @_;
   my($align,$i,$j,$k,$n);
   my($label);

   $align = {};

   $align->{ids}  = ();
   $align->{list} = ();
   $align->{nseq} = 0;

   my($name_count) = 0;
   my($title_next) = 0;
   for($n=0; $n<=$#data; ++$n) {
        $_ = $data[$n];
	chop;
	if(/^>P1;/) {
		$label = substr($_,4);
		$label =~ s/ .*//;
#		print "Found label $label\n";
		$postbumpf = $_; $postbumpf =~ s/>[^ ]* //;
	  	if(defined($align->{ids}{$label})) {
			$label = $label . "-" . $name_count;
	  	}
	  	$align->{list}[$name_count]=$label;
		$align->{ids}{$label}{seq} = "";
		$name_count ++;
		$in_seq=1;
		$title_next = 1;
	} elsif(($in_seq==1) && ($title_next==0)) {
		if(defined($align->{ids}{$label})) {
			$_ =~ s/-/ /g;
			$align->{ids}{$label}{seq} .= $_;
		}
		if(/\*/) {
			$in_seq=0;
			$align->{ids}{$label}{seq} =~ s/\*//;
		}
	} elsif($title_next==1) {
		$align->{ids}{$label}{postbumpf} = $_;
		$title_next = 0;
	}
   }
   $align->{nseq} = $name_count;
   $align->{alen} = length($align->{ids}{$align->{list}[0]}{seq});
   $align->{id_prs} = get_print_string($align);
   return $align;
}

sub write_apir {

	my($i,$j,$k);
	my($align) = $_[0];
	my($nogap);
	my($id);
  	my($tmp) = \*OUT;

	if(defined($_[2])) { $nogap = $_[2]; }
	else { $nogap = 0; }


	my($outfile) = $_[1];
	open(OUT,">$outfile") || die "Error opening output file $outfile\n";

	foreach $id (keys %{$align->{ids}}) {
		if(defined($align->{ids}{$id}{start})) {
			$align->{ids}{$id}{newid} = $id . "/" . ($align->{ids}{$id}{start}+1) . "-" . ($align->{ids}{$id}{end}+1);
			if(defined($align->{ids}{$id}{ranges})) {
				$align->{ids}{$id}{newid} .= $align->{ids}{$id}{ranges};
			}
		} else {
			$align->{ids}{$id}{newid} = $id;
		}
	}
	for($i=0; $i<$align->{nseq}; ++$i) {
		$id = $align->{list}[$i];
		print OUT ">P1;",$id,"\n";
                
                print OUT "structureX:",$id,":1: :";
                $leni=get_len($align,$id);
                printf OUT " %d ",$leni;
                print OUT " : : ";
		if(defined($align->{ids}{$id}{prebumpf})) { print OUT $align->{ids}{$id}{prebumpf}," "; }
                print OUT " : : -1.00 :-1.00";
		print OUT "\n";
		for($j=0; $j<length($align->{ids}{$id}{seq}); ++$j) {
                        $_=substr($align->{ids}{$id}{seq},$j,1);
                        tr/[a-z]/[A-Z]/;
                        if ($_ eq " ") {$_="-";}
			print OUT $_;
			if((($j+1)%60)==0) { print OUT "\n"; }
		}
		print OUT "*\n";
	}

	close(OUT);
}

sub get_amenta {

# Get aligned MENTA data (i.e. afasta = aligned MENTA)
#
# Assumes that only "-" characters are gaps (treats as spaces internally)
#  (actually, of course this means spaces will be gaps, or maybe they should be ignored
#  strange dilemma that I have got myself into.  For the minute they will just be kept as is
#  (spaces are gaps as well).  This could lead to some problems.

   my(@data) = @_;
   my($align,$i,$j,$k,$n);
   my($label);

   $align = {};

   $align->{ids}  = ();
   $align->{list} = ();
   $align->{nseq} = 0;

   my($name_count) = 0;
   my($title_next) = 0;
   for($n=0; $n<=$#data; ++$n) {
        $_ = $data[$n];
	chop;
	if(/^>/) {
		$label = substr($_,4);
		$label =~ s/ .*//;
#		print "Found label $label\n";
		$postbumpf = $_; $postbumpf =~ s/>[^ ]* //;
	  	if(defined($align->{ids}{$label})) {
			$label = $label . "-" . $name_count;
	  	}
	  	$align->{list}[$name_count]=$label;
		$align->{ids}{$label}{seq} = "";
		$name_count ++;
		$in_seq=1;
		$title_next = 0;
	} elsif(($in_seq==1) && ($title_next==0)) {
		if(defined($align->{ids}{$label})) {
			$_ =~ s/-/ /g;
			$align->{ids}{$label}{seq} .= $_;
		}
		if(/\*/) {
			$in_seq=0;
			$align->{ids}{$label}{seq} =~ s/\*//;
		}
	} elsif($title_next==1) {
		$align->{ids}{$label}{postbumpf} = $_;
		$title_next = 0;
	}
   }
   $align->{nseq} = $name_count;
   $align->{alen} = length($align->{ids}{$align->{list}[0]}{seq});
   $align->{id_prs} = get_print_string($align);
   return $align;
}

sub write_amenta {

	my($i,$j,$k);
	my($align) = $_[0];
	my($nogap);
	my($id);
  	my($tmp) = \*OUT;

	if(defined($_[2])) { $nogap = $_[2]; }
	else { $nogap = 0; }


	my($outfile) = $_[1];
	open(OUT,">$outfile") || die "Error opening output file $outfile\n";

	foreach $id (keys %{$align->{ids}}) {
		if(defined($align->{ids}{$id}{start})) {
			$align->{ids}{$id}{newid} = $id . "/" . ($align->{ids}{$id}{start}+1) . "-" . ($align->{ids}{$id}{end}+1);
			if(defined($align->{ids}{$id}{ranges})) {
				$align->{ids}{$id}{newid} .= $align->{ids}{$id}{ranges};
			}
		} else {
			$align->{ids}{$id}{newid} = $id;
		}
	}
	for($i=0; $i<$align->{nseq}; ++$i) {
		$id = $align->{list}[$i];
		print OUT ">M1 ",$id,"\n";
		for($j=0; $j<length($align->{ids}{$id}{seq}); ++$j) {
                        $_=substr($align->{ids}{$id}{seq},$j,1);
                        tr/[a-z]/[A-Z]/;
                        if ($_ eq " ") {$_="-";}
			print OUT $_;
			if((($j+1)%60)==0) { print OUT "\n"; }
		}
		print OUT "\n";
	}

	close(OUT);
}

sub get_pfam {

   my(@data) = @_;
   my($align,$i,$j,$k,$n);

   $align = {};

   $align->{ids}  = ();
   $align->{list} = ();
   $align->{nseq} = 0;

   my($name_count) = 0;
   my($label,$seq);
   for($n=0; $n<=$#data; ++$n) {
      $_ = $data[$n];
      chop;
      if((/^[^ ]/) && (length($_)>1) && (!/^[A-Z][A-Z]   /) && (!/\/\//)) {
        ($label,$seq,$accession) = (split(/ +/));
	$seq =~ s/[\.\-]/ /g;
#	print "Calling $_ sequence data.  Length so far is ",$align->{alen},"\n";
        if(defined($align->{ids}{$label})) {
                $align->{ids}{$label}{seq} .= $seq;
        } else {
                $align->{list}[$name_count]=$label;
                $align->{ids}{$label}{seq} = $seq;
		$align->{ids}{$label}{prebumpf} = $accession;
                $name_count++;
        }
      }
   }
   $align->{nseq} = $name_count;
   $align->{alen} = length($align->{ids}{$align->{list}[0]}{seq});

   $align->{id_prs} = get_print_string($align);
   return $align;
}
sub get_hmmalign {

# Currently this just ignores the "#" lines and the other junk at the start


   my(@data) = @_;
   my($align,$i,$j,$k,$n);
   my($name_count) = 0;
   my($label,$seq);

   $align = {};

   $align->{ids}  = ();
   $align->{list} = ();
   $align->{nseq} = 0;


   for($n=0; $n<=$#data; ++$n) {
      $_ = $data[$n];
      @f = split(/[ \t]+/,$data[$n]);
      chop;
      if((!/^#/) && (!/HMM /) && (!/HMMER /) && (!/- - - -/) && (!/Sequence/) && (!/Copyright/) && (/^[^ ]/) && (length($_)>1) && (!/^[A-Z][A-Z]   /) && (!/\/\//)  && (!/Freely/) ) {
        ($label,$seq) = (split(/ +/));
	$seq =~ s/[\.\-]/ /g;
#	print "Calling $_ sequence data.  Length so far is ",$align->{alen},"\n";
        if(defined($align->{ids}{$label})) {
                $align->{ids}{$label}{seq} .= $seq;
        } else {
                $align->{list}[$name_count]=$label;
                $align->{ids}{$label}{seq} = $seq;
                $name_count++;
        }
      }
   }
   $align->{nseq} = $name_count;
   $align->{alen} = length($align->{ids}{$align->{list}[0]}{seq});

   $align->{id_prs} = get_print_string($align);
   return $align;
}

sub merge_align { # arguments: align1 align2 (linker)

	# Merges alignments

	my($i,$j,$k);
	my($align1) = $_[0];
	my($align2) = $_[1];

	my($align3) = {};
	my($linker);

	my($debug) = 0;
 	my($id,$id1,$id2);


	if($debug==1) { print "Trying to link ",$align1->{nseq}," len=",$align1->{alen}," and ",$align2->{nseq}," len=",$align2->{alen},"\n"; }

	if(defined($_[2])) { # identifier to link with has been specified
		$linker = $_[2];
		if((!defined($align1->{ids}{$linker})) || (!defined($align2->{ids}{$linker}))) {
			die "Error: defined linker $linker missing from one or both of the alignments\n";
		}
	} else {
		# Just use the first identifier found in common
		if($debug==1) { print "Sizes are $align1->{nseq} and $align2->{nseq}\n"; } 
		foreach $id1 (keys %{$align1->{ids}}) {
			foreach $id2 (keys %{$align2->{ids}}) {
				if($id1 eq $id2) { $linker = $id1; last; }
			}
			last if(defined($linker));
		}
		if(!defined($linker)) {
			die "Error: no link found\n";
		}
	}
	if($debug==1) { 
		print "\nWill link on $linker\n"; 
		print "seq1: |",$align1->{ids}{$linker}{seq},"|\nseq2: |",$align2->{ids}{$linker}{seq},"|\n\n";
	}
	$align3->{nseq} = 0;
	$align3->{list} = ();
     	# Create a blank entry for each entry in the old alignments
	if($debug==1) { print "Adding entries corresponding to the first alignment total is ",$align1->{nseq},"\n"; }
	for($i=0; $i<$align1->{nseq}; ++$i) {
		$id = $align1->{list}[$i];
		$align1->{ids}{$id}{seq} =~ s/[-\.]/-/g;
		$align3->{ids}{$id}{seq} = "";
		$align3->{list}[$align3->{nseq}] = $id;
		if(defined($align1->{ids}{$id}{prebumpf})) {
			$align3->{ids}{$id}{prebumpf} = $align1->{ids}{$id}{prebumpf};
		}
		if(defined($align1->{ids}{$id}{postbumpf})) {
			$align3->{ids}{$id}{postbumpf} = $align1->{ids}{$id}{postbumpf};
		}
		$align3->{nseq}++;
#		if($debug==1) { print "Added $id\n"; }
	}
#	if($debug==1) { print "Adding entries corresponding to the second alignment\n"; }
	for($i=0; $i<$align2->{nseq}; ++$i) {
		$id = $align2->{list}[$i];
		$align2->{ids}{$id}{seq} =~ s/[-\.]/-/g;
		if((!defined($align3->{ids}{$id})) && ($id ne $linker)) {
			$align3->{list}[$align3->{nseq}] = $id;
			$align3->{ids}{$id}{seq} = "";
			if(defined($align2->{ids}{$id}{prebumpf})) {
				$align3->{ids}{$id}{prebumpf} = $align2->{ids}{$id}{prebumpf};
			}
			if(defined($align2->{ids}{$id}{postbumpf})) {
				$align3->{ids}{$id}{postbumpf} = $align2->{ids}{$id}{postbumpf};
			}
			$align3->{nseq}++;
		}
#		if($debug==1) { print "Added $id\n"; }
	}

	# Now just go through the alignments and index them to the common sequence
  	#  bomb out if there is a mis-match between the sequences

	$pointer1 = 0; $counter1 = 0;
	$pointer2 = 0; $counter2 = 0;
	$align3->{alen} = 0;
#	print $align1->{alen}," vs ",length($align1->{ids}{$align1->{list}[0]}{seq})," ";
#	print $align2->{alen}," vs ",length($align2->{ids}{$align2->{list}[0]}{seq}),"\n";
	if($debug==1) { print "Doing the merging...\n"; } 
	while(($counter1<$align1->{alen}) || ($counter2<$align2->{alen})) {


#		print "Current values: $counter1 / $align1->{alen} $counter2 / $align2->{alen} ";
#		print substr($align1->{ids}{$linker}{seq},$counter1,1)," ",substr($align2->{ids}{$linker}{seq},$counter2,1),"\n";
#		print "Here1\n";
		while(substr($align1->{ids}{$linker}{seq},$counter1,1) eq " ") { # if blank, copy the rest and put blanks in the other alignment
			last if(($counter1>=$align1->{alen}) && ($counter2>=$align2->{alen}));
#			printf("1 %4d %4d %4d %4d - \n",$pointer1,$counter1,$pointer2,$counter2);
			foreach $id1 (keys %{$align1->{ids}}) {
				$align3->{ids}{$id1}{seq} .= substr($align1->{ids}{$id1}{seq},$counter1,1);
#				print  substr($align1->{ids}{$id1}{seq},$counter1,1);
			}
			foreach $id2 (keys %{$align2->{ids}}) {
				if(!defined($align1->{ids}{$id2})) {
					$align3->{ids}{$id2}{seq} .=  " ";
#					print  "-";
				}
			}
			if(defined($align1->{p}[$counter1]{pre})) {
				$align3->{p}[$align3->{alen}]{pre} = $align1->{p}[$counter1]{pre};
			}
			if(defined($align1->{p}[$counter1]{post})) {
				$align3->{p}[$align3->{alen}]{post} = $align1->{p}[$counter1]{post};
			}
			$counter1++;
			$align3->{alen}++;
#			print "\n";
		}
#		print "Here2\n";
		last if(($counter1>=$align1->{alen}) && ($counter2>=$align2->{alen}));
		while(substr($align2->{ids}{$linker}{seq},$counter2,1) eq " ") { # if blank, copy the rest and put blanks in the other alignment
			last if(($counter1>=$align1->{alen}) && ($counter2>=$align2->{alen}));
#			printf("2 %4d %4d %4d %4d - \n",$pointer1,$counter1,$pointer2,$counter2);
			foreach $id1 (keys %{$align1->{ids}}) {
				$align3->{ids}{$id1}{seq} .= " ";
#				print  ".";
			}
			foreach $id2 (keys %{$align2->{ids}}) {
				if(!defined($align1->{ids}{$id2})) {
					$align3->{ids}{$id2}{seq} .= substr($align2->{ids}{$id2}{seq},$counter2,1);
#					print substr($align2->{ids}{$id2}{seq},$counter2,1);
				}
			}
			if(defined($align2->{p}[$counter2]{pre})) {
				$align3->{p}[$align3->{alen}]{pre} = $align2->{p}[$counter2]{pre};
			}
			if(defined($align2->{p}[$counter2]{post})) {
				$align3->{p}[$align3->{alen}]{post} = $align2->{p}[$counter2]{post};
			}
			$counter2++;
			$align3->{alen}++;
#			print "\n";
		}
		last if(($counter1>=$align1->{alen}) && ($counter2>=$align2->{alen}));
		if(uc(substr($align1->{ids}{$linker}{seq},$counter1,1)) ne uc(substr($align2->{ids}{$linker}{seq},$counter2,1))) {
			print STDERR "Error: linking with sequence ",$linker," mismatch found at alignment positions $counter1 and $counter2: ",
				substr($align1->{ids}{$linker}{seq},$counter1,1),
				" versus ",
				substr($align2->{ids}{$linker}{seq},$counter2,1),
				"\nSeq1 = $align1->{ids}{$linker}{seq}\n",
				"Seq2 = $align2->{ids}{$linker}{seq}\n";
			print STDERR "       ";
			for($k=0; $k<$counter2; ++$k) {
				print STDERR ".";
			}
			print STDERR "\n";
			die;
		}
		while((uc(substr($align1->{ids}{$linker}{seq},$counter1,1)) eq uc(substr($align2->{ids}{$linker}{seq},$counter2,1))) &&
		      (substr($align1->{ids}{$linker}{seq},$counter1,1) ne " ")) {
			last if(($counter1>=$align1->{alen}) && ($counter2>=$align2->{alen}));
#			printf("3 %4d %4d %4d %4d - \n",$pointer1,$counter1,$pointer2,$counter2);
			foreach $id1 (keys %{$align1->{ids}}) {
                                $align3->{ids}{$id1}{seq} .= substr($align1->{ids}{$id1}{seq},$counter1,1);
#				print substr($align1->{ids}{$id1}{seq},$counter1,1);
                        }
			foreach $id2 (keys %{$align2->{ids}}) {
                                if(!defined($align1->{ids}{$id2})) {
                                        $align3->{ids}{$id2}{seq} .= substr($align2->{ids}{$id2}{seq},$counter2,1);
#					print substr($align2->{ids}{$id2}{seq},$counter2,1);
                                }
                        }
			$counter1++; $counter2++; $pointer1++; $pointer2++;
			if(defined($align1->{p}[$counter1]{pre})) {
				$align3->{p}[$align3->{alen}]{pre} = $align1->{p}[$counter1]{pre};
			}
			if(defined($align1->{p}[$counter1]{post})) {
				$align3->{p}[$align3->{alen}]{post} = $align1->{p}[$counter1]{post};
			}
			$align3->{alen}++;
#			print "\n";
		}

	}
	

        $align3->{id_prs} = get_print_string($align3);
	return $align3;
}

sub do_clustalw_align {
	
	# Given a reference containing sequences, this runs
	#  Clustalw and returns the result
	# Needs tayloring to system as there is an executable below

	$clustalw = "/apps/ATG/bin/clustalw";

	my($seqs) = $_[0];
	my($command,$pirfile);
	my($align);
	my(@data);

	$prefix = "/tmp/do_clustalw_align." . $$ ;

	$pirfile = $prefix . ".pir";
	$alnfile = $prefix . ".aln";

	write_apir($seqs,$pirfile);

	system($command);
	$command = $clustalw . " " . $pirfile . " 1> /dev/null 2> /dev/null";
	system($command);

	@data = read_align($alnfile);
	$align = get_clustal(@data);
	@data=();

#	$command = "/bin/rm " . $prefix . "*";
#	system($command);

	return $align;
}

sub extract_align {

	my($i,$j);	
	my($align) = $_[0];
	my($strings) = $_[1];
	my($invert) = 0;
	if(defined($_[2])) { $invert = $_[2]; }
	my(@extract);

	my($align2) = {};
	my($nseq_new);
	my($new_strings);
	my($ignore);

	# Given a reference containing sequences and an array of
	#  strings, this creates a new reference containing
	#  only those that match one or more of the strings

	$name_count = 0;

	$align2->{ids}  = ();
	$align2->{list} = ();
	$align2->{nseq} = 0;
	$align2->{alen} = 0;
	$align2->{p} = $align->{p};

	if($invert == 1 ) { # This is like a grep -v 
	    @extract = (split(/ +/,$strings));
	    $new_strings = "";
            for($i=0; $i<$align->{nseq}; ++$i) {
                $id = $align->{list}[$i];
		$ignore = 0;
		for($j=0; $j<=$#extract; ++$j) {
		   if($id =~ /$extract[$j]/) {
			$ignore = 1;
		   }
		}
		if($ignore==0) {
		    $new_strings .= $id . " ";
	 	}
	    }
	    $strings = $new_strings;
	}
	@extract = (split(/ +/,$strings));
#	print "Extract string is $extract\n";

	for($j=0; $j<=$#extract; ++$j) {
       	   for($i=0; $i<$align->{nseq}; ++$i) {
                $id = $align->{list}[$i];
			if($id =~ /$extract[$j]/) {
				$align2->{ids}{$id}{seq} = $align->{ids}{$id}{seq};
				$align2->{list}[$name_count] = $id;
				$align2->{nseq}++;
				$name_count++;
				if(defined($align->{ids}{$id}{prebumpf})) {
					$align2->{ids}{$id}{prebumpf} = $align->{ids}{$id}{prebumpf};
				}
				if(defined($align->{ids}{$id}{postbumpf})) {
					$align2->{ids}{$id}{postbumpf} = $align->{ids}{$id}{postbumpf};
				}
				last;
			}
	   }
	}
	$align2->{alen} = length($align->{ids}{$align->{list}[0]}{seq});

   	$align2->{id_prs} = get_print_string($align2);
	return $align2;
}
sub get_phylip {

   my(@data) = @_;
   my($align,$i,$j,$k,$n);

   $align = {};

   $align->{ids}  = ();
   $align->{list} = ();
   $align->{nseq} = 0;

   my($in_first_block) = 1;
   my($name_count) = 0;
   my($label,$seq);
   $name_count = 0;
   for($n=0; $n<=$#data; ++$n) {
      $_ = $data[$n];
      chop;
      if($n==0) { 
         $_=~ s/^ *//;
         ($nseq,$alen) = (split(/ +/));
#	 print "Nseq $nseq, length $alen\n";
      } else {
        if(((/[^ ]/) && (!/[^ ]/)) || (length($_)<1)) { # blank line
           $in_first_block = 0;
	   if($name_count>0) { 
		# This permits extra blank lines at end of file
		$align->{nseq} = $name_count; 
	   } 
           $name_count = 0;
        } elsif($in_first_block) {
           s/^ *//;
           ($id,$seq) = (split(/ +/,$_,2));
#	   print "Read $id, $seq\n";
	   $seq =~ s/ //g;
           $align->{list}[$name_count] = $id;
           $align->{ids}{$id}{seq} = $seq;
	   $name_count++;
        } else {
	   $seq = $_;
	   $seq =~ s/ //g;
           $id = $align->{list}[$name_count];
           if(!defined($align->{ids}{$id})) { 
               die "Error in get_phylip reading and id called $id\n";
 	   } 
           $align->{ids}{$id}{seq} .= $seq;
	   $name_count++;
        }
      }
   }
   $align->{alen} = length($align->{ids}{$align->{list}[0]}{seq});
#   print "Align nseq = $align->{nseq} length $align->{alen}\n";

   $align->{id_prs} = get_print_string($align);
   return $align;
}
sub write_phylip {

	my($align) = $_[0];
	my($i,$j,$k,$id);

	my($outfile) = $_[1];
        open(OUT,">$outfile") || die "Error opening output file $outfile\n";



	printf(OUT "%4d %4d\n\n",$align->{nseq},$align->{alen});

	foreach $id (keys %{$align->{ids}}) {
		if(defined($align->{ids}{$id}{start})) {
			$align->{ids}{$id}{newid} = $id . "/" . ($align->{ids}{$id}{start}+1) . "-" . ($align->{ids}{$id}{end}+1);
			if(defined($align->{ids}{$id}{ranges})) {
				$align->{ids}{$id}{newid} .= $align->{ids}{$id}{ranges};
			}
		} else {
			$align->{ids}{$id}{newid} = $id;
		}
	}
	$i=0;
	while($i<$align->{alen}) {	
	 	for($k=0; $k<$align->{nseq}; ++$k) {
	 		$id = $align->{list}[$k];
			if($i<60)  { printf(OUT "%-30s ",$align->{ids}{$id}{newid}); } 
			for($j=0; $j<60; ++$j) {
				last if(($i+$j)>=$align->{alen});
				print OUT substr($align->{ids}{$id}{seq},($i+$j),1);
			}
			printf(OUT "\n");
		}
		$i+=60;
		printf(OUT "\n");
	}
	close(OUT);
}

sub get_print_string {
	my($align) = $_[0];
	my($max_len);
	my($i);

	$max_len = 0;
	for($i=0; $i<$align->{nseq}; ++$i) {
		$this_len = length($align->{list}[$i]);
		if($this_len > $max_len) { $max_len = $this_len; }
	}
	return $max_len;
}

sub read_exchange_matrix {

##  Matrix made by matblas from blosum62.iij
##  * column uses minimum score
##  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
##  Blocks Database = /data/blocks_5.0/blocks.dat
##  Cluster Percentage: >= 62
##  Entropy =   0.6979, Expected =  -0.5209
#   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
#A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
#R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 

	my($i,$j,$n);
	my($f) = $_[0];
	my($mat) = {};
	
	open(M,$f) || die "Error opening file $f\n";
#	print "Reading file from $f\n";
	$n = 0;
	while(<M>) {
		chop;
		if(!/^\#/) {
		   if($n==0) {
			$_ =~ s/ //g;
			$mat->{aa}  = $_;
			$mat->{naa} = length($_);
		   } else {
			$x = substr($_,0,1);
			for($i=0; $i<$mat->{naa}; ++$i) {
				$y = substr($mat->{aa},$i,1);
				$mat->{"m"}{$x}{$y} = substr($_,1+3*$i,3);
#				printf("%3s -> %3s = %4d \n",$x,$y,$mat->{"m"}{$x}{$y});
			}
		  }
		  $n++;
	        }
	}
	return $mat;
}

sub prettify {


	my($hydrophobic)= "ACFILMVWYHTacfilmvwyht";
	my($polar)      = "DEGHKNPQRSTYdeghknpqrsty";
	my($small)      = "ACGSTPacgstp";
	my($ss_string)  = "HGEBSTICXhgebsticx-.";
	my($helix_def)  = "HGhg";
	my($strand_def) = "Ee";
	
	my($align) = $_[0];
	my($outfile) = $_[1];
	my($set) = $_[2];
	my($min_frac) = $_[3];
#	print "passed arguments are @_\n";
	$stamp_char = $_[4];  
	$stamp_cut  = $_[5];  
	$stamp_win  = $_[6];  
#	if(defined($stamp_char)) { print "IN PRETTIFY stamp char is $stamp_char\n"; }

	my($nseqs);
	my($i,$j,$k);
	my($aa,$ss);

        open(OUT,">>$outfile") || die "Error opening output file $outfile\n";

	# Crude:
	# Cysteines, Conserved small (ATGPS), Conserved polar (RKDEQNSTYH), Conserved hydrophobic (ILMVFWY)

	$nseqs=0;

	$start = -1;
	$end = -1;
#	print "Sequences to use are $set\n";
	for($i=0; $i<$align->{nseq}; ++$i) {
		$id=$align->{list}[$i];
		if($set =~ /$id/) {
			$use{$id} = 1;
#			print "Using $id\n";
			if(($id!~ /_dssp/) && ($id !~ /Sec_/) && ($id !~ /^space/)) { 
				if($start==-1) { $start = $i; }
				$end = $i;
				$nseqs++;
			}
		} else {
			$use{$id} = 0;
		}
	}
	$min_aa = $min_frac * $nseqs; 
	for($i=0; $i<$align->{alen}; ++$i) {
		$n_cysteine     = 0;
		$n_hydrophobic  = 0;
		$n_polar        = 0;
		$n_small        = 0;
		$n_id           = 0;
		for($j=0; $j<$align->{nseq}; ++$j) {
			$id = $align->{list}[$j];
			$aa = uc(substr($align->{ids}{$id}{seq},$i,1));
			if($j==0) { $first_aa = $aa; }
			if(($aa =~ /[A-Za-z]/) && ($use{$id}==1) && 
				($id !~ /_dssp/) && ($id !~ /Sec_/)) {
				if($aa eq "C") { $n_cysteine++; }
				if($small =~ /$aa/) { $n_small++; }
				if($hydrophobic =~ /$aa/) { $n_hydrophobic++; }
				if($polar =~ /$aa/) { $n_polar++; }
				if($aa eq $first_aa) { $n_id++; }
			}
		}	
#		print "Pos ",$i," cys: ",$n_cysteine," small: ",$n_small," hydrophobic: ",$n_hydrophobic," polar ",$n_polar," min: ",$min_aa,"\n";
		if($n_id==$nseqs) {
		        print OUT "BOX_REGION   ",$i+1," ",$start+1," ",$i+1," ",$end+1," 0.0 \n";
		}
		if($n_cysteine>=$min_aa) {
			$r = get_ranges_prop($align,$i,$start,$end,"C");
			for($j=0; $j<=$#{$r->{range}}; ++$j) {
			   print OUT "COLOUR_REGION   ",$i+1," ",$r->{range}[$j]{st}," ",$i+1," ",$r->{range}[$j]{en}," 2 \n"; 
			   print OUT "INVERSE_CHARS C ",$i+1," ",$r->{range}[$j]{st}," ",$i+1," ",$r->{range}[$j]{en},"\n"; 
			}
		} elsif($n_small>=$min_aa) {
			$r = get_ranges_prop($align,$i,$start,$end,$small);
			for($j=0; $j<=$#{$r->{range}}; ++$j) {
		           print OUT "COLOUR_REGION   ",$i+1," ",$r->{range}[$j]{st}," ",$i+1," ",$r->{range}[$j]{en}," 13\n"; 
			}
		} elsif($n_hydrophobic>=$min_aa) { 
			$r = get_ranges_prop($align,$i,$start,$end,$hydrophobic);
			for($j=0; $j<=$#{$r->{range}}; ++$j) {
		           print OUT "COLOUR_REGION   ",$i+1," ",$r->{range}[$j]{st}," ",$i+1," ",$r->{range}[$j]{en}," 8\n"; 
			}
		} elsif($n_polar>=$min_aa) { 
			$r = get_ranges_prop($align,$i,$start,$end,$polar);
			for($j=0; $j<=$#{$r->{range}}; ++$j) {
		           print OUT "CCOL_CHARS  ",$polar," ",$i+1," ",$r->{range}[$j]{st}," ",$i+1," ",$r->{range}[$j]{en}," 2\n"; 
			}
		}
	}
	# Do secondary structures
	for($i=0; $i<$align->{nseq}; ++$i) {
		$id = $align->{list}[$i];
#		print "Considering $id\n";
		if(($use{$id}==1) && (($id =~ /Sec_/) || ($id =~ /_dssp/))) { # This is a secondary structure
			print "$id is a secondary structure\n";
			for($j=0; $j<length($ss_string); ++$j) {
				print OUT "SUB_CHARS 1 ",$i+1," ",$align->{alen}," ",$i+1," ",substr($ss_string,$j,1)," SPACE\n";
			}
			print OUT "# Sec strucs for $id\n";
			for($j=0; $j<$align->{alen}; ++$j) {
		 		$ss = substr($align->{ids}{$id}{seq},$j,1);
				$ss =~ s/\?/ /;
#				print "SS is $ss\n";
				
				if($helix_def =~ /$ss/) { # Helix
					$sec_start = $j+1;
					while((($helix_def =~ /$ss/) || ($ss eq " ")) && ($j<$align->{alen})) {
						if($ss ne " ") { $sec_end = $j; }
						$ss = substr($align->{ids}{$id}{seq},$j,1); $ss =~ s/\?/ /;
						$j++;
					}
					print OUT "HELIX ",$sec_start," ",$i+1," ",$sec_end,"\n";
					print OUT "COLOUR_TEXT_REGION ",$sec_start," ",$i+1," ",$sec_end," ",$i+1," 9 \n";
				}
				if($strand_def =~ /$ss/) { # Strand
					$sec_start = $j+1;
					while((($strand_def =~ /$ss/) || ($ss eq " ")) && ($j<$align->{alen})) {
						if($ss ne " ") { $sec_end = $j; }
						$ss = substr($align->{ids}{$id}{seq},$j,1); $ss =~ s/\?/ /;
						$j++;
					}
					print OUT "STRAND ",$sec_start," ",$i+1," ",$sec_end,"\n";
					print OUT "COLOUR_TEXT_REGION ",$sec_start," ",$i+1," ",$sec_end," ",$i+1," 6\n";
				}
			  }
		}
	}
	# Do boxing for stamp equivalences
	if((defined($stamp_char)) && (defined($stamp_cut)) && (defined($stamp_win))) {
            @equiv = get_stamp_equiv($align,$stamp_char,$stamp_cut,$stamp_win); 
#	    print "Read equiv for ",$stamp_char,$stamp_cut,$stamp_win,"\n";
	    if(@equiv) { # give up if no STAMP (or wrong field) found
               for($i=0; $i<$align->{alen}; ++$i) {
                   if($equiv[$i] ==1) {
                       print OUT "BOX_REGION ",$i+1," ",$start+1," ";
                       while(($equiv[$i]==1) && ($i<$align->{alen})) { $i++; }
                       print OUT $i," ",$end+1,"\n";
                  } 
               } 
            } else { 
		die "Error reading equivalences\n";
	    }
	} 
}

sub get_stamp_equiv {
        my($i,$j);
	my($align)      = $_[0];
	my($stamp_char) = $_[1];
	my($stamp_win)  = $_[2];
	my($stamp_cut)  = $_[3];
	my($dat);
	my($stretch);

	if(defined($align->{stamp}{$stamp_char})) {
#	   print "Here size of alignment is $align->{alen}\n";
	   my(@rel,$val);

           # Store rel 1/0 array
           for($i=0; $i<$align->{alen}; ++$i) {
             if(defined($align->{p}[$i]{post})) {
		$dat = $align->{p}[$i]{post};
		$dat =~ s/^ *//;
                $val=(split(/ +/,$dat))[$align->{stamp}{$stamp_char}];
                if(defined($val)) {
#		   print "$i value is $val\n";
		   if((($stamp_char ne "A") && ($val>=$stamp_cut)) || (($stamp_char eq "A") && ($val<=$stamp_cut))) {
		      push(@rel,1);
		   } else {
		      push(@rel,0);
	           }
                 } else {
                      push(@rel,0);
                 }
               } else {
                      push(@rel,0);
               }
	   }
           # Smooth rel 1/0 array by window size
           for($i=0; $i<$align->{alen}; ++$i) {
               if($rel[$i] == 1) { # Check neighbours
		  $stretch = 1;
                  $j=$i-1; while($rel[$j]==1) { $stretch++; $j--; }
                  $j=$i+1; while($rel[$j]==1) { $stretch++; $j++; }
		  if($stretch<$stamp_win) { $rel[$i]=0; }
               }
#	       print "Position $i rel[i] is $rel[$i]\n";
           }
           return @rel;
        } 
}

sub get_ranges_prop {

# This delineates a position in an alignment into
#  a set of ranges that match what is taken to be a conserved property

   my($i,$j);
   my($align) = $_[0];
   my($pos)   = $_[1];
   my($start) = $_[2];
   my($end)   = $_[3];
   my($set)   = $_[4];
   my($aa);
   my($r) = {};
   my($id);
   my($nr);

   $nr=0;

   for($i=$start; $i<=$end; ++$i) {
      $id = $align->{list}[$i];
      $aa = substr($align->{ids}{$id}{seq},$pos,1);
      if($set =~ /$aa/) { # This will be a start of a range
	 $r->{range}[$nr]{st} = $i+1;
         while(($set =~ /$aa/) && ($i<=$end)) {
            $id = $align->{list}[$i];
            $aa = substr($align->{ids}{$id}{seq},$pos,1);
	    if($set =~ /$aa/) { $r->{range}[$nr]{en} = $i+1 }
	    $i++;
         }
	 $i--;
	 $nr++;
      }
   }
   return $r;
}

sub write_one_fasta {
	my($outfile,$id,$title,$seq) = @_;

        open(OUT,">$outfile") || die "Error opening output file $outfile\n";
	print OUT ">",$id," ";
	if(defined($title)) { print OUT $title; }
	print OUT "\n";
	for($j=0; $j<length($seq); ++$j) {
		$aa = substr($seq,$j,1);
		if($aa eq " ") { print OUT "-"; }
		else { print OUT $aa; }
#		print OUT substr($seq,$j,1);
		if((($j+1)%60)==0) { print OUT "\n"; }
	}
	print OUT "\n";
	close(OUT);
}

sub get_len {
    my($align) = $_[0];
    my($name) = $_[1];

    my($i,$j,$len,$alen,$aa);
  
    $len=0;
    for($j=0; $j<length($align->{ids}{$name}{seq}); ++$j) {
       $aa=substr($align->{ids}{$id}{seq},$j,1);
       $aa = substr($align->{ids}{$name}{seq},$j,1);
       if($aa=~/[A-Z]/ || $aa=~/[a-z]/) { $len++; }
    }

    return $len;
}

