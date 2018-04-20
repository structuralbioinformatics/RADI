#!/usr/bin/perl

 $sol={};
 $sol={ A => 115.0, R => 263.0, N => 184.0, D => 170.0, C => 149.0, 
        Q => 208.0, E => 207.0,
        G => 86.0, H => 206.0, I => 187.0, L => 192.0, K => 222.0, 
        M => 210.0, F => 230.0,
        P => 140.0 , S => 140.0, T => 164.0, W => 269.0, Y => 257.0, V => 161.0};

 $bur={};
 $bur={ 9 => 'W', 8 => 'B', 7 => 'C', 6 => 'D', 5 => 'E', 
        4 => 'F', 3 => 'G', 2 => 'H', 1 => 'I', 0 => 'A'};

 open (INP,"<$ARGV[0]");

  while (<INP>){
    if (/  #  RESIDUE AA STRUCTURE BP1 BP2  ACC   N-H-->O  O-->H-N  N-H-->O  O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA/){$begin=1; next;}
    if ($begin==1){
      $aa=substr($_,13,1);
      if ($aa eq "!"){next;}
      push @aa,$aa;
      $s=substr($_,16,1); 
      if ($s eq " " || $s eq "-"){
        push @ss,"-";
      }else{
        push @ss,$s; 
      }
      $b=substr($_,35,3);
      $bb=int(10*$b/$sol->{$aa});
      if ($bb < 0){$bb=0;}
      if ($bb >= 10) {$bb= 9;}
      push @bb,$bb;
    }
  }

  close (INP);

  print ">P1;$ARGV[0]Seq\n";
  print "$ARGV[0]Seq \n";
  $n=0;
  for $i (0..$#aa){
   print "$aa[$i]";
   $n++;
   if ($n >79 ){print "\n";$n=0;}
  }
  print "*\n";

  print ">P1;$ARGV[0]SS\n";
  print "$ARGV[0]SS \n";
  $n=0;
  for $i (0..$#ss){
   if ($ss[$i] eq "-") { $new_ss="C";}else{$new_ss=$ss[$i];}
   print "$new_ss";
   $n++;
   if ($n >79 ){print "\n";$n=0;}
  }
  print "*\n";
 
