use strict;
use Bio::Root::Root;
use Bio::SeqIO;
use vars qw( @ISA $CODONS $synsites %synchanges );
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::Seq::EncodedSeq;
use Bio::AlignIO;
#use Bio::Tools::Run::Alignment::Clustalw;
#use SeqStatistics;
use Bio::Align::DNAStatistics;

#### here, the quartets are true quartets, i.e., thiplets and doublets are NOT included

my $BOOTNUM = 100;

@ISA = qw( Bio::Align::DNAStatistics );

my @t = split '', "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
$CODONS = $Bio::Align::DNAStatistics::CODONS;
$synsites = $Bio::Align::DNAStatistics::synsites;
%synchanges = %Bio::Align::DNAStatistics::synchanges;


# perl program is to cat the two fasta files into one
# according to the gene colinearity in japonica genome

# the colinearity infomation is contained in the directory ~/jpduplication/block/

# the input is the japonica and indica orthologues in dimerfasta
# the output is the 4 homologues in two subspecies which is positionally related

# perl thisname block.*.txt


print "\nUsage: perl perlname  blockfile outfile\n\n";

#### only for checking whether it is right of aligned (rundom) sequences
#my $tmpfile = "check.aln.txt";
#my $os_check = Bio::AlignIO -> new(-file=>">".$tmpfile, -format=>"CLUSTALW");

#my $input = "/home/wangxy/sorghum1.3/data/RAP_Sbi1.4_Anchors/os_sb_quartets.txt";
my $input = $ARGV[0]; #### quartet file;
open(IN, $input) or die "cannot open block files $input due to $!.\n";

my @arr = split(/\//, $input);
my $output = "mh_oj.merge.quartet.for.wholly.P.CV.PaPs.txt";
open(CV, ">".$output) or die "cannot open outfile $output due to $!.\n";
print CV "quartet\tPs(para1)\tPa(para2)\tPs(para2)\tPa(ortho1)\tPs(ortho1)\tPa(ortho1)\tPs(ortho2)\tKa(para1)\tKs(para1)\tKa(para2)\tKs(para2)\tKa(ortho1)\tKs(ortho1)\tKa(ortho1)\tKs(ortho2)\tcv1\tcv2\tpv(Ps1Pa1Ps2Pa2)\n";
my $lineno = 0;
while(<IN>)
{
   if($lineno eq 0){$lineno++; next;}
   
   $_ =~ s/[\n\r]//g;
   my @temarr = split(/\s+/, $_);
#   my $quartetno = $temarr[0];
   my $osid1 = $temarr[1];
   my $osid2 = $temarr[3];
   my $sbid1 = $temarr[2];
   my $sbid2 = $temarr[4];

#   if($osid1 !~ /^Oj11/){next;}
   
   my $quartetid = $osid1."_".$osid2."_".$sbid1."_".$sbid2;
print $quartetid."-----------------------\n";
   my @dnaobj1 = ();
   my %dna_hash1;
   my $protos1 = Bio::SeqIO -> new(-file=> ">prot1.fasta", -format=>"fasta");

   my $seqno = 0;
   for(my $i=1; $i<=$#temarr; $i++)
   {
      if($temarr[$i] !~ /^\w/){next;}
      if($temarr[$i] =~ /^Mh/)
      {
         my $fasta = "/Users/wangjinpeng/Desktop/Wang-Jianyu/5.25/CV/cds/Mh/".$temarr[$i].".fasta";
         my $fastaio = new Bio::SeqIO(-format => "FASTA", -file => $fasta);
         $dnaobj1[$seqno] = $fastaio -> next_seq;
         $dnaobj1[$seqno] -> display_id($dnaobj1[$seqno] -> display_id());
         $dnaobj1[$seqno] -> seq(uc($dnaobj1[$seqno] -> seq));
         $dna_hash1{$dnaobj1[$seqno]->display_id()} = $dnaobj1[$seqno];
         $protos1 -> write_seq($dnaobj1[$seqno] -> translate());
      }
      elsif($temarr[$i] =~ /^Oj/)
      {
         my $fasta = "/Users/wangjinpeng/Desktop/Wang-Jianyu/5.25/CV/cds/Oj/".$temarr[$i].".fasta";
         my $fastaio = new Bio::SeqIO(-format => "FASTA", -file => $fasta);
         $dnaobj1[$seqno] = $fastaio -> next_seq;
         $dnaobj1[$seqno] -> display_id($dnaobj1[$seqno] -> display_id());
         $dnaobj1[$seqno] -> seq(uc($dnaobj1[$seqno] -> seq));
         $dna_hash1{$dnaobj1[$seqno]->display_id()} = $dnaobj1[$seqno];
         $protos1 -> write_seq($dnaobj1[$seqno] -> translate());
      }
      $seqno ++;
   }

###### without outgroup alignment
   system("clustalw2 -infile=prot1.fasta");

   my $is_prot_aln1 = Bio::AlignIO -> new(-file=>"prot1.aln", -format=>"CLUSTALW");
   my $prot_aln1 = $is_prot_aln1 -> next_aln();

   #### length of alignment
   my $len  = $prot_aln1 -> length;

   my @id;
   my @protseq;
   my $seq_no = $prot_aln1 -> no_sequences;

   #### read the id and sequence
   for(my $i=0; $i<$seq_no; $i++)
   {
       $id[$i] = $prot_aln1 -> get_seq_by_pos($i+1) -> display_id();
       my @tmpseq = split(//, $prot_aln1 -> get_seq_by_pos($i+1) -> seq());
       $protseq[$i] = \@tmpseq;
   }

   #### minimum aa identity among genes
   my $mini_identity = 1;
   my $max_gap_level = 0;
   my $max_gap_num   = 0;
   for(my $k=0; $k<$seq_no; $k++)
   {
      for(my $j=$k+1; $j<$seq_no; $j++)
      {
          my $gapnum = 0;
          my $idennum= 0;
          for(my $i=0; $i<$len; $i++)
          {
             if($protseq[$k][$i] eq "-" || $protseq[$j][$i] eq "-"){$gapnum ++;}
             if($protseq[$k][$i] eq $protseq[$j][$i]){$idennum ++;}
          }
          
          my $identity;
          if($len - $gapnum ne 0)
          {
            $identity = $idennum/($len-$gapnum);
          }
          else
          {
            $identity = 0;
          }
          my $gap_level = $gapnum/$len;
          if($identity < $mini_identity){$mini_identity = $identity;}
          if($gap_level > $max_gap_level){$max_gap_level = $gap_level;}
          if($gapnum > $max_gap_num){$max_gap_num = $gapnum;}
      }
   }

   #### a restriction to sequence similarity     $max_gap_level > 0.20 ||
   #if($mini_identity < 0.40 ||  $len - $max_gap_num < 50){next;}

system("rm prot1.fasta prot1.aln prot1.dnd");

   my $dna_aln1 = &aa_to_dna_aln($prot_aln1, \%dna_hash1);

   my $os_dna_aln1 = Bio::AlignIO -> new(-file=>">/Users/wangjinpeng/Desktop/Wang-Jianyu/5.25/CV/wholly/out.cv.mh/".$quartetid.".dna1.aln", -format=>"CLUSTALW");
   $os_dna_aln1 -> write_aln($dna_aln1);
  
   my @newidseq  = remove_gap_new($dna_aln1);
 
   if(length($newidseq[1]) eq 0){next;}
  
   my @newseqobj;
   my %dna_hash;
   my @prot_arr;
   for(my $i=0; $i<=$#newidseq; $i=$i+2)
   {
      my $newseqobj = Bio::Seq -> new( -id => $newidseq[$i], -seq  => $newidseq[$i+1]);
      $newseqobj[$i/2] = $newseqobj;
      $dna_hash{$newidseq[$i]} = $newseqobj;
#print CV $newidseq[$i]." ".$newseqobj->seq."\n";
      push @prot_arr, $newseqobj -> translate();
   }

   my $protos = Bio::SeqIO -> new(-file=> ">prot.fasta", -format=>"fasta");
   for(my $i=0; $i<=$#prot_arr; $i++)
   {$protos -> write_seq($prot_arr[$i]);}
   system("clustalw2 -infile=prot.fasta");

   my $is_prot_aln = Bio::AlignIO -> new(-file=>"prot.aln", -format=>"CLUSTALW");
   my $prot_aln = $is_prot_aln -> next_aln();
  
   system("rm prot.fasta prot.aln prot.dnd");

   my $os_prot_aln = Bio::AlignIO -> new(-file=>">/Users/wangjinpeng/Desktop/Wang-Jianyu/5.25/CV/wholly/protein.aln.mh/".$quartetid.".prot.aln", -format=>"CLUSTALW");
   $os_prot_aln -> write_aln($prot_aln);

   my $dna_aln = &aa_to_dna_aln($prot_aln, \%dna_hash);
   my $os_dna_aln = Bio::AlignIO -> new(-file=>">/Users/wangjinpeng/Desktop/Wang-Jianyu/5.25/CV/wholly/dna.aln.mh/".$quartetid.".dna.aln", -format=>"CLUSTALW");
   $os_dna_aln -> write_aln($dna_aln);

#   my $cdsalnfile = "alignment/pamldna/".$quartetid.".NUC";
#   write_aln_paml($cdsalnfile, $dna_aln);

#   my $seqfile = "alignment/pamldna/".$quartetid.".NUC";
#   my $outfile = "KV.yn00";
#   set_ctl_parameters("seqfile", $seqfile,
#                      "outfile", $outfile);
#   system("yn00");

   my $os_para = $osid1."-".$osid2;
   my $sb_para = $sbid1."-".$sbid2;
   my $os_sb_ortho1 = $osid1."-".$sbid1;
   my $os_sb_ortho2 = $osid2."-".$sbid2;

   my (%PaPs) = read_NG($os_para, $sb_para, $os_sb_ortho1, $os_sb_ortho2, $dna_aln);

   my @PaPs = %PaPs;
#   print CV "PaPs @PaPs[0..$#PaPs]\n";

#print CV $os_para."\n";
#print CV $sb_para."\n";
#print CV $os_sb_ortho1."\n";
#print CV $os_sb_ortho2."\n";

    my (%Pa, %Ps);
    foreach my $pair(keys(%PaPs))
    {
#print "pair is $pair....$PaPs{$pair}\n";
       my @pairvalue = split(/ /, $PaPs{$pair});
       $Pa{$pair} = $pairvalue[0];
       $Ps{$pair} = $pairvalue[1];
    }
    print  CV $quartetid."\t";
    print  CV $Pa{$os_para}."\t".$Ps{$os_para}."\t".$Pa{$sb_para}."\t".$Ps{$sb_para}."\t";
    print  CV $Pa{$os_sb_ortho1}."\t".$Ps{$os_sb_ortho1}."\t".$Pa{$os_sb_ortho2}."\t".$Ps{$os_sb_ortho2}."\t";

    print  CV jk($Pa{$os_para})."\t".jk($Ps{$os_para})."\t".jk($Pa{$sb_para})."\t".jk($Ps{$sb_para})."\t";
    print  CV jk($Pa{$os_sb_ortho1})."\t".jk($Ps{$os_sb_ortho1})."\t".jk($Pa{$os_sb_ortho2})."\t".jk($Ps{$os_sb_ortho2})."\t";

   my $is_conversion_os=0;
   my $is_conversion_sb=0;

   if($Ps{$os_para} <= $Ps{$os_sb_ortho1} && $Ps{$os_para} <= $Ps{$os_sb_ortho2}
    &&$Pa{$os_para} <= $Pa{$os_sb_ortho1} && $Pa{$os_para} <= $Pa{$os_sb_ortho2}
   )
   {
      $is_conversion_os = 1; print CV "R\t";
   }
   else
   {  print CV "N\t";
   }

   if($Ps{$sb_para} <= $Ps{$os_sb_ortho1} && $Ps{$sb_para} <= $Ps{$os_sb_ortho2}
   && $Pa{$sb_para} <= $Pa{$os_sb_ortho1} && $Pa{$sb_para} <= $Pa{$os_sb_ortho2}
   )
   {
      $is_conversion_sb = 1; print CV "S\t";
   }
   else
   {  print CV "N\t";
   }

   my @bootprob;
   if($is_conversion_os +  $is_conversion_sb ne 0)
   {
      @bootprob = bootstrap_global($os_para, $sb_para, $os_sb_ortho1, $os_sb_ortho2, $dna_aln);
   }
   if($is_conversion_os eq 1 &&  $is_conversion_sb eq 1)
   {
     print CV "@bootprob[0..$#bootprob]\n";
   }
   elsif($is_conversion_os eq 1)
   {
     print CV "@bootprob[0..1] - -\n";
   }
   elsif($is_conversion_sb eq 1)
   {
     print CV "- - @bootprob[2..3]\n";
   }
   else
   {
     print CV "\n";
   }
   
   $lineno ++;
}
close($ARGV[0]);

sub read_NG
{
   my $os_para = $_[0];
   my $sb_para = $_[1];
   my $os_sb_ortho1 = $_[2];
   my $os_sb_ortho2 = $_[3];
   my $dna_aln = $_[4];
#print "os $os_para\n";
#print "sb $sb_para\n";
#print "ot $os_sb_ortho1\n";
#print "ot $os_sb_ortho2\n";

   my @seq = ();
   my @id  = ();
   my $seqno = 0;
foreach my $seq ($dna_aln -> each_seq())
{
   $seq[$seqno] = $seq->seq();
   $id[$seqno] = $seq ->display_id();
   $seqno ++;
   #print CV $seq -> display_id().":\n".$seq->seq()."\n";
}

   my %PaPs;

   for(my $i=0; $i<=$#seq; $i++)
   {
      for(my $j=$i+1; $j<=$#seq; $j++)
      {
         my ($syn_count, $non_syn_count, $aa_count, $r_p_count, $gap_cnt) = analyse_mutations($seq[$i], $seq[$j]);
         my $av_s_site = count_av_syn_sites($seq[$i], $seq[$j], $synsites);
         my $av_ns_syn_site = length($seq[0]) - $gap_cnt- $av_s_site;
         my $syn_prop = $syn_count / $av_s_site;
         my $nc_prop = $non_syn_count / $av_ns_syn_site;
         $PaPs{$id[$i]."-".$id[$j]} = $nc_prop." ".$syn_prop;
         $PaPs{$id[$j]."-".$id[$i]} = $nc_prop." ".$syn_prop;
      }
   }

   return ($os_para, $PaPs{$os_para}, $sb_para, $PaPs{$sb_para}, $os_sb_ortho1, $PaPs{$os_sb_ortho1}, $os_sb_ortho2, $PaPs{$os_sb_ortho2});
}

sub bootstrap_global()
{
#   print "in bootstraping ................\n";
   my $os_para = $_[0];
   my $sb_para = $_[1];
   my $os_sb_ortho1 = $_[2];
   my $os_sb_ortho2 = $_[3];
   my $dna_aln = $_[4];

#   $os_check -> write_aln($dna_aln);

   my $dna_aln_rand = $dna_aln;

   my $align_len = $dna_aln -> length();
   my %align_seq;
   my @seq;
   my %align_seq_rand;
   my @tmpid;
   my $i = 0;
   foreach my $seq_obj ($dna_aln -> each_seq())
   {
      $align_seq{$seq_obj -> display_id()} = $seq_obj -> seq();
      my @arr = split(//, $seq_obj -> seq());
      $seq[$i] = \@arr;
#      print "\n############\n".$seq_obj -> seq()."\n";
      $align_seq_rand{$seq_obj -> display_id()} = "";
      $tmpid[$i] = $seq_obj -> display_id();
      $i = $i + 1;
   }
   my @id = sort(@tmpid);

   #### type 1 bootstrapping will be unvalid when the sequences are quite similar
   #### say, if only 5% mutation among homologs, it will be ~5% of the mutation-involved
   #### codon columns to be selected for bootstrap, resulting the identidcal random seqs
   #### then, conversion bootstrap values will <~ 5%

   #### for type 2 bootstrapping approach, the percentages of the identical columns, mutation-
   #### involved columns and gap-involved columns in the random sequences will be equal to that
   #### in the actual sequences, which will be a prerequisition to construct random sequences.
   my (@diffpos, @idenpos, @gappos);
   my ($k1, $k2, $k3) = (0, 0, 0);
   for(my $i=0; $i<$align_len-3; $i=$i+3)
   {
      my $isgap=0;
      my $isiden=0;
      my $codon0 = $seq[0][$i].$seq[0][$i+1].$seq[0][$i+2];
      if($codon0 =~ /-/){$isgap = 1;}
      if($isgap eq 0)
      {
         for(my $j=1; $j<=$#id; $j++)
         {
            my $codon1 = $seq[$j][$i].$seq[$j][$i+1].$seq[$j][$i+2];
            if($codon1 =~ /-/){$isgap = 1; last;}
            if($codon1 eq $codon0){$isiden ++;}
         }
      }
      if($isgap eq 1){$gappos[$k1] = int($i/3); $k1++;}
      elsif($isiden eq 3){$idenpos[$k2] = int($i/3);  $k2++;}
      else{$diffpos[$k3] = int($i/3); $k3++;}
   }

   my ($count_conversion_os_Ks, $count_conversion_os_Ka) = (0, 0);
   my ($count_conversion_sb_Ks, $count_conversion_sb_Ka) = (0, 0);

   for(my $r=1; $r<=$BOOTNUM; $r++)
   {
#print "\n----------->>>".$r."\n";

       my @random_num = create_rand($align_len/3, $align_len/3);

       #my @random_num1 = create_rand2(@diffpos);
       #my @random_num2 = create_rand2(@idenpos);
       #my @random_num3 = create_rand2(@gappos);
       #my @random_num = (@random_num1, @random_num2, @random_num3);

       @random_num = permut(@random_num);

       for(my $i=0; $i<=$#random_num; $i++)
       {
          if($i eq 0)
          {
             foreach my $so ($dna_aln -> each_seq())
             {$align_seq_rand{$so -> display_id()} = "";}
          }
          for(my $j=0; $j<=$#id; $j++)
          {
             my $codon = substr($align_seq{$id[$j]}, ($random_num[$i]-1) * 3, 3);

             ##### purge off gaps after selecting the codon columns in the random sequences
             if($codon !~ /-/){$align_seq_rand{$id[$j]} .= $codon;}

#             print $random_num[$i]." ".$j." ".substr($align_seq{$id[$j]}, ($random_num[$i]-1) * 3, 3)." ";
          }
       }

       my %dna_hash;
       my $protos = Bio::SeqIO -> new(-file=> ">prot.fasta", -format=>"fasta");

       foreach my $so ($dna_aln_rand -> each_seq())
       {
          my $len = length($align_seq_rand{$so -> display_id()});
 #         print "length ...".$len."\n";
          $so -> seq($align_seq_rand{$so -> display_id()});
          $dna_hash{$so->display_id()} = $so;
          $protos -> write_seq($so -> translate());
 #         print "\n===========\n".$align_seq_rand{$so -> display_id()}."\n";
       }
       ###### without outgroup alignment
      system("clustalw2 -infile=prot.fasta");

      if(!(-e "prot.aln")){$r= $r-1; next;}
      if(-z "prot.aln"){print "prot.aln is zero-sized\n"; $r = $r-1; next;}

      my $is_prot_aln = Bio::AlignIO -> new(-file=>"prot.aln", -format=>"CLUSTALW");
      my $prot_aln = $is_prot_aln -> next_aln();

      system("rm prot.fasta prot.aln prot.dnd");

      my $dna_aln = &aa_to_dna_aln($prot_aln, \%dna_hash);

#      $os_check -> write_aln($dna_aln);

      ##### in essense, the following Ka Ks are Pa Ps
      my (%KaKs) = read_NG($os_para, $sb_para, $os_sb_ortho1, $os_sb_ortho2, $dna_aln);

      my (%Ka, %Ks);
      foreach my $pair(keys(%KaKs))
      {
         print $pair." ".$KaKs{$pair}."\n";
         my @pairvalue = split(/ /, $KaKs{$pair});
         $Ka{$pair} = $pairvalue[0];
         $Ks{$pair} = $pairvalue[1];
      }

      if($Ks{$os_para} <= $Ks{$os_sb_ortho1} && $Ks{$os_para} <= $Ks{$os_sb_ortho2})
      {
         $count_conversion_os_Ks ++;
      }

      if($Ks{$sb_para} <= $Ks{$os_sb_ortho1} && $Ks{$sb_para} <= $Ks{$os_sb_ortho2})
      {
         $count_conversion_sb_Ks ++;
      }
      
      if($Ka{$os_para} <= $Ka{$os_sb_ortho1} && $Ka{$os_para} <= $Ka{$os_sb_ortho2})
      {
         $count_conversion_os_Ka ++;
      }

      if($Ka{$sb_para} <= $Ka{$os_sb_ortho1} && $Ka{$sb_para} <= $Ka{$os_sb_ortho2})
      {
         $count_conversion_sb_Ka ++;
      }
   }
   return ($count_conversion_os_Ks/$BOOTNUM, $count_conversion_os_Ka/$BOOTNUM, 
           $count_conversion_sb_Ks/$BOOTNUM, $count_conversion_sb_Ka/$BOOTNUM);
}

sub create_rand2()
{
   my @arr = @_;
   my $num = $#arr+1;
   my $max = $#arr+1;

   my @random_num = create_rand($num, $max);
   my @random_pos;
   for(my $i=0; $i<$num; $i++)
   {
      $random_pos[$i] = $arr[$random_num[$i]];
   }
   return @random_pos;
}

sub permut
{
   my @arr = @_;
   my $num = $#arr+1;
   my $max = $#arr+1;

   my @random_num = create_rand($num, $max);
   my @arr_permuted;
   for(my $i=0; $i<$num; $i++)
   {
      $arr_permuted[$i] = $arr[$random_num[$i]];
   }
   return @arr_permuted;
}

sub create_rand()
{
   my $num = $_[0];
   my $max = $_[1];
   my @random_num;
   for(my $i=0; $i<$num; $i++)
   {
      $random_num[$i] = int(rand($max));
   }
   return @random_num;
}

sub remove_gap_new()
{
   my $aln = $_[0];
   my $aln_len = $aln -> length();
   my @oldidseq = ();
   my @newidseq;
   my $i = 0;

   foreach my $seq ( $aln -> each_alphabetically() )
   {
       print $i." ".$seq -> display_id."\n";
       my @t_arr = split(//, $seq -> seq);
       $oldidseq[$i+1] = \@t_arr;
       $oldidseq[$i] = $seq -> display_id;
       $i = $i+2;
   }

   for($i=0; $i<=$#oldidseq; $i++)
   {
      if($i%2 eq 0){$newidseq[$i] = $oldidseq[$i];}else{$newidseq[$i] = "";}
   }

   for(my $k = 0; $k < $aln_len; $k ++)
   {
      my $isgap = 0;

      for($i=1; $i<=$#oldidseq; $i=$i+2)
      {
         if($oldidseq[$i][$k] eq "-"){$isgap = 1; last;}
      }

      if($isgap eq 1)
      {; }
      else
      {
         for($i=1; $i<=$#oldidseq; $i=$i+2)
         {
            $newidseq[$i] .= $oldidseq[$i][$k];
         }
      }
   }

#print "@newidseq[0..$#newidseq]\n";
   return(@newidseq);
}
sub analyse_mutations {
#  print ">>> in SeqStatistics::analysie_mutations\n";
  #compares 2 sequences to find the number of synonymous/non synonymous
  # mutations between them
  my ($seq1, $seq2) = @_;
  my %mutator = (2=> {0=>[[1,2], #codon positions to be altered depend on which is the same
			  [2,1]],
		      1=>[[0,2],
			  [2,0]],
		      2=>[[0,1],
			  [1,0]],	},
		 3=> [		#all need to be altered 
		      [0,1,2],
		      [1,0,2],
		      [0,2,1],
		      [1,2,0],
		      [2,0,1],
		      [2,1,0] ],
		);
  my $TOTAL = 0;		# total synonymous changes
  my $TOTAL_n = 0;      	# total non-synonymous changes
  my $TOTAL_a = 0;		# total amino-acid changes
  my @TOTAL_p = ( 0, 0, 0 );	# total uncleotide acid changes in 3 positions of codons
  my $gap_cnt = 0;

  my %input;
  my $seqlen = length($seq1);
  for (my $j=0; $j< $seqlen; $j+=3) {
    $input{'cod1'} = substr($seq1, $j,3);
    $input{'cod2'} = substr($seq2, $j,3);

    #ignore codon if beeing compared with gaps! 
    if ($input{'cod1'} =~ /\-/ || $input{'cod2'} =~ /\-/) {
      $gap_cnt += 3; #just increments once if there is a apair of gaps
      next;
    }

    $TOTAL_p[0]++ unless substr($input{'cod1'}, 0, 1) eq substr($input{'cod2'}, 0, 1);
    $TOTAL_p[1]++ unless substr($input{'cod1'}, 1, 1) eq substr($input{'cod2'}, 1, 1);
    $TOTAL_p[2]++ unless substr($input{'cod1'}, 2, 1) eq substr($input{'cod2'}, 2, 1);

    # different animo acids
    $TOTAL_a++ unless _is_same_aa(\%input);

    my ($diff_cnt, $same) = Bio::Align::DNAStatistics::count_diffs(\%input);

    #ignore if codons are identical
    next if $diff_cnt == 0 ;
    if ($diff_cnt == 1) {
      $TOTAL += $synchanges{$input{'cod1'}}{$input{'cod2'}};
      $TOTAL_n += 1 - $synchanges{$input{'cod1'}}{$input{'cod2'}};
      #print " \nfordiff is 1 , total now $TOTAL, total n now $TOTAL_n\n\n"
    } elsif ($diff_cnt ==2) {
      my $s_cnt = 0;
      my $n_cnt = 0;
      my $tot_muts = 4;
      #will stay 4 unless there are stop codons at intervening point
    OUTER:for my $perm (@{$mutator{'2'}{$same}}) {
	my $altered = $input{'cod1'};
	my $prev= $altered;
	#		print "$prev -> (", $t[$CODONS->{$altered}], ")";
	for my $mut_i (@$perm) { #index of codon mutated
	  substr($altered, $mut_i,1) = substr($input{'cod2'}, $mut_i, 1);
	  if ($t[$CODONS->{$altered}] eq '*') {
	    $tot_muts -=2;
	    #print "changes to stop codon!!\n";
	    next OUTER;
	  } else {
	    $s_cnt += $synchanges{$prev}{$altered};
	    #					print "$altered ->(", $t[$CODONS->{$altered}], ") ";
	  }
	  $prev = $altered;
	}
	#		print "\n";
      }
      if ($tot_muts != 0) {
	$TOTAL += ($s_cnt/($tot_muts/2));
	$TOTAL_n += ($tot_muts - $s_cnt)/ ($tot_muts / 2);
      }
 
    } elsif ($diff_cnt ==3 ) {
      my $s_cnt = 0;
      my $n_cnt = 0;
      my $tot_muts = 18;	#potential number  of mutations
    OUTER: for my $perm (@{$mutator{'3'}}) {
	my $altered = $input{'cod1'};
	my $prev= $altered;
	#	print "$prev -> (", $t[$CODONS->{$altered}], ")";
	for my $mut_i (@$perm) { #index of codon mutated
	  substr($altered, $mut_i,1) = substr($input{'cod2'}, $mut_i, 1);
	  if ($t[$CODONS->{$altered}] eq '*') {
	    $tot_muts -=3;
	    #	print "changes to stop codon!!\n";
	    next OUTER;
						
	  } else {
	    $s_cnt += $synchanges{$prev}{$altered};
	    #			print "$altered ->(", $t[$CODONS->{$altered}], ") ";
	  }
	  $prev = $altered;
	}
	#	print "\n";
			 
      }				#end OUTER loop
      #calculate number of synonymous/non synonymous mutations for that codon
      # and add to total
      if ($tot_muts != 0) {
	$TOTAL += ($s_cnt / ($tot_muts /3));
	$TOTAL_n += 3 - ($s_cnt / ($tot_muts /3));
      }
    }				#endif $diffcnt = 3
  }				#end of sequencetraversal
#  print " there are $TOTAL syn mutations and $TOTAL_n non-syn  mutations $TOTAL_a aa\n";
#  print "<<< get out of SeqStatistics::analysie_mutations\n";
  return ($TOTAL, $TOTAL_n, $TOTAL_a, \@TOTAL_p, $gap_cnt);
}

sub _is_same_aa {
  my ($r_input) = @_;
  my $cod1 = uc $r_input->{'cod1'};
  my $cod2 = uc $r_input->{'cod2'};
  return ( $t[$CODONS->{$cod1}] eq $t[$CODONS->{$cod2}] );
}

sub pc {
  my ($self, $p_a) = @_;
  return -1 * log(1 - $p_a);
}

sub pc_var {
  my ($p, $n) = @_;
  return $p / log((1 - $p) * $n);
}

sub jk {
	my ($p) = $_[0];
	if ($p >= 0.75) {
#		print " Jukes Cantor won't  work -too divergent!";
		return -1;
		}
	return -1 * (3/4) * (log(1 - (4/3) * $p));
}

# Override for debug
sub count_syn_sites {
    #counts the number of possible synonymous changes for sequence
    my ($seq, $synsite) = @_;
    die "not integral number of codons" if length($seq) % 3 != 0;
    my $S = 0;
    for (my $i = 0; $i< length($seq); $i+=3) {
	my $cod = substr($seq, $i, 3);
	next if $cod =~ /\-/;	#deal with alignment gaps
	$S +=  $synsite->{$cod}{'s'};
    }
    #print "S is $S\n";
    return $S;
}

sub count_av_syn_sites {
  # counts the average number of possible synonymous changes
  # for a pair of sequences
  my ($seq1, $seq2, $synsite) = @_;
  die "not integral number of codons" if length($seq1) % 3;
  my $S = 0;
  for (my $i = 0; $i < length($seq1); $i += 3) {
    my $cod1 = substr($seq1, $i, 3);
    my $cod2 = substr($seq2, $i, 3);
    next if ( $cod1 =~ /\-/ || $cod2 =~ /\-/ );
    $S += $synsite->{$cod1}{'s'} + $synsite->{$cod2}{'s'};
  }
  return $S / 2;
}
