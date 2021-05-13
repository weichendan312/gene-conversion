use strict;
use Bio::SeqIO;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::Seq::EncodedSeq;
use Bio::AlignIO;

my $input = "mh_oj.merge.quartet";
open(IN, $input) or die "cannot open block files $input due to $!.\n";
open(CV,">mh_oj.merge.quartet1.cv.out") or die "can not open cv out file due to $!.\n";
open(OUT, ">mh_oj.merge.quartet1.newcv2.crossover") or die "cannot open output file due to $!.\n";
open(SIM, ">mh_oj.merge.quartet1.crossover") or die "cannot open output file due to $!.\n";

my $lineno = 0;
while(<IN>)
{
   if($lineno < 0){$lineno ++; next;}
   $_ =~ s/[\n\r]//g;
   my $quartetid = join("_", 
       split(/\s/, $_));
   print CV "$lineno $quartetid\n";
   $lineno++;
   main($quartetid);
}

sub main()
{
   my $quartetid = $_[0];
   my @temarr = split(/_/, $quartetid);
   my ($A, $C, $At, $Ct) = @temarr;
#print OUT $quartetid."\n";

   my @dnaobj1 = ();
   my %dna_hash1;
   my $protos1 = Bio::SeqIO -> new(-file=> ">prot1.fasta", -format=>"fasta");

   my $seqno = 0;
   my $isallfileexist = 1;

   for(my $i=0; $i<=$#temarr; $i++)
   {
      if($temarr[$i] !~ /^\w/){next;}
      if($temarr[$i] =~ /^Mh/)
      {
          #my $fasta = "/home/wangjp/myfile/paper/Rice4/CV/CDS/OB/".$temarr[$i].".fasta";
         my $fasta = "/Users/wangjinpeng/Desktop/Wang-Jianyu/5.25/CV/cds/Mh/".$temarr[$i].".fasta";
         if(!(-e $fasta)){next; $isallfileexist = 0;}

         my $fastaio = new Bio::SeqIO(-format => "FASTA", -file => $fasta);
         $dnaobj1[$seqno] = $fastaio -> next_seq;
         $dnaobj1[$seqno] -> display_id($dnaobj1[$seqno] -> display_id());
         $dna_hash1{$dnaobj1[$seqno]->display_id()} = $dnaobj1[$seqno];
         $protos1 -> write_seq($dnaobj1[$seqno] -> translate());
      }
      elsif($temarr[$i] =~ /^Oj/)
      {
         my $fasta = "/Users/wangjinpeng/Desktop/Wang-Jianyu/5.25/CV/cds/Oj/".$temarr[$i].".fasta";

         if(!(-e $fasta)){next; $isallfileexist = 0;}

         my $fastaio = new Bio::SeqIO(-format => "FASTA", -file => $fasta);
         $dnaobj1[$seqno] = $fastaio -> next_seq;
         $dnaobj1[$seqno] -> display_id($dnaobj1[$seqno] -> display_id());
         $dna_hash1{$dnaobj1[$seqno]->display_id()} = $dnaobj1[$seqno];
         $protos1 -> write_seq($dnaobj1[$seqno] -> translate());
      }
#      elsif($temarr[$i] =~ /^GSBRNA/)
#      {
#         my $fasta = "/home/wangjp/myfile/paper/teacher/Bnapus/data/Bn.cds/".$temarr[$i].".fasta";
#         my $fastaio = new Bio::SeqIO(-format => "FASTA", -file => $fasta);

#         if(!(-e $fasta)){next; $isallfileexist = 0;}

#         $dnaobj1[$seqno] = $fastaio -> next_seq;
#print $dnaobj1[$seqno] -> display_id()."\n";
#print $dnaobj1[$seqno] -> seq()."\n";
#         $dnaobj1[$seqno] -> display_id($dnaobj1[$seqno] -> display_id());
#         $dna_hash1{$dnaobj1[$seqno]->display_id()} = $dnaobj1[$seqno];
#         $protos1 -> write_seq($dnaobj1[$seqno] -> translate());
#      }
      $seqno ++;
   }

   if($isallfileexist eq 0){print OUT "file num $quartetid\n"; next;}

###### without outgroup alignment
   #system("clustalw -infile=prot1.fasta");
   system("clustalw2 -infile=prot1.fasta");
######3clustalw2
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
          my $alignlen = 0;
          for(my $i=0; $i<$len; $i++)
          {
             if($protseq[$k][$i] eq "-" || $protseq[$j][$i] eq "-"){$gapnum ++;}
             elsif($protseq[$k][$i] eq $protseq[$j][$i]){$idennum ++;}

	     if($protseq[$k][$i] =~ /[a-zA-Z]/ &&  $protseq[$j][$i] =~ /[a-zA-Z]/){$alignlen ++;}
          }
          
          my $identity;
	  my $gap_level;
          if($alignlen ne 0)
          {
            $identity = $idennum/$alignlen;
	    $gap_level = $gapnum/$alignlen;
          }
          else
          {
            $identity = -1; $gap_level = -1;
          }

          if($identity < $mini_identity){$mini_identity = $identity;}
          if($gap_level > $max_gap_level){$max_gap_level = $gap_level;}
          if($gapnum > $max_gap_num){$max_gap_num = $gapnum;}
          print SIM $k."\t".$j."\t".$id[$k]."\t".$id[$j]."\t".$identity."\n";
      }
   }

#   return 0;

   #### a restriction to sequence similarity     $max_gap_level > 0.20 ||
   if($mini_identity < 0.40 ||  $len - $max_gap_num < 50){print OUT "low similarity $quartetid\n"; next;}

   # system("rm prot1.fasta prot1.aln prot1.dnd");
   my $os_prot_aln = Bio::AlignIO -> new(-file=>">/Users/wangjinpeng/Desktop/Wang-Jianyu/5.25/CV/recently/clustprot/".$quartetid.".prot.aln", -format=>"CLUSTALW");
   ##C:\\Users\\Administrator\\Desktop\\legume\\legume_old_new\\
   $os_prot_aln -> write_aln($prot_aln1);

   my $dna_aln = &aa_to_dna_aln($prot_aln1, \%dna_hash1);
   my $os_dna_aln = Bio::AlignIO -> new(-file=>">/Users/wangjinpeng/Desktop/Wang-Jianyu/5.25/CV/recently/clustdna/".$quartetid.".dna.aln", -format=>"CLUSTALW");
   $os_dna_aln -> write_aln($dna_aln);

   my $len  =$dna_aln -> length();

   my @id;
   my @seq;
   my $seq_no = $dna_aln->no_sequences;

   for(my $i=1; $i<=$seq_no; $i++)
   {

      my $id = $dna_aln -> get_seq_by_pos($i) -> display_id();
      my @tmpseq = split(//, $dna_aln -> get_seq_by_pos($i) -> seq());
#print OUT $i." ".$id."\n";
      if($id =~ /$A/)
      {
#print OUT "0\n";
         $id[0] = $id; $seq[0] = \@tmpseq;
      }
      elsif($id =~ /$C/)
      {
#print OUT "1\n";
         $id[1] = $id; $seq[1] = \@tmpseq;
      }
      elsif($id =~ /$At/)
      {
#print OUT "2\n";
         $id[2] = $id; $seq[2] = \@tmpseq;
      }
      elsif($id =~ /$Ct/)
      {
#print OUT "3\n";
         $id[3] = $id; $seq[3] = \@tmpseq;
      }
   }
   ###################caculate CV
   my $matchlen = 0;
   my $para_A = 0;
   my $para_B = 0;
   my $ortho_A1 = 0;
   my $ortho_A2 = 0;

   for(my $s=0; $s<$len; $s++)
   {
	   if($seq[0][$s] eq "-" || $seq[1][$s] eq "-" || $seq[2][$s] eq "-" || $seq[3][$s] eq "-"){next;}
	   $matchlen ++;
	   for(my $p=0; $p<=3; $p++)
	   {
		   for(my $q=$p+1; $q<=3; $q++)
		   {
			   my $linkedid = $id[$p].$id[$q];
			   if($linkedid =~/$A/ && $linkedid =~/$At/)
			   {
				   if($seq[$p][$s] eq $seq[$q][$s]){$para_A++;}
			   }
			   if($linkedid =~/$C/ && $linkedid =~/$Ct/)
			   {
				   if($seq[$p][$s] eq $seq[$q][$s]){$para_B++;}
			   }
			   if($linkedid =~/$A/ && $linkedid =~/$C/)
			   {
				   if($seq[$p][$s] eq $seq[$q][$s]){$ortho_A1++;}
			   }
			   if($linkedid =~/$At/ && $linkedid =~/$Ct/)
			   {
				   if($seq[$p][$s] eq $seq[$q][$s]){$ortho_A2++;}
			   }
		   }
	   }
   }

   print CV "siimilarity pA: ".$para_A." pB: ".$para_B." oA1: ".$ortho_A1." oA2: ".$ortho_A2."\n";
   my $para_Ar = $para_A/$matchlen;
   my $para_Br = $para_B/$matchlen;
   my $ortho_A1r = $ortho_A1/$matchlen;
   my $ortho_A2r = $ortho_A2/$matchlen;
   print CV "siimilarity pAr: ".$para_Ar." pBr: ".$para_Br." oA1r: ".$ortho_A1r." oA2r: ".$ortho_A2r."\n";

   if($para_Ar > $ortho_A1r && $para_Ar > $ortho_A2r)
   {
	   print CV "$quartetid: likely gene conversion occurred in A_species between $A and $At\t";
	   if($ortho_A1r > $ortho_A2r){print CV "likely converted from $At to $A\n";}
	   elsif($ortho_A1r < $ortho_A2r){print CV "likely converted from $A to $At\n";}
	   else{print CV "likely converted from NO_KNOWN to NO_KNOWN\n";}
   }
   if($para_Br > $ortho_A1r && $para_Br > $ortho_A2r)
   {
	   print CV "$quartetid: likely gene conversion occurred in C_species between $C and $Ct\t";
	   if($ortho_A1r > $ortho_A2r){print CV "likely converted from $Ct to $C\n";}
	   elsif($ortho_A1r < $ortho_A2r){print CV "likely converted from $C to $Ct\n";}
	   else{print CV "likely converted from NO_KNOWN to NO_KNOWN\n";}
   }
    print CV "\n";
   #
   ######################################
#   for(my $s=0; $s<$len; $s++)
#   {
#      if($seq[0][$s] eq "-" || $seq[1][$s] eq "-" || $seq[2][$s] eq "-" || $seq[3][$s] eq "-"){next;}
#      if($seq[0][$s] eq $seq[1][$s]){next;}
#      if($seq[0][$s] eq $seq[2][$s] && $seq[0][$s] eq $seq[3][$s])
#      {
#         print CV "@id[0..3] ".$s." ".$seq[0][$s]." ".$seq[1][$s]." ".$seq[2][$s]." ".$seq[3][$s]." C2A\n";
#      }
#      elsif($seq[1][$s] eq $seq[2][$s] && $seq[1][$s] eq $seq[3][$s])
#      {
#         print CV "@id[0..3] ".$s." ".$seq[0][$s]." ".$seq[1][$s]." ".$seq[2][$s]." ".$seq[3][$s]." A2C\n";
#      }
#      elsif($seq[0][$s] eq $seq[2][$s] && $seq[1][$s] eq $seq[3][$s])
#      {
#         print CV "@id[0..3] ".$s." ".$seq[0][$s]." ".$seq[1][$s]." ".$seq[2][$s]." ".$seq[3][$s]." no_conversion\n";
#      }
#      elsif($seq[0][$s] eq $seq[2][$s] && $seq[1][$s] ne $seq[3][$s])
#      {
#         print CV "@id[0..3] ".$s." ".$seq[0][$s]." ".$seq[1][$s]." ".$seq[2][$s]." ".$seq[3][$s]." C/Ct mutation\n";
#      }
#      elsif($seq[0][$s] ne $seq[2][$s] && $seq[1][$s] eq $seq[3][$s])
#      {
#         print CV "@id[0..3] ".$s." ".$seq[0][$s]." ".$seq[1][$s]." ".$seq[2][$s]." ".$seq[3][$s]." A/At mutation\n";
#      }
#      elsif($seq[0][$s] eq $seq[3][$s] && $seq[1][$s] eq $seq[2][$s] && $seq[0][$s] ne $seq[1][$s])
#      {
#	      print CV "@id[0..3] ".$s." ".$seq[0][$s]." ".$seq[1][$s]." ".$seq[2][$s]." ".$seq[3][$s]." A<->C crossover\n";
#      }

#      elsif($seq[0][$s] ne $seq[2][$s] && $seq[1][$s] ne $seq[3][$s])
#      {
#         print CV "@id[0..3] ".$s." ".$seq[0][$s]." ".$seq[1][$s]." ".$seq[2][$s]." ".$seq[3][$s]." Any mutation\n";
#      }
#   }
#  print CV "\n\n"; 

}

