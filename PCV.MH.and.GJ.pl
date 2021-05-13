#use strict;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Seq::EncodedSeq;

#### So when the sequences are not very similar, adopt Clustalw
$ENV{CLUSTALDIR} = '/Users/wangjinpeng/Downloads/clustalw-2.1-macosx/clustalw-2.1-macosx/clustalw2';
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Tools::Run::Alignment::TCoffee;

use SeqStatistics;
use Bio::Align::Utilities qw(aa_to_dna_aln);

printTime();

#### set values of parameters
my $gappenalty = -2;
my $BOOTNUM = 100;

#### create codon hash table
my %codonhash;
getcodon();

#### OUTPUT files
my $outinfo = $ARGV[0].".0.4.local.out";
unlink($outinfo);
open(INFO, ">>".$outinfo) or die "cannot open output file ".$outinfo." due to $!.";

#### the following input lines after this part will be changed to read the sequences separately
#    according to some matching information
#### containing the tree topology info of the homologs
#
#    though the following lines contain the words "para" and "ortholog"
#    this perl script is not only for general cases of paralogs and orthlogs
#
#    literally, "para" denote two homologs that should be a little farther from one other,
#               "ortho" denote two homologs that should be a little nearer 
#
#               and
#               
#               "species" denote groups of homologs that should be nearer
#
#    Therefore, for duplicated colinear genes in rice japonica and indica genomes,
#    the literal definations are just as what they appearingly are
#    (
#       This is Model 1: ((I1, J1), (I2, J2)) ===> ((A_1, B_1), (A_2, B_2))
#
#       ------I1-----
#       ------J1-----
#       ------I2-----
#       ------J2-----
#    )
#
#    When the tandem paralogs produced before japonica-indica divergence are considerred,
#    "para" denotes two tandem paralogs in each subspecies, while "ortho" denotes two corresponding orthologs
#    in a subspecies.
#    (
#       This is Model 2: ((I1, J1), (I2, J2)) ====> ((A_1, B_1), (A_2, B_2)) 
#
#       ---- I1 --- I2 ---
#       ---- J1 --- J2 --- 
#    )
#    
#    When the tandem paralogs produced before large duplication events are considerred,
#    "para" denotes two tandem paralogs produced before large duplication,while
#    "ortho" denotes two homoelogs in two duplicated segments
#    {
#      This is Model 3: ((J1, J2), (J3, J4)) ====> ((A_1, B_1), (A_2, B_2))
#
#      ---- J1 --- J3 ---
#      ---- J2 --- J4 ---
#    }
#
#    The Model 3 is not quite realistic for the scarcity of colinear duplicates in rice genome.
#    But may be used in poplar genome, where colinear genes are considerably preserved
#  
my $colinearno;  
open(HOMO, $ARGV[0]) or die "cannot open homolog infomation file due to $!\n";
#my $conversionid = $ARGV[1];
#open (CV,">",$conversionid) or die "can not open conversion id file $conversionid due to $!.\n";

while(<HOMO>)
{
   $_ =~ s/[\n\r]//g;
   my @temarr = split(/\s+/, $_);

   my @osid1 = split(/[;\/\|]/, $temarr[0]);
   for(my $i1=0; $i1<=$#osid1; $i1++)
   {
      my @osid2 = split(/[;\/\|]/, $temarr[2]);
      for(my $i2=0; $i2<=$#osid2; $i2++)
      {
         my @sbid1 = split(/[;\/\|]/, $temarr[1]);
         for(my $i3=0; $i3<=$#sbid1; $i3++)
         {
            my @sbid2 = split(/[;\/\|]/, $temarr[3]);
            for(my $i4=0; $i4<=$#sbid2; $i4++)
            {
               print CV "groupsize ".$#osid1." ".$#osid2." ".$#sbid1." ".$#sbid2."\n";
               my $quartetid = $osid1[$i1]."_".$osid2[$i2]."_".$sbid1[$i3]."_".$sbid2[$i4];
#               print CV $quartetid."\n";
               main($quartetid);
            }
         }
      }
   }
}


sub main()
{
   my @tmparr = split(/_/, $_[0]);

   my %id_hash;
   my %id_hash_r;
   my %dna_hash;
   my @prot_arr;
   my @species = (A, B);

   my $quartet = $_[0];

      my $expected_top = 0;
       
      $id_hash{"A_1"} = $tmparr[0];
      $id_hash_r{$tmparr[0]} = "A_1";
      $id_hash{"A_2"} = $tmparr[1];
      $id_hash_r{$tmparr[1]} = "A_2";
      $id_hash{"B_1"} = $tmparr[2];
      $id_hash_r{$tmparr[2]} = "B_1";
      $id_hash{"B_2"} = $tmparr[3];
      $id_hash_r{$tmparr[3]} = "B_2";


      my $species_no = $#species + 1;

      print "the involved species are "."@species[0..$#species]"."\n";

      my $dna_aln_file = "/Users/wangjinpeng/Desktop/Wang-Jianyu/5.25/CV/wholly/dna.aln.mh/".$quartet.".dna.aln";     
     print INFO $dna_aln_file."\n";

      if(!(-e $dna_aln_file)){next;}

      my $is_aln = Bio::AlignIO -> new(-file=>$dna_aln_file, -format=>"CLUSTALW");
      my $dna_aln = $is_aln -> next_aln();

### change seq names in alignment
#   I have to redo the alignment
      foreach my $seqobj ($dna_aln -> each_seq())
      {
#         print "\n".$seqobj -> display_id()."\n".$seqobj -> seq()."\n";
         $seqobj -> display_id($id_hash_r{$seqobj -> display_id()});
         my $seq = $seqobj -> seq();
         $seq =~ s/-//g;
         $seqobj -> seq($seq);
#         print "\n".$seqobj -> display_id()."\n".$seqobj -> seq()."\n";
         $dna_hash{$seqobj -> display_id()} = $seqobj;
         my $prot_obj = $seqobj -> translate();
         push(@prot_arr, $prot_obj); 
      }      
      @params = ('ktuple' => 2, 'matrix' => 'BLOSUM', 'quiet' => 1, 'pwgapopen' => 10);
      
      my $factory = Bio::Tools::Run::Alignment::Clustalw -> new (@params);

      my $prot_aln = $factory -> align(\@prot_arr);

      my $dna_aln = &aa_to_dna_aln($prot_aln, \%dna_hash);


   #### length of alignment
   my $len  = $prot_aln -> length;

   my @id;
   my @protseq;
   my $seq_no = $prot_aln -> no_sequences;


   #### read the id and sequence
   for(my $i=0; $i<$seq_no; $i++)
   {
       $id[$i] = $prot_aln -> get_seq_by_pos($i+1) -> display_id();
       my @tmpseq = split(//, $prot_aln -> get_seq_by_pos($i+1) -> seq());
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
   if($mini_identity < 0.40 ||  $len - $max_gap_num < 50){next;}
   #

      #### global tree topology patterns based on Ks
#      my @comment_G = global_conv(\@species, $dna_aln);
     
      #### partial gene conversion detection
      
      my $is_global_conversion = 0;
      my $is_local_conversion = 0;
      for(my $i=0; $i<$species_no; $i++)
      {
#         print INFO "specise=".$species[$i]." ".$comment_G[$i]."\n";
        
#         if($comment_G[$i] =~ /no global conversion/ && $comment_G[$i] !~ /singleton/)
#         {
            my $real_max = local_conv($i, \@species, $dna_aln);
            if($real_max >= 1){$is_local_conversion = 1;} 
#         }
#         else
#         {
#            $is_global_conversion = 1;
#         }
      }
#      if($is_global_conversion + $is_local_conversion >= 1)
#      {
         #$dnaaln_os2 -> write_aln($dna_aln);
#      }
      

}

printTime();

sub local_conv()
{
   #### Lin et al. 2006 PNAS, their method is used here to find local converted regions in paralogs as compared
   #    with the orthologs from other species
   #
   #    the method involves a simple dynamic programming algorithm
   #    and binormal test
   #     
   my $thissp = $_[0]; # check the paralogs in this species for possible conversion
   my $species= $_[1]; # this is a reference
   my $dna_aln= $_[2];
   #print "in local_conv ...\n";

   my $species_no = $#$species + 1;

   my %dna_aln;
   foreach my $seqobj($dna_aln -> each_seq())
   {
      #print $seqobj -> display_id()."\n".$seqobj -> seq()."\n";
      $dna_aln{$seqobj -> display_id()} = $seqobj;
   }               
   my @seq_id = keys(%dna_aln);

   my $para11 = $$species[$thissp]."_1";
   my $para12 = $$species[$thissp]."_2";
   if(isexist($para11, @seq_id) * isexist($para12, @seq_id) eq 0){return "missing paralog in ".$$species[$thissp];}
 
   #### check the ratio of gaps in the alignment
   #    if too much, discard the task to reveal conversion
   #    
   my $gapratio = check_gap($dna_aln{$para11}, $dna_aln{$para12});
   if($gapratio >= 0.5){return "too many gaps";}

   #### the distance per site (0 or 1) between the paralogs
   #
   my @parasd = calculate_site_distance($dna_aln{$para11}, $dna_aln{$para12});
   #print "parasd\n@parasd[0..$#parasd]\n";

#print "this is species is ".$thissp."\n";

   #### infer conversion in paralogs of species $thissp
   #    as against the orthologs of species $i
   #    
   for(my $i=0; $i<$species_no; $i++)
   {
#print "i is ".$i." thissp is ".$thissp."\n";

       if($i eq $thissp){next;}
       my $para21=$$species[$i]."_1";
       my $para22=$$species[$i]."_2";
       #print "orthologs are ".$para21." and ".$para22."\n";
 
       #### the distance per site between the orthologs
      #
       my @ortho_obj = ("", "");
       if(isexist($para21, @seq_id)){$ortho_obj[0] = $dna_aln{$para21};}
       if(isexist($para22, @seq_id)){$ortho_obj[1] = $dna_aln{$para22};}
 
       my (@orthosd, @orthosd1, @orthosd2);
       print "ortholog nunmber is ".$#ortho_obj."\n";
       #print $ortho_obj[0] -> display_id()."  .... \n";
       #print $ortho_obj[1] -> display_id()."  ......\n";

       if($#ortho_obj eq 0)
       {
          @orthosd = calculate_site_distance($dna_aln{$para11}, $dna_aln{$para21});
       } 
       if($#ortho_obj eq 1)
       {
          @orthosd1 =  calculate_site_distance($dna_aln{$para11}, $dna_aln{$para21});
          @orthosd2 =  calculate_site_distance($dna_aln{$para12}, $dna_aln{$para22});
        #  print "orthosd1\n@orthosd1[0..$#orthosd1]\n";
        #  print "orthosd2\n@orthosd2[0..$#orthosd2]\n";

	#  print $dna_aln{$para11}->seq()."\n";
          @orthosd = merge_distance(\@orthosd1, \@orthosd2);
       }
       #print "orthosd\n@orthosd[0..$#orthosd]\n";

       #### the difference between the othosd and parasd
       #
       my @distdiff = distance_difference(\@orthosd, \@parasd);
       #print "distdiff\n@distdiff[0..$#distdiff]\n";

       #### detect the maximal score of distdiff according to Lin et al. 2006
       #
       
       my @score_pointer = scoring_rigor(@distdiff);
       my @score = @score_pointer[0..int($#score_pointer/2)];
       my @pointer=@score_pointer[int($#score_pointer/2)+1..$#score_pointer];
       #print "score\n@score[0..$#score]\n";
       #print "pointer\n@pointer[0..$#pointer]\n";
       #### retrieve the maximal value and its position
       #print "retrieving the maximal score ...\n";
       my @max_value_pos = my_max(@score);
       my $max = $max_value_pos[0];
       my $real_max = $max;
       my $end = $max_value_pos[1];

     while($max > 0)
     {
       print "the max is ".$max." at ".$end."\n";
       #### backtracing procedure to reveal the longest similar segments
       my $start = $end;
       my $pr = $pointer[$end];
       for(my $j=$end - 1; $j >= 0; $j--)
       {
          #print $j." ".$pr."\n";
          if($pr <= 0){last;}
          if($pr > 0)
          {
             $start = $pr;
             $pr = $pointer[$pr];
          }
       }
       my $conv_len = $end - $start + 1;
       if($conv_len < 10){$max = -1; next;}
       #print "supposed converted segment length ".$conv_len."\n";  

#print INFO $para11."|".$id_hash_r{$para11}."|".substr($dna_aln{$para11}->seq, $start, $end-$start)."\n";
#print INFO $para12."|".$id_hash_r{$para12}."|".substr($dna_aln{$para12}->seq, $start, $end-$start)."\n";
#print INFO $para21."|".$id_hash_r{$para21}."|".substr($dna_aln{$para21}->seq, $start, $end-$start)."\n";
#print INFO $para22."|".$id_hash_r{$para22}."|".substr($dna_aln{$para22}->seq, $start, $end-$start)."\n";


       #### calculating the binormal probability according to Lin et al. PNAS. 2006
       my ($para_dis, $ortho_dis) = calculate_DB(\@parasd, \@orthosd, $start, $end); ### D or B in Lin et al.
       my $p_value = binorm_p($conv_len, $para_dis, $ortho_dis/$conv_len);
     
       #### bootstrapping to evaluate the possibility of emerging D < B in the aligned homologs
       #my @b_value = bootstrap_local_1(\@parasd, \@orthosd, $end-$start+1, $para_dis, $ortho_dis);

       my $b_value = bootstrap_local_2(\@distdiff, $end-$start+1, $max);

       #### deciding who is donor, $para11, $para12
       #    as compared with $para21, $para22
       my $donor = "";
       my $fourexist = isexist($para11, @seq_id) * isexist($para12, @seq_id) *  isexist($para21, @seq_id) *  isexist($para22, @seq_id);
       if($fourexist eq 1)
       {
           my $psd = sum(@parasd[$start..$end]);
           my $osd1= sum(@orthosd1[$start..$end]);
           my $osd2= sum(@orthosd2[$start..$end]);

           #print "psd is ".$psd."osd1 is ".$osd1." osd2 is ".$osd2."\n";
           if($osd1 < $osd2 && $psd <= $osd1)
           {
               $donor = $para11;
               #print "possible donor is ".$para11."\n";
           }
           elsif($osd2 < $osd1 && $psd <= $osd2)
           {
               $donor = $para12
               #print "possible donor is ".$para12."\n";
           }
           else
           {
               $donor = "NOKNOWN";  
           }
       }

       #print INFO "species=".$$species[$thissp]." local conversion "." from=".$start." to=".$end." D=".$para_dis." B=".$ortho_dis." p-value=".$p_value." boot-value1=".$b_value[0]." boot-value2=".$b_value[1]." donor=".$donor."\n";
       #print INFO "species=".$$species[$thissp]." local conversion "." from=".$start." to=".$end." D=".$para_dis." B=".$ortho_dis." boot-value1=".$b_value[0]." boot-value2=".$b_value[1]." donor=".$donor."\n";
       if($b_value > 0.8)
       {
         print INFO "species=".$$species[$thissp]." local conversion "." from=".$start." to=".$end." D=".$para_dis." B=".$ortho_dis." boot-value2=".$b_value." donor=".$donor."\n";
       } 
       #### masking the supposed converted region
       #print "masking ...\n";
       my $j = $start +1;
       $score[$start] = -1;
       while($pointer[$j] > 0)
       {
          $pointer[$j] = -1;
          $score[$j] = -1;
          $j += 1;
       }

       #### retrieve the maximal value and its position
       #print "retrieving the secondary maximal score ...\n";
       @max_value_pos = my_max(@score);
       $max = $max_value_pos[0];
       $end = $max_value_pos[1];

     } #### ending while 
   }
   return $real_max;
}

sub bootstrap_local_2()
{
   my $rdistdiff = $_[0];
   my $conv_len = $_[1];
   my $max = $_[2];
   my $seq_len = $#$rdistdiff;

   my $random_conv = 0;
   
   for(my $i=0; $i<$BOOTNUM; $i++)
   {
       my @random_num=create_rand($seq_len, $seq_len);

       my @distdiffrand = ();
       for(my $j=0; $j<$seq_len; $j++)
       {$distdiff_rand[$j] = $$rdistdiff[$random_num[$j]];}

       my @score_pointer = scoring_rigor(@distdiff_rand);
       my @score = @score_pointer[0..int($#score_pointer/2)];
       my @pointer=@score_pointer[int($#score_pointer/2)+1..$#score_pointer];

       my @max_value_pos = my_max(@score);
       my $max_rand = $max_value_pos[0];
       my $end = $max_value_pos[1];

       my $start = $end;
       my $pr = $pointer[$end];
       for(my $j=$end - 1; $j >= 0; $j--)
       {
          if($pr <= 0){last;}
          if($pr > 0)
          {
             $start = $pr;
             $pr = $pointer[$pr];
          }
       }
       my $conv_len_rand = $end - $start + 1;     
      # print "max ".$max." max_rand ".$max_rand." conv_len_rand ".$conv_len_rand." conv_len ".$conv_len."\n";
 
       if($max_rand < $max && $conv_len_rand < $conv_len)
       {$random_conv ++;}
   }
   return ($random_conv/$BOOTNUM);
}

sub bootstrap_local_1()
{
   my $parasd = $_[0];
   my $orthosd= $_[1];
   my $conv_len = $_[2];
   my $paradis = $_[3];
   my $orthodis= $_[4];
   my $seq_len = $#$parasd;

   my $random_conv1 = 0;
   my $random_conv2 = 0;
   for(my $i=0; $i<$BOOTNUM; $i++)
   {
       my @random_num=create_rand($conv_len, $seq_len);   
       my ($D, $B) = (0,0);
       for(my $j=0; $j<$conv_len; $j++)
       {
          if($$parasd[$random_num[$j]] eq -1 || $$orthosd[$random_num[$j]] eq -1)
          {next;}
         
          $D += $$parasd[$random_num[$j]];
          $B += $$orthosd[$random_num[$j]]; 
          #print $random_num[$j]." ".$$parasd[$random_num[$j]]." ".$$orthosd[$random_num[$j]]." ".$D." ".$B." ";
       }
       #print "D B = ".$D." ".$B."\n";
       if($D < $B){$random_conv1 ++; if($B - $D >= $orthodis - $paradis){$random_conv2 ++;}}
   }
   return (1-$random_conv1/$BOOTNUM, 1-$random_conv2/$BOOTNUM);
}

sub calculate_DB()
{
       my $parasd = $_[0];
       my $orthosd= $_[1];
       my $start  = $_[2];
       my $end    = $_[3];

       my $para_dis = 0;
       my $ortho_dis= 0;
       for(my $i=$start; $i<=$end; $i++)
       {
          $para_dis += $$parasd[$i];
          $ortho_dis+= $$orthosd[$i];
       }
       return($para_dis, $ortho_dis);
}
sub scoring_rigor()
{
       my @distdiff = @_;
       my (@score, @pointer);
       $score[0] = $distdiff[0];
       $pointer[0]=-1;

       for(my $i=1; $i<=$#distdiff; $i++)
       {
          if($distdiff[$i] eq -2)
          {
             if($score[$i-1] > 0)
             {$score[$i] = $score[$i-1] + $gappenalty; $pointer[$i] = $i - 1;}
             else
             {$score[$i] = $gappenalty; $pointer[$i] = -1;}
          }
          else
          {
             if($score[$i-1] > 0)
             {$score[$i] = $score[$i-1] + $distdiff[$i]; $pointer[$i] = $i - 1;}
             else
             {$score[$i] = $distdiff[$i]; $pointer[$i] = -1;}
          }
       }

       #print "score\n@score[0..$#score]\n";
       return(@score, @pointer);
}

sub distance_difference()
{
   my $orthosd = $_[0];
   my $parasd= $_[1];
   
   my @distdiff; ## = orthosd - parasd
   for(my $i=0; $i<=$#$orthosd; $i++)
   {
      if($$orthosd[$i] eq -1 || $$parasd[$i] eq -1){$distdiff[$i] = -2;}
      else{$distdiff[$i] = $$orthosd[$i] - $$parasd[$i];}
   }
   return @distdiff;
}

sub calculate_site_distance()
{
   my $seqobj1 = $_[0];
   my $seqobj2 = $_[1];
   my @seq1 = split(//, $seqobj1 -> seq());
   my @seq2 = split(//, $seqobj2 -> seq());

   my @distance;
   for(my $i=0; $i<=$#seq1; $i++)
   {
      if($seq1[$i] eq $seq2[$i]){$distance[$i] = 0;}
      elsif($seq1[$i] eq "-" || $seq2[$i] eq "-"){$distance[$i] = -1;}
      else{$distance[$i] = 1;}
   }
   return @distance;
}

sub merge_distance()
{
   my $rdis1 = $_[0];
   my $rdis2 = $_[1];
   my @distance;
   for(my $i=0; $i<=$#$rdis1; $i++)
   {
      if($$rdis1[$i] eq 0)
      {
          if($$rdis2[$i] eq 0)
         {
            $distance[$i] = 0;
         }
         elsif($$rdis2[$i] eq 1)
         {
            $distance[$i] = 1/2;
         }
          else
         {
            $distance[$i] = -1;
         }
      }
      elsif($$rdis1[$i] eq 1)
      {
         if($$rdis2[$i] eq 0)
         {
            $distance[$i] = 1/2;
         }
         elsif($$rdis2[$i] eq 1)
         {
            $distance[$i] = 1;
         }
          else
         {
            $distance[$i] = -1;
         }
      }
      else
      {
         $distance[$i] = -1;
      }
   }
   return @distance;
}
sub global_conv()
{
   my $species = $_[0];
   my $dna_aln = $_[1];

   my %Distance = get_distance2($dna_aln);  
foreach my $key(sort(keys(%Distance)))
{
   print $key." ".$Distance{$key}."\n";
}
 
   my @pair_id = keys(%Distance);

   my @comment = ();
   my $species_no = $#$species + 1;
   
   for(my $i=0; $i<$species_no; $i++)
   {
      #### get the paralogs in the species $i
      #    if not found, indicating a missing paralog
      my $para11=$$species[$i]."_1";
      my $para12=$$species[$i]."_2";
      my $para = join("-", $para11, $para12);
print "para is $para\n";

      if(isexist($para, @pair_id) eq 0){$comment[$i] = "singleton in species ".$$species[$i]; next;}
      if($Distance{$para} eq ""){$comment[$i] = "too divergent paralogs in species ".$$species[$i]; next;} 

      #### if($isexist eq 1)
      $comment[$i] = "";
      for(my $j=0; $j<$species_no; $j++)
      {
         if($i eq $j){next;}
         my $para21=$$species[$j]."_1";
         my $para22=$$species[$j]."_2";
         my $ortho1=join("-", sort($para11, $para21));
         my $ortho2=join("-", sort($para12, $para22));
         if(isexist($ortho1, @pair_id) eq 1 && isexist($ortho2, @pair_id) eq 1)
         {
            if($Distance{$para} < $Distance{$ortho1} && $Distance{$ortho1} < $Distance{$ortho2} )
            {
               my $p_value = bootstrap_global($para, $ortho1, $dna_aln);
               $comment[$i] .= "global conversion=TYPE1 ".$para."=".$Distance{$para}." ".$ortho1."=".$Distance{$ortho1}." ".$ortho2."=".$Distance{$ortho2}." boot=".$p_value." donor=".$para11;
            } 
            elsif($Distance{$para} < $Distance{$ortho2} && $Distance{$ortho2} < $Distance{$ortho1} )
            {
               my $p_value = bootstrap_global($para, $ortho2, $dna_aln);
               $comment[$i] .= "global conversion=TYPE2 ".$para."=".$Distance{$para}." ".$ortho1."=".$Distance{$ortho1}." ".$ortho2."=".$Distance{$ortho2}." boot=".$p_value." donor=".$para12;
            }
            elsif($Distance{$para} < $Distance{$ortho1} && $Distance{$ortho2} == $Distance{$ortho1} )
            {
               my $p_value = bootstrap_global($para, $ortho1, $dna_aln);
               $comment[$i] .= "global conversion=TYPE3 ".$para."=".$Distance{$para}." ".$ortho1."=".$Distance{$ortho1}." ".$ortho2."=".$Distance{$ortho2}." boot=".$p_value;
            }
         }
         elsif(isexist($ortho1, @pair_id) eq 1) 
         {
            if($Distance{$para} < $Distance{$ortho1})
            {
               my $p_value = bootstrap_global($para, $ortho1, $dna_aln);
               $comment[$i] .= "global conversion=TYPE4 ".$para."=".$Distance{$para}." ".$ortho1."=".$Distance{$ortho1}." ".$ortho2."=".$Distance{$ortho2}." boot=".$p_value;
            }
         }
         elsif(isexist($ortho2, @pair_id) eq 1)
         {
            if($Distance{$para} < $Distance{$ortho2})
            {
               my $p_value = bootstrap_global($para, $ortho2, $dna_aln);
               $comment[$i] .= "global conversion=TYPE5 ".$para."=".$Distance{$para}." ".$ortho1."=".$Distance{$ortho1}." ".$ortho2."=".$Distance{$ortho2}." boot=".$p_value;
            }
         }
      }
      if($comment[$i] eq ""){$comment[$i] = "no global conversion";}
   }
   
   return @comment;
}

sub bootstrap_global()
{
   my $para = $_[0];
   my $ortho= $_[1];
   my $dna_aln = $_[2];
   my $dna_aln_rand = $dna_aln;

   my $align_len = $dna_aln -> length();
   my %align_seq;
   my %align_seq_rand;
   my @tmpid;
   my $i = 0;
   foreach my $seq_obj ($dna_aln -> each_seq())
   {
      $align_seq{$seq_obj -> display_id()} = $seq_obj -> seq();
#      print "\n############\n".$seq_obj -> seq()."\n";
      $align_seq_rand{$seq_obj -> display_id()} = "";
      $tmpid[$i] = $seq_obj -> display_id();
      $i = $i + 1;
   } 
   my @id = sort(@tmpid);   

   my $count = 0;
   for(my $r=1; $r<=$BOOTNUM; $r++)
   {
#print "\n----------->>>".$r."\n";

      my @random_num = create_rand($align_len/3, $align_len/3);
       
       for(my $i=0; $i<=$#random_num; $i++)
       {
          if($i eq 0)
          {
             foreach my $so ($dna_aln -> each_seq())
             {$align_seq_rand{$so -> display_id()} = "";}
          } 
          for(my $j=0; $j<=$#id; $j++)
          {
             $align_seq_rand{$id[$j]} .= substr($align_seq{$id[$j]}, ($random_num[$i]-1) * 3, 3);
 #            print $random_num[$i]." ".substr($align_seq{$id[$j]}, ($random_num[$i]-1) * 3, 3)." ";
          }
       }
       foreach my $so ($dna_aln_rand -> each_seq())
       {
          my $len = length($align_seq_rand{$so -> display_id()});
 #         print "length ...".$len."\n";
          $so -> seq($align_seq_rand{$so -> display_id()});
 #         print "\n===========\n".$align_seq_rand{$so -> display_id()}."\n";
       }
 #      print "count ".$count."\n";
       my $aln_len = $dna_aln_rand -> length();
 #      print "the alignment length is ".$aln_len."\n";
       
 #      my %Ks = get_Ks($dna_aln_rand);
       my %Distance = get_distance2($dna_aln_rand);
 
 #      if($Ks{$para} < $Ks{$ortho}){$count++;}
       if($Distance{$para} < $Distance{$ortho}){$count++;}
   }
   return $count/$BOOTNUM;
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

sub isexist()
{
   my $id = $_[0];
   my @id = $_[1..$#_];
   my $isexist = 0;
   for(my $j=0; $j<=$#id; $j++)
   {
      if($id eq $id[$j]){$isexist = 1;}
   }
   return $isexist;
}
sub check_gap()
{
   my $seqobj1 = $_[0];
   my $seqobj2 = $_[1];
   my $seq1 = $seqobj1 -> seq();
   my $seq2 = $seqobj2 -> seq();

   my @seq1 = split(//, $seq1);
   my @seq2 = split(//, $seq2);

   my $num = 0;
   for(my $i = 0 ; $i <= $#seq1; $i ++)
   {
      if($seq1[$i] eq "-"){$num ++;}
      if($seq2[$i] eq "-"){$num ++;} 
   }
   my $gapratio = $num/(2*$#seq1 + 2);
   return $gapratio;
}

sub binorm_p()
{
   my $n = $_[0];
   my $x = $_[1];
   my $p = $_[2];
   if($p eq 1){$p = 0.999999;}
   if($p eq 0){$p = 0.000001;}
   #print "\n".$n." ".$x." ".$p."\n";
   my $pv = 0;
   for(my $k = 0; $k <= $x; $k ++)
   {
      my $value = exp( logsum($n)-logsum($k)-logsum($n-$k) + log($p)*$k + log(1-$p)*($n-$k));
      $pv += $value;
    #  print $k." ".$pv."\n";
   }
   return $pv;
}

sub logsum()
{
   my $n = $_[0];
   my $sum = 0;
   
   if($n > 1)
   {
      for(my $i = 1; $i <= $n; $i++ )
      {
         $sum += log($i);
      }
   }
   return $sum;
}

sub my_max()
{
   my $max_value = 0;
   my $position;
   for(my $i=0; $i<=$#_; $i++)
   {
       if($_[$i] > $max_value)
       {$max_value = $_[$i]; $position = $i;}
   }
   return($max_value, $position);
}

sub min2()
{
   if($_[0] < $_[1]){return $_[0];}else{return $_[1];}
}
###### not used at present
sub get_Ks()
{
   my $stats = new SeqStatistics;
   my $result = $stats->calc_all_KaKs_pairs($_[0]); #### the input is dna align object: $dna_aln
   my (%Ks);
   for my $an (@$result)
   {
#       print "comparing ". $an->{'Seq1'}." and ". $an->{'Seq2'}. " \n";
       my $order = 0;
       my @tmp = sort($an->{'Seq1'}, $an->{'Seq2'});
       my $pair_id = join("-", @tmp);
  
       for (sort keys %$an )
       {
          next if /Seq/;
#          printf("%-9s %.4f \n",$_ , $an->{$_});
          if($_ eq "D_s")
          {
            $Ks{$pair_id} = $an->{$_};
          }
       }
#       print INFO $pair_id." ".$Ks{$pair_id}."\n";
   }
   return(%Ks);
}

sub get_distance()
{
   my $stats = new SeqStatistics;
   my $result = $stats->calc_all_KaKs_pairs($_[0]); #### the input is dna align object: $dna_aln
   my (%distance);
   for my $an (@$result)
   {
       print "comparing ". $an->{'Seq1'}." and ". $an->{'Seq2'}. " \n";
       my $order = 0;
       my @tmp = sort($an->{'Seq1'}, $an->{'Seq2'});
       my $pair_id = join("-", @tmp);
       my ($S_d, $N_d, $S, $N) = (0, 0, 0, 0); 

       for (sort keys %$an )
       {
          next if /Seq/;
          printf("%-9s %.4f \n",$_ , $an->{$_});
          if($_ eq "S_d"){$S_d = $an->{$_};}
          if($_ eq "N_d"){$N_d = $an->{$_};}
          if($_ eq "S"){$S = $an->{$_};}
          if($_ eq "N"){$N = $an->{$_};}
       }
       $distance{$pair_id} = ($S_d + $N_d)/($S + $N);
       #print INFO $pair_id." ".$distance{$pair_id}."\n";
   }
   return(%distance);
}

sub get_distance2()
{
   my $aln = $_[0]; #### the input is dna align object: $dna_aln
   my (%distance);
   
   my @seq = ();
   my $seqno = 0;
   foreach my $seq_obj ($aln -> each_seq())
   {
      $seq[$seqno] = $seq_obj;
      $seqno ++;
   }

   for(my $i=0; $i<=$#seq; $i++)
   {
      for(my $j=$i+1; $j<=$#seq; $j++)
      {
         my @seq1 = split(//, $seq[$i] -> seq());
         my @seq2 = split(//, $seq[$j] -> seq());
         my $id1 = $seq[$i] -> display_id();
         my $id2 = $seq[$j] -> display_id();

         my $ndiff = 0;
         my $npos  = 0;
         for(my $k=0; $k<=$#seq1; $k++)
         {
            if($seq1[$k] eq "-" || $seq2[$k] eq "-"){next;}
            $npos ++;
            if($seq1[$k] ne $seq2[$k]){$ndiff ++;}
         }
         if($npos >0){$distance{$id1."-".$id2} = $ndiff/$npos; $distance{$id2."-".$id1} = $ndiff/$npos;}
         else{$distance{$id1."-".$id2} = -1; $distance{$id2."-".$id1} = -1;}

print "distance $id1 $id2".$distance{$id1."-".$id2}."\n";
      }
   }

   return(%distance);
}


sub isSynonymous
{
    #print "are the two codons synonymous\n";
    my($codon1)=$_[0];
    my($codon2)=$_[1];
    
    if($codonhash{$codon1} eq $codonhash{$codon2})
    {return 1;}
    else
    {return 0;}
}

sub is4FoldDegenerate
{
    #print "are the two codons 4 fold degenerate\n";
    my($codon1)=$_[0];
    my($codon2)=$_[1];
    
    if($codonhash{$codon1} eq $codonhash{$codon2})
    {
       my $aminoacid = $codonhash{$codon1};
       if($aminoacid eq "Val" || $aminoacid eq "Pro" || $aminoacid eq "Thr" ||$aminoacid eq "Ala"
        ||$aminoacid eq "Gly")
       {return 1;}
       elsif($aminoacid eq "Leu" && $codon1 ne "TTA" && $codon1 ne "TTG" && $codon2 ne "TTA" && $codon2 ne "TTG")
       {return 1;}
       elsif($aminoacid eq "Ser" && $codon1 ne "AGT" && $codon1 ne "AGC" && $codon2 ne "AGT" && $codon2 ne "AGC")
       {return 1;}
       elsif($aminoacid eq "Arg" && $codon1 ne "AGA" && $codon1 ne "AGG" && $codon2 ne "AGA" && $codon2
 ne "AGG")
       {return 1;}
       else
       {return 0;}
    }
    else
    {return 0;}
}

sub printArrInFile
{	
        my $fh = $_[0];
	my @array = @_;
	my $i;
	
	for($i=1; $i<=$#array; $i++)
	{
		print $fh $array[$i]." ";
	}
}

sub sumArr
{
	my @array = @_;
	my $i;
        my $sum=0;
	for($i=0; $i<=$#array; $i++)
	{
		$sum+=$array[$i];
	}
	return $sum;
}

sub printArr
{
	my @array = @_;
	my $i;
	for($i=0; $i<=$#array; $i++)
	{
		print $array[$i]." ";
	}
	print "\n";
}

sub chisq_value
{
        my (@array1,@array2);

        $array1[0] = $_[0];
        $array1[1] = $_[1];
        $array2[0] = $_[2];
        $array2[1] = $_[3];

        my @arr2D;
        $arr2D[0] = \@array1;
        $arr2D[1] = \@array2;
        
        my ($i, $j);
        my $item;
        my $chisq_value=0;
        for($i=0; $i<=1; $i++)
        {
           for($j=0; $j<=1; $j++)
           {  #print "arr2D ".$arr2D[$i][$j]." ".$item." ";
              $item = sum(@array1,@array2)*($arr2D[$i][$j]-sum($arr2D[$i][0],$arr2D[$i][1])*sum($arr2D[0][$j],$arr2D[1][$j])/sum(@array1,@array2))**2/(sum($arr2D[$i][0],$arr2D[$i][1])*sum($arr2D[0][$j],$arr2D[1][$j]));
              $chisq_value = $chisq_value + $item;
           }
        }
#        print " sums ".sum(@array1,@array2)." ".sum(@array2)." ".sum(@array1)." ";
#        print "chisq_value ".$chisq_value."\n";
        return $chisq_value;
}

sub chisq_pvalue
{
    my $r_shell = $_[0];
    my @chisq_value = @_;
    #print "$chisq_value[0], @chisq_value[1..$#chisq_value]";
    open R_SHELL, ">$r_shell" or die "Cannot create r_shell file\n";
    print R_SHELL "#!/bin/bash","\n";
    print R_SHELL "cat <<EOF | R --vanilla --slave", "\n";
    #print R_SHELL "pv <- pchisq(".$chisq_value.",1,lower.tail = F);","\n";
    print R_SHELL join("", "pv <- pchisq(", join("", "c(",join(",", @chisq_value[1..$#chisq_value]), ")"), ",1,lower.tail = F);"), "\n";
    print R_SHELL 'cat(pv,"\n");',"\n";
    print R_SHELL 'EOF',"\n";
    close R_SHELL;
    
    my $result = `bash $r_shell`;
    return $result;
}

sub significantregions
{
    my $rchisqvalue = $_[0];
    my $winNum = $_[1];
 
    my %significantregions;

    my $i;
    my $signstart = -1;

    for($i = 0; $i < $winNum; $i ++)
    {
        if($$rchisqvalue[$i] > 3.84)
        {
            if($signstart eq -1)
            {$signstart = $i;}
        }
        else
        {
            if($signstart ne -1)
            {
               my $key = "win";
               my $chisqv = 0;
               for(my $j = $signstart; $j <= $i-1; $j ++)
               {
                    $key .= "-".$j;
                    if($$rchisqvalue[$j] > $chisqv){$chisqv = $$rchisqvalue[$j];}
               } ### fetch the largest the chisqvalue as the extended window value
               $significantregions{$key} = $chisqv;
               $signstart = -1;
               my @arr = %significantregions;
               #printArrInFile($fh, "\n", $#arr, @arr);
            }
        }
    }
    return %significantregions;
}

sub getcodon()
{	
	$codonhash{"TTT"}="Phe";$codonhash{"TCT"}="Ser";$codonhash{"TAT"}="Tyr";$codonhash{"TGT"}="Cys";
	$codonhash{"TTC"}="Phe";$codonhash{"TCC"}="Ser";$codonhash{"TAC"}="Tyr";$codonhash{"TGC"}="Cys";
	$codonhash{"TTA"}="Leu";$codonhash{"TCA"}="Ser";$codonhash{"TAA"}="Ter";$codonhash{"TGA"}="Ter";
	$codonhash{"TTG"}="Leu";$codonhash{"TCG"}="Ser";$codonhash{"TAG"}="Ter";$codonhash{"TGG"}="Trp";

	$codonhash{"CTT"}="Leu";$codonhash{"CCT"}="Pro";$codonhash{"CAT"}="His";$codonhash{"CGT"}="Arg";
	$codonhash{"CTC"}="Leu";$codonhash{"CCC"}="Pro";$codonhash{"CAC"}="His";$codonhash{"CGC"}="Arg";
	$codonhash{"CTA"}="Leu";$codonhash{"CCA"}="Pro";$codonhash{"CAA"}="Gln";$codonhash{"CGA"}="Arg";
	$codonhash{"CTG"}="Leu";$codonhash{"CCG"}="Pro";$codonhash{"CAG"}="Gln";$codonhash{"CGG"}="Arg";

	$codonhash{"ATT"}="Ile";$codonhash{"ACT"}="Thr";$codonhash{"AAT"}="Asn";$codonhash{"AGT"}="Ser";
	$codonhash{"ATC"}="Ile";$codonhash{"ACC"}="Thr";$codonhash{"AAC"}="Asn";$codonhash{"AGC"}="Ser";
	$codonhash{"ATA"}="Ile";$codonhash{"ACA"}="Thr";$codonhash{"AAA"}="Lys";$codonhash{"AGA"}="Arg";
	$codonhash{"ATG"}="Met";$codonhash{"ACG"}="Thr";$codonhash{"AAG"}="Lys";$codonhash{"AGG"}="Arg";

	$codonhash{"GTT"}="Val";$codonhash{"GCT"}="Ala";$codonhash{"GAT"}="Asp";$codonhash{"GGT"}="Gly";
	$codonhash{"GTC"}="Val";$codonhash{"GCC"}="Ala";$codonhash{"GAC"}="Asp";$codonhash{"GGC"}="Gly";
	$codonhash{"GTA"}="Val";$codonhash{"GCA"}="Ala";$codonhash{"GAA"}="Glu";$codonhash{"GGA"}="Gly";
	$codonhash{"GTG"}="Val";$codonhash{"GCG"}="Ala";$codonhash{"GAG"}="Glu";$codonhash{"GGG"}="Gly";
}

sub sum
{
      my $sum = 0;
      for(my $i=0; $i<=$#_; $i++)
      {                          
          $sum = $sum + $_[$i];
      }           
      return $sum;
}


sub printTime()
{	
	my $t = time();

       ($sec,$min,$hour,$dom,$mon,$year,$wday,$yday,$isdst)=localtime($t);

	printf("time : %2.2d:%2.2d:%2.2d\n",$hour,$min,$sec);
}	
		
sub max()
{ 
   my $max = 0;
   for(my $i = 0; $i <= $#_; $i ++)
   {
      if($_[$i] > $max)
      {$max = $_[$i];}
   }
   return $max;
}    




