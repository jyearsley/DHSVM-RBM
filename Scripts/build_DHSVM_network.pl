# Perl script to build the network file neeeded by RBM
# Also creates a file <Project>.map which is a mapping of
# the segment numbers used by RBM to the segment numbers
# output by DHSVM.  This can be useful when plotting output
# from the stream temperature model.
#
#!/usr/bin/perl

use Date::Calc qw(:all);
#
# Default value of ndelta, the number of segments in each stream cell
#
$n_default = 2;
#
print "NOTE: This is an updated version of build_DHSVM_network\n";
print "      that builds the network file for the version of RBM\n";
print "      that incorporates variable Mohseni and Leopold parameters.\n";
print "\n";
print "      The number of reach segments, ndelta is set at 2. Specific changes,\n";
print "      where necessary, must be made in the *.net file\n";
print "\n";
print "      Input ProjectName for topology file: <Project>\n";
print "      This script will build a network file: <Project>.net and a <Project>.segmap file\n";
chomp($project=<STDIN>);
#
$dhsvm=$project.'.dir';
print "$dhsvm\n";
open DHSVM, "$dhsvm" or die "Cannot open $dhsvm\n"; 
$dhsvm_net=$project.'.net';
open NET, ">$dhsvm_net" or die "Cannot open NET file\n";
$forcing_file=$project.'.forcing';
#
# Prepare the header for the nework file:<ProjectName>.net";
#
print NET "Temperature simulation for the project: $project\n";
# 
print NET "$forcing_file\n";
#
#
$node_map=$project.'.segmap';
open MAP, ">$node_map" or die "Cannot open MAP file\n";
#
# Read the DHSVM topology file
#
while ($zz=<DHSVM>) {
  $n_seg++;
  ($seg_id,$ds_seg,$x_seg,$elev,$avg_az,@tribs)=split/\s+/,$zz;
  $x_meters[$seg_id]=$x_seg;
  $dist_feet[$seg_id]=$x_seg*3.2808;
  $seg_no[$n_seg]=$seg_id;
  $ds_seg_no[$seg_id]=$ds_seg;
  $elev[$seg_id]=$elev;
  $n_tribs=$#tribs+1;
  $n_us_tribs[$seg_id]=$n_tribs;
  for $n (1..$n_tribs) {
    $us_tribs[$seg_id][$n]=$tribs[$n-1];
  }
#
# If the segment ID is less than 0 this is the outlet segment
#
  if ($ds_seg_no[$seg_id] < 0) {
     $ds_seg_no[$seg_id]=0;
     $outlet_seg=0;
#
# Add a dummy segment to handle the possibility of a tributary to the
# actual outlet segment 
#
     $avg_az[0]=0;
     $us_tribs[0][1]=$seg_id;
     $ds_seg_no[0]=-1;
     $n_us_tribs[0]=1;
  }
#
# Test for headwater cell
#
   if ($tribs[0] eq undef && $seg_id > 0) {
     $n_head++;
     $no_seg_head[$n_head]=$seg_id;
     $ds_head_seg[$n_head]=$ds_seg_no[$seg_id];      
   }
}
#
# Print the number of headwaters to DHSVM
#
$total_head=$n_head;
$total_segs=$n_seg;
print NET "$n_head Headwaters $n_seg Stream Segments in Project: $project\n";
print MAP "$n_head $n_seg\n";
#
# First find the Main Stem
#
$max_segs=-99999;
#
for $nh (1..$n_head) {
  $n_segs=0;
  $start_seg=$no_seg_head[$nh];
  while ($start_seg > 0) {
    $n_hw_seg[$start_seg]++;
    $hw_seg[$start_seg][$n_hw_seg[$start_seg]]=$nh;
    $n_segs++;
    $start_seg=$ds_seg_no[$start_seg];
    if ($n_segs > $max_segs) {
      $last_seg=$start_seg;
      $max_segs=$n_segs;
      $main_stem=$nh;
      $n_seg_hw[$nh]=$n_segs
    }
  }
}
print "MAX SEGS $max_segs\n";
$dmmy_seg=$no_seg_head[$main_stem];
for $ns (1..$max_segs) {
  print "$ns\n";
  $branch_seg[1][1][$ns]=$dmmy_seg;
  $path_seg[$main_stem][$ns]=$dmmy_seg;
  $dmmy_seg=$ds_seg_no[$dmmy_seg];
}
#
# Stream order of the Main stem is = 1.  Stream order increases as stream size decreases
# the reverse of the Strahler classification
#
$stream_order=1;
$stream_order[$main_stem]=$stream_order;
#
# Main Stem has only one (1) branch
#
$no_branches[$stream_order]=1;
$branch_headwaters[$no_branches[$stream_order]][$stream_order]=$main_stem;
$conf_seg[$main_stem]=$outlet_seg;
$n_seg[$main_stem][1]=$max_segs;
#
# Starting segment
#
$n_head_rem=$n_head-1;
#
while ($n_head_rem != 0 && $stream_order < 10) {
  $stream_order++;
  $str_m=$stream_order-1;
BRANCH:  for $nb (1..$no_branches[$str_m]) {
    $nh=$branch_headwaters[$nb][$str_m];
    $n_seghw=$n_seg_hw[$nh];
    for $n_s (1..$n_seghw) {
      $nn_s = $n_seghw-$n_s+1;
      $nn_sm=$nn_s-1;
      $start_seg=$branch_seg[$str_m][$nb][$nn_s];
      $us_start_seg=$branch_seg[$str_m][$nb][$nn_sm];
#
#  Following the main stem upstream.  At each node, check to see how many
#  streams are contributing at this point
#
      for $n_us (1..$n_us_tribs[$start_seg]) {
        $trib_seg=$us_tribs[$start_seg][$n_us];
        if ($trib_seg ne $us_start_seg) {
          $no_branches[$stream_order]++;
          $n_head=$n_hw_seg[$trib_seg];
#
# Find the "main stems" of the next higher order (smaller) streams
# 
          $max_segs=-99999.;
          for $n_hd (1..$n_head) {
            $no_hw=$hw_seg[$trib_seg][$n_hd];
#
#  Going downstream to find maximum number of cells in this branch
#  
            $test_seg=$no_seg_head[$no_hw];
            $n_segs=0;        
            while ($test_seg != $start_seg) {
              $n_segs++;
              $test_seg=$ds_seg_no[$test_seg]; 
              if ($n_segs > $max_segs) {
                $main_stem=$no_hw;
                $max_segs=$n_segs;
              }
            }
# End of examining headwaters loop     
          }    
          $n_seg_hw[$main_stem]=$max_segs;
          $branch_headwaters[$no_branches[$stream_order]][$stream_order]=$main_stem;
#
#  Identify the segments in this branch
#
          $dmmy_seg=$no_seg_head[$main_stem];
          for $ns (1..$max_segs) {
            $branch_seg[$stream_order][$no_branches[$stream_order]][$ns]=$dmmy_seg;
            $path_seg[$main_stem][$ns]=$dmmy_seg;
            $dmmy_seg=$ds_seg_no[$dmmy_seg];
          }
# End of IF statement excluding the upstream segment
        }
# End of upstream tributaries loop
      }
# End of number of segments loop for a specific branch 
    }
# End of branch loop
  }
# End of main while loop
}
#
#  Counting the total number of segments after adding a segment
#  for those branches that had only one (1)
#
$seg_seq=0;
for $n_stordr (1..$stream_order) {
  $nn_st_m=$stream_order-$n_stordr+1;
  for $nb (1..$no_branches[$nn_st_m]) {
    $nh=$branch_headwaters[$nb][$nn_st_m];
    $nn_seg=$n_seg_hw[$nh];
    for $ns (1..$nn_seg) {
      $seg_seq++;
      $path=$path_seg[$nh][$ns];
      $trib_ndx[$ds_seg_no[$path]]=$seg_seq+1;
    }
  }
}
$seg_seq=0;
for $n_stordr (1..$stream_order) {
  $nn_st_m=$stream_order-$n_stordr+1;
  for $nb (1..$no_branches[$nn_st_m]) {
    $nh=$branch_headwaters[$nb][$nn_st_m];
    $nn_seg=$n_seg_hw[$nh];
    $r_feet0[$nn_seg]=0.0;
    for $ns (1..$nn_seg) {
      $ni=$nn_seg-$ns+1;
      $r_feet1[$ni]=$r_feet0[$ni]+$dist_feet[$path_seg[$nh][$ni]];
      $r_feet0[$ni-1]=$r_feet1[$ni];
    }
    for $ns (1..$nn_seg) {
      $seg_seq++;
      if ($ns ==1) {
        $path=$path_seg[$nh][$nn_seg];
        $path2=$ds_seg_no[$path];
        printf NET "#_Segments %5d Headwaters %5d TribCell   %6d %6d %6d\n"
                   ,$n_seg_hw[$nh],$nh,$trib_ndx[$ds_seg_no[$path]],$path,$path2;
      }  
      printf NET "Seq   %5d Path %5d X_0  %10.3f X_1  %10.3f Elevation %8.2f %5d\n",$seg_seq, $path_seg[$nh][$ns]
                 ,$r_feet0[$ns],$r_feet1[$ns],$elev[$path_seg[$nh][$ns]],$n_default;
      printf MAP " Sequence %5d Path %5d\n",$seg_seq, $path_seg[$nh][$ns]; 
    }
    $nsm=$n_seg_hw[$nh];
  }
}
