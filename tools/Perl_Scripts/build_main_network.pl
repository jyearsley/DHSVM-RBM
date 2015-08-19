#
# Perl script to build the network file neeeded by RBM
# for the Connectcut River main stem (a special case)
#
# Also creates a file <Project>.map which is a mapping of
# the segment numbers used by RBM to the segment numbers
# output by DHSVM.  This can be useful when plotting output
# from the stream temperature model.
#
#!/usr/bin/perl

use Date::Calc qw(:all);
#
print "Input ProjectName \n";
print "Requires files: <ProjectName>,Segments - file defining segments\n";
print "              : <ProjectName>.Main - mapping RBM/DHSVM segments for main stem\n";
print "              : <ProjectName>.Trib - mapping RBM/DHSVM segments for tributaries\n";
print "This script will build a network file: <ProjectName>.net\n";
chomp($project=<STDIN>);
#
$segment_file=$project.'.Segments';
open SEGMENTS, "$segment_file" or die "Cannot open $segment_file\n";
$main_file=$project.'.Main_Stem';
open MAIN, "$main_file" or die "Cannot open $main_file\n"; 
$trib_file=$project.'.Trib';
open TRIB, "$trib_file" or die "Cannot open $trib_file\n"; 
$conn_only=$project.'.Only_Table';
open ONLY, "$conn_only" or die "Cannot open $conn_only\n"; 
$dhsvm_net=$project.'.net';
open NET, ">$dhsvm_net" or die "Cannot open NET file\n";
#
# Open SOURCE file
#
$source_net=$project.'.Source';
open SOURCE, ">$source_net" or die "Cannot open SOURCE file\n";
#
# Prepare the header for the nework file:<ProjectName>.net";
#
print "Input parameters for initial (headwaters) temperatures \n";
print "and hydraulic parameters, depth and stream speed \n";
print "Input values are separated by commas\n";
print "    \n";
print "Input parameter, <smooth>, for smoothing daily air temperatures \n";
#chomp($smooth=<STDIN>);
#
print "Input Mohseni nonlinear regregression parameters \n"; 
print "<alpha>,<beta>,<gamma>,<mu>\n";
#chomp($_=<STDIN>);
#($alpha,$beta,$gamma,$mu)=split/,/;
#
print "Leopold coefficients,<U_a>, <U_B> + u_min for stream speed, u \n";
print "where u = <U_a>*Q**<U_b> and u_min is a threshold speed (English units \n";
#chomp($_=<STDIN>);
#($U_a,$U_b,$u_min)=split/,/;
#
print "Leopold coefficients,<D_a>, <D_b> + d_min for stream depth, D \n";
print "where D = <D_a>*Q**<D_b> and D_min is a threshold depth (English units) \n";
#chomp($_=<STDIN>);
#($D_a,$D_b,$D_min)=split/,/;
#
print NET "Temperature simulation for the project: $project\n";
#
$forcing_file=$project.'forcing'; 
#
print NET "$forcing_file\n";
#
# Dummy Mohseni parameters and Leopold coefficients
#
$alpha=20; $beta=16; $gamma=0.2; $mu=0.1; $smooth=0.1;
#
printf NET "%5.1f %4.1f %4.2f %4.1f %4.2f Mohseni parameters\n",
           $alpha,$beta,$gamma,$mu,$smooth;
#
$D_a=0.4; $D_b=0.4; $D_min==0.5;
#
printf NET "%5.2f %4.2f %4.1f           Leopold coefficients for depth\n",
            $D_a,$D_b,$D_min;
#
$U_a=0.4; $U_b=0.4; $U_min=0.5;
#
printf NET "%5.2f %4.2f %4.1f           Leopold coefficients for stream speed\n",
            $U_a,$U_b,$u_min; 
#
$node_map=$project.'.segmap';
open MAP, ">$node_map" or die "Cannot open MAP file\n";
#
# Read the MAIN file
#
<MAIN>;
$total_main=0;
while ($_=<MAIN>) {
  $total_main++;
# Low_Conn 309 319 1 2545.58  
  ($Basin,$RBM_seg,$DHSVM_seg,$new_seg,$dist_meters)=split/\s+/;
  $xref_RBM[$new_seg]=$RBM_seg;
  $xref_DHSVM[$new_seg]=$DHSVM_seg;
  $dist_feet[$new_seg]=3.2808*$dist_meters;
}
#
# Read the TRIB file
#
<TRIB>;
#
while ($_=<TRIB>) {
#Upper_Conn 373 361 181 356 199 249 
  ($Basin,$RBM_seg,$DHSVM_seg,$RBM_trib,$DHSVM_trib,$new_main_seg,$new_trib_seg)=split/\s+/; 
  $Basin[$new_main_seg]=$Basin;
  $no_sources[$new_main_seg]++;
  $source[$new_main_seg][$no_sources[$new_main_seg]]=$new_trib_seg;
  $trib_seg_ref[$new_main_seg][$no_sources[$new_main_seg]]=$new_trib_seg;
}
#
# Read the ONLY Table
# 
<ONLY>;
$_=<ONLY>;
(@Only_map)=split/\s+/;
$total_tribs=$#Only_map;
print SOURCE "$total_tribs\n";
#
# Read the SEGMENTS file
#
$no_segs=0;
<SEGMENTS>;
while ($_=<SEGMENTS>) {
  (@main_seg)=split/\s+/;
  $no_in_this_row=$#main_seg;
  $no_segs=$no_segs + $no_in_this_row;
}
#
# First and last segments
#
$first_seg = $main_seg[0];
$last_seg  = $main_seg[$no_segs];
print "Number of Segments $no_segs $first_seg $last_seg\n";
#
$net_index=0;
$ntrb=0;
for $ns (0..$no_segs) {
  $main_seg=$main_seg[$ns];
  $net_index++;
  $net_seg[$main_seg]=$net_index;
  $seg_xref[$net_index]=$main_seg;
  $no_net_sources[$net_index]=$no_sources[$main_seg];
  if ($no_sources[$main_seg] > 0) {
    for $nsrc (1..$no_sources[$main_seg]) {
      $trib_seg_xref[$net_index][$nsrc]=$trib_seg_ref[$main_seg][$nsrc];
      $net_sources[$net_index][$nsrc]=$source[$main_seg][$nsrc];
      printf SOURCE "%12s %5d %5d ", $Basin[$main_seg],$net_index, $main_seg;
      printf SOURCE " $nsrc $source[$main_seg][$nsrc]";
      $nntrb=$ntrb+1;
      print SOURCE " $nntrb $Only_map[$ntrb]";
      print SOURCE "\n";
      $ntrb++;
    }
  $total_dist=$total_dist+$dist_feet[$main_seg];
  }
}
print "Trib map - $conn_only\n"; 
for $ntrb (0..$total_tribs) {

}
$miles=$total_dist/5280;
print "Total Distance $total_dist (feet) $miles\n";
print "Total number of segments - $no_segs\n";
$no_segs++;
print "Total number of segments - $no_segs\n";
#
# Print the number of headwaters to DHSVM
#
$n_head = 1;
#
print NET "$n_head Headwaters $no_segs Stream Segments in Project: $project\n";
print MAP "$n_head $no_segs\n";
#
printf NET "#_Segments %5d Headwaters %5d TribCell   %6d %6d %6d\n"
                   ,$no_segs,$n_head;
$dist=$total_dist;
$X_1=$total_dist;
for $ns (1..$no_segs) {
#
  $X_0=$X_1 - $dist_feet[$seg_xref[$ns]];
#   
  printf NET "Seq   %5d Path %5d X_0  %10.3f X_1  %10.3f Elevation %8.2f \n"
               ,$ns,$seg_xref[$ns],$X_0,$X_1,500;
  $X_1=$X_0;
}
#      printf MAP " Sequence %5d Path %5d\n",$seg_seq, $path_seg[$nh][$ns]; 
#    }
#    $nsm=$n_seg_hw[$nh];
#  }
#}
