#!/usr/bin/perl
use 5.8.1; use strict; use warnings;

my ($cmd, $path_exe, $path_in, $path_out, $genome, $fasta);
my ($motif_base, $motif_ext, $n_motifs, $pwm_ext, $dist_ext, $bin_ext, $pwm_log);
my ($i, $spacer_min, $spacer_max, $genome_prom, $pvalue_ce, $log_pvalue_bonf, $asym_rat);

if(scalar(@ARGV)==0){ die "Wrong arguments!";}

$path_exe=           $ARGV[0]; # path to executable
$path_in=            $ARGV[1]; # input path, fasta, motifs
$path_out=           $ARGV[2]; # output path
$motif_base=         $ARGV[3]; # motif base name
$motif_ext=          $ARGV[4]; # motif extention name
$n_motifs=           $ARGV[5]; # number of motifs
$genome_prom=        $ARGV[6]; # genome promoters fasta
$fasta      =        $ARGV[7]; # forground fasta 
$spacer_min    =     $ARGV[8]; # spacer min
$spacer_max    =     $ARGV[9]; # spacer max
$pvalue_ce     =     $ARGV[10]; # pvalue CE
$log_pvalue_bonf =   $ARGV[11]; # pvalue bomferroni
$asym_rat      =     $ARGV[12]; # threshold asymmetry ratio


$pwm_ext = ".pwm";
$dist_ext = ".dist";
$bin_ext = ".binary";
$pwm_log ="pwm.log";

if ( -d "$path_out"){
    print "Directory already exist.\n";
}
else{
   mkdir($path_out)
   or die("Can't create directory \"$path_out\": $!\n");
}

for($i=1;$i<=$n_motifs;$i++)
{
$cmd= "$path_exe/pfm_to_pwm/pfm_to_pwm ${path_in}${motif_base}${i}${motif_ext} $path_out/${motif_base}${i}${pwm_ext}";
print "$cmd\n";
system $cmd;

$cmd= "$path_exe/pwm_thr_err/pwm_thr_err ${path_in}${motif_base}${i}${motif_ext} ${path_out}${motif_base}${i}${pwm_ext} ${genome_prom} ${path_out}${motif_base}${i}${dist_ext} ${path_out}${motif_base}${bin_ext} 0.002 0.0000005 ${pwm_log} 0.00002 0";
print "$cmd\n";
system $cmd;
}

$cmd= "$path_exe/denovo/denovo ${path_in}${fasta} ${path_out}${motif_base}${bin_ext} ${n_motifs} ${spacer_min} ${spacer_max} ${pvalue_ce} ${log_pvalue_bonf} ${asym_rat}";
print "$cmd\n";
system $cmd;

