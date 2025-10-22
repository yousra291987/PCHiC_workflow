#!/usr/bin/env perl
#use List::Util qw(sum);
#use List::MoreUtils qw(uniq);
#files generated: readnames(random)\tstrand1(0:1)\tchr1\tpos1(coord)\tfrag1(fend)\tstrand2\tchr2\tpos2\tfrag2

$mats=$ARGV[0];
$fends=$ARGV[1];
$outfile=$ARGV[2];

open(MATS, $mats) || die $mats;

open(FENDS,$fends)|| die $fends;

open(OUT,'>', $outfile);
%fend_table;
%param_table;
while(<FENDS>){
	chomp;
	@tab=split(/\t/);
	#if( $tab[4] eq "+") {
	#	$tab[4]=0;
	#}else { $tab[4]=1;
	#}
	$fend_table{$tab[3]}=join("\t",$tab[0],$tab[1]);
	
}

#while (<PARAM>) {
#	chomp;
#	@t=split(/\t/);
#	$param_table{$t[6]}=$t[0];
#}

while (<MATS>) {
	chomp;
	@f=split(/\t/);
	for (my $i =1; $i <= $f[2]; $i++) {
		print OUT int(rand(1000)),"\t",0,"\t",$fend_table{$f[0]},"\t",$f[0],"\t",1,"\t",$fend_table{$f[1]},"\t",$f[1],"\t","40","\t","40","\n";
	}
}

