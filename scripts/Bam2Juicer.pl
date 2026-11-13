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
	s/\r$//;   # tolerate Windows (CRLF) line endings so the id key has no trailing \r
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

$skipped = 0;
$written = 0;
while (<MATS>) {
	chomp;
	s/\r$//;   # tolerate Windows (CRLF) line endings
	@f=split(/\t/);
	next if (@f < 3);
	($frag1, $frag2, $count) = ($f[0], $f[1], $f[2]);
	# Only emit rows whose BOTH fragments are known in the rmap and whose count
	# is a plain integer. Each fend_table value expands to "chr\tstart" (two
	# fields), so skipping unknown fragments guarantees every printed line has
	# the full 11 columns juicer's `pre` expects (instead of silently writing a
	# short, malformed line that breaks downstream).
	unless (defined $fend_table{$frag1} && defined $fend_table{$frag2} && $count =~ /^\d+$/) {
		$skipped++;
		next;
	}
	for (my $i = 1; $i <= $count; $i++) {
		print OUT int(rand(1000)),"\t",0,"\t",$fend_table{$frag1},"\t",$frag1,"\t",1,"\t",$fend_table{$frag2},"\t",$frag2,"\t","40","\t","40","\n";
		$written++;
	}
}
if ($skipped > 0) {
	print STDERR "[Bam2Juicer] skipped $skipped MAT record(s) with an unknown fragment id or non-integer count; wrote $written row(s).\n";
}

