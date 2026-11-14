#!/usr/bin/env perl.exe
use strict;
use diagnostics;

if ($#ARGV != 7) {
	print STDERR "usage: $0 <bedpe> <baitmap> <rmap file> <output> <read1 chr col> <read1 coord col> <read 2 chr col> <read 2 coord col>\n";
	exit 1;
}

my $bedpe = $ARGV[0];
my $bmap = $ARGV[1];
my $rmap = $ARGV[2];
my $ofn = $ARGV[3];
my $leftchrcol = $ARGV[4]-1;
my $leftcoordcol = $ARGV[5]-1;
my $rightchrcol = $ARGV[6]-1;
my $rightcoordcol = $ARGV[7]-1;

# Captured (bait) fragments, keyed by fragment id. Used to flag each end of an
# interaction as a probe (1) or not (0) in the .mat's ProbeF1/ProbeF2 columns.
# The capture-efficiency step (00_CaptureHiC_CaptureEfficiency_Control.r) reads
# those flags from columns 2 and 4; a standard CHiCAGO .rmap has only 4 columns
# (chr start end fragID) and no probe flag, so taking it from the rmap left those
# columns empty -> P0/P1/P2 all 0 and an empty CaptureEfficiency plot.
print "Traversing bait map...\n";
my %baits;
open IN,$bmap or die $!;
while (my $line = <IN>) {
	chomp $line;
	$line =~ s/\r$//;   # tolerate Windows (CRLF) line endings
	my @f = split(/\t/,$line);
	next unless (defined $f[3] && $f[3] ne "");
	$baits{$f[3]} = 1;
}
close IN or die $!;

#Fragment info
print "Traversing fragment map...\n";
my %frags;
my %coord2index;
open IN,$rmap or die $!;
while (my $line = <IN>) {
	chomp $line;
	$line =~ s/\r$//;   # tolerate Windows (CRLF) line endings
	my @f = split(/\t/,$line);
	my $chr = $f[0];                 # keep the rmap chromosome name exactly as given
	(my $key = $chr) =~ s/^chr//i;   # normalized key so "chr1" (bed) and "1" (rmap) match
	$frags{$f[3]} = {chrom=>$chr, start=>$f[1], end=>$f[2], probe=>$f[4]};

	$coord2index{$key}->{$f[3]} = {start=>$f[1], end=>$f[2], probe=>$f[4]};
}

close IN or die $!;
foreach my $chrom(keys %coord2index) {
	my @sorted = sort {$a <=> $b} keys %{$coord2index{$chrom}};
	$coord2index{$chrom}->{sorted} = \@sorted;
}

print "Traversing bedpe file...\n";
my $reads = 0;
my $kept = 0;
my %ints;
my %intrep;
open IN,$bedpe or die $!;
while (my $line = <IN>) {
	++$reads;
	print $reads,"...\n" if ($reads % 1000000 == 0);
	chomp $line;
	$line =~ s/\r$//;   # tolerate Windows (CRLF) line endings
	my @f = split(/\t/,$line);
	#next unless ($f[$leftchrcol] eq $f[$rightchrcol]);
	my $leftfrag = coordtofrag($f[$leftchrcol],$f[$leftcoordcol]);
	next if ($leftfrag == -1);
	my $rightfrag = coordtofrag($f[$rightchrcol],$f[$rightcoordcol]);
	next if ($rightfrag == -1);
	++$kept;
	my $bcount = 0;
	if ((exists $ints{$leftfrag}->{$rightfrag}) ||(exists $ints{$rightfrag}->{$leftfrag}) ) {
		++$ints{$leftfrag}->{$rightfrag};
		
	}
	else {
		$ints{$leftfrag}->{$rightfrag}=1;
		#push @{$intrep{$leftfrag}},$rightfrag;
	}
}

close IN or die $!;
print "$kept reads out of $reads processed\n";

my $cis;
print "Writing output file...\n";
open OUT,">",$ofn or die $!;
print OUT "Fragment1\tProbeF1\tFragment2\tProbeF2\tN\n";

my @baitsort = sort {$a <=> $b} keys %ints;
foreach my $bait(@baitsort) {
	my @oesort = sort {$a <=> $b} keys %{$ints{$bait}};
	$cis += scalar @oesort;
	foreach my $oe(@oesort) {
		print OUT $bait,"\t",(exists $baits{$bait} ? 1 : 0);
		print OUT "\t",$oe,"\t",(exists $baits{$oe} ? 1 : 0);
		print OUT "\t",$ints{$bait}->{$oe},"\n";
	}
}
close OUT or die $!;
print "$cis different interactions output\n";



sub binary_search {
	my $arr = shift;
	my $value = shift;
	my $left = 0;
	my $right = $#$arr;
	while ($left <= $right) {
		my $mid = ($right+$left) >> 1;
		my $ind = $arr->[$mid];
		my $mstart = $frags{$ind}->{start};
		my $mend = $frags{$ind}->{end};
		if ($value >= $mstart and $value <= $mend) {
			return $mid;
		}
		elsif ($value < $mstart) {
			$right = $mid-1;
		}
		else {
			$left = $mid+1;
		}
	}
	$left = -1 if ($left>$#$arr);
	$right = -1 if ($right<0);
	return $left-1;
}

sub coordtofrag {
	my ($chr,$coord) = (@_);
	(my $key = $chr) =~ s/^chr//i;   # match the normalized rmap keys
	my $sorted = $coord2index{$key}->{sorted};
	# Chromosome not present in the rmap (e.g. a scaffold or a naming mismatch):
	# skip the read quietly instead of flooding the log with "uninitialized value"
	# warnings for every read, which previously produced multi-GB log files.
	return (-1) unless (defined $sorted && @$sorted);
	my $index_p = binary_search($sorted,$coord);
	return (-1) unless (defined $index_p && $index_p != -1);
	my $frag = $sorted->[$index_p];
	return (defined $frag ? $frag : -1);
}
