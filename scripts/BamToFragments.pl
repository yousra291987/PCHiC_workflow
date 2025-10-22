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

#Lookup between baitID and baitname
#print "Traversing bait map...\n";
#my %baits;
#open IN,$bmap or die $!;
#while (my $line = <IN>) {
#	chomp $line;
#	my @f = split(/\t/,$line);
#	$baits{$f[3]} = $f[4];
#}
#close IN or die $!;

#Fragment info
print "Traversing fragment map...\n";
my %frags;
my %coord2index;
open IN,$rmap or die $!;
while (my $line = <IN>) {
	chomp $line;
	my @f = split(/\t/,$line);
	my $chr = $f[0];
	$chr = "chr".$f[0] unless ($f[0] =~ /chr/);
	$frags{$f[3]} = {chrom=>$chr, start=>$f[1], end=>$f[2], probe=>$f[4]};
	
	$coord2index{$chr}->{$f[3]} = {start=>$f[1], end=>$f[2], probe=>$f[4]};
}

close IN or die $!;
foreach my $chrom(keys %coord2index) {
	my @sorted = sort {$a <=> $b} keys %{$coord2index{$chrom}};
	$coord2index{$chrom}->{sorted} = \@sorted;
}

print "Traversing bedpe file...\n";
my $reads;
my $kept;
my %ints;
my %intrep;
open IN,$bedpe or die $!;
while (my $line = <IN>) {
	++$reads;
	print $reads,"...\n" if ($reads % 1000000 == 0);
	chomp $line;
	my @f = split(/\t/,$line);
	#next unless ($f[$leftchrcol] eq $f[$rightchrcol]);
	my $leftfrag = coordtofrag($f[$leftchrcol],$f[$leftcoordcol]);
	next if ($leftfrag == -1);
	my $rightfrag = coordtofrag($f[$rightchrcol],$f[$rightcoordcol]);
	next if ($rightfrag == -1);
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
		print OUT $bait,"\t",$frags{$bait}->{probe};
		print OUT "\t",$oe,"\t",$frags{$oe}->{probe};
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
	my $index_p = binary_search($coord2index{$chr}->{sorted},$coord);
	return (-1) if ($index_p == -1);
	my $frag = $coord2index{$chr}->{sorted}[$index_p];
	return($frag);
}
