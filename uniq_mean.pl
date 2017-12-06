
#!perl -w
use strict;

open IN, shift or die $!;
my (%sum, %count, %hash);
while(<IN>){
	chomp;
	my @F = split / /;
	if (exists $sum{$F[5]}){
		$sum{$F[5]} += $F[4];
		$count{$F[5]} ++;
	}else{
		$hash{$F[0]}{$F[1]} = $F[5];
		$sum{$F[5]} = $F[4];
		$count{$F[5]} = 1;
	}
}

for my $k1(1..22,"X","Y"){
	if (exists $hash{"chr$k1"}){
		for my $k2(sort {$a<=>$b} keys %{$hash{"chr$k1"}}){
			my $value = $sum{$hash{"chr$k1"}{$k2}} / $count{$hash{"chr$k1"}{$k2}};
			print "chr$k1\t$k2\t$hash{\"chr$k1\"}{$k2}\t$value\n";
		}
	}
}

