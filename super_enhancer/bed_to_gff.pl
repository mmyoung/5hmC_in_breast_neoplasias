#!/bin/perl
use strict;
use warnings;


my($filebed,$fileout)=@ARGV;
open BED, "<$filebed" or die $!;
open OUT, ">$fileout" or die $!;
our %bed;our %input;

while (my $peak=<BED>)
{
	chomp $peak;
	my @arr=split /\t/,$peak;
	my $start=int($arr[1])+1;
	my $end=int($arr[2])+1;
	print OUT "$arr[0]\t$arr[3]\t\t$start\t$end\t\t$arr[5]\t\t$arr[3]\n";
}
