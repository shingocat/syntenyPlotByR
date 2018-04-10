#!/usr/bin/perl
use strict;
use warnings;

open (IN1,"$ARGV[0]") || die "can't open file $ARGV[0]\n";

my %hash1;my %hash2;

while(<IN1>){
chomp;
my @a=split;
if($hash1{$a[1]}){
	$hash1{$a[1]}+=$a[8] if $a[9] eq "+";
	$hash1{$a[1]}-=$a[8] if $a[9] eq "-";
	}else{
	$hash1{$a[1]}=0;
	$hash1{$a[1]}+=$a[8] if $a[9] eq "+";
	$hash1{$a[1]}-=$a[8] if $a[9] eq "-";
	}
}

close IN1;

foreach (keys %hash1){
if($hash1{$_}>=0){
	$hash1{$_}="+";
	}else{
	$hash1{$_}="-";
	}
}

open (IN1,"$ARGV[0]") || die "can't open file $ARGV[0]\n";
my %hash3;
while(<IN1>){
chomp;
my @a=split;
if($a[9] eq $hash1{$a[1]}){
	$hash3{$a[1]}=0 if !$hash2{$a[1]};
	if($hash3{$a[1]}<$a[8]){$hash2{$a[1]}=$a[4];$hash3{$a[1]}=$a[8];}
	}
}

close IN1;
open (IN1,"$ARGV[0]") || die "can't open file $ARGV[0]\n";

while(<IN1>){
chomp;
my $all=$_;
my @a=split;
print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$hash2{$a[1]}\t$a[0]\n";
}
