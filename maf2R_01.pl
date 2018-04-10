#!/usr/bin/perl
use strict;
use warnings;

open (IN1,"$ARGV[0]") || die "can't open file $ARGV[0]\n";

my $on;my ($ref,$ref_len,$ref_start,$ref_end,$ref_between,$query,$query_len,$query_start,$query_end,$query_order);

while(<IN1>){
chomp;
my @a=split;
next if @a==0;
if($a[0] eq "a"){
	$on=1;
	}elsif($a[0]=~/^\s*$/){
	$on=0;
	}elsif($a[0] eq "s" && $on==1){
	$ref=$a[1];
	$ref_len=$a[5];
	$ref_start=$a[2];
	$ref_end=$a[2]+$a[3];
	$ref_between=$a[3];
	$on++;
	}elsif($a[0] eq "s" && $on==2){
	$query=$a[1];
	$query_len=$a[5];
	$query_start=$a[2];
	$query_end=$a[2]+$a[3];
	$query_order=$a[4];
	if($query_order eq "+"){
		print "$ref\t$query\t$ref_len\t$query_len\t$ref_start\t$ref_end\t$query_start\t$query_end\t$ref_between\t$query_order\n";
		}elsif($query_order eq "-"){
		print "$ref\t$query\t$ref_len\t$query_len\t$ref_start\t$ref_end\t$query_end\t$query_start\t$ref_between\t$query_order\n";
		}
	}
}
