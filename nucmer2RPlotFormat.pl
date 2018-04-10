#!/usr/bin/perl
#
# Author: qinmao
# Email: mqin@ymail.com
#
use 5.010; # using perl version 5.010;
use warnings; #using warnings to show message;

my $file = shift(@ARGV);
die("could not open $file!\n") unless(open(IN, $file));
chomp(my @data = <IN>);
close($file);
my $records_num = scalar(@data);
my $index = 2; # the real alignment start
my $line = "";
my $refid = "";
my $qryid = "";
my $reflen = 0;
my $qrylen = 0;
#my @reformat_dat;
my %reformat_dat = ();
my $reformat_index = 0;
while(1)
{
	last if($index == $records_num);
	$line = $data[$index];
	if($line =~ /^>/)
	{
		my @temp = split(/\s+/, $line);
		$refid = substr($temp[0],1);
		$qryid = $temp[1];
		$reflen = $temp[2];
		$qrylen = $temp[3];
	} else
	{
		my @temp = split(/\s+/, $line);
		if(scalar(@temp) != 1){
			if(exists($reformat_dat{$qryid}))
			{
				my $arrs = $reformat_dat{$qryid};
				push(@$arrs, ([$refid, $qryid, $reflen, $qrylen, $temp[0], $temp[1], 
				$temp[2], $temp[3], $temp[1] - $temp[0], 0, "NA"]));
			} else
			{
				my @arrs = ([$refid, $qryid, $reflen, $qrylen, $temp[0], $temp[1], 
					$temp[2], $temp[3], $temp[1] - $temp[0], 0, "NA"]);
				$reformat_dat{$qryid} = \@arrs;
			}
#$reformat_dat[$reformat_index] = ([$refid, $qryid, $reflen, $qrylen, $temp[0], $temp[1], 
#				$temp[2], $temp[3], $temp[1] - $temp[0], 0, "NA"]);
#		$reformat_index = $reformat_index + 1;
		}
	}
	$index = $index + 1;
}
#print(@reformat_dat);
#say($reformat_index." records");
foreach $key(keys(%reformat_dat))
{
	my $arrs = $reformat_dat{$key};
	my %ref_hash = ();
	my $max = 0;
	my $max_ref = "";
# compute belong ref
	foreach $record (@$arrs)
	{
#say(@$record);
		my @temp = @$record;
		if(exists($ref_hash{$temp[0]}))
		{
			my $value = $ref_hash{$temp[0]}[0] + $temp[8];
			$ref_hash{$temp[0]}[0] = $value;
			my $hash_arrs = $ref_hash{$temp[0]};
			push(@$hash_arrs, $record);
			if($value > $max)
			{
				$max = $value;
				$max_ref = $temp[0];
			}
		} else
		{
			my @hash_arrs = ($temp[8], $record);
			$ref_hash{$temp[0]} = \@hash_arrs;
			if($temp[8] > $max)
			{
				$max = $temp[8];
				$max_ref = $temp[0];
			}
		}
	}
# comopute ref start;
	 my $ref_start = 0;
	 my $ref_belong_address = $ref_hash{$max_ref};
	 my $ref_ollen = 0;
	 my @belong_ref = @$ref_belong_address;
	for($i = 1; $i < scalar(@belong_ref); $i = $i + 1)
	{
		my $record_address = $belong_ref[$i];
		my @temp_record = @$record_address;
#say(@temp_record);
		if($temp_record[8] > $ref_ollen)
		{
			$ref_ollen = $temp_record[8];
			$ref_start = $temp_record[4];
		}
	}
	foreach $record (@$arrs)
	{
		print(@$record[0]."\t"); # refid
		print(@$record[1]."\t"); # qryid
		print(@$record[2]."\t"); # reflen
		print(@$record[3]."\t"); # qrylen
		print(@$record[4]."\t"); # refstart
		print(@$record[5]."\t"); # refend
		print(@$record[6]."\t"); # qrystart
		print(@$record[7]."\t"); # qryend
		print(@$record[8]."\t"); # refollen
		print($ref_start."\t"); #belong ref start
		print($max_ref."\n"); # belong ref
	}
}
