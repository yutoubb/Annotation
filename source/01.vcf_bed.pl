#!perl -w
use strict;
use FindBin qw($Bin);
use Getopt::Long;
my($in,$out,$bed,$help);

=head1 #Parameter
perl vcf_bed.pl
	-i 		the input file. The format of input file should be .vcf or .vcf.gz/zip
	-b 		bed file.default for BRCA1/2 genes' bed region 
	-h|help 	die usage
=cut
=head1 #Example
perl vcf_bed.pl -i input.vcf

=cut

GetOptions(
	'i=s'=>\$in,
	'b:s'=>\$bed,
	'h|help'=>\$help,
);



die `pod2text $0` unless($in && !defined $help);
#==================open files============================#
$bed = join "/",$Bin,'BRCA.bed' unless($bed);
open IN1, "$bed" || die "#$bed $!\n";
if($in =~ /vcf.[gzip]+/){
	open IN2, "gunzip -dc $in|" || die "#IN2$!\n";;
}
elsif($in =~ /vcf$/){
	open IN2, "$in" || die "#IN2$!\n";
}
else{
	die "The format of input file should be .vcf or .vcf.gz/zip";
}
($out) = $in =~ /^(.*).vcf/;
open OUT, ">$out.BRCA.vcf" || die "#OUT$!\n";
#==================files open end========================#

#=============hash for BED region filter=================#
my(@a,$pos,%hash);
while(<IN1>){
	chomp;
	@a = split /\t/;
	foreach my $i($a[1]..$a[2]){
		$pos = join ":",$a[0],$i;
		$hash{$pos} = 1;
	}
}
close IN1;
#=====================hash end==========================#
while(<IN2>){
	chomp;
	if(/^#/){
		print OUT "$_\n";
	}
	else{
		@a = split /\t/;
		$pos = join ":",$a[0],$a[1];
		if(exists $hash{$pos}){
			print OUT "$_\n";
		}
	}
}
close IN2;
close OUT;
