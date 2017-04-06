#!perl -w
use strict;
use FindBin qw($Bin);
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;

=head1 #Parameter

=head1 perl variant_effect_predictor.BRCA_onlineV1.pl
	-i		*.vcf file
	-ens		[options]ensembl anno file,[default means the ensembl anno file shares the same path with vcf file]
	-o		[options]output file,[default means the output file shares the same path with vcf file]
	-h|help         die usage

=cut
=head1 #Example

perl variant_effect_predictor.BRCA_onlineV1.pl -i *.vcf -o output.file

=cut

my($vcf,$ensembl,$out,$help);
GetOptions(
	'i=s'=>\$vcf,
	'ens:s'=>\$ensembl,
	'o:s'=>\$out,
	'h|help'=>\$help,
);
die `pod2text $0` unless($vcf && !defined $help);
#=============================================open files=========================================================#
my($in_name,$in_path,$in_suffix) = fileparse($vcf);

my ($prefix) = $in_name =~ m/(.*).vcf/;
$ensembl = join "",$in_path,$prefix unless($ensembl);
my $out_1 = join ".",$prefix,'T.anno';
$out = join "",$in_path,$out_1 unless($out);
open IN1, "$vcf" or die "#$vcf $!\n";
open IN2, "$ensembl" or die "#$ensembl $!\n#The author warnning you that you need to use -ens and -o parameter to defined the path of ensembl annotation file and output file.";
open OUT, ">$out" or die "#$out $!\n#The author warnning you that you need to use -ens and -o parameter to defined the path of ensembl
 annotation file and output file.";
#============================================files open end======================================================#

#================================================================================================================#
#=================================read vcf file for a.Gtype b.AD (c.ratio,based on b.AD)=========================#

my (@a,$variation,$Gt,%hash_Gt,%hash_Refalt,%hash_Refalt1,%hash_AD,$start,@ref,@alt,@alt1,@alt2,$ref,$alt,$alt1,$alt2,%hash_start,$len,%hash_forad3,$ref_a,$i,$j,$alt_a,$intron,@exon,$AD_sum);
while(<IN1>){
	chomp;
	unless(m/^#/){
		@a = split /\t+/;
		if ($a[2] ne '.'){#situations that mutations have rsID in vcf file
			$variation = $a[2];
			($hash_AD{$variation}) = $a[-1] =~ m/\d\/\d:([0-9,]+):/;
			($Gt) = $a[7] =~ m/AC=(\d)/;
			if($Gt){
			$hash_Gt{$variation} = 'Het' if $Gt eq 1;
			$hash_Gt{$variation} = 'Hom' if $Gt eq 2;
			}
			else{
			$hash_Gt{$variation} = '-';
			}
			$hash_Refalt{$variation} = join "\t",$a[3],$a[4];
			$hash_Refalt1{$variation} = join ">",$a[3],$a[4];
			$hash_Refalt1{$variation} =~ s/,/\//g;
			$hash_start{$variation} = $a[1];
			$hash_forad3{$variation} = join "/",$a[3],$a[4];#build hash to store AD value
			$hash_forad3{$variation} =~ s/,/\//g;
		}
		else{#situations that mutations have no rsID in vcf file
			$len = length($a[3]) + length ($a[4]);
			if($len > 2){#not point to point mutations
				$start = $a[1] + 1;
				@ref = split //,$a[3];
				$ref_a = shift @ref;#This is why $start need to plus 1
				if (@ref){
					$ref = join "",@ref;
				}
				else{
					$ref = '-';
				}
				if($a[4] =~ m/,/){#multiple mutations at the same position
					@alt = split /,/,$a[4];
					$i = @alt;
					for($j=0,$j<$i,$j++){
						@alt1 = split //,$alt[$j];
						shift @alt1;
						if(@alt1){
							$alt1 = join "",@alt1;
							push(@alt2,$alt1);
						}
						else{
							$alt1 = '-';
							push(@alt2,$alt1);
						}
					}
					$alt = join "/",@alt2;
					@alt2=();
				}
				else{#normal mutations, not multiple mutations at the same position
					@alt1 = split //,$a[4];
					$alt_a = shift @alt1;
					if(@alt1){
						$alt = join "",@alt1;
					}
					else{
						$alt = '-';
					}
				}
				if($ref_a eq $alt_a){
				$variation = join "_",$a[0],$start,$ref;
				$variation = join "/",$variation,$alt;
				}
				else{
				$variation = join "_",$a[0],$a[1],$a[3];
				$variation = join "/",$variation,$a[4];	
				}
				$variation =~ s/,/\//g;
				$alt = undef;
			}	
			else{#point to point mutations
			$variation = join "_",$a[0],$a[1],$a[3];
			$variation = join "/",$variation,$a[4];	
			}
			$variation =~ s/chr//g;
			($hash_AD{$variation}) = $a[-1] =~ m/\d\/\d:([0-9,]+):/;
			($Gt) = $a[7] =~ m/AC=(\d)/;
			if($Gt){
			$hash_Gt{$variation} = 'Het' if $Gt eq 1;
			$hash_Gt{$variation} = 'Hom' if $Gt eq 2;
			}
			else{
			$hash_Gt{$variation} = '-';
			}
			$hash_Refalt{$variation} = join "\t",$a[3],$a[4];
			$hash_Refalt1{$variation} = join ">",$a[3],$a[4];
			$hash_start{$variation} = $a[1];
			$hash_forad3{$variation} = join "/",$a[3],$a[4];
			$hash_forad3{$variation} =~ s/,/\//g;
		}
#		print "*vcf*$variation\t*$hash_AD{$variation}\n";
	}
}	
close IN1;
#==================================================read vcf file end=============================================#
#================================================================================================================#

#================================================================================================================#
#========================read ensembl anno file to export basic annotation information===========================#

my (%hash_line,$chr,@AD,$ratio,$ratio1,$ratio2,$Mtype,$gene,$fuction,$exon,$exon_all,$trans,$cHGVS,$pHGVS,@b,$rs_id,@c,$end,$c_len,@variation,$line);
print OUT "#chr\tstart\tend\tGene\tTrans\tcHGVS\tpHGVS\tExon\tMtype\tGtype\tFuction\tg.HGVS\tRef\tAlt\tEnsambl_ID\trs_ID\tAD\tratio\n";
while(<IN2>){
	chomp;
	if(m/^#/){#filter the description of ensembl annotation
	}
	else{
		@a = split /\t/;
		$line = join "-",$a[0],$a[4];#remove duplication of same annotation record
		if (exists $hash_line{$line}){
		}
		else{
			$hash_line{$line} = 1;
			($chr) = $a[1] =~ m/([0-9XY]+):\d+/;
			if(exists $hash_start{$a[0]}){
				$start = join "",$hash_start{$a[0]},$hash_Refalt1{$a[0]};#for example:31000243T>C,vcf position & ref>alt
#print "**1$a[0]*$start**\n";
			}
			else{
				$start = '-';
			}
			$start = '-' unless($start);
			$trans = $a[4];
			$Mtype = $a[17];
			$gene = $a[18];
			$fuction = $a[6];
			$exon = $a[34];
			$intron = $a[35];
			if($a[34] ne '-' || $a[35] ne '-'){
#print "**2$a[0]*$start**\n";
				if($a[34] ne '-'){
					@exon = split /\//,$a[34];
					if($exon[0] eq $exon[1]){
						$exon = "EX$exon[0]\E";
					}
					else{
						$exon = "EX$exon[0]";
					}
				}
				else{
					@exon = split /\//,$a[35];
					$exon = "IVS$exon[0]";
				}
			}
			else{
				$exon = '-';
			}
			$rs_id = $a[0] if $a[0] =~ m/rs/;
			unless ($rs_id) {
				$rs_id = '-';
			}
#			($cHGVS) = $a[37] =~ m/^[A-Z0-9\._]+:(c\.[0-9a-zATCG\_\>\+\-\*]+)/;
			($cHGVS) = $a[37] =~ m/^[A-Z0-9\._]+:(c\.\w+)/;
			($pHGVS) = $a[38] =~ m/^[A-Z0-9\._]+:(p\.\w+)/;
			unless($cHGVS){
				$cHGVS = '-';
			}
			unless($pHGVS){
				if($a[6] =~ m/^synonymous/){
					$pHGVS = 'p.[=]';
				}
				else{
					$pHGVS = '-';
				}
			}
			$pHGVS =~ s/Ter/*/g;
			if(exists $hash_AD{$a[0]}){
#print "**3$a[0]*$hash_AD{$a[0]}**\n";
				@AD = split /,/,$hash_AD{$a[0]};
				my $AD = @AD;
				if($AD == 2){
#print "**4$a[0]*$hash_AD{$a[0]}**\n";
					$ratio = sprintf ("%.2f",$AD[1] / ($AD[0] + $AD[1]) * 100);
					if(($gene eq 'BRCA1' && $trans eq 'NM_007294.3') || ($gene eq 'BRCA2' && $trans eq 'NM_000059.3')){
#print "**5$a[0]*$hash_AD{$a[0]}**\n";
					@b = split /_/,$a[0];
					@c = split /\//,$b[-1];
					$c_len = length($c[0]);
					$end = $b[1] + $c_len - 1;
					print OUT "chr$b[0]\t$b[1]\t$end\t$gene\t$trans\t$cHGVS\t$pHGVS\t$exon\t$Mtype\t$hash_Gt{$a[0]}\t$fuction\tg.$start\t$hash_Refalt{$a[0]}\t$a[0]\t$rs_id\t$hash_AD{$a[0]}\t$ratio\%\n";
					}
				}
				else{
#print "**4-1$a[0]*$hash_AD{$a[0]}**\n";	
					if(($gene eq 'BRCA1' && $trans eq 'NM_007294.3') || ($gene eq 'BRCA2' && $trans eq 'NM_000059.3')){
#print "**5-1$a[0]*$hash_AD{$a[0]}**\n";
						@b = split /_/,$a[0];
						@c = split /\//,$b[-1];
						$c_len = length($c[0]);
						$end = $b[1] + $c_len - 1;
						@variation =  split "/",$hash_forad3{$variation};#this hash stores reference/alt
						foreach $i (0..$AD){#denominator of ratio
						$AD_sum+=$AD[$i];
						}
						foreach $i (1..$AD){
							$ratio1 = sprintf ("%.2f",$AD[$i] / $AD_sum * 100);
							if($a[2] == $variation[$i]){
								print OUT "chr$b[0]\t$b[1]\t$end\t$gene\t$trans\t$cHGVS\t$pHGVS\t$exon\t$Mtype\t$hash_Gt{$a[0]}\t$fuction\tg.$start\t$hash_Refalt{$a[0]}\t$a[0]\t$rs_id\t$hash_AD{$a[0]}\t$ratio1\%\n";
							}
						}
					}
				}
				$cHGVS = undef;
				$pHGVS = undef;
			}
			else{
#print "**6$a[0]*$gene*$trans**\n";
				if(($gene eq 'BRCA1' && $trans eq 'NM_007294.3') || ($gene eq 'BRCA2' && $trans eq 'NM_000059.3')){
#print "**6-1$a[0]**\n";
					@b = split /:/,$a[1];
					@c = split /-/,$b[-1];
					$start = $c[0] - 1;
					$end = $c[1];
					if($a[0] =~ /deletion/){
#print "**7$a[0]*$start*$exon**\n";
						print OUT "chr$b[0]\t$start\t$end\t$gene\t$trans\t$exon deletion\t$exon deletion\t$exon\tCNV\t-\t$fuction\t$exon deletion\t-\t-\t$a[0]\t$rs_id\t-\t-\n";
					}
					elsif($a[0] =~ /duplication/){
#print "**8$a[0]*$start*$exon**\n";						
						print OUT "chr$b[0]\t$start\t$end\t$gene\t$trans\t$exon duplication\t$exon duplication\t$exon\tCNV\t-\t$fuction\t$exon duplication\t-\t-\t$a[0]\t$rs_id\t-\t-\n";
					}
				}
			}
		}
	}
}			
close IN2;
close OUT;

