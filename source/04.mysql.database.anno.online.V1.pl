#!/use/bin/perl -w
use strict;
use DBI;
my $anno_file=shift;#注释程序得到的文件
my $db=&connect_mysql;
my $utf8=$db->prepare("set names utf8");
$utf8->execute();
###表头###
print "gene\ttrans\tchr\tstart\tend\tExIn_ID\tref\tmut\tAD\tratio\tMtype\tGtype\tg.HGVS\tc.HGVS\tp.HGVS\tfunction\trs_id\tClinical_Classification_BGI\tfunction_predict\thapmap_AF\tdbsnp_maf\tnifty_maf\t1000_AF\t1000_EAS_AF\t1000_AMR_AF\t1000_AFR_AF\t1000_EUR_AF\t1000_SAS_AF\tESP6500_AA_AF\tESP6500_EA_AF\tExAC_AF\tExAC_Adj_AF\tExAC_AFR_AF\tExAC_AMR_AF\tExAC_EAS_AF\tExAC_FIN_AF\tExAC_NFE_AF\tExAC_SAS_AF\tExAC_nonTCGA_Adj_AF\tExAC_nonTCGA_AFR_AF\tExAC_nonTCGA_AMR_AF\tExAC_nonTCGA_EAS_AF\tExAC_nonTCGA_FIN_AF\tExAC_nonTCGA_NFE_AF\tExAC_nonTCGA_SAS_AF\tExAC_nonpsych_AF\tExAC_nonpsych_Adj_AF\tExAC_nonpsych_AFR_AF\tExAC_nonpsych_AMR_AF\tExAC_nonpsych_EAS_AF\tExAC_nonpsych_FIN_AF\tExAC_nonpsych_NFE_AF\tExAC_nonpsych_SAS_AF\tSIFT_pred\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_pred\tLRT_pred\tMutationAssessor_pred\tMetaSVM_pred\tMetaLR_pred\tBRCA_inhouse_clinical_category\tbic_Clinical_Classification_Category\tumd_biological_significance\tClinvar_Clinical_Significance\tlovd_FunctionalAnalysis_or_Result\tdbscsnv_ada_score\tdbscsnv_rf_score\n";
###
open IN,"<$anno_file" or die $!;
while(<IN>)
{
	chomp;
	next if(/^#/);
#	print "#not next\n";
	my($chr,$start,$end,$ref_hg19,$Mtype,$AD,$GType,$Genotype,$trans,$gene,$exon,$function,$cHGVS,$pHGVS,$Ratio)=(split "\t",$_)[0,1,2,12,8,16,9,14,4,3,7,10,5,6,17];
#	print "chr=$chr,start=$start,end=$end,ref=$ref_hg19,Mtype=$Mtype,AD=$AD,GType=$GType,Genotype=$Genotype,trans=$trans,gene=$gene,exon=$exon,function=$function,cHGVS=$cHGVS,pHGVS=$pHGVS,Ratio=$Ratio\n";
	next unless $gene=~/BRCA1|BRCA2/;
#	print "BRCAnot next\n";
	my $pos=$start; #注释程序结果中发现的规律
#	my ($HGVS_c1,$HGVS_c2)=(split ";",$cHGVS)[0,1];
	my @Genotype = split /_/,$Genotype;
	my @Genotype_split=(split "/",$Genotype[-1]);
#	my $HGVS_c=(($Genotype_split[0] ne $ref_hg19)and($Genotype_split[1] ne $ref_hg19)and($Genotype_split[0] ne $Genotype_split[1]))?$HGVS_c1.";c.".$HGVS_c2:$HGVS_c1;
#	$HGVS_c1=$cHGVS;
#	my $HGVS_c=$cHGVS;
#	my ($HGVS_p1,$HGVS_p2,$HGVS_p)=(".") x 3;
#	if($pHGVS ne ".")
#	{
#		($HGVS_p1,$HGVS_p2)=(split ";",$pHGVS)[0,1];
#		$HGVS_p=$HGVS_c=~/\;/?$HGVS_p1.";p.".$HGVS_p2:$HGVS_p1;
#	}
#	$HGVS_c1=~s/\[//g; $HGVS_c1=~s/\]//g;
#	$HGVS_c2=~s/\[//g; $HGVS_c2=~s/\]//g;
#	$HGVS_c=~s/\[//g; $HGVS_c=~s/\]//g;
#	my $HGVS_p1=$pHGVS;
#	 my $HGVS_p=$pHGVS;
#	$HGVS_p1=~s/\[//g; $HGVS_p1=~s/\]//g;
#	$HGVS_p2=~s/\[//g; $HGVS_p2=~s/\]//g;
#	$HGVS_p=~s/\[//g; $HGVS_p=~s/\]//g;
#	my $Ratio=$HGVS_c=~/\;/?(split ",",$AD)[0]/((split ",",$AD)[0]+(split ",",$AD)[1]).";".(split ",",$AD)[1]/((split ",",$AD)[0]+(split ",",$AD)[1]):(split ",",$AD)[1]/((split ",",$AD)[0]+(split ",",$AD)[1]);
	my ($mut1,$mut2)=($Genotype_split[0],$Genotype_split[1]);
	unless($mut2){
		$mut1 = $Mtype;
		$mut2 = $Mtype;
	}
	my $mut=$cHGVS=~/\;/?$mut1.",".$mut2:$mut2;
	$mut = $Mtype unless($mut);
	my $HGVS_g=$cHGVS=~/\;/?"g.".$pos.$ref_hg19.">".$mut1.";"."g.".$pos.$ref_hg19.">".$mut2:"g.".$pos.$ref_hg19.">".$mut2;

#print "\n*$gene\t$trans\t$chr\t$pos\t$end\t$exon\t$ref_hg19\t$mut\t$AD\t$Ratio\t$Mtype\t$GType\t$HGVS_g\t$cHGVS\t$pHGVS\t$function\n";
	#####调整HGVS_g的表示方式####
	my $ref_hg19_reverse=$ref_hg19;
	$ref_hg19_reverse=~tr/ATCG/TAGC/;
	my @mut_all=($mut1,$mut2);
	my $mut_reverse1=$mut1;
	my $mut_reverse2=$mut2;
	$mut_reverse1=~tr/ATCG/TAGC/;
	$mut_reverse2=~tr/ATCG/TAGC/;
	my @mut_reverse=($mut_reverse1,$mut_reverse2);
	
	my @HGVS_g_bic_lovd;
	my @HGVS_g_bic_lovd_reverse;
	foreach my $i(0..1)
	{
		if(length($ref_hg19)==length($mut_all[$i])) #snv
		{
			push(@HGVS_g_bic_lovd,"g.".$pos.$ref_hg19.">".$mut_all[$i]);
			push(@HGVS_g_bic_lovd_reverse,"g.".$pos.$ref_hg19_reverse.">".$mut_reverse[$i]) if($gene eq "BRCA1");
		}
		else
		{
			if(length($ref_hg19)>length($mut_all[$i]))
			{
				if(length($mut_all[$i])>1) #delins
				{ 
					push(@HGVS_g_bic_lovd,"g.".($pos+1)."_".$end."delins".substr($mut_all[$i],1));
					push(@HGVS_g_bic_lovd_reverse,"g.".($pos+1)."_".$end."delins".substr($mut_reverse[$i],1)) if($gene eq "BRCA1");
				}
				else	#del
				{
					if(length($ref_hg19)==2){push(@HGVS_g_bic_lovd,"g.".$end."del");} #不写出具体的缺失的碱基，因为位点一样的话缺失的碱基是一样的
					else{push(@HGVS_g_bic_lovd,"g.".($pos+1)."_".$end."del");}
				}
			}
			elsif(length($ref_hg19)<length($mut_all[$i])){ #ins
				push(@HGVS_g_bic_lovd,"g.".$pos."_".($pos+1)."ins".substr($mut_all[$i],1));
				push(@HGVS_g_bic_lovd_reverse,"g.".$pos."_".($pos+1)."ins".substr($mut_reverse[$i],1)) if($gene eq "BRCA1");
			}
		}
	}
	my $HGVS_g_bic_lovd1=$HGVS_g_bic_lovd[0]; #bic,lovd数据库HGVS_g的表示方法,snv,ins,del,delins,这里只需要考虑ref到mut2的突变形式，因为在这两个数据库中每一行只有一个点突变的记录
	my $HGVS_g_bic_lovd2=$HGVS_g_bic_lovd[1];
	my $HGVS_g_bic_lovd_reverse1="NA";
	my $HGVS_g_bic_lovd_reverse2="NA";
	if(scalar(@HGVS_g_bic_lovd_reverse)==2)
	{
		$HGVS_g_bic_lovd_reverse1=$HGVS_g_bic_lovd_reverse[0];
		$HGVS_g_bic_lovd_reverse2=$HGVS_g_bic_lovd_reverse[1];
	}


	#umd数据库HGVS_g的表示方法,snv,ins,del,delins
	#umd genomic可能的表示方法，一个点的缺失：g.41276136del,g.41276136_41276136del;多个点的缺失:g.41276136del,g.41276136_41276138del,不论多少个位点的缺失都会有两种表示方法。有时以hg19的为参考，有时又取其互补碱基,因为BRCA1是负链基因，BRCA2是正链基因
	#ins更是乱七八糟，很任性的胡乱大小写，有时插入的是具体的碱基，有时插入的又是数字，甚至有时插入的是氨基酸.但是好在表示方式一致，g.***ins	
	#delins还算正常
	
	my @HGVS_g_umd1;
	my @HGVS_g_umd2;
	my @HGVS_g_umd_reverse;
	foreach my $i(0..1)
	{
		if(length($ref_hg19)==length($mut_all[$i])){
			push(@HGVS_g_umd1,"g.".$pos.$ref_hg19.">".$mut_all[$i]);
			push(@HGVS_g_umd_reverse,"g.".$pos.$ref_hg19_reverse.">".$mut_reverse[$i]) if($gene eq "BRCA1");
		} #考虑正负链表示方法的不同,因为BRCA1是负链基因，有些数据库中的表示方法不太一样
		else
		{
			if(length($ref_hg19)>length($mut_all[$i]))
			{
				if(length($mut_all[$i])>1){
					push(@HGVS_g_umd1,"g.".($pos+1)."_".$end."delins".substr($mut_all[$i],1));
					push(@HGVS_g_umd_reverse,"g.".($pos+1)."_".$end."delins".substr($mut_reverse[$i],1)) if($gene eq "BRCA1");
				}
				else{push(@HGVS_g_umd1,"g.".($pos+1)."del");}
			}
			elsif(length($ref_hg19)<length($mut_all[$i])){
				push(@HGVS_g_umd1,"g.".$pos."ins".substr($mut_all[$i],1));
				push(@HGVS_g_umd_reverse,"g.".$pos."ins".substr($mut_reverse[$i],1)) if($gene eq "BRCA1");
			}
		}

		if(length($ref_hg19)==length($mut_all[$i])){push(@HGVS_g_umd2,"g.".$pos.$ref_hg19.">".$mut_all[$i]);}
		else
		{
			if(length($ref_hg19)>length($mut_all[$i]))
			{
				if(length($mut_all[$i])>1){push(@HGVS_g_umd2,"g.".($pos+1)."_".$end."delins".substr($mut_all[$i],1));}
				else{push(@HGVS_g_umd2,"g.".($pos+1)."_".$end."del");}
			}
			elsif(length($ref_hg19)<length($mut_all[$i])){push(@HGVS_g_umd2,"g.".$pos."ins".substr($mut_all[$i],1));}
		}
	}
	my ($HGVS_g_umd1_1,$HGVS_g_umd1_2)=($HGVS_g_umd1[0],$HGVS_g_umd1[1]);
	my ($HGVS_g_umd2_1,$HGVS_g_umd2_2)=($HGVS_g_umd2[0],$HGVS_g_umd2[1]);
	my $HGVS_g_umd_reverse1="NA";
	my $HGVS_g_umd_reverse2="NA";
	if(scalar(@HGVS_g_umd_reverse)==2)
	{
		$HGVS_g_umd_reverse1=$HGVS_g_umd_reverse[0];
		$HGVS_g_umd_reverse2=$HGVS_g_umd_reverse[1];
	}
	###########
	
	###dbsnp:找到对应突变位点的rs编号###
	my $rs_id="NA";
	my $rs_id_output;
	$chr =~ s/chr//g;
	my $dbsnp=$db->prepare(qq{select * from dbsnp where chrom="$chr" and pos="$pos" and ref="$ref_hg19" and (alt like "%$mut1%" or alt like "%$mut2%")});
	$dbsnp->execute();
	my $i=0;
	while(my $ref=$dbsnp->fetchrow_hashref())
	{
		my $mut_db=$ref->{"alt"};
		my @mut_db=split /,/,$mut_db;
		foreach $mut_db (@mut_db){
			if($mut_db eq $mut1 || $mut_db eq $mut2){
				$rs_id=$ref->{"id"};
				$rs_id_output=$rs_id;
			}
		}
		$i++;
	}
	if($i==0){$rs_id="NA";$rs_id_output="-";}
	
	#print $rs_id,"\n";

	###BIC:临床分类信息(Clinical_Classification_Category)###
	my $bic_Clinical_Classification_Category;
	my $bic_Clinical_Classification_Evidence;
	my @Category;
	my @Evidence;
	#my $bic=$db->prepare(qq{select * from bic where HGVS_Genomic like "%$HGVS_g_bic_lovd1%" or HGVS_Genomic like "%$HGVS_g_bic_lovd2%" or dbSnp="$rs_id"});
	my $bic=$db->prepare(qq{select * from bic where ((HGVS_Genomic like "%$HGVS_g_bic_lovd1%" or HGVS_Genomic like "%$HGVS_g_bic_lovd2%") and dbSnp="$rs_id") or ((HGVS_Genomic="" or HGVS_Genomic="-") and dbSnp="$rs_id") or ((HGVS_Genomic like "%$HGVS_g_bic_lovd1%" or HGVS_Genomic like "%$HGVS_g_bic_lovd2%") and (dbSnp="" or dbSnp="-"))});  #根据HGVS_genomic或rs_id进行筛选,因为有些点并不存在rs编号
	$bic->execute();
	#print $bic->rows(),"\n";
	#print $HGVS_g_bic_lovd1,"\t",$HGVS_g_bic_lovd2,"\n";
	if($bic->rows()==0){$bic_Clinical_Classification_Category="-";$bic_Clinical_Classification_Evidence="-";}
	elsif($bic->rows()==1)
	{
		while (my $ref=$bic->fetchrow_hashref())
		{
			$bic_Clinical_Classification_Category=$ref->{"Clinical_Classification_Category"};
			$bic_Clinical_Classification_Evidence=$ref->{"Clinical_Classification_Evidence"};
		}
	}
	else
	{
		while (my $ref=$bic->fetchrow_hashref())
		{
			push(@Category,$ref->{"Clinical_Classification_Category"});
			push(@Evidence,$ref->{"Clinical_Classification_Evidence"});
		}
		$bic_Clinical_Classification_Category=join(";",@Category);
		$bic_Clinical_Classification_Evidence=join(";",@Evidence);
	}
	
	###ncbi_snp:MAF###
	my $ncbi_snp_maf;
	my $ncbi_snp=$db->prepare(qq{select * from ncbi_snp where rs_id="$rs_id"}); #ncbi_snp数据库中每条记录都有rs编号，这是最方便的查询方法
	$ncbi_snp->execute();
	while(my $ref=$ncbi_snp->fetchrow_hashref())
	{
		$ncbi_snp_maf=$ref->{"GLOBAL_MAF"} eq ""?"-":$ref->{"GLOBAL_MAF"};
		unless($ncbi_snp_maf eq '-'){
			$ncbi_snp_maf = $ncbi_snp_maf =~ /=(.*)\//;
		}
	}
	if($ncbi_snp->rows()==0){$ncbi_snp_maf="-";}


	###hapmap:AF频率###
	my $hapmap_AF;
	my $hapmap=$db->prepare(qq{select * from hapmap where dbsnp="$rs_id"}); #rs编号筛选
	$hapmap->execute();
	while (my $ref=$hapmap->fetchrow_hashref())
	{
		my $frequence=$ref->{"frequence"};
		if($frequence=~/AF/){my $AF=(split "=",(split ";",$frequence)[1])[1];$hapmap_AF=$AF;}
		else{$hapmap_AF="-";}
	}
	if($hapmap->rows()==0){$hapmap_AF="-";}

	###1000_hg19:AF频率，包括对不同地区的分类###
	my($hg19_1000_AF,$hg19_1000_EAS_AF,$hg19_1000_AMR_AF,$hg19_1000_AFR_AF,$hg19_1000_EUR_AF,$hg19_1000_SAS_AF);
	my $hg19_1000=$db->prepare(qq{select * from 1000_hg19 where pos="$pos" and ((ref="$ref_hg19" and (mut like "%$mut1%" or mut like "%$mut2%")) or (ref="$ref_hg19_reverse" and (mut like "%$mut_reverse1%" or mut like "%$mut_reverse2%")))});#1000人的数据库中存在一些mut不单一的情况，例如“A,C”等等
	$hg19_1000->execute();
	while(my $ref=$hg19_1000->fetchrow_hashref())
	{
		my $info=$ref->{"info"};
		$info=~/;AF=([\d\.]+)/;$hg19_1000_AF=$1;
		$info=~/;EAS_AF=([\d\.]+)/;$hg19_1000_EAS_AF=$1;
		$info=~/;AMR_AF=([\d\.]+)/;$hg19_1000_AMR_AF=$1;
		$info=~/;AFR_AF=([\d\.]+)/;$hg19_1000_AFR_AF=$1;
		$info=~/;EUR_AF=([\d\.]+)/;$hg19_1000_EUR_AF=$1;
		$info=~/;SAS_AF=([\d\.]+)/;$hg19_1000_SAS_AF=$1;
	}
	if($hg19_1000->rows()==0){$hg19_1000_AF="-";$hg19_1000_EAS_AF="-";$hg19_1000_AMR_AF="-";$hg19_1000_AFR_AF="-";$hg19_1000_EUR_AF="-";$hg19_1000_SAS_AF="-";}

	###NIFTY###
	my $nifty_maf;	#nifty数据库只输出对应的maf
	my $ref_hg19_mut1=$ref_hg19."/".$mut1;
	my $ref_hg19_mut2=$ref_hg19."/".$mut2;
	my $ref_hg19_mut_reverse1=$ref_hg19_reverse."/".$mut_reverse1;
	my $ref_hg19_mut_reverse2=$ref_hg19_reverse."/".$mut_reverse2;
	my $nifty=$db->prepare(qq{select * from nifty where chr="$chr" and pos="$pos" and (gene_type="$ref_hg19_mut1" or gene_type="$ref_hg19_mut2" or gene_type="$ref_hg19_mut_reverse1" or gene_type="$ref_hg19_mut_reverse2")});
	$nifty->execute();
	while(my $ref=$nifty->fetchrow_hashref())
	{
		$nifty_maf=$ref->{"maf"};
	}
	if($nifty->rows()==0){$nifty_maf="-";}
	#print $nifty_maf,"\n";
	
	###dbnsfp:ESP6500,SIFT,Polyphen2,EXAC######
	my($ESP6500_AA_AF,$ESP6500_EA_AF,$ExAC_AF,$ExAC_Adj_AF,$ExAC_AFR_AF,$ExAC_AMR_AF,$ExAC_EAS_AF,$ExAC_FIN_AF,$ExAC_NFE_AF,$ExAC_SAS_AF,$ExAC_nonTCGA_Adj_AF,$ExAC_nonTCGA_AFR_AF,$ExAC_nonTCGA_AMR_AF,$ExAC_nonTCGA_EAS_AF,$ExAC_nonTCGA_FIN_AF,$ExAC_nonTCGA_NFE_AF,$ExAC_nonTCGA_SAS_AF,$ExAC_nonpsych_AF,$ExAC_nonpsych_Adj_AF,$ExAC_nonpsych_AFR_AF,$ExAC_nonpsych_AMR_AF,$ExAC_nonpsych_EAS_AF,$ExAC_nonpsych_FIN_AF,$ExAC_nonpsych_NFE_AF,$ExAC_nonpsych_SAS_AF,$SIFT_pred,$Polyphen2_HDIV_pred,$Polyphen2_HVAR_pred,$LRT_pred,$MutationAssessor_pred,$MetaSVM_pred,$MetaLR_pred);
	my $dbnsfp=$db->prepare(qq{select * from dbnsfp where hg19_pos_1_based="$pos" and ((ref="$ref_hg19" and (alt="$mut1" or alt="$mut2")) or (ref="$ref_hg19_reverse" and (alt="$mut_reverse1" or alt="$mut_reverse2")))});#这个数据库记录的都是点突变，在搜索的时候可以考虑不需要使用rs_id
	$dbnsfp->execute();
	while(my $ref=$dbnsfp->fetchrow_hashref())
	{
		$ESP6500_AA_AF=$ref->{"ESP6500_AA_AF"}eq"."?"-":$ref->{"ESP6500_AA_AF"};
		$ESP6500_EA_AF=$ref->{"ESP6500_EA_AF"}eq"."?"-":$ref->{"ESP6500_EA_AF"};
		$ExAC_AF=$ref->{"ExAC_AF"}eq"."?"-":$ref->{"ExAC_AF"};
		$ExAC_Adj_AF=$ref->{"ExAC_Adj_AF"}eq"."?"-":$ref->{"ExAC_Adj_AF"};
		$ExAC_AFR_AF=$ref->{"ExAC_AFR_AF"}eq"."?"-":$ref->{"ExAC_AFR_AF"};
		$ExAC_AMR_AF=$ref->{"ExAC_AMR_AF"}eq"."?"-":$ref->{"ExAC_AMR_AF"};
		$ExAC_EAS_AF=$ref->{"ExAC_EAS_AF"}eq"."?"-":$ref->{"ExAC_EAS_AF"};
		$ExAC_FIN_AF=$ref->{"ExAC_FIN_AF"}eq"."?"-":$ref->{"ExAC_FIN_AF"};
		$ExAC_NFE_AF=$ref->{"ExAC_NFE_AF"}eq"."?"-":$ref->{"ExAC_NFE_AF"};
		$ExAC_SAS_AF=$ref->{"ExAC_SAS_AF"}eq"."?"-":$ref->{"ExAC_SAS_AF"};
		$ExAC_nonTCGA_Adj_AF=$ref->{"ExAC_nonTCGA_Adj_AF"}eq"."?"-":$ref->{"ExAC_nonTCGA_Adj_AF"};
		$ExAC_nonTCGA_AFR_AF=$ref->{"ExAC_nonTCGA_AFR_AF"}eq"."?"-":$ref->{"ExAC_nonTCGA_AFR_AF"};
		$ExAC_nonTCGA_AMR_AF=$ref->{"ExAC_nonTCGA_AMR_AF"}eq"."?"-":$ref->{"ExAC_nonTCGA_AMR_AF"};
		$ExAC_nonTCGA_EAS_AF=$ref->{"ExAC_nonTCGA_EAS_AF"}eq"."?"-":$ref->{"ExAC_nonTCGA_EAS_AF"};
		$ExAC_nonTCGA_FIN_AF=$ref->{"ExAC_nonTCGA_FIN_AF"}eq"."?"-":$ref->{"ExAC_nonTCGA_FIN_AF"};
		$ExAC_nonTCGA_NFE_AF=$ref->{"ExAC_nonTCGA_NFE_AF"}eq"."?"-":$ref->{"ExAC_nonTCGA_NFE_AF"};
		$ExAC_nonTCGA_SAS_AF=$ref->{"ExAC_nonTCGA_SAS_AF"}eq"."?"-":$ref->{"ExAC_nonTCGA_SAS_AF"};
		$ExAC_nonpsych_AF=$ref->{"ExAC_nonpsych_AF"}eq"."?"-":$ref->{"ExAC_nonpsych_AF"};
		$ExAC_nonpsych_Adj_AF=$ref->{"ExAC_nonpsych_Adj_AF"}eq"."?"-":$ref->{"ExAC_nonpsych_Adj_AF"};
		$ExAC_nonpsych_AFR_AF=$ref->{"ExAC_nonpsych_AFR_AF"}eq"."?"-":$ref->{"ExAC_nonpsych_AFR_AF"};
		$ExAC_nonpsych_AMR_AF=$ref->{"ExAC_nonpsych_AMR_AF"}eq"."?"-":$ref->{"ExAC_nonpsych_AMR_AF"};
		$ExAC_nonpsych_EAS_AF=$ref->{"ExAC_nonpsych_EAS_AF"}eq"."?"-":$ref->{"ExAC_nonpsych_EAS_AF"};
		$ExAC_nonpsych_FIN_AF=$ref->{"ExAC_nonpsych_FIN_AF"}eq"."?"-":$ref->{"ExAC_nonpsych_FIN_AF"};
		$ExAC_nonpsych_NFE_AF=$ref->{"ExAC_nonpsych_NFE_AF"}eq"."?"-":$ref->{"ExAC_nonpsych_NFE_AF"};
		$ExAC_nonpsych_SAS_AF=$ref->{"ExAC_nonpsych_SAS_AF"}eq"."?"-":$ref->{"ExAC_nonpsych_SAS_AF"};
		$SIFT_pred=$ref->{"SIFT_pred"}eq"."?"-":$ref->{"SIFT_pred"};
		$Polyphen2_HDIV_pred=$ref->{"Polyphen2_HDIV_pred"}eq"."?"-":$ref->{"Polyphen2_HDIV_pred"};
		$Polyphen2_HVAR_pred=$ref->{"Polyphen2_HVAR_pred"}eq"."?"-":$ref->{"Polyphen2_HVAR_pred"};
		$LRT_pred=$ref->{"LRT_pred"}eq"."?"-":$ref->{"LRT_pred"};
		$MutationAssessor_pred=$ref->{"MutationAssessor_pred"}eq"."?"-":$ref->{"MutationAssessor_pred"};
		$MetaSVM_pred=$ref->{"MetaSVM_pred"}eq"."?"-":$ref->{"MetaSVM_pred"};
		$MetaLR_pred=$ref->{"MetaLR_pred"}eq"."?"-":$ref->{"MetaLR_pred"};
	}
	if($dbnsfp->rows()==0)
	{
		($ESP6500_AA_AF,$ESP6500_EA_AF,$ExAC_AF,$ExAC_Adj_AF,$ExAC_AFR_AF,$ExAC_AMR_AF,$ExAC_EAS_AF,$ExAC_FIN_AF,$ExAC_NFE_AF,$ExAC_SAS_AF,$ExAC_nonTCGA_Adj_AF,$ExAC_nonTCGA_AFR_AF,$ExAC_nonTCGA_AMR_AF,$ExAC_nonTCGA_EAS_AF,$ExAC_nonTCGA_FIN_AF,$ExAC_nonTCGA_NFE_AF,$ExAC_nonTCGA_SAS_AF,$ExAC_nonpsych_AF,$ExAC_nonpsych_Adj_AF,$ExAC_nonpsych_AFR_AF,$ExAC_nonpsych_AMR_AF,$ExAC_nonpsych_EAS_AF,$ExAC_nonpsych_FIN_AF,$ExAC_nonpsych_NFE_AF,$ExAC_nonpsych_SAS_AF,$SIFT_pred,$Polyphen2_HDIV_pred,$Polyphen2_HVAR_pred,$LRT_pred,$MutationAssessor_pred,$MetaSVM_pred,$MetaLR_pred)=("-") x 32;
		
	}

	###clinvar:Clinical Significance###
	my $Clinical_Significance;
	my $clinvar=$db->prepare(qq{select * from ncbi_clinvar where rsid="$rs_id"});#每一行都有rs编号，并且不会存在rs编号重复的情况
	$clinvar->execute();
	while(my $ref=$clinvar->fetchrow_hashref())
	{
		$Clinical_Significance=$ref->{"Clinical_Significance"}eq""?"-":$ref->{"Clinical_Significance"};
	}
	if($clinvar->rows()==0){$Clinical_Significance="-";}

	###UMD###
	#print $HGVS_g,"\n";
	my $umd_biological_significance="-";
	my $umd=$db->prepare(qq{select * from umd where genomic like "%$HGVS_g_umd1_1%" or genomic like "%$HGVS_g_umd1_2%" or genomic like "%$HGVS_g_umd2_1%" or genomic like "%$HGVS_g_umd2_2%" or genomic like "%$HGVS_g_umd_reverse1%" or genomic like "%$HGVS_g_umd_reverse2%"}); #该数据库没有rs编号的信息，只能通过genomic的变化去查找，不过需要对genomic进行形式上的调整
	$umd->execute();
	while(my $ref=$umd->fetchrow_hashref())
	{
		if($ref->{"genomic"}=~/\>/ or $ref->{"genomic"}=~/del$/)
		{
			$umd_biological_significance=$ref->{"Biological_significance"};
		}
		elsif($ref->{"genomic"}=~/del([ATCG]+)/)
		{
			$umd_biological_significance=$ref->{"Biological_significance"} if(length($1)==(length($ref_hg19)-1));
		}
		elsif($ref->{"genomic"}=~/ins([ATCG]+)/)
		{
			$umd_biological_significance=$ref->{"Biological_significance"} if(substr($mut1,1) eq $1 or substr($mut_reverse1,1) eq $1 or substr($mut2,1) eq $1 or substr($mut_reverse2,1) eq $1);
		}
		elsif($ref->{"genomic"}=~/del(\d+)/ or $ref->{"genomic"}=~/ins(\d+)/)
		{
			$umd_biological_significance=$ref->{"Biological_significance"} if($1==(length($ref_hg19)-1));
		}
		$umd_biological_significance=~s/\s//g;
	}
	if($umd->rows()==0){$umd_biological_significance="-";}
	#print $umd_biological_significance,"\n";

	###BRCA_inhouse###
	#print $HGVS_c,"\n";
	my $clinical_category;
#	my $Description;
#	my $Suggestion;
#	my $Publication;
	#my $utf8=$db->prepare("set names utf8");
	#$utf8->execute();
	my $brca_inhouse=$db->prepare(qq{select * from brca_in_house where c_HGVS="$cHGVS" and gene="$gene"}); #该数据库既没有genomic的变化又没有rs编号，只能通过cds去匹配
=h
	if($cHGVS eq "")
	{
		$brca_inhouse=$db->prepare(qq{select * from brca_in_house where c_HGVS="$cHGVS" and gene="$gene"}); #该数据库既没有genomic的变化又没有rs编号，只能通过cds去匹配
	}
	else
	{
		$brca_inhouse=$db->prepare(qq{select * from brca_in_house where c_HGVS like "%$cHGVS%" and gene="$gene"}); #该数据库既没有genomic的变化又没有rs编号，只能通过cds去匹配
		
	}
=cut
	$brca_inhouse->execute();
	while(my $ref=$brca_inhouse->fetchrow_hashref())
	{
		$clinical_category=$ref->{"clinical_category"};
#		$Description=$ref->{"Description"};
#		$Suggestion=$ref->{"Suggestion"};
#		$Publication=$ref->{"Publication"};
	}
	if($brca_inhouse->rows()==0){$clinical_category="-";}
	#print $clinical_category,"\n";


	###LOVD:需要的信息，不同的方法会得到不同的结果，可以整合到一起进行输出###
	my $lovd=$db->prepare(qq{select * from lovd where ((DNA_change_genomic_hg19 like "%$HGVS_g_bic_lovd1%" or DNA_change_genomic_hg19  like "%$HGVS_g_bic_lovd2%") and dbSNP_IN="$rs_id") or ((DNA_change_genomic_hg19 like "%$HGVS_g_bic_lovd1%" or DNA_change_genomic_hg19 like "%$HGVS_g_bic_lovd2%") and (dbSNP_IN="-" or dbSNP_IN="")) or ((DNA_change_genomic_hg19="" or DNA_change_genomic_hg19="-") and dbSNP_IN="$rs_id")}); #该数据库每一行都存在genomic的记录,部分存在rs编号
	$lovd->execute();
	my @FunctionalAnalysis_or_Technique;
	my @FunctionalAnalysis_or_Result;
	my $lovd_result;
	my @lovd_result_arry;
	while(my $ref=$lovd->fetchrow_hashref())
	{
		my $technique=$ref->{"FunctionalAnalysis_or_Technique"};
		my $result=$ref->{"FunctionalAnalysis_or_Result"};
		if($technique ne "-" and $technique ne "?" and $result ne "-" and $result ne "?")
		{
			push(@FunctionalAnalysis_or_Technique,$technique);
			push(@FunctionalAnalysis_or_Result,$result);
		}
	}
	if(scalar(@FunctionalAnalysis_or_Technique)==0){$lovd_result="-";}
	else
	{
		foreach my $j(0..(scalar(@FunctionalAnalysis_or_Technique)-1))
		{
			push(@lovd_result_arry,$FunctionalAnalysis_or_Technique[$j].":".$FunctionalAnalysis_or_Result[$j]);
		}
		$lovd_result=join(";",@lovd_result_arry);
	}
	#print $lovd_result,"\n";
	#print scalar(@FunctionalAnalysis_or_Technique),"\n";

	###dbscsnv###
	my $dbscsnv_ada_score;
	my $dbscsnv_rf_score;
	my $dbscsnv=$db->prepare(qq{select * from dbscsnv where pos="$pos" and ((ref="$ref_hg19" and (alt="$mut1" or alt="$mut2")) or (ref="$ref_hg19_reverse" and (alt="$mut_reverse1" or alt="$mut_reverse2")))});#这个数据库记录的都是点突变，在搜索的时候可以考虑不使用rs_id
	$dbscsnv->execute();
	while(my $ref=$dbscsnv->fetchrow_hashref())
	{
		$dbscsnv_ada_score=$ref->{"ada_score"} eq "."?"-":$ref->{"ada_score"}; 
		$dbscsnv_rf_score=$ref->{"rf_score"} eq "."?"-":$ref->{"rf_score"};
	}
	if($dbscsnv->rows()==0){$dbscsnv_ada_score="-";$dbscsnv_rf_score="-";}
my $final_class;
if ($clinical_category ne '-'){
$final_class = $clinical_category;
}
elsif($bic_Clinical_Classification_Category ne '-'){
$final_class = $bic_Clinical_Classification_Category;
}
elsif($umd_biological_significance ne '-'){
$final_class = $umd_biological_significance;
}
else{
if($Mtype eq 'CNV'){
$final_class = 'Class4-Likely_pathogenic';
}
elsif ($function =~ /Missense/i){
$final_class = 'Class3-VUS';
}
elsif ($function =~ /Nonsense/i || $function =~ /stop_gain/i){
$final_class = 'Class5-Pathogenic';
}
else{
$final_class = 'Class2-Likely_non_pathogenic ';
}}

if($function =~ /^synonymous_variant/)   ### add by lishiyong
{
	$final_class = 'Class1-Non_pathogenic';
}
 
my ($final_pred,@final_pred,%hash_pred);
$final_pred = join "",$SIFT_pred,$Polyphen2_HDIV_pred,$Polyphen2_HVAR_pred,$LRT_pred,$MutationAssessor_pred,$MetaSVM_pred,$MetaLR_pred;

$final_pred =~ s/\.//g;
$final_pred =~ s/\;//g;
$final_pred =~ s/\s//g;
$final_pred =~ s/\-//g;
unless($final_pred){
	$final_pred = '-';
}
@final_pred = split //,$final_pred;
foreach $final_pred (@final_pred){
	$hash_pred{$final_pred}+=1;
}
foreach my $i (sort {$hash_pred{$b}<=>$hash_pred{$a}} keys %hash_pred){
	$final_pred = $i;
	last;
}
if($final_pred eq '-'){
	$final_pred = 'A';
}
	print "$gene\t$trans\tchr$chr\t$pos\t$end\t$exon\t$ref_hg19\t$mut\t$AD\t$Ratio\t$Mtype\t$GType\t$HGVS_g\t$cHGVS\t$pHGVS\t$function\t$rs_id_output\t$final_class\t$final_pred\t$hapmap_AF\t$ncbi_snp_maf\t$nifty_maf\t$hg19_1000_AF\t$hg19_1000_EAS_AF\t$hg19_1000_AMR_AF\t$hg19_1000_AFR_AF\t$hg19_1000_EUR_AF\t$hg19_1000_SAS_AF\t";
	print "$ESP6500_AA_AF\t$ESP6500_EA_AF\t$ExAC_AF\t$ExAC_Adj_AF\t$ExAC_AFR_AF\t$ExAC_AMR_AF\t$ExAC_EAS_AF\t$ExAC_FIN_AF\t$ExAC_NFE_AF\t$ExAC_SAS_AF\t$ExAC_nonTCGA_Adj_AF\t$ExAC_nonTCGA_AFR_AF\t$ExAC_nonTCGA_AMR_AF\t$ExAC_nonTCGA_EAS_AF\t$ExAC_nonTCGA_FIN_AF\t$ExAC_nonTCGA_NFE_AF\t$ExAC_nonTCGA_SAS_AF\t$ExAC_nonpsych_AF\t$ExAC_nonpsych_Adj_AF\t$ExAC_nonpsych_AFR_AF\t$ExAC_nonpsych_AMR_AF\t$ExAC_nonpsych_EAS_AF\t$ExAC_nonpsych_FIN_AF\t$ExAC_nonpsych_NFE_AF\t$ExAC_nonpsych_SAS_AF\t$SIFT_pred\t$Polyphen2_HDIV_pred\t$Polyphen2_HVAR_pred\t$LRT_pred\t$MutationAssessor_pred\t$MetaSVM_pred\t$MetaLR_pred\t";
	print "$clinical_category\t$bic_Clinical_Classification_Category\t$umd_biological_significance\t$Clinical_Significance\t$lovd_result\t$dbscsnv_ada_score\t$dbscsnv_rf_score\n";
#	print "$Description\t$Suggestion\t$Publication\n";
	$final_pred = undef;
	
	######################################################################################

}
$db->disconnect; #断开数据库的连接


#连接数据库
sub connect_mysql   
{
	my $dsn = "DBI:mysql:brcav2:192.168.32.150:9906";
	my $user = "ionadmin";
	my $pass = "ionadmin";
	my $dbh = DBI->connect($dsn, $user, $pass, {RaiseError=>0, PrintError => 1}) or die "Could not connect to mysql server: $DBI::err($DBI::errstr)\n";
}
