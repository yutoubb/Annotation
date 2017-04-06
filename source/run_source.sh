#!/bin/bash
if [ ! -n "$1" ] ;then  
    echo "Please specify the name of the input vcf file!" 
else
    VCF=$1
    Sample=${VCF%%.*}
    perl /opt/bin/01.vcf_bed.pl -i /mnt/work/job/input/$VCF && \
    perl /opt/bin/ensembl-tools-release-87/scripts/variant_effect_predictor/variant_effect_predictor.pl --refseq --offline --dir_cache /opt/data/VEP --fasta /opt/data/hg19/hg19.fasta --everything --force_overwrite --tab -i /mnt/work/job/input/$Sample\.BRCA.vcf -o /mnt/work/job/output/$Sample && \
    perl /opt/bin/03.variant_effect_predictor.BRCA.pl -i /mnt/work/job/input/$Sample\.BRCA.vcf -ens /mnt/work/job/output/$Sample -o /mnt/work/job/output/$Sample\.anno && \
    perl /opt/bin/liuzu.read_mysql.T.pl /mnt/work/job/output/$Sample\.anno > /mnt/work/job/output/$Sample\.anno.db && \
    echo "Annotation finished! Result saved in file "$Sample".anno.db!"
fi
