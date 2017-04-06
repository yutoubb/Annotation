#!/bin/bash
if [ ! -n "$1" ] ;then  
    echo "Please specify the name of the input vcf file!" 
else
    VCF=$1
    Sample=${VCF%%.*}
    /opt/bin/vcf_bed.exe -i /mnt/work/job/input/$VCF && \
    perl /opt/bin/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl --refseq --offline --dir_cache /opt/data/VEP --fasta /opt/data/hg19/hg19.fasta --everything --force_overwrite --tab -i /mnt/work/job/input/$Sample\.BRCA.vcf -o /mnt/work/job/output/$Sample && \
    /opt/bin/vep_BRCA.exe -i /mnt/work/job/input/$Sample\.BRCA.vcf -ens /mnt/work/job/output/$Sample -o /mnt/work/job/output/$Sample\.anno && \
    /opt/bin/read_mysql.exe /mnt/work/job/output/$Sample\.anno > /mnt/work/job/output/$Sample\.anno.db && \
    echo "Annotation finished! Result saved in file "$Sample".anno.db!"
fi
