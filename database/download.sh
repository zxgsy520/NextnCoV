#wget -C ftp://download.big.ac.cn/GVM/Coronavirus/vcf/2019-nCoV_total.vcf
/nextomics/Pipeline/NextnCoV/v1.0.0/scripts/stat_group_vcf.py --fasta reference.fasta --vcf 2019-nCoV_total.vcf >mutation_rate.tsv

