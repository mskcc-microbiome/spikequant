#/bin/bash
set -eoux
outdir=$1

mkdir -p ${outdir}

curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/013/045/GCA_000013045.1_ASM1304v1/GCA_000013045.1_ASM1304v1_genomic.fna.gz | gunzip > ${outdir}/Salinibacter_ruber.fa
curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/167/675/GCF_000167675.1_v2.0/GCF_000167675.1_v2.0_genomic.fna.gz | gunzip  > ${outdir}/Trichoderma_reesei.fa
curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/223/905/GCF_000223905.1_ASM22390v1/GCF_000223905.1_ASM22390v1_genomic.fna.gz | gunzip  > ${outdir}/Haloarcula_hispanica.fa
