#/bin/bash
set -eoux
outdir=$1

mkdir -p ${outdir}

curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/013/045/GCA_000013045.1_ASM1304v1/GCA_000013045.1_ASM1304v1_genomic.fna.gz | gunzip > ${outdir}/Salinibacter_ruber.fa
# getting Trichoderma reesei  QM6a instead of  Trichoderma reesei ATCC 13631 .  Might be fine license-wise but not sure
curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/167/675/GCF_000167675.1_v2.0/GCF_000167675.1_v2.0_genomic.fna.gz | gunzip  > ${outdir}/Trichoderma_reesei.fa
# getting CBA1121 instead of Haloarcula hispanica ATCC 33960  for the same reason.
curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/729/095/GCF_008729095.1_ASM872909v1/GCF_008729095.1_ASM872909v1_genomic.fna.gz | gunzip  > ${outdir}/Haloarcula_hispanica.fa
