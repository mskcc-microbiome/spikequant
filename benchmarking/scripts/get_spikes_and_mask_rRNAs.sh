#/bin/bash
set -eux
set -o pipefail
outdir=$1

mkdir -p ${outdir}

curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/013/045/GCA_000013045.1_ASM1304v1/GCA_000013045.1_ASM1304v1_genomic.gff.gz | gunzip | grep  $'\trRNA\t' > ${outdir}/Salinibacter_ruber.bed
bedtools maskfasta -fi $(dirname $outdir)/split/Salinibacter_ruber_primary.fa -bed ${outdir}/Salinibacter_ruber.bed -fo ${outdir}/Salinibacter_ruber_primary_masked.fa
bedtools maskfasta -fi $(dirname $outdir)/split/Salinibacter_ruber_other.fa -bed ${outdir}/Salinibacter_ruber.bed -fo ${outdir}/Salinibacter_ruber_other_masked.fa

# Trichoderma reesei is left alone as it doesn't have any annotated 16S - we could consider masking ITS etc
cp $(dirname $outdir)/split/Trichoderma_reesei_*.fa $outdir
#curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/167/675/GCF_000167675.1_v2.0/GCF_000167675.1_v2.0_genomic.fna.gz | gunzip > ${outdir}/Trichoderma_reesei.fa

# getting CBA1121 instead of Haloarcula hispanica ATCC 33960  for the same reason.
curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/729/095/GCF_008729095.1_ASM872909v1/GCF_008729095.1_ASM872909v1_genomic.gff.gz | gunzip | grep  $'\trRNA\t' > ${outdir}/Haloarcula_hispanica.bed
bedtools maskfasta -fi $(dirname $outdir)/split/Haloarcula_hispanica_primary.fa -bed ${outdir}/Haloarcula_hispanica.bed -fo ${outdir}/Haloarcula_hispanica_primary_masked.fa
bedtools maskfasta -fi $(dirname $outdir)/split/Haloarcula_hispanica_other.fa -bed ${outdir}/Haloarcula_hispanica.bed -fo ${outdir}/Haloarcula_hispanica_other_masked.fa
#cp $(dirname $outdir)/split/Haloarcula_hispanica_other.fa $outdir
