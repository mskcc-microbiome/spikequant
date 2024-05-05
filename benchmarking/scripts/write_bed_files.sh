#/bin/bash
set -eux
set -o pipefail
outdir=$1

mkdir -p ${outdir}

curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/013/045/GCA_000013045.1_ASM1304v1/GCA_000013045.1_ASM1304v1_genomic.gff.gz | gunzip | grep  $'\trRNA\t' > ${outdir}/Salinibacter_ruber_rrna.bed
# create genome file
samtools faidx benchmarking/spikes/raw/Salinibacter_ruber.fa
# create full genome bed
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' benchmarking/spikes/raw/Salinibacter_ruber.fa.fai > ${outdir}/Salinibacter_ruber_full.bed
# create first chromosome bed
head -n 1 ${outdir}/Salinibacter_ruber_full.bed > ${outdir}/Salinibacter_ruber_primary.bed
# create non-rRNA bed
bedtools complement -i ${outdir}/Salinibacter_ruber_rrna.bed -g benchmarking/spikes/raw/Salinibacter_ruber.fa.fai > ${outdir}/Salinibacter_ruber_nonrrna.bed
# create marker bed
awk '{{ if ($3 >= $2) print  }}' benchmarking/spikes/markers/Salinibacter_ruber.bed > hh.tmp
    # extract and flip th incorrectly oriented regions
    awk '{{ if ($3 < $2) print  }}' benchmarking/spikes/markers/Salinibacter_ruber.bed  | awk '{{ print $1 "\t" $3 "\t" $2 }}' >> hh.tmp
    # and sort
    sort -k 1,1 -k2,2n hh.tmp > ${outdir}/Salinibacter_ruber_markers.bed
    rm hh.tmp





# Trichoderma reesei is left alone as it doesn't have any annotated 16S - we could consider masking ITS etc
samtools faidx benchmarking/spikes/raw/Trichoderma_reesei.fa
# create full genome bed
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' benchmarking/spikes/raw/Trichoderma_reesei.fa.fai > ${outdir}/Trichoderma_reesei_full.bed
# create first chromosome bed
head -n 1 ${outdir}/Trichoderma_reesei_full.bed > ${outdir}/Trichoderma_reesei_primary.bed




# getting CBA1121 instead of Haloarcula hispanica ATCC 33960  for the same reason.
curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/729/095/GCF_008729095.1_ASM872909v1/GCF_008729095.1_ASM872909v1_genomic.gff.gz | gunzip | grep  $'\trRNA\t' > ${outdir}/Haloarcula_hispanica_rrna.bed
# create genome file
samtools faidx benchmarking/spikes/raw/Haloarcula_hispanica.fa
# create full genome bed
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' benchmarking/spikes/raw/Haloarcula_hispanica.fa.fai > ${outdir}/Haloarcula_hispanica_full.bed
# create first chromosome bed
head -n 1 ${outdir}/Haloarcula_hispanica_full.bed > ${outdir}/Haloarcula_hispanica_primary.bed
# create non-rRNA bed
bedtools complement -i ${outdir}/Haloarcula_hispanica_rrna.bed -g benchmarking/spikes/raw/Haloarcula_hispanica.fa.fai > ${outdir}/Haloarcula_hispanica_nonrrna.bed
# create marker bed
awk '{{ if ($3 >= $2) print  }}' benchmarking/spikes/markers/Haloarcula_hispanica.bed > hh.tmp
    # extract and flip th incorrectly oriented regions
    awk '{{ if ($3 < $2) print  }}' benchmarking/spikes/markers/Haloarcula_hispanica.bed  | awk '{{ print $1 "\t" $3 "\t" $2 }}' >> hh.tmp
    # and sort
    sort -k 1,1 -k2,2n hh.tmp > ${outdir}/Haloarcula_hispanica_markers.bed
    rm hh.tmp
