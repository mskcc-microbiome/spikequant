#/bin/bash
set -eux
set -o pipefail
outdir=$1

for prog in blastn barrnap seqkit bowtie2-inspect;
do
    which  $prog
done

mkdir -p ${outdir}

# get markers
if [ ! -f "${outdir}/choco_2024.fa" ] ; then
bowtie2-inspect   /data/brinkvd/resources/dbs/metaphlan/mpa_vJun23_CHOCOPhlAnSGB_202403/mpa_vJun23_CHOCOPhlAnSGB_202403 > ${outdir}/choco_2024.fa
# sruber is SGB15515
# hispanica is SGB225
seqkit grep   --by-name --use-regexp --pattern "SGB15515$" ${outdir}/choco_2024.fa > ${outdir}/sruber_markers.fa
seqkit grep   --by-name --use-regexp --pattern "SGB225$" ${outdir}/choco_2024.fa > ${outdir}/hhispanica_markers.fa
fi


blastn -subject benchmarking/spikes/raw/Salinibacter_ruber.fa -query ${outdir}/sruber_markers.fa -out ${outdir}/tmp_sruber.bed -outfmt "6 sacc sstart send"   -max_target_seqs 1
blastn -subject benchmarking/spikes/raw/Salinibacter_ruber.fa -query ${outdir}/hhispanica_markers.fa -out ${outdir}/tmp_hhispanica.bed -outfmt "6 sacc sstart send"   -max_target_seqs 1


clean_bed(){
    bed=$1
    newbed=$2
    # get the correct oriented ones
    awk '{{ if ($3 >= $2) print  }}' $bed  > ${newbed}_tmp
    # extract and flip th incorrectly oriented regions
    awk '{{ if ($3 < $2) print  }}' $bed  | awk '{{ print $1 "\t" $3 "\t" $2 }}' >> ${newbed}_tmp
    # and sort
    sort -k 1,1 -k2,2n ${newbed}_tmp > ${newbed}
    rm ${newbed}_tmp
}

clean_bed ${outdir}/tmp_sruber.bed ${outdir}/Salinibacter_ruber_markers.bed
clean_bed ${outdir}/tmp_hhispanica.bed ${outdir}/Haloarcula_hispanica_markers.bed

make_full_primary_and_no_rrna_beds(){
    fa=$1
    outpre=$2
    kingdom=$3
    samtools faidx $fa
    # this is specifically to fix the h hisp where the accessions are out of order
    sort ${fa}.fai  > ${fa}.fai.sorted
    barrnap --kingdom $kingdom $fa | cut -f 1,4,5 | grep -v gff > ${outpre}_rrna.bed.tmp
    # create full genome bed
    awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${fa}.fai.sorted | sort -k 1,1 -k2,2n > ${outpre}_full.bed
    # create first chromosome bed; in the case of h hispanica the sort order is misleading, so we use the .fai order
#    head -n 1 ${outpre}_full.bed > ${outpre}_primary.bed
    awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${fa}.fai | head -n 1 > ${outpre}_primary.bed
    # create non-rRNA bed
    bedtools complement -i ${outpre}_rrna.bed.tmp -g ${fa}.fai.sorted > ${outpre}_nonrrna.bed

}

make_full_primary_and_no_rrna_beds  benchmarking/spikes/raw/Salinibacter_ruber.fa ${outdir}/Salinibacter_ruber bac
make_full_primary_and_no_rrna_beds  benchmarking/spikes/raw/Trichoderma_reesei.fa ${outdir}/Trichoderma_reesei euk
make_full_primary_and_no_rrna_beds  benchmarking/spikes/raw/Haloarcula_hispanica.fa ${outdir}/Haloarcula_hispanica  arc





rm ${outdir}/*.tmp
