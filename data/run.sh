set -x -e
scratch="/sc/arion/scratch/arayan01/projects/chipseq_tut"

function rawdata {
    rm -rf ${scratch}/data/samples
    rm -rf ${scratch}/data/refgen
    mkdir -p ${scratch}/data/samples
    mkdir -p ${scratch}/data/refgen
    
    sampledir="${scratch}/data/samples"
    refgendir="${scratch}/data/refgen"

    mkdir -p "${sampledir}" "${refgendir}"

    while read -r sample url; do
        if [[ "$sample" == "mm10" ]]; then
            wget -c -O "${refgendir}/${sample}.zip" "${url}"
        else
            wget -c -O "${sampledir}/${sample}.fastq.gz" "${url}"
        fi
    done < url_links.txt
}

function chipseq_data {
    rm -rf ${scratch}/data/chipseqdata
    mkdir -p ${scratch}/data/chipseqdata
    chipseq="${scratch}/data/chipseqdata"

    while read -r sample url; do
        wget -c -O "${chipseq}/${sample}.bam" "${url}"
    done < chipseq_samples.txt
}

#rawdata
chipseq_data