#!/bin/bash
set -x -e


scratch="/sc/arion/scratch/arayan01/projects/chipseq_tut"
sampledir="${scratch}/data/samples"
adapter_path="/sc/arion/work/arayan01/test-env/envs/atacseq/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa"

read1="/sc/arion/scratch/arayan01/projects/chipseq_tut/data/samples/wt_H3K4me3_read1.fastq.gz"
read2="/sc/arion/scratch/arayan01/projects/chipseq_tut/data/samples/wt_H3K4me3_read2.fastq.gz"


function fastqc {
    rm -rf "${sampledir}/fastqc"
    mkdir -p "${sampledir}/fastqc"
    fastqc -t 2 "${read1}" -o "${sampledir}/fastqc"
    fastqc -t 2 "${read2}" -o "${sampledir}/fastqc"
}

function trimming {
    rm -rf "${sampledir}/trimming"
    mkdir -p "${sampledir}/trimming"
    out="${sampledir}/trimming"
    trim_galore --quality 15 --stringency 3 --report --paired -o ${out} ${read1} ${read2}
    #quality - Trim low-quality ends from reads in addition to adapter removal, default is 20
    #stringency - Overlap with adapter sequence required to trim a sequence
    #Defaults to a very stringent setting of 1, i.e. even a single base pair of overlapping sequence will be trimmed of the 3' end of any read.
    #report - Generate a report file
} 
function mapping {
    rm -rf ${scratch}/results/mapping
    mkdir "${scratch}/results/mapping"
    mapdir="${scratch}/results/mapping"
    refgen="${scratch}/data/refgen/mm10"
    output="${scratch}/results/mapping"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "mapping_job.txt" \
        bowtie2 -x ${refgen} \
            -1 ${read1} \
            -2 ${read2} \
            -S ${mapdir}/output.sam \
            --met-file ${mapdir}/bowtie2_mapping_stats.txt
    #-x is our refrence genome in this case (Mouse)
    #-1 and -2 are our reads from our trimming output
    #--met-file is our statistics of our mapping 
}

function chipqc {
    rm -rf ${scratch}/results/chipqc
    mkdir -p ${scratch}/results/chipqc
    
    output="${scratch}/results/chipqc"
    data="${scratch}/data/chipseqdata"

    #make sure these files are indexed for input into multibamsum.

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "chipqc_job.txt" \
        multiBamSummary bins -b ${data}/*.bam \
        -o ${output}/output_matrix.npz \
        -bs 1000 \
        -n 500 \
        -r "chrX" 
}
function corplot {
    rm -rf "${scratch}/results/plot/correlationplot"
    mkdir -p "${scratch}/results/plot/correlationplot"

    data="${scratch}/results/chipqc/output_matrix.npz"
    output="${scratch}/results/plot/correlationplot.png"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "corplot_job.txt" \
        plotCorrelation -in ${data} \
        -c pearson \
        -p heatmap \
        -o ${output}
}

function ipplot {
    data="${scratch}/data/chipseqdata"
    output="${scratch}/results/plot/ip_plot.png"

    experiment="${data}/wt_input_rep1.bam"
    control="${data}/wt_H3K4me3_rep1.bam"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "ipplot_job.txt" \
        plotFingerprint -b ${experiment} ${control} \
        -plot ${output} \
        -r "chrX" \
        -v \
        -bs 10000
}

function norm_stats {
    rm -rf "${scratch}/results/norm"
    mkdir -p "${scratch}/results/norm"

    data="${scratch}/data/chipseqdata"
    
    experiment="${data}/wt_input_rep1.bam"
    control="${data}/wt_H3K4me3_rep1.bam"
    outdir="${scratch}/results/norm"
    
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "norm_job1.txt" \
        "samtools idxstats "$experiment" > ${outdir}/input_rep.idxstats.txt"
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "norm_job2.txt" \
        "samtools idxstats "$control" > ${outdir}/4me3_rep1.idxstats.txt"
}
function norm_filecov {
    data="${scratch}/data/chipseqdata"

    experiment="${data}/wt_input_rep1.bam"
    control="${data}/wt_H3K4me3_rep1"
    outdir="${scratch}/results/norm"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "norm_filecov_bedgraph_job1.txt" \
        bamCoverage --bam ${control}.bam -o ${outdir}/wt_H3K4me3_rep1.bw --binSize 25 --normalizeUsing RPGC --effectiveGenomeSize 2308125349 -r chrX 
}
#fastqc
#trimming
#mapping
#chipqc
#corplot
#ipplot
#norm_stats
norm_filecov
