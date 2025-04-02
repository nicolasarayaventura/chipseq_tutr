#!/bin/bash
set -x -e


scratch="/sc/arion/scratch/arayan01/projects/chipseqtut"
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

    control="${data}/wt_input_rep1.bam"
    experiment="${data}/wt_H3K4me3_rep1.bam"

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
    
    control="${data}/wt_input_rep1.bam"
    experiment="${data}/wt_H3K4me3_rep1.bam"
    outdir="${scratch}/results/norm"
    
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "norm_job1.txt" \
        "samtools idxstats "$control" > ${outdir}/input_rep.idxstats.txt"
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "norm_job2.txt" \
        "samtools idxstats "$experiment" > ${outdir}/4me3_rep1.idxstats.txt"
}
function bamcov {
    data="${scratch}/data/chipseqdata"

    control="${data}/wt_input_rep1.bam"
    experiment="${data}/wt_H3K4me3_rep1.bam"
    outdir="${scratch}/results/norm"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "norm_filecov_bedgraph_job1.txt" \
        "bamCoverage -b ${experiment} -o ${outdir}/wt_H3K4me3_rep1.bw \
        --outFileFormat bigwig \
        --binSize 25 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2308125349 \
        --region chrX \
        --verbose"
}

function bamcomp {
    data="${scratch}/data/chipseqdata"

    control="${data}/wt_input_rep1.bam"
    experiment="${data}/wt_H3K4me3_rep1.bam"
    outdir="${scratch}/results/norm"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "norm_comp_job1.txt" \
        bamCompare -b1 ${experiment} -b2 ${control} -o ${outdir}/comp.bw \
        -bs 50 \
        -of bigwig \
        -r chrX \
        --operation log2
}

function peakcalling {
    rm -rf ${scratch}/results/peakcalling
    mkdir "${scratch}/results/peakcalling"
    
    data="${scratch}/data/chipseqdata"

    control="${data}/wt_input_rep1.bam"
    experiment="${data}/wt_H3K4me3_rep1.bam"
    outdir="${scratch}/results/norm"
    
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "peakcalling_job1.txt" \
        macs2 callpeak -t ${experiment} \
        -c ${control} \
        -f BAMPE \
        --g mm \
        --name wt_H3K4me3_rep1_peak \
        --format BAM \
        --tsize 75 \
        --outdir ${scratch}/results/peakcalling

}

function prep1 {
    rm -rf ${scratch}/results/signalcomp_plots
    mkdir "${scratch}/results/signalcomp_plots"

    data="${scratch}/data/chipseqdata"
    experiment="${data}/wt_CTCF_rep1.bam"
    control="${data}/wt_input_rep1.bam"
    
    outdir="${scratch}/results/signalcomp_plots"

    # Run bamCompare
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${outdir}/prep1_job.txt" \
        bamCompare -b1 ${experiment} -b2 ${control} -o ${outdir}/wt_CTCF_rep1.bw \
            -bs 50 \
            -of bigwig \
            -r chrX \
            --operation log2
}

function prep2 {
    data="${scratch}/data/chipseqdata"
    experiment="${data}/wt_CTCF_rep1.bam"
    control="${data}/wt_input_rep1.bam"
    
    outdir="${scratch}/results/signalcomp_plots"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${outdir}/prep2a_job.txt" \
        macs2 callpeak -t ${experiment} -c ${control} -f BAMPE --g mm --format BAM --outdir ${outdir} -n wt_CTCF_rep1
}

function prep3 {
    data="${scratch}/data/chipseqdata"
    experiment="${data}/wt_CTCF_rep1.bam"
    experiment2="${scratch}/results/peakcalling/wt_H3K4me3_rep1_peak_peaks.narrowPeak"
    control="${data}/wt_input_rep1.bam"

    outdir="${scratch}/results/signalcomp_plots"

    #Concatenate
    cat ${outdir}/wt_CTCF_rep1_peaks.narrowPeak ${experiment2} > ${outdir}/concatenated_peaks.bed
    concatenated="${outdir}/concatenated_peaks.bed"

    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${outdir}/prep3a_job.txt" \
        "bedtools sort -i ${concatenated} > ${outdir}/sorted_peaks.bed"
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${outdir}/prep3b_job.txt" \
        "bedtools merge -i ${outdir}/sorted_peaks.bed > ${outdir}/merged_peaks.bed"
}
function heatmap {
    data="${scratch}/data/chipseqdata"
    outdir="${scratch}/results/signalcomp_plots"

    bw1="${outdir}/wt_CTCF_rep1.bw"
    bw2="${scratch}/results/peakcalling/wt_H3K4me3_rep1_peak_peaks.narrowPeak"
    
    merged_peaks="${outdir}/merged_peaks.bed"
#BOOK MARK 4/2/2025
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${outdir}/computeMatrix_job.txt" \
        "computeMatrix scale-regions -S ${bw1} ${bw2} \
            -R ${merged_peaks} \
            --referencePoint center \
            --upstream 3000 \
            --downstream 3000 \
            -o ${outdir}/matrix.gz \
            --skipZeros"
    bsub -P acc_oscarlr -q premium -n 2 -W 24:00 -R "rusage[mem=8000]" -o "${outdir}/plotHeatmap_job.txt" \
        plotHeatmap -m ${outdir}/matrix.gz \
            --heatmapWidth 8 \
            --heatmapHeight 8 \
            --colorMap RdYlBu \
            --outFileName ${outdir}/heatmap.png \
            --showAdvancedOptions \
            --refPointLabel "Center of Region" \
            --kmeans 2
}


#fastqc
#trimming
#mapping
#chipqc
#corplot
#ipplot
#norm_stats
#bamcov
#bamcomp
#peakcalling
#prep1
#prep2
#prep3
heatmap