import pandas as pd
from collections import defaultdict
import re
import glob

sample_file = config["sample_file"]
sample_df=pd.read_csv(sample_file, header=0, sep = "\t")
library=sample_df['SampleName'].tolist()
chip_lib=[x for x in library if not x.startswith(('Input', 'ATAC'))]

ASSAYS = set()
library_pe = set()
library_se = set()
LIBRARIES = defaultdict(lambda : defaultdict(set))
for lib in library:
    assay, tissue, replicate = lib.split("_")
    LIBRARIES[assay][tissue].add(replicate)
    ASSAYS.add(assay)
    num=len(glob.glob("Raw_Reads/"+ lib + "*.fq.gz"))
    if num == 1:
        library_se.add(lib)
    elif num == 2:
        library_pe.add(lib)
    else:
        raise Exception('Inputs are not either single or paired end, even or file names are not right,  please check input files!!!\n\n')
if "Input" in ASSAYS:
    ASSAYS.remove('Input')
PEAK_ASSAYS=list(ASSAYS)
library_pe = list(library_pe)
library_se = list(library_se)
bams=expand('Aligned_Reads/{library}.bam', library=library)
tagfile=expand('Aligned_Reads/{library}.tagAlign.gz', library=library)
peak_libs=filter(lambda x: x.startswith(tuple(PEAK_ASSAYS)), library)
spp=expand('Metrics/{library}.spp_stats.txt', library=peak_libs)
peakfiles=expand('Peak_Calls/{library}_Peaks.bed', library=peak_libs)
frip=expand('Metrics/{library}_FRiP.txt', library=peak_libs)
jsd=expand('Metrics/{library}_DeepTools_Metrics.txt', library=peak_libs)
AlgnState=expand('Metrics/{library}_Alignment_Stats.json',library=library)
AlignTXT = expand('Tables/{assay}_Alignment_Summary.txt',assay=PEAK_ASSAYS)
AlignCSV = expand('Tables/{assay}_Alignment_Summary.csv', assay=PEAK_ASSAYS)
mtxTXT =expand('Tables/{assay}_Quality_Metrics.txt', assay=[x for x in PEAK_ASSAYS if not x.startswith("ATAC")])
mtxCSV = expand('Tables/{assay}_Quality_Metrics.csv', assay=[x for x in PEAK_ASSAYS if not x.startswith("ATAC")])
PeakSum = expand('Tables/{assay}_Peak_Summary.txt',assay=PEAK_ASSAYS)
Peakcsv = expand('Tables/{assay}_Peak_Summary.csv', assay=PEAK_ASSAYS)
bw=expand('Track_Hub/{library}_Coverage.bw', library=library)
QCall="Tables/ChIP_assays_summary.txt"
trackdb = 'Track_Hub/trackDb.txt'

rule all:
    input: bams, spp, tagfile, peakfiles, frip, jsd, AlgnState, AlignTXT, AlignCSV, mtxTXT, mtxCSV, PeakSum, Peakcsv, bw, QCall, trackdb

rule bwa_index:
    input: 
        config["genome"]
    output: 
        config["genome"] + ".amb",
        config["genome"] + ".ann",
        config["genome"] + ".bwt",
        config["genome"] + ".pac",
        config["genome"] + ".sa"
    params:
        prefix = config["genome"],
        algorithm = "bwtsw"
    conda:
        'Envs/bwa.yaml'
    threads: 24
    shell:
        'bwa index -a {params.algorithm} {input}'

ruleorder: trim_paired_reads > trim_reads
rule trim_reads:
    input:
        'Raw_Reads/{library}.fq.gz'
    output:
        reads = 'Trimmed_Reads/{library}_trimmed.fq.gz',
        report = 'Trimmed_Reads/{library}.fq.gz_trimming_report.txt'
    conda:
        'Envs/trimgalore.yaml'
    threads: 12
    shell: 
        'trim_galore -q {config[mapq]} --cores {threads} {input} -o Trimmed_Reads'

rule trim_paired_reads:
    input:
        'Raw_Reads/{library}_R1.fq.gz',
        'Raw_Reads/{library}_R2.fq.gz'
    output: 
        reads1 = 'Trimmed_Reads/{library}_R1_val_1.fq.gz',
        reads2 = 'Trimmed_Reads/{library}_R2_val_2.fq.gz',
        report = 'Trimmed_Reads/{library}_R1.fq.gz_trimming_report.txt',
        report2 = 'Trimmed_Reads/{library}_R2.fq.gz_trimming_report.txt'
    conda:
        'Envs/trimgalore.yaml'
    threads: 12
    shell: 
        'trim_galore -q {config[mapq]} --cores {threads} --paired {input} -o Trimmed_Reads'

rule bwa_mem:
    input: 
        index = rules.bwa_index.output,
        reads = rules.trim_reads.output.reads
    output: 
        bam = config['tempdir'] + '/{library}.aligned.bam'
    threads: 12
    conda:
        'Envs/bwa.yaml'
    shell: 
        "bwa mem -M -t {threads} -R '@RG\\tID:{wildcards.library}\\tSM:{wildcards.library}' {config[genome]} {input.reads} | samtools view -bS - > {output}"

rule bwa_mem_paired:
    input: 
        rules.bwa_index.output,
        reads1 = rules.trim_paired_reads.output.reads1,
        reads2 = rules.trim_paired_reads.output.reads2 
    output: 
        bam = config['tempdir'] + '/{library}.aligned.bam'
    threads: 12
    conda:
        'Envs/bwa.yaml'
    shell: 
        "bwa mem -M -t {threads} -R '@RG\\tID:{wildcards.library}\\tSM:{wildcards.library}' {config[genome]} {input.reads1} {input.reads2} | samtools view -bS - > {output}"

rule filter_alignments:
    input: 
        config['tempdir'] + '/{library}.aligned.bam'
    output: 
        bam = config['tempdir'] + '/{library}.filtered.bam'
    conda:
        'Envs/samtools.yaml'
    shell: 
        'samtools view -h -q {config[mapq]} {input} | samtools view -S -b - > {output}'

rule namesort_bam:
    input:
        rules.filter_alignments.output
    output:
        temp(config['tempdir'] + '/{library}.namesorted.bam')
    threads: 12
    conda:
        'Envs/samtools.yaml'
    shell:
        'samtools sort -@ {threads} -n -o {output} {input}'

rule fixmate_bam:
    input:
        rules.namesort_bam.output
    output:
        temp(config['tempdir'] + '/{library}.fixmate.bam')
    threads: 12
    conda:
        'Envs/samtools.yaml'
    shell:
        'samtools fixmate -m {input} {output}'

rule sort_bam:
    input: 
        rules.fixmate_bam.output
    output: 
        temp(config['tempdir'] + '/{library}.sorted.bam')
    threads: 12
    conda:
        'Envs/samtools.yaml'
    shell:
        'samtools sort -@ {threads} -o {output} {input}'

rule remove_duplicates:
    input:
        rules.sort_bam.output
    output:
        bam = 'Aligned_Reads/{library}.bam'
    threads: 12
    conda:
        'Envs/samtools.yaml'
    shell:
        'samtools markdup -r {input} {output.bam}'

rule idx_bam:
    input:
        rules.remove_duplicates.output
    output:
        bai='Aligned_Reads/{library}.bam.bai'
    conda:
        'Envs/samtools.yaml'
    shell:
        'samtools index {input}'

rule bam_to_tagalign:
    input:
        rules.remove_duplicates.output# 'Aligned_Reads/{library}.bam'
    output:
        'Aligned_Reads/{library}.tagAlign.gz'
    conda:
        'Envs/bedtools.yaml'
    threads: 1
    resources:
         mem_mb=2000
    shell:
        """bedtools bamtobed -i {input} | awk 'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}' | gzip -nc > {output}"""

rule spp_stats:
    input:
        rules.bam_to_tagalign.output #'Aligned_Reads/{library}.tagAlign.gz'
    output:
        stats = 'Metrics/{library}.spp_stats.txt',
        figure = 'Metrics/{library}_Cross_Correlation.pdf'
    threads: 4
    #wildcard_constraints:
    #    library=chip_lib
    conda:
        'Envs/r.yaml'
    shell:
        'Rscript {config[workdir]}/Scripts/run_spp.R -c={input} -rf -out={output.stats} -p={threads} -s=0:2:400 -savp={output.figure} -tmpdir={config[tempdir]}'

def get_assay_type(wildcards):
    if wildcards.library in library_pe:
        return 'Paired'
    elif wildcards.library in library_se:
        return 'Single'

def call_peak_inputs(wildcards):
    inputs = {}
    assay, tissue, rep = wildcards.library.split('_')
    inputs['bam']=["Aligned_Reads/{}.bam".format(wildcards.library)]
    inputs['bai']=["Aligned_Reads/{}.bam.bai".format(wildcards.library)]
    inputs['chip']=["Aligned_Reads/{}.tagAlign.gz".format(wildcards.library)]
    if assay != "ATAC":
        inputs['control']=["Aligned_Reads/Input_" + tissue +"_"+ rep + ".bam"]
        inputs['cbai']=["Aligned_Reads/Input_" + tissue +"_"+ rep + ".bam.bai"]
        inputs['ctl']=["Aligned_Reads/Input_" + tissue +"_"+ rep + ".tagAlign.gz"]
    return inputs

rule call_peaks:
    input: 
        unpack(call_peak_inputs)
    output:    
        peaks = 'Peak_Calls/{library}_Peaks.bed'
    params:
        assaytype=get_assay_type
    conda:
        'Envs/macs.yaml'
    #wildcard_constraints:
    #    library=chip_lib
    threads: 4
    script: 
        'Scripts/CallPeaks.py'

rule individual_frip:
    input:
        bam='Aligned_Reads/{library}.bam',
        bai='Aligned_Reads/{library}.bam.bai',
        peaks='Peak_Calls/{library}_Peaks.bed'
    output:
        metrics = 'Metrics/{library}_FRiP.txt',
        figure = temp('Metrics/{library}_FRiP.png')
    conda:
        'Envs/deeptools.yaml'
    #wildcard_constraints:
    #    library=chip_lib
    threads: 4
    script:
        'Scripts/Calculate_FRiP.py'

def deeptools_jsd_inputs(wildcards):
    assay, tissue, rep = wildcards.library.split('_')
    inputs = {'bam': ['Aligned_Reads/{}.bam'.format(wildcards.library)]}
    inputs['bai']="Aligned_Reads/{}.bam.bai".format(wildcards.library)
    inputs['control']="Aligned_Reads/Input_" + tissue +"_"+ rep + ".bam"
    inputs['cbai']="Aligned_Reads/Input_" + tissue +"_"+ rep + ".bam.bai"
    return inputs

rule deeptools_jsd:
    input:
        unpack(deeptools_jsd_inputs)
    output:
        metrics = 'Metrics/{library}_DeepTools_Metrics.txt',
        png = temp('Metrics/{library}_Temp.png')
    conda:
        'Envs/deeptools.yaml'
    #wildcard_constraints:
    #    library=chip_lib
    threads: 4
    shell:
        'plotFingerprint -b {input.bam} {input.control} -p={threads} -plot {output.png} --extendReads=200 --skipZeros --outQualityMetrics={output.metrics} --JSDsample={input.control}'

def alignment_stats_inputs(wildcards):
    inputs = {}
    if wildcards.library in library_pe:
        inputs['trimmed_fq'] = 'Trimmed_Reads/{library}_R1_val_1.fq.gz'.format(library=wildcards.library)
        inputs['trim_report'] = 'Trimmed_Reads/{library}_R1.fq.gz_trimming_report.txt'.format(library=wildcards.library)
    else:
        inputs['trimmed_fq'] = 'Trimmed_Reads/{library}_trimmed.fq.gz'.format(library=wildcards.library)
        inputs['trim_report'] ='Trimmed_Reads/{library}.fq.gz_trimming_report.txt'.format(library=wildcards.library)
    inputs['aligned_bam'] = config['tempdir'] + '/{library}.aligned.bam'.format(library=wildcards.library)
    inputs['filtered_bam'] = config['tempdir'] + '/{library}.filtered.bam'.format(library=wildcards.library)
    inputs['deduped_bam'] = 'Aligned_Reads/{library}.bam'.format(library=wildcards.library)
    if ( not wildcards.library.startswith(('ATAC', 'Input'))):
        inputs['spp_stats'] = 'Metrics/{library}.spp_stats.txt'.format(library=wildcards.library)
    return inputs

rule get_alignment_stats:
    input:
        unpack(alignment_stats_inputs)
    output:
        json = 'Metrics/{library}_Alignment_Stats.json'
    params:
        assay_type = get_assay_type
    script:
        'Scripts/Get_Alignment_Stats.py'


def libraries(assay=None, tissue=None, rep=None, sample=None, has_input=None, skip_input=True, omit_training=False, peak_assay_only=False):
    """ Returns a list of libraries matching the given parameters. """
    if sample:
        tissue, rep = sample.split("_")
    for _assay in sorted(LIBRARIES):
        if skip_input and _assay in config['inputs']:
            continue
        if has_input and _assay in config['no_input']:
            continue
        if has_input == False and _assay not in config['no_input']:
            continue
        if peak_assay_only and _assay not in PEAK_ASSAYS:
            continue
        for _tissue in sorted(LIBRARIES[_assay]):
            for _rep in sorted(LIBRARIES[_assay][_tissue]):
                if (not assay or assay == _assay) and (not tissue or tissue == _tissue) and (not rep or rep == _rep):
                    lib = '{assay}_{tissue}_{rep}'.format(assay=_assay, tissue=_tissue, rep=_rep)
                    if omit_training and 'ChromHMM_training_omit' in config and lib in config['ChromHMM_training_omit']:
                        continue
                    yield lib

rule make_alignment_report:
    input:
        stats = lambda wildcards: expand('Metrics/{library}_Alignment_Stats.json', library=libraries(wildcards.PEAK_ASSAYS, skip_input=False)), 
        bams = lambda wildcards: expand('Aligned_Reads/{library}.bam', library=libraries(wildcards.PEAK_ASSAYS, skip_input=False))
    output:
        txt = 'Tables/{PEAK_ASSAYS}_Alignment_Summary.txt',
        csv = 'Tables/{PEAK_ASSAYS}_Alignment_Summary.csv'
    params:
        libraries = lambda wildcards: expand('{library}', library=libraries(assay=wildcards.PEAK_ASSAYS, skip_input=False))
    script:
        'Scripts/Make_Alignment_Summary.py'

rule make_quality_metrics_report:
    input:
        stats = lambda wildcards: expand('Metrics/{library}_Alignment_Stats.json', library=libraries(assay=wildcards.PEAK_ASSAYS, skip_input=False)),
        deeptools = lambda wildcards: expand('Metrics/{library}_DeepTools_Metrics.txt', library=libraries(assay=wildcards.PEAK_ASSAYS, skip_input=False))
    output:
        txt = 'Tables/{PEAK_ASSAYS}_Quality_Metrics.txt',
        csv = 'Tables/{PEAK_ASSAYS}_Quality_Metrics.csv'
    params:
        libraries = lambda wildcards: expand('{library}', library=libraries(assay=wildcards.PEAK_ASSAYS, skip_input=False)),
        assay = "{PEAK_ASSAYS}"
    script: 'Scripts/Make_Quality_Metrics.py'

def replicates_for_assay(assay):
    """ Given an assay, returns a list of all replicates with libraries of that assay """
    reps = set()
    for tissue in LIBRARIES[assay]:
        reps.update(LIBRARIES[assay][tissue])
    return sorted(list(reps))

rule make_peak_summary_report:
    input: 
        lambda wildcards: expand('Peak_Calls/{library}_Peaks.bed', library=libraries(assay=wildcards.PEAK_ASSAYS, skip_input=True)),
        lambda wildcards: expand('Metrics/{library}_FRiP.txt', library=libraries(assay=wildcards.PEAK_ASSAYS, skip_input=True))
    output: 
        txt = 'Tables/{PEAK_ASSAYS}_Peak_Summary.txt', 
        csv = 'Tables/{PEAK_ASSAYS}_Peak_Summary.csv'
    params: 
        tissues = lambda wildcards: list(LIBRARIES[wildcards.PEAK_ASSAYS].keys()), 
        reps = lambda wildcards: replicates_for_assay(wildcards.PEAK_ASSAYS),
        library=lambda wildcards: expand("{library}", library=libraries(assay=wildcards.PEAK_ASSAYS, skip_input=True))
    script: 
        'Scripts/Peak_Call_Summary.py'

rule bam_coverage:
    input: 
        bam = 'Aligned_Reads/{library}.bam', 
        bai = 'Aligned_Reads/{library}.bam.bai'
    output: 
        'DeepTools/{library}.bw'
    conda: 'Envs/deeptools.yaml'
    threads: 4
    shell: 
        'bamCoverage -p={threads} -b {input.bam} -o {output} -of=bigwig --normalizeUsing RPKM --effectiveGenomeSize {config[genomesize]} --ignoreDuplicates --extendReads=200 -bs 20'

rule BamCoverage_BigWig:
    input:
        bdg = 'DeepTools/{library}.bw',
    output:
        'Track_Hub/{library}_Coverage.bw'
    shell:
        'ln -sr {input.bdg} {output}'

rule Create_TrackDB:
    input:
        fes = expand('Track_Hub/{library}_Coverage.bw', library=library),
    output:
        trackdb = 'Track_Hub/trackDb.txt'
    script:
        'Scripts/Make_TrackDB.py'

rule Aggregate_QC:
    input:
        lambda wildcards: expand('Tables/{assays}_Alignment_Summary.csv', assays=PEAK_ASSAYS),
        lambda wildcards: expand('Tables/{assays}_Peak_Summary.csv', assays=PEAK_ASSAYS),
        lambda wildcards: expand('Tables/{assays}_Peak_Summary.csv', assays=[x for x in PEAK_ASSAYS if not x.startswith('ATAC')])
    output:
        "Tables/ChIP_assays_summary.txt"
    params:
       assays = ' '.join(PEAK_ASSAYS)
    conda:
        'Envs/r.yaml'
    shell:
        'Rscript {config[workdir]}/Scripts/QC_Summary.R {params}'
