"""
Author: Y. Ahmed-Braimah
--- RNA-seq snakemake workflow for running on CBSU cluster
"""

import json
import os
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output

##--------------------------------------------------------------------------------------##
## Functions
##--------------------------------------------------------------------------------------##

# To print process messages
def message(x):
  print()

# To remove suffix from a string
def rstrip(text, suffix):
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

## define environment variables

##--------------------------------------------------------------------------------------##
## Global config files: 
##--------------------------------------------------------------------------------------##

configfile: 'config.yml'

# Full path to an uncompressed FASTA file with all chromosome sequences.
DNA = config['DNA']

# Full path to an uncompressed GTF file with known gene annotations.
GTF = config['GTF']

# Full path to a folder where output files will be created.
OUT_DIR = config['OUT_DIR']

## set the working directory for each job
# USER = os.environ.get('USER')
# JOB_ID = os.environ.get('JOB_ID')
# WORK_DIR = "/workdir/"

# message("The current working directory is " + WORK_DIR)

# Samples and their corresponding filenames.
# paired-end:
peFILES = json.load(open(config['PE_SAMPLES_JSON'])) # The "samples.json" file can be created with
peSAMPLES = sorted(peFILES.keys())             # the "make_json_SE/PE_samples.py" script
# single-end
seFILES = json.load(open(config['SE_SAMPLES_JSON'])) # The "samples.json" file can be created with
seSAMPLES = sorted(seFILES.keys())                  # the "make_json_SE/PE_samples.py" script
# read both
#FILES = json.load(open(config['SAMPLES_JSON']))
combinedSam = [peSAMPLES, seSAMPLES]
SAMPLES = [y for x in combinedSam for y in x]  

##--------------------------------------------------------------------------------------##
## RULES
##--------------------------------------------------------------------------------------##

## Final expected output(s)
rule all: 
    input: 
      expand(join(OUT_DIR, 'Tophat', '{sample}',  'accepted_hits.bam'), sample = SAMPLES),
      expand(join(OUT_DIR, 'Cufflinks', '{sample}', 'transcripts.gtf'), sample = SAMPLES),
      expand(join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html'), sample = SAMPLES),
      expand(join(OUT_DIR, 'fastQC', '{sample}' + '.R2_fastqc.html'), sample = peSAMPLES),
      join(OUT_DIR, 'MultiQC', 'multiqc_report.html'),
      join(OUT_DIR, 'Cuffmerge', 'merged.gtf'),
      # expand(join(OUT_DIR, 'Cuffquant', '{sample}',  'abundances.cxb'), sample = SAMPLES),
      # join(OUT_DIR, 'Cuffnorm', 'expression_data', 'run.info'),
      expand(join(OUT_DIR, 'eXpress', '{sample}', 'results.xprs'), sample = SAMPLES),
      join(OUT_DIR, 'eXpress', 'genes.TMM.EXPR.matrix')
        



## Rule to generate bowtie2 genome index file
rule index:
    input:
        dna = DNA
    output:
        index = join(dirname(DNA), rstrip(DNA, '.fa') + '.rev.1.bt2'),
        bt2i = join(dirname(DNA), rstrip(DNA, '.fa') + '.ok')
    log:
        join(dirname(DNA), 'bt2.index.log')
    benchmark:
        join(dirname(DNA), 'bt2.index.benchmark.tsv')
    message: 
        """--- Building bowtie2 genome index """
    run:
        # shell('mkdir -p ')
        shell('samtools faidx {input.dna}')
        shell('bowtie2-build {input.dna} ' + join(dirname(DNA), rstrip(DNA, '.fa')) + ' > {log} 2>&1')
        shell('touch ' + join(dirname(DNA), rstrip(DNA, '.fa') + '.ok'))

## Rule to check quality of trimmed reads
rule fastqc_pe:
    input:
        r1 = lambda wildcards: peFILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: peFILES[wildcards.sample]['R2']
    output:
        r1 = join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html'),
        r2 = join(OUT_DIR, 'fastQC', '{sample}' + '.R2_fastqc.html')
    log:
        join(OUT_DIR, 'fastQC', '{sample}', 'fastQC_pe.log')
    benchmark:
        join(OUT_DIR, 'fastQC', '{sample}', 'fastQC_pe.benchmark.tsv')
    message: 
        """--- Checking trimmed read quality of sample "{wildcards.sample}" with FastQC """
    run:
        if not os.path.exists(join(OUT_DIR, 'fastQC')):
            os.makedirs(join(OUT_DIR, 'fastQC'))
        shell('fastqc'
                ' -o ' + join(OUT_DIR, 'fastQC') + 
                ' {input.r1} {input.r2} > {log} 2>&1')

## Rule to check quality of trimmed reads
rule fastqc_se:
    input:
        r1 = lambda wildcards: seFILES[wildcards.sample]['R1']
    output:
        r1 = join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html')
    log:
        join(OUT_DIR, 'fastQC', '{sample}', 'fastQC_pe.log')
    benchmark:
        join(OUT_DIR, 'fastQC', '{sample}', 'fastQC_pe.benchmark.tsv')
    message: 
        """--- Checking trimmed read quality of sample "{wildcards.sample}" with FastQC """
    run:
        if not os.path.exists(join(OUT_DIR, 'fastQC')):
            os.makedirs(join(OUT_DIR, 'fastQC'))
        shell('fastqc'
                ' -o ' + join(OUT_DIR, 'fastQC') + 
                ' {input.r1} > {log} 2>&1')

## Rule for mapping reads to the genome with Tophat
rule tophat_pe:
    input:
        r1 = lambda wildcards: peFILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: peFILES[wildcards.sample]['R2'],
        idx = rules.index.output.bt2i
    output: 
        bam = join(OUT_DIR, 'Tophat', '{sample}', 'accepted_hits.bam')
    params: 
        gtf = GTF
    log:
        join(OUT_DIR, 'Tophat', '{sample}', 'tophat.map.log')
    benchmark:
        join(OUT_DIR, 'Tophat', '{sample}', 'tophat.map.benchmark.tsv')
    message: 
        """--- Mapping sample "{wildcards.sample}" with Tophat."""
    run: 
          shell('tophat'                                     
                ' -o ' + join(OUT_DIR, 'Tophat', '{wildcards.sample}') + '/'    
                ' -G {params.gtf}'                                 
                ' -p 8 ' + join(dirname(DNA), rstrip(DNA, '.fa')) +                           
                ' {input.r1} {input.r2} > {log} 2>&1')

# Rule for mapping reads to the genome with Tophat
rule tophat_se:
    input:
        r1 = lambda wildcards: seFILES[wildcards.sample]['R1'],
        idx = rules.index.output.bt2i
    output: 
        bam = join(OUT_DIR, 'Tophat', '{sample}', 'accepted_hits.bam')
    params: 
        gtf = GTF
    log:
        join(OUT_DIR, 'Tophat', '{sample}',  'tophat.map.log')
    benchmark:
        join(OUT_DIR, 'Tophat', '{sample}', 'tophat.map.benchmark.tsv')
    message: 
        """--- Mapping sample "{wildcards.sample}" with Tophat."""
    run: 
          shell('tophat'                                     
                ' -o ' + join(OUT_DIR, 'Tophat', '{wildcards.sample}') + '/'    
                ' -G {params.gtf}'                                 
                ' -p 8 ' + join(dirname(DNA), rstrip(DNA, '.fa')) +                           
                ' {input.r1} > {log} 2>&1')

## Rule to collate fastQC outputs with multiQC
rule multiQC:
    input:
        expand(join(OUT_DIR, 'Tophat', '{sample}', 'accepted_hits.bam'), sample = SAMPLES),
        expand(join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html'), sample = SAMPLES),
        expand(join(OUT_DIR, 'fastQC', '{sample}' + '.R2_fastqc.html'), sample = peSAMPLES),
        expand(join(OUT_DIR, 'eXpress', '{sample}', 'results.xprs'), sample = SAMPLES)

    output:
        file = join(OUT_DIR, 'MultiQC', 'multiqc_report.html')
    log:
        join(OUT_DIR, 'MultiQC', 'multiQC.log')
    benchmark:
        join(OUT_DIR, 'MultiQC', 'multiQC.benchmark.tsv')
    message: 
        """--- Running MultiQC """
    run:
        shell('ls -1 ' + join(OUT_DIR) + '/Tophat/*/align_summary.txt > ' + join(OUT_DIR, 'summary_files.txt'))
        shell('ls -1 ' + join(OUT_DIR) + '/fastQC/*fastqc.zip >> ' + join(OUT_DIR, 'summary_files.txt'))
        shell('ls -1 ' + join(OUT_DIR) + '/eXpress/*/eXpress.log >> ' + join(OUT_DIR, 'summary_files.txt'))
        shell('multiqc'
                ' -f'
                ' -o ' + join(OUT_DIR, 'MultiQC') + ' -d -dd 2 -l ' + join(OUT_DIR, 'summary_files.txt') +
                ' > {log} 2>&1')

## Rule for assembling rtansfrags with Cufflinks
rule cufflinks:
    input: 
        expand(join(OUT_DIR, 'Tophat', '{sample}', 'accepted_hits.bam'), sample = SAMPLES)
    output: 
        gtf = join(OUT_DIR, 'Cufflinks', '{sample}', 'transcripts.gtf')
    params:  
        gtf=GTF
    log:
        join(OUT_DIR, 'Cufflinks', '{sample}', 'cufflinks.log')
    benchmark:
        join(OUT_DIR, 'Cufflinks', '{sample}', 'cufflinks.benchmark.tsv')
    message: 
        """--- Assembling "{wildcards.sample}" transcripts with cufflinks."""
    run: 
        shell('cufflinks'
              ' -g {params.gtf}'
              ' -p 8'
              ' -o ' +  join(OUT_DIR, 'Cufflinks', '{wildcards.sample}') + '/'
              ' {input}'
              ' &> {log}')


## Rule for merging Cufflinks assembled transcripts
rule cuffmerge:
    input:
        asmblys = expand(join(OUT_DIR, 'Cufflinks', '{sample}', 'transcripts.gtf'), sample = SAMPLES)
    output: 
        merged = join(OUT_DIR, 'Cuffmerge', 'merged.gtf')
    params: 
        gtf =GTF, 
        fa = DNA
    log:
        join(OUT_DIR, 'Cuffmerge', 'cuffmerge.log')
    benchmark:
        join(OUT_DIR, 'Cuffmerge', 'cuffmerge.benchmark.tsv')
    message: 
        "--- Comparing transcripts to the reference and outputting merged gtf file."
    run: 
        # generate the assemblies text file
        shell('ls -1 ' + join(OUT_DIR) + '/Cufflinks/*/transcripts.gtf > ' + join(OUT_DIR, 'Cuffmerge', 'assemblies.txt'))
        ## an alternative way to create the "assemblies.txt" file:
        ## with open (output.txt, 'w') as out:
        ##  print(*input, sep="\n", file=out)

        # run cuffmerge
        shell('cuffmerge'
              ' -o ' + join(OUT_DIR, 'Cuffmerge') +
              ' -g {params.gtf}'
              ' --keep-tmp'
              ' -s {params.fa}' 
              ' -p 2 ' + join(OUT_DIR, 'Cuffmerge', 'assemblies.txt') + 
              ' &> {log}')

# ## Rule for quantifying abundance with Cuffquant
# rule cuffquant:
#     input:
#         bam = expand(join(OUT_DIR, 'Tophat', '{sample}', 'accepted_hits.bam'), sample = SAMPLES),
#         gtf = rules.cuffmerge.output.merged,
#         dna = DNA
#     output: 
#         join(OUT_DIR, 'Cuffquant', '{sample}', 'abundances.cxb')
#     log:
#         join(OUT_DIR, 'Cuffquant', '{sample}', 'cffqnt.map.log')
#     benchmark:
#         join(OUT_DIR, 'Cuffquant', '{sample}', 'cffqnt.benchmark.tsv')
#     message: 
#         """--- Quantifying abundances with Cuffquant for sample "{wildcards.sample}"."""
#     run: 
#         # run cuffmerge
#         shell('cuffquant'
#               ' -o ' + join(OUT_DIR, 'Cuffquant', '{wildcards.sample}') +
#               # ' -p 8'
#               ' -b {input.dna}'
#               ' -u' 
#               ' {input.gtf}'
#               ' {input.bam}'
#               ' &> {log}')

# ## Rule for quantifying abundance with Cuffquant
# rule cuffnorm:
#     input:
#         cxb = expand(join(OUT_DIR, 'Cuffquant', '{sample}', 'abundances.cxb'), sample=SAMPLES),
#         gtf = rules.cuffmerge.output.merged
#     output: 
#         join(OUT_DIR, 'Cuffnorm', 'expression_data', 'run.info')
#     log:
#         join(OUT_DIR, 'Cuffnorm', 'cffnrm.map.log')
#     benchmark:
#         join(OUT_DIR, 'Cuffnorm', 'cffnrm.benchmark.tsv')
#     message: 
#         """--- Merge Cuffquant abundances with Cuffnorm."""
#     run:
#         # create sample sheet
#         shell('echo -e "sample_name\tgroup" > ' + join(OUT_DIR, 'Cuffnorm', 'sample_sheet.txt'))
#         shell('ls -1 ' + join(OUT_DIR, 'Cuffquant', '*', 'abundances.cxb') + ' > ' + join(OUT_DIR, 'Cuffnorm', 'cq_abndces.txt'))
#         shell('cat ' + join(OUT_DIR, 'Cuffnorm', 'cq_abndces.txt') + ' | sed "s/\/Cuffquant.*//g" | sed "s/.*\///g" | paste -d"\t" ' + join(OUT_DIR, 'Cuffnorm', 'cq_abndces.txt') +' - >> ' + join(OUT_DIR, 'Cuffnorm', 'sample_sheet.txt'))
#         # run cuffnorm
#         shell('cuffnorm'
#               ' --use-sample-sheet'
#               ' -o ' + join(OUT_DIR, 'Cuffnorm', 'expression_data') +
#               # ' -p 8'
#               ' {input.gtf} ' + join(OUT_DIR, 'Cuffnorm', 'sample_sheet.txt') +
#               ' &> {log}')

# Rule for making a cDNA file from the new GTF file and the DNA, then prep index for bowtie2.
rule make_cdna:
    input:
        dna = DNA,
        gtf = rules.cuffmerge.output.merged
    output:
        cdna = join(OUT_DIR, 'transcriptome', 'gffread_transcripts.fa'),
        geneTrans = join(OUT_DIR, 'transcriptome', 'gffread_transcripts.gene_trans_map'),
        bt2_trans_indx = join(OUT_DIR, 'transcriptome', 'gffread_transcripts.fa.bowtie2.ok')
    log:
        gffread = join(OUT_DIR, 'transcriptome', 'logs', 'gffread.log'),
        trans_bt2 = join(OUT_DIR, 'transcriptome', 'logs', 'trans_bt2.log')
    benchmark:
        join(OUT_DIR, 'transcriptome', 'logs', 'gffread.benchmark.tsv')
    message: 
        "--- Building bowtie2 transcriptome index for gffread transcripts."
    run:
        # Extract a sequence for each transcript in the GTF file.
        shell('gffread -F -w {output.cdna} -g {input.dna} {input.gtf} > {log.gffread}')
        # Extract the FASTA header from the cDNA file and make into
        # trans_map file.
        shell('grep ">" {output.cdna} | sed "s/>//g" | sed "s/gene=//g" | awk \'{{print $2"\t"$1}}\' | sort -u > {output.geneTrans}')
        # And finally make the index files.
        shell('align_and_estimate_abundance.pl' 
              ' --transcripts {output.cdna}'
              ' --gene_trans_map {output.geneTrans}'
              ' --est_method eXpress'
              ' --aln_method bowtie2'
              ' --prep_reference'
              ' --output_dir ' + join(OUT_DIR, 'transcriptome') +
              ' > {log.trans_bt2} 2>&1')

# Rule for mapping reads to the new transcriptome file with bowtie2 and quantifying abundance with eXpress
rule express_pe:
    input:
        read1 = lambda wildcards: peFILES[wildcards.sample]['R1'],
        read2 = lambda wildcards: peFILES[wildcards.sample]['R2'],
        cdna = rules.make_cdna.output.cdna,
        geneTrans = rules.make_cdna.output.geneTrans
    output:
        results = join(OUT_DIR, 'eXpress', '{sample}', 'results.xprs')
    log:
        join(OUT_DIR, 'eXpress', '{sample}', 'eXpress.log')
    benchmark:
        join(OUT_DIR, 'eXpress', '{sample}', 'eXpress.benchmark.tsv')
    message: 
        """--- Mapping "{wildcards.sample}" reads to transcriptome with bowtie2 and quantifying abundance with eXpress."""
    run:
        shell('align_and_estimate_abundance.pl' 
              ' --transcripts {input.cdna}'
              ' --seqType fq'
              ' --left {input.read1}'
              ' --right {input.read2}'
              ' --gene_trans_map {input.geneTrans}'
              ' --thread_count 8'  
              ' --est_method eXpress'
              ' --aln_method bowtie2'
              ' --output_dir ' + join(OUT_DIR, 'eXpress', '{wildcards.sample}') +  
              ' > {log} 2>&1')

# Rule for mapping reads to the new transcriptome file with bowtie2 and quantifying abundance with eXpress
rule express_se:
    input:
        read1 = lambda wildcards: seFILES[wildcards.sample]['R1'],
        cdna = rules.make_cdna.output.cdna,
        geneTrans = rules.make_cdna.output.geneTrans
    output:
        results = join(OUT_DIR, 'eXpress', '{sample}', 'results.xprs')
    log:
        join(OUT_DIR, 'eXpress', '{sample}', 'eXpress.log')
    benchmark:
        join(OUT_DIR, 'eXpress', '{sample}', 'eXpress.benchmark.tsv')
    message: 
        """--- Mapping "{wildcards.sample}" reads to transcriptome with bowtie2 and quantifying abundance with eXpress."""
    run:
        shell('align_and_estimate_abundance.pl' 
              ' --transcripts {input.cdna}'
              ' --seqType fq'
              ' --single {input.read1}'
              ' --gene_trans_map {input.geneTrans}'
              ' --thread_count 8'  
              ' --est_method eXpress'
              ' --aln_method bowtie2'
              ' --output_dir ' + join(OUT_DIR, 'eXpress', '{wildcards.sample}') +  
              ' > {log} 2>&1')

rule merge_abundance:
    input:
        quants = expand(join(OUT_DIR, 'eXpress', '{sample}', 'results.xprs'), sample = SAMPLES)
    output:
        abundances = join(OUT_DIR, 'eXpress', 'genes.TMM.EXPR.matrix'),
        samplesList = join(OUT_DIR, 'eXpress', 'genes.samples.list')
    log:
        join(OUT_DIR, 'eXpress', 'abnd_merge.log')
    benchmark:
        join(OUT_DIR, 'eXpress', 'abnd_merge.benchmark.tsv')        
    message: 
        "--- Merging eXpress outputs from all samples"
    run:
        shell('ls -1 ' + join(OUT_DIR, 'eXpress', '*', 'results.xprs.genes') + ' > ' + join(OUT_DIR, 'eXpress', 'genes.samples.list'))
        shell('ls -1 ' + join(OUT_DIR, 'eXpress', '*', 'results.xprs') + ' > ' + join(OUT_DIR, 'eXpress', 'isoforms.samples.list'))
        shell('cd ' + join(OUT_DIR, 'eXpress') +
                ' && abundance_estimates_to_matrix.pl'
                ' --est_method eXpress'
                ' --name_sample_by_basedir'
                ' --out_prefix genes'
                ' ' + join(OUT_DIR, 'eXpress', 'genes.samples.list') +
                ' > {log} 2>&1')
        shell('cd ' + join(OUT_DIR, 'eXpress') +
                ' && abundance_estimates_to_matrix.pl'
                ' --est_method eXpress'
                ' --name_sample_by_basedir'
                ' --out_prefix isoforms'
                ' ' + join(OUT_DIR, 'eXpress', 'isoforms.samples.list') +
                ' > {log} 2>&1')
