
if config['options']['paired']:

    rule fastqc_raw:
        input:
            path.join(raw_dir, '{sample}{lane}_{pair}.fastq.gz')
        output:
            html = path.join(qc_dir, 'fastqc_raw', '{sample}{lane}_{pair}_fastqc.html'),
            zip = path.join(qc_dir, 'fastqc_raw', '{sample}{lane}_{pair}_fastqc.zip')
        params:
            output_dir = path.join(qc_dir, 'fastqc_raw')
        shell:
            'mkdir -p {params.output_dir} && '
            'fastqc --quiet --outdir {params.output_dir} {input}'

    rule fastqc_trim:
        input:
            path.join(raw_dir, '{sample}{lane}_{pair}_trimmed.fastq.gz')
        output:
            html = path.join(qc_dir, 'fastqc_trim', '{sample}{lane}_{pair}_trimmed_fastqc.html'),
            zip = path.join(qc_dir, 'fastqc_trim', '{sample}{lane}_{pair}_trimmed_fastqc.zip')
        params:
            output_dir = path.join(qc_dir, 'fastqc_trim')
        shell:
            'mkdir -p {params.output_dir} && '
            'fastqc --quiet --outdir {params.output_dir} {input}'


else:

    rule fastqc_raw:
        input:
            path.join(raw_dir, '{sample}{lane}.fastq.gz')
        output:
            html = path.join(qc_dir, 'fastqc_raw', '{sample}{lane}_fastqc.html'),
            zip = path.join(qc_dir, 'fastqc_raw', '{sample}{lane}_fastqc.zip')
        params:
            output_dir=path.join(qc_dir, 'fastqc_raw')
        shell:
            'mkdir -p {params.output_dir} && '
            'fastqc --quiet --outdir {params.output_dir} {input}'


    rule fastqc_trim:
        input:
            path.join(raw_dir, '{sample}_{lane}_trimmed.fastq.gz')
        output:
            html = path.join(qc_dir, 'fastqc_trim', '{sample}{lane}_trimmed_fastqc.html'),
            zip = path.join(qc_dir, 'fastqc_trim', '{sample}{lane}_trimmed_fastqc.zip')
        params:
            output_dir=path.join(qc_dir, 'fastqc_trim')
        shell:
            'mkdir -p {params.output_dir} && '
            'fastqc --quiet --outdir {params.output_dir} {input}'


def multiqc_prerequisite(wildcards):

# bam indexing
    all_inputs = expand(path.join(bam_dir, '{sample}', 'Aligned.sort.bam.bai'),
                            sample = all_samples)
# genome unloading
    if config['options']['analysis_type'] != 'chimeric':
        all_inputs = all_inputs + [path.join(bam_dir, 'star_unload_2nd.done')]

# fastqc for raw file
    fastqc_temp_raw = expand(path.join(qc_dir, 'fastqc_raw', '{sample}{lane}{{pair}}_fastqc.zip'), 
        sample = all_samples, lane = all_lanes)
    pairs = ["_R1", "_R2"] if config['options']['paired'] else [""]
    fastqc_raw = expand(fastqc_temp_raw, pair=pairs)

# fastqc for trimmed file 
    fastqc_temp_trim = expand(path.join(qc_dir, 'fastqc_trim', '{sample}{lane}{{pair}}_trimmed_fastqc.zip'), 
        sample = all_samples, lane = all_lanes)
    pairs = ["_R1", "_R2"] if config['options']['paired'] else [""]
    fastqc_trim = expand(fastqc_temp_trim, pair=pairs)

# merge list for prerequisite
    #all_inputs = all_inputs + fastqc_raw + fastqc_trim
    all_inputs = all_inputs + fastqc_raw
    return all_inputs


rule multiqc:
    input:
        multiqc_prerequisite
    output:
        path.join(qc_dir, 'multiqc_report.html')
    params:
        qc_dir = qc_dir,
        out_dir = qc_dir
    shell:
        'multiqc --outdir {params.out_dir} {params.qc_dir}'


