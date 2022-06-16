#####################################################
# preparing STAR alignment 
#####################################################

# input fastq files
def align_inputs(wildcards):
    fastq_path = path.join(raw_dir, "{sample}{{lane}}{{pair}}.fastq.gz".format(
        sample=wildcards.sample))
    pairs = ["_R1", "_R2"] if config['options']['paired'] else [""]
    lanes = ["_L00"+str(i) for i in range(1,config['options']['n_lane']+1)] if config['options']['n_lane']>1 else [""]
    list_input = expand(fastq_path, pair = pairs, lane = lanes)

    if config['options']['paired']:
        R1 = ','.join(list(filter(lambda x: '_R1' in x, list_input)))
        R2 = ','.join(list(filter(lambda x: '_R2' in x, list_input)))
        return [R1, R2]
    else:
        return list_input

#####################################################
# STAR alignment 
#####################################################

if config['options']['analysis_type'] == 'chimeric':
   
    # STAR 2-pass alignment for fusion analysis
    # In this mode, shared momeory mode cannot be used
    rule star_align_fusion:
#        input:
#            align_inputs
        output:
            bam = temp(path.join(bam_dir, '{sample}', 'Aligned.out.bam')),
            chim = path.join(bam_dir, '{sample}', 'Chimeric.out.junction')
        params:
            input = align_inputs,
            index = config['star_align']['index'],
            options = format_options(config['star_align']['options_fusion']),
            out_prefix = path.join(bam_dir, '{sample}/'),
            qc_prefix = path.join(qc_dir, 'star_align','{sample}'),
            qc_dir = path.join(qc_dir, 'star_align')
        conda:
            '../envs/star.yaml'
        threads:
            config['star_align']['threads']
        log:
            path.join(log_dir, 'star_align', '{sample}_align_fusion.log')
        shell:
            'STAR {params.options} --genomeDir {params.index} '
            '--outFileNamePrefix {params.out_prefix} --runThreadN {threads} '
            '--readFilesCommand gunzip -c --outSAMtype BAM Unsorted '
            '--outReadsUnmapped None --twopassMode Basic '
            '--outSAMstrandField intronMotif --outSAMunmapped Within '
            '--readFilesIn {params.input} 2> {log} '
            '&& mkdir -p {params.qc_dir} '
            '&& mv {params.out_prefix}Log.final.out {params.qc_prefix}.Log.final.out'

    if config['options']['star_fusion']:
        rule star_fusion:
            input:
                path.join(bam_dir, '{sample}', 'Chimeric.out.junction')
            output:
                path.join(fusion_dir, '{sample}', 'star-fusion.fusion_predictions.abridged.coding_effect.tsv')
            params:
                lib_dir = config['star_fusion']['reflib'],
                out_dir = path.join(fusion_dir, '{sample}'),
                options = format_options(config['star_fusion']['options'])
            conda:
                '../envs/star_fusion.yaml'
            threads:
                config['star_fusion']['threads']
            log:
                path.join(log_dir, 'star_fusion', '{sample}_star_fusion.log')
            shell:
                'STAR-Fusion {params.options} --genome_lib_dir {params.lib_dir} '
                '-J {input} --output_dir {params.out_dir} '
                '--CPU {threads} 2> {log}'
                #'&& rm -rf {params.out_dir}/_starF_checkpoints '
                #'&& rm -rf {params.out_dir}/star-fusion.preliminary'


# 2-pass alignment based on multi sample-driven splicing junction database
elif config['options']['analysis_type'] == 'junction':

    # Using shared memory with genomeLoad
    rule star_preload:
        params:
            config['star_align']['index']
        output:
            touch(path.join(bam_dir, 'star_preload.done'))
        conda:
            '../envs/star.yaml'
        log:
            path.join(log_dir, 'star_align', 'genome_preload.log')
        shell:
            'STAR --genomeLoad LoadAndExit --genomeDir {params}'
 
    rule star_align_1st:
        input:
            fq = align_inputs,
            tmp = path.join(bam_dir, 'star_preload.done')
        output:
            path.join(bam_dir, 'star_1', '{sample}', 'SJ.out.tab')
        params:
            index = config['star_align']['index'],
            options = format_options(config['star_align']['options']),
            out_prefix = path.join(bam_dir, 'star_1', '{sample}/')
        conda:
            '../envs/star.yaml'
        threads:
            config['star_align']['threads']
        log:
            path.join(log_dir, 'star_align', '{sample}_1.log')
        shell:
            'STAR {params.options} --genomeDir {params.index} '
            '--outFileNamePrefix {params.out_prefix} --runThreadN {threads} '
            '--readFilesCommand gunzip -c --outSAMtype BAM Unsorted '
            '--genomeLoad LoadAndKeep --readFilesIn {input.fq} 2> {log}'

    # filtering out the low confident junctionis
    # remove non-canonical junctions ($5>0)
    # remove junctions with # of uniquely mapped reads crossing the junction <=2 ($7>2)
    # remove annotated junctions because they will be included later stage ($6==0)
    rule filter_SJ:
        input:
            sj = expand(path.join(bam_dir, 'star_1', '{sample}', 'SJ.out.tab'), sample = all_samples)
        output:
            path.join(bam_dir, 'SJ_db', 'SJ.out.comb.tab')
        params:
            path.join(bam_dir, 'star_1')
        shell:
            "cat {input.sj} | awk '($5>0 && $7>2 && $6==0)' | cut -f1-6 | sort | uniq > {output} "
            "&& rm -rf {params}"

    rule genome_unload_1st:
        input:
            path.join(bam_dir, 'SJ_db', 'SJ.out.comb.tab')
        output:
            touch(path.join(bam_dir, 'star_unload_1st.done'))
        params:
            config['star_align']['index']
        conda:
            '../envs/star.yaml'
        shell:
            'STAR --genomeLoad Remove --genomeDir {params}'
        
    # indexing genome with the junctions obtained from all samples
    rule star_SJ_indexing:
        input:
            path.join(bam_dir, 'SJ_db', 'SJ.out.comb.tab')
        output:
            path.join(bam_dir, 'SJ_db', 'sjdbList.out.tab')
        params:
            sjdb = path.join(bam_dir, 'SJ_db'),
            fa = config['star_align']['fasta'],
            gtf = config['star_align']['gtf']
        conda:
            '../envs/star.yaml'
        threads:
            config['star_align']['threads_indexing']
        log:
            path.join(log_dir, 'star_align', 'star_SJ_indexing.log')
        shell:
            'STAR --runMode genomeGenerate --genomeDir {params.sjdb} '
            '--genomeFastaFiles {params.fa} --sjdbGTFfile {params.gtf} '
            '--runThreadN {threads} --sjdbFileChrStartEnd {input} 2>{log}'

    # Using shared memory for newly indexed genome
    rule star_preload_2nd:
        input:
            pre = path.join(bam_dir, 'star_unload_1st.done'),
            sj = path.join(bam_dir, 'SJ_db', 'sjdbList.out.tab')
        output:
            touch(path.join(bam_dir, 'star_preload_2nd.done'))
        params:
            path.join(bam_dir, 'SJ_db')
        conda:
            '../envs/star.yaml'
        log:
            path.join(log_dir, 'star_align', 'star_preload_2nd.log')
        shell:
            'STAR --genomeLoad LoadAndExit --genomeDir {params} 2> {log}'

    # STAR-2nd
    rule star_align_2nd:
        input:
            fq = align_inputs,
            tmp = path.join(bam_dir, 'star_preload_2nd.done'),
            sj = path.join(bam_dir, 'SJ_db', 'sjdbList.out.tab')
        output:
            temp(path.join(bam_dir, '{sample}','Aligned.out.bam'))
        params:
            index = path.join(bam_dir, 'SJ_db'),
            options = format_options(config['star_align']['options']),
            out_prefix = path.join(bam_dir, '{sample}/'),
            qc_prefix = path.join(qc_dir, 'star_align','{sample}'),
            qc_dir = path.join(qc_dir, 'star_align')
        conda:
            '../envs/star.yaml'
        threads:
            config['star_align']['threads']
        log:
            path.join(log_dir, 'star_align', '{sample}_2nd.log')
        shell:
            'STAR {params.options} --genomeDir {params.index} '
            '--outFileNamePrefix {params.out_prefix} --runThreadN {threads} '
            '--readFilesCommand gunzip -c --outSAMtype BAM Unsorted '
            '--genomeLoad LoadAndKeep --readFilesIn {input.fq} 2> {log} '
            '&& mkdir -p {params.qc_dir} '
            '&& mv {params.out_prefix}Log.final.out {params.qc_prefix}.Log.final.out'

    # unloading the loaded genomes
    rule genome_unload_2nd:
        input:
            expand(path.join(bam_dir, '{sample}','Aligned.out.bam'),
                sample = all_samples)
        output:
            touch(path.join(bam_dir, 'star_unload_2nd.done'))
        params:
            sj = path.join(bam_dir, 'SJ_db'),
            rm = path.join(bam_dir, 'star_1') 
        conda:
            '../envs/star.yaml'
        shell:
            'STAR --genomeLoad Remove --genomeDir {params.sj} '
            '&& rm -rf {params.rm}'


# alignment for just quantification
else:
    # Using shared memory with genomeLoad
    rule star_preload:
        params:
            config['star_align']['index']
        output:
            touch(path.join(bam_dir, 'star_preload.done'))
        conda:
            '../envs/star.yaml'
        log:
            path.join(log_dir, 'star_align', 'genome_preload.log')
        shell:
            'STAR --genomeLoad LoadAndExit --genomeDir {params}'
 
    rule star_align_normal:
        input:
            fq = align_inputs,
            tmp = path.join(bam_dir, 'star_preload.done')
        output:
            temp(path.join(bam_dir, '{sample}', 'Aligned.out.bam'))
        params:
            index = config['star_align']['index'],
            options = format_options(config['star_align']['options']),
            out_prefix = path.join(bam_dir, 'star', '{sample}/'),
            qc_prefix = path.join(qc_dir, 'star_align','{sample}'),
            qc_dir = path.join(qc_dir, 'star_align')
        conda:
            '../envs/star.yaml'
        threads:
            config['star_align']['threads']
        log:
            path.join(log_dir, 'star_align', '{sample}_normal.log')
        shell:
            'STAR {params.options} --genomeDir {params.index} '
            '--outFileNamePrefix {params.out_prefix} --runThreadN {threads} '
            '--readFilesCommand gunzip -c --outSAMtype BAM Unsorted '
            '--genomeLoad LoadAndKeep --readFilesIn {input.fq} 2> {log} '
            '&& mkdir -p {params.qc_dir} '
            '&& mv {params.out_prefix}Log.final.out {params.qc_prefix}.Log.final.out'


