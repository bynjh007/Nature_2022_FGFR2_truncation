rule star_fusion:
    input:
        path.join(chim_dir, '{sample}', 'Chimeric.out.junction')
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