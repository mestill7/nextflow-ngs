process RUN_CLEANUP {
    publishDir "${params.project_dir}/output", mode: 'copy'

    input:
    tuple val(pair_id), path(bam_file)

    output:
    path("removal_summary.txt"), emit: report


    script:
    """
    if [[ ${params.experiment_type} == 'rnaseq' ]]; then
       rm ${params.project_dir}/output/bam/*.sorted.bam ${params.project_dir}/output/bam/*.sorted.bam.bai ${params.project_dir}/output/bam/*.unique.sorted.rmdup.bam
       echo "File cleanup completed" > removal_summary.txt
    elif [[ ${params.experiment_type} == 'cutrun' ]]; then
       rm ${params.project_dir}/output/bam/*_paired.bam ${params.project_dir}/output/bam/*_dis.bam ${params.project_dir}/output/bam/*_dis_sorted.bam ${params.project_dir}/output/bam/*.sorted.bam ${params.project_dir}/output/bam/*.sorted.bam.bai ${params.project_dir}/output/bam/*.unique.sorted.rmdup.bam
       echo "File cleanup completed" > removal_summary.txt
    else
       rm ${params.project_dir}/output/bam/*.sorted.bam ${params.project_dir}/output/bam/*.sorted.bam.bai ${params.project_dir}/output/bam/*.unique.sorted.rmdup.bam
       echo "File cleanup completed" > removal_summary.txt
    fi
    """
}
