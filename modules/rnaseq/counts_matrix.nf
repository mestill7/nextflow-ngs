process COUNTS_MATRIX {
    
    publishDir "${params.project_dir}/output/counts", mode: 'copy'

    input:
    path("output/counts/*txt")

    output:
    path("featureCounts_genic.txt"), emit: count_matrix_genic
    path("featureCounts_exonic.txt"), emit: count_matrix_exonic

    script:
    """
    #!/usr/bin/env python

    import sys
    import os 
    import pandas as pd

    
    gene_counts = {}
    exon_counts = {}
    
    for files in os.listdir("${params.project_dir}output/counts/"):
        if files.endswith("gene.txt"):
            counts_dict = gene_counts
        elif files.endswith("exon.txt"):
            counts_dict = exon_counts
        else:
            continue
            
        sample = files.split(".")[0] 
        counts_dict[sample] = {}
        with open("${params.project_dir}output/counts/"+files, "r") as infile:
            next(infile)
            next(infile)
            for lines in infile:
                lines = lines.strip().split("\t")
                counts_dict[sample][lines[0]] = int(float(lines[-1]))
    
    gene_dataframe = pd.DataFrame(gene_counts)
    exon_dataframe = pd.DataFrame(exon_counts)
    gene_dataframe.to_csv("featureCounts_genic.txt", sep="\t")
    exon_dataframe.to_csv("featureCounts_exonic.txt", sep="\t")
    """


}
