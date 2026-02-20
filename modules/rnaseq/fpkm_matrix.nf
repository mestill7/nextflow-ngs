process FPKM_MATRIX {

    publishDir "${params.project_dir}/output/counts", mode: 'copy'

    input:
    path("output/counts/*txt")

    output:
    path("FPKM_genic_matrix.txt"), emit: fpkm_matrix_genic
    path("FPKM_exonic_matrix.txt"), emit: fpkm_matrix_exonic

    script:
    """
    #!/usr/bin/env python

    import sys
    import os
    import pandas as pd

    a=os.listdir(path="${params.project_dir}output/counts/")
    gene_files=[x for x in a if x.endswith(".gene.txt")]
    exon_files=[x for x in a if x.endswith(".exon.txt")]

    ## Treat each count file individually
    def myfun(x):
        x_counts = pd.read_csv("${params.project_dir}output/counts/"+x, sep='\t', comment='#')
        x_counts.sort_values('Geneid',axis=0,inplace=True)
        gl = x_counts.Length
        colsum = x_counts.iloc[:,6].sum()/(1e6)
        rpkm = x_counts.iloc[:,6]/(gl/1000)/colsum
        x_counts['FPKM']=rpkm
        return x_counts

    gene_fpkm_list = [myfun(x) for x in gene_files]
    exon_fpkm_list = [myfun(x) for x in exon_files]

    def save_fpkm(mylist,output_file):
        if len(set([x.shape for x in mylist]))==1 :
        ## Collect sample names
            sample_labels = [list(x.columns)[-2].replace('.unique.sorted.rmdup.bam','') for x in mylist]
            ## Merge gene info and fpkm for all samples
            gene_info = mylist[0].iloc[:,:6]
            cur_fpkm = pd.DataFrame({sample_labels[i]:mylist[i].FPKM for i in range(len(mylist))})
            cur_fpkm = pd.concat([gene_info,cur_fpkm],ignore_index=False,axis=1)
            cur_fpkm.to_csv(output_file, index=False, sep='\t')
            return cur_fpkm

    save_fpkm(exon_fpkm_list,"FPKM_exonic_matrix.txt")
    save_fpkm(gene_fpkm_list,"FPKM_genic_matrix.txt")
    """
}
