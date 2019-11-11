Channel
    .fromFilePairs("${params.input_dir}/raw/*.{bed,fam,bim}",size:3)
    .ifEmpty {error "No files in  ${params.input_dir}."}
    .filter  {key, filename  -> key in params.input_pat}
    .into {bfile_ch1; bfile_ch2}

process make_snpblocks{

    publishDir "${params.output_dir}/snpblocks/",
                pattern:"*.var.ranges",
                overwrite:true,
                mode:'copy'

    
    input:
     set bfile, file(bfileNames) from bfile_ch1

    output:
    file("*.var.ranges") into snpblock_files_ch

    script:
    """
    plink --bfile ${bfile} --write-var-ranges ${params.nblocks} --out ${bfile}
    """

}


snpblock_files_ch
    .map{file ->
            def key  = file.name.tokenize(".").get(0)
            def snps = file.splitCsv(header:true, sep: "\t")
            return tuple(key, snps)
    }
    .transpose()
    .set{snpblocks_ch}


process compute_maxT{
   
    input:
        each snp from snpblocks_ch
        tuple bfile, file(bfileList) from bfile_ch2

    output:
        file("*.assoc.mperm") into maxTfile_ch
        file("*.assoc") into assocfile_ch

    when:
        snp[0]==bfile

    script:
    """
    plink --bfile ${bfile} --assoc mperm=${params.mperm} --seed 567489 --mperm-save --out ${bfile} --snps ${snp[1].FIRST}-${snp[1].LAST}
    """
}


maxTfile_ch
    .collectFile(keepHeader:true, storeDir:"${params.output_dir}/data")
    .set{maxT_ch}

assocfile_ch
    .collectFile(keepHeader:true, storeDir:"${params.output_dir}/data")
    .set{assoc_ch}


process plot_graphs{

    echo true
    publishDir "${params.output_dir}/",
                 pattern:"*.png",
                 overwrite:true,
                 mode:'copy' 

    input:
        file(assoc) from assoc_ch
        file(maxT) from maxT_ch

    output:
        file("*.png")
    script:
    """
    python3 ${params.script_dir}/graphs.py --assoc ${assoc} --maxT ${maxT}
    """
} 

