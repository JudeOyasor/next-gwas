Channel
    .fromFilePairs("${params.input_dir}/raw/*.{bed,fam,bim}",size:3)
    .ifEmpty {error "No files in  ${params.input_dir}."}
    .filter  {key, filename  -> key in params.input_pat}
    .into {bfile_ch1; bfile_ch2}

process make_snpblocks{

    echo true
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


process maxT{

    echo true
    publishDir "${params.output_dir}/maxT/",
                pattern:"*.model.best.mperm",
                overwrite: true,
                mode:'copy'

    input:
        each snp from snpblocks_ch
        tuple bfile, file(bfileList) from bfile_ch2

    output:
        file("*.model.best.mperm") into analysis_ch

    when:
        snp[0]==bfile

    script:
    """
    plink --bfile ${bfile}  --model mperm=${params.mperm} --mperm-save --out ${bfile} --snps ${snp[1].FIRST}-${snp[1].LAST}
    """
}


process merge_files{

    echo true
    
    input:
        file(maxT) from analysis_ch

    output:

    script:
    """
        echo $maxT
    """
}

