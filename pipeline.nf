Channel
    .fromFilePairs("${params.input_dir}/raw/*.{bed,fam,bim}",size:3)
    .ifEmpty {error "No files in  ${params.input_dir}."}
    .filter  {key, filename  -> key in params.input_pat}
    .into {bfile_ch1; bfile_ch2}

process list_snps{

    echo true
    publishDir "${params.output_dir}/snplist/",
                pattern:"*.snplist",
                overwrite:true,
                mode:'copy'

    
    input:
     set bfile, file(bfileNames) from bfile_ch1

    output:
    file("*.snplist") into snpfile_ch

    script:
    """
    plink --bfile ${bfile} --write-snplist --out ${bfile}
    """

}


snpfile_ch
    .map{file ->
            def key  = file.name.tokenize(".").get(0)
            def snps = file.splitText()
            return tuple(key, snps)
    }
    .transpose()
    .set{snplist_ch}


process maxT_per_snp{

    echo true
    publishDir "${params.output_dir}/maxT/",
                pattern:"*.mperm.*",
                overwrite: true,
                mode:'copy'

    input:
        each snp from snplist_ch
        tuple bfile, file(bfileList) from bfile_ch2

    output:
        file("*.mperm.*") into analysis_ch

    when:
        snp[0]==bfile

    script:
    """
    plink --bfile ${bfile}  --model mperm=${params.mperm} --mperm-save --out ${bfile}-\$( echo '${snp[1]}'|cut -d':' -f 2) --snp ${snp[1]}
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

