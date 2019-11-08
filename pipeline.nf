Channel
    .fromFilePairs("${params.input_dir}/raw/*.{bed,fam,bim}",size:3)
    .ifEmpty {error "No files in  ${params.input_dir}."}
    .filter  {key, filename  -> key in params.input_pat}
    .into {bfile_ch1; bfile_ch2}

process list_snps{

    cpus 2 
    memory "4GB"
    maxForks 2

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
    plink --noweb --bfile ${bfile} --silent --write-snplist --out ${bfile}
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


process analysis_per_snp{

    echo true
    publishDir "${params.output_dir}/maxT/",
                pattern:"*.mperm.*",
                overwrite: true,
                mode:'copy'
    
    publishDir "${params.output_dir}/assoc/",
                pattern: "*.assoc",
                overwrite:true,
                mode:"copy"

    cpus 3
    memory "6GB"
    maxForks 10

    input:
        each snp from snplist_ch
        tuple bfile, file(bfileList) from bfile_ch2

    output:
    tuple file("*.assoc"), file("*.mperm.*") into analysis_ch

    when:
        snp[0]==bfile

    script:
    """
    plink --noweb --bfile ${bfile} --assoc --mperm ${params.mperm} --mperm-save --out ${bfile}-\$( echo '${snp[1]}'|cut -d':' -f 2) --silent --snp ${snp[1]}
    """
}

process merge_files{

    echo true
    
    input:
      set file(assoc), file(maxT) from analysis_ch

    output:

    script:
    """
        echo $assoc
        echo $maxT
    """
}

