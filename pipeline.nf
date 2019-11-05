Channel
    .fromFilePairs("${params.input_dir}/raw/*.{bed,fam,bim}",size:3)
    .ifEmpty {error "No files in  ${params.input_dir}."}
    .filter  {key, filename  -> key in params.input_pat}
    .into {bfile_ch1; bfile_ch2; bfile_ch3}

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
    .into {snplist_ch1; snplist_ch2}


process chi2_per_snp{
    echo true
    publishDir "${params.output_dir}/assoc/",
                pattern: "*.assoc",
                overwrite:true,
                mode:"copy"

    cpus 3
    memory "6GB"
    maxForks 5

    input:
        each snp from snplist_ch1
        tuple bfile, file(bfileNames) from bfile_ch2

    output:
        file("*.assoc") into snp_chi2_ch

    when:
        snp[0]==bfile

    script:
    """
    plink --noweb --bfile ${bfile} --assoc --out ${bfile}-\$( echo '${snp[1]}'| cut -d':' -f 2) --silent --snp ${snp[1]}
    """

}

process maxT_per_snp{

    echo true
    publishDir "${params.output_dir}/maxT/",
                pattern:"*mperm",
                overwrite: true,
                mode:'copy'

    cpus 3
    memory "6GB"
    maxForks 10

    input:
        each snp from snplist_ch2
        tuple bfile, file(bfileList) from bfile_ch3

    output:
        file("*mperm") into snp_maxT_ch

    when:
        snp[0]==bfile

    script:
    """
    echo plink --noweb --bfile ${bfile} --assoc --mperm ${params.mperm} --out ${bfile}-\$( echo '${snp[1]}'| cut -d':' -f 2) --silent --snp ${snp[1]}
    """
}

