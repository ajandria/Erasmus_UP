#!usr/bin/env  nextflow

/*
 * Define output folders
 */ 


fastq_in = '/home/ajan/fastq_data'
fastqc_raw_out = '/home/ajan/shared_folder/nextflow_liam/fastqc_raw'
fastp_out = '/home/ajan/shared_folder/nextflow_liam/fastp'
fastqc_fastp_out = '/home/ajan/shared_folder/nextflow_liam/fastqc_fastp_out'
salmon_index_in = '/home/ajan/shared_folder/salmon_index/index'
salmon_quant_out = '/home/ajan/shared_folder/nextflow_liam/salmon_quant_out'
gtf_in = '/home/ajan/shared_folder/star_hg19_run/chr_primary_assembly/Homo_sapiens.GRCh37.87.chr.gtf'
multiqc_fastqc_raw_out = '/home/ajan/shared_folder/nextflow_liam/multiqc/multiqc_fastqc_raw'
multiqc_fastqc_fastp_out = '/home/ajan/shared_folder/nextflow_liam/multiqc/multiqc_fastqc_fastp'
multiqc_fastp_out = '/home/ajan/shared_folder/nextflow_liam/multiqc/multiqc_fastp'
multiqc_salmon_quant_out = '/home/ajan/shared_folder/nextflow_liam/multiqc/multiqc_salmon_quant'

log.info """\


===========================================================================================================
 Transcriptomics 0.2a                                                                                     
==========================================================================================================
|                                                                                                          
| fastqRaw path                          : $fastq_in                                                      
| Salmon index dir                       : $salmon_index_in                                               
| GTF path                               : $gtf_in
| Head output dir                        : /home/ajan/shared_folder/nextflow_liam                                             
|                                                                                
=========================================================================================================


"""

/*
 *  Parse the input parameters
 */
Channel
.fromFilePairs("$fastq_in/*_{1,2}.fastq.gz", flat: true)
.into{ raw_fastq_fastqc; raw_fastq_fastp}

/*
 * Process 1: Fastqc on raw fastq
 */

process fastqc_raw {

    conda '/home/ajan/.conda/envs/salmon_hpc'

    publishDir "$fastqc_raw_out", mode:'copy'
   
    input:
        set val(id), file(fastq_1), file(fastq_2) from raw_fastq_fastqc
    
    output:
        file("*.zip") into fastqc_multiqc
        file "*.html"

    script:
        """
        fastqc -t ${task.cpus} ${fastq_1} ${fastq_2}
	"""
}

/*
 * Process 2: Fastp on raw fastq
 */
 
process fastp {
    
    conda '/home/ajan/.conda/envs/salmon_hpc'

    publishDir "$fastp_out", mode:'copy', pattern: '*.fastq.gz'
    publishDir "$fastp_out", mode:'copy', pattern: '*fastp.json'
    publishDir "$fastp_out", mode:'copy', pattern: '*fastp.html'

    input:
        set val(id), file(fastq_1), file(fastq_2) from raw_fastq_fastp
    
    output:
        set val(id), file("${id}_1_fastp.fastq.gz"), file("${id}_2_fastp.fastq.gz") into (fastp_fastqc, fastp_salmon)
        file("${id}_fastp.json") into fastp_multiqc
        file "${id}_fastp.html"

    script:
        """
	fastp -i ${fastq_1} -I ${fastq_2} -o ${id}_1_fastp.fastq.gz -O ${id}_2_fastp.fastq.gz -j ${id}_fastp.json -h ${id}_fastp.html -q 30 -e 25 -n 3 -l 40 -c -x 20 -w ${task.cpus}

	"""
}  

 /*
  * Process 3: Fastqc on filtered with fastp fastq
  */

 process fastqc_fastp_out {

     conda '/home/ajan/.conda/envs/salmon_hpc'

     publishDir "$fastqc_fastp_out", mode:'copy'
 
     input:
         set val(id), file(fastq_1), file(fastq_2) from fastp_fastqc
  
     output:
         file("*.zip") into fastqc_fastp_multiqc
         file "*.html"

     script:
         """
         fastqc -t ${task.cpus} ${fastq_1} ${fastq_2}
 	"""
 }

 /*
  * Process 4: Run salmon 
  */

 process salmon_quant {

    conda '/home/ajan/.conda/envs/salmon'

     publishDir "$salmon_quant_out", mode:'copy'
   
     input:
         set val(id), file(fastq_1), file(fastq_2) from fastp_salmon
    
     output:
         path id into salmonResults_dir

     script:
         """
         salmon quant -i $salmon_index_in -l ISR -1 ${fastq_1} -2 ${fastq_2} -p ${task.cpus} -o $id -g $gtf_in
 	"""
 }

/*
 * Process 5: Run multiqc on fastp filtered fastq
 */

process multiqc_fastp {

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_fastp_out", mode:'copy'
   
    input:
        file("*") from fastp_multiqc.collect()
    
    output:
        file("multiqc_fastp.html")

    script:
        """
        multiqc *.json -n multiqc_fastp
	"""
}

/*
 * Process 6: Run multiqc on raw fastq
 */

process multiqc_fastqc_raw {

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_fastqc_raw_out", mode:'copy'
   
    input:
        file("*") from fastqc_multiqc.collect()
    
    output:
        file("multiqc_fastqc_raw.html")

    script:
        """
        multiqc *.zip

        mv multiqc_report.html multiqc_fastqc_raw.html
	"""
}

/*
 * Process 7: Run multiqc on fastp filtered FastQC
 */

process multiqc_fastqc_fastp {

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_fastqc_fastp_out", mode:'copy'
   
    input:
        file("*") from fastqc_fastp_multiqc.collect()
    
    output:
        file("multiqc_fastqc_fastp.html")

    script:
        """
        multiqc *.zip

        mv multiqc_report.html multiqc_fastqc_fastp.html
	"""
}

/*
 * Process 8: Run multiqc on salmon
 */

process multiqc_salmon_quant {

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_salmon_quant_out", mode:'copy'
   
    input:
        file("*") from salmonResults_dir.collect()
    
    output:
        file("multiqc_salmon_quant.html")

    script:
        """
        multiqc .

        mv multiqc_report.html multiqc_salmon_quant.html
	"""
}

workflow.onComplete {
log.info ( workflow.success ? "\n The workflow was complete!" : "Oops .. something went wrong" )
}
