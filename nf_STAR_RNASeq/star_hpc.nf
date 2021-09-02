
/*
 * Define output folders
 */ 

fastq_in = '/home/ajan/fastq_data'
fastqc_raw_out = '/home/ajan/shared_folder/nextflow_liam_star/fastqc_raw'
fastp_out = '/home/ajan/shared_folder/nextflow_liam_star/fastp'
fastqc_fastp_out = '/home/ajan/shared_folder/nextflow_liam_star/fastqc_fastp_out'
star_index_in = '/home/ajan/shared_folder/star_hg19_run/index'
star_alignments_out = '/home/ajan/shared_folder/nextflow_liam_star/star_alignments'
gtf_in = '/home/ajan/shared_folder/star_hg19_run/chr_primary_assembly/Homo_sapiens.GRCh37.87.chr.gtf'
qualimap_out = '/home/ajan/shared_folder/nextflow_liam_star/qualimap'
picard_out = '/home/ajan/shared_folder/nextflow_liam_star/picard'
multiqc_fastqc_raw_out = '/home/ajan/shared_folder/nextflow_liam_star/multiqc/fastqc_raw'
multiqc_fastqc_fastp_out = '/home/ajan/shared_folder/nextflow_liam_star/multiqc/fastqc_fastp'
multiqc_fastp_out = '/home/ajan/shared_folder/nextflow_liam_star/multiqc/fastp'
multiqc_star_alignments_out = '/home/ajan/shared_folder/nextflow_liam_star/multiqc/star_alignments'
multiqc_qualimap_out = '/home/ajan/shared_folder/nextflow_liam_star/multiqc/qualimap'
multiqc_picard_out = '/home/ajan/shared_folder/nextflow_liam_star/multiqc/picard'
ref_flat_in = '/home/ajan/shared_folder/star_hg19_run/chr_primary_assembly/noChr_refFlat.txt'
samtools_index_out = star_alignments_out
samtools_flagstat_out = '/home/ajan/shared_folder/nextflow_liam_star/samtools_flagstat'
multiqc_samtools_flagstat_out = '/home/ajan/shared_folder/nextflow_liam_star/multiqc/samtools_flagstat'
multiqc_featureCounts_out = '/home/ajan/shared_folder/nextflow_liam_star/multiqc/featureCounts'
featureCounts_out = '/home/ajan/shared_folder/nextflow_liam_star/featureCounts'

log.info """\


====================================================================================================================

 Transcriptomics 0.2a    

=====================================================================================================================
|                                                                                                                   
| .fastq files                           : $fastq_in                                                                
| STAR index                             : $star_index_in                                                           
| GTF                                    : $gtf_in
| refFlat                                : $ref_flat_in 
| Head output dir                        : /home/ajan/shared_folder/nextflow_liam_star                              
|                                                                                                                     
=====================================================================================================================


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
        file ("*.zip") into fastqc_multiqc
        file "*.html"

    script:
        """
        fastqc -t ${task.cpus} ${fastq_1} ${fastq_2}
	"""
}

/*
 * Process 2: Filtering
 */
 
process fastp {

    conda '/home/ajan/.conda/envs/salmon_hpc'
    
    publishDir "$fastp_out", mode:'copy', pattern: '*.fastq.gz'
    publishDir "$fastp_out", mode:'copy', pattern: '*fastp.json'
    publishDir "$fastp_out", mode:'copy', pattern: '*fastp.html'

    input:
        set val(id), file(fastq_1), file(fastq_2) from raw_fastq_fastp
    
    output:
        set val(id), file("${id}_1_fastp.fastq.gz"), file("${id}_2_fastp.fastq.gz") into (fastp_fastqc, fastp_star)
        file("${id}_fastp.json") into fastp_multiqc
        file "${id}_fastp.html"

    script:
        """
	fastp -i ${fastq_1} -I ${fastq_2} -o ${id}_1_fastp.fastq.gz -O ${id}_2_fastp.fastq.gz -j ${id}_fastp.json -h ${id}_fastp.html -q 30 -e 25 -n 3 -l 40 -c -x 20 -w ${task.cpus}

	"""
}  

 /*
  * Process 3: Fastqc on filtered
  */

 process fastqc_fastp {

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
  * Process 4: Run STAR 
  */

 process star_align {

     conda '/home/ajan/.conda/envs/star'

     publishDir "$star_alignments_out", mode:'copy'
   
     input:
         set val(id), file(fastq_1), file(fastq_2) from fastp_star
    
     output:
         set val(id), file ("${id}*.bam") into (star_qualimap, star_samtools_index, star_samtools_flagstat, star_picard, star_subread, star_featureCounts)
         path ("${id}*") into star_alignments_path

     script:
         """
         STAR --runThreadN ${task.cpus} --runMode alignReads --genomeDir $star_index_in --readFilesIn ${fastq_1} ${fastq_2} --readFilesCommand zcat --readStrand Reverse --outSAMtype BAM SortedByCoordinate --sjdbGTFfile $gtf_in --sjdbOverhang 50 --outFileNamePrefix ${id}
 	"""
 }

 /*
  * Process 5: Run Qualimap 
  *
  * Process *5*: Run samtools sort by name for qualimap to omit // this is possible but to tidy the code up option:
  *                                                           // -Djava.io.tmpdir=<my/temp/path>
  *                                                           // https://www.biostars.org/p/42613/
  *                                                           // will be used to set temp directory for bam files for the time of sorting
  */

 process qualimap {

    conda '/home/ajan/.conda/envs/qualimap'

     publishDir "$qualimap_out", mode:'copy'
   
     input:
         set val(id), file(bam) from star_qualimap
    
     output:
         path id into qualimap_out_path

     script:
         """
         export JAVA_OPTS="-Djava.io.tmpdir=/home/ajan/shared_folder/nextflow_liam_star/trash"

         qualimap rnaseq -bam ${bam} -gtf ${gtf_in} -outdir ${id}/ -pe -p strand-specific-reverse -outformat PDF:HTML --java-mem-size=16G
	"""
 }

 process picard_matrix {

    conda '/home/ajan/.conda/envs/picard'

     publishDir "$picard_out", mode:'copy'
   
     input:
         set val(id), file(bam) from star_picard
    
     output:
         path id into picard_out_path

     script:
         """
         picard CollectRnaSeqMetrics -I ${bam} -O ${id} --REF_FLAT ${ref_flat_in} --STRAND SECOND_READ_TRANSCRIPTION_STRAND
	"""
 }

 process samtools_index {

    conda '/home/ajan/.conda/envs/samtools'

     publishDir "$samtools_index_out", mode:'copy'
   
     input:
         set val(id), file(bam) from star_samtools_index
    
     output:
         file("*")

     script:
         """
         samtools index -@ 4 ${bam}
	"""
 }

 process samtools_flagstat {

    conda '/home/ajan/.conda/envs/samtools'

     publishDir "$samtools_flagstat_out", mode:'copy'
   
     input:
         set val(id), file(bam) from star_samtools_flagstat
    
     output:
         path ("${id}.txt") into samtools_flagstat_out_path

     script:
         """
         samtools flagstat ${bam} -@ 4 > ${id}.txt
	"""
 }

 process featureCounts {

    conda '/home/ajan/.conda/envs/subread'

     publishDir "$featureCounts_out", mode:'copy'
   
     input:
         set val(id), file(bam) from star_featureCounts
    
     output:
         path ("*.summary") into featureCounts_multiqc
         file ("*")

     script:
         """
         featureCounts -t gene -s 2 -T ${task.cpus} --verbose -p -a ${gtf_in} -o ${id}_featureCounts_matrix.txt ${bam}
	"""
 }

/*
 * Process 5: Run multiqc on qualimap
 */

process multiqc_featureCounts {

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_featureCounts_out", mode:'copy'
   
    input:
        file ("*") from featureCounts_multiqc.collect()
    
    output:
        file("multiqc_featureCounts.html")

    script:
        """
        multiqc . -n multiqc_featureCounts
	"""
}

process multiqc_samtools_flagstat {

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_samtools_flagstat_out", mode:'copy'
   
    input:
        file("*") from samtools_flagstat_out_path.collect()
    
    output:
        file("multiqc_samtools_flagstat.html")

    script:
        """
        multiqc . -n multiqc_samtools_flagstat
	"""
}

/*
 * Process 5: Run multiqc on qualimap
 */

process multiqc_qualimap {

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_qualimap_out", mode:'copy'
   
    input:
        file("*") from qualimap_out_path.collect()
    
    output:
        file("multiqc_qualimap.html")

    script:
        """
        multiqc . -n multiqc_qualimap
	"""
}

/*
 * Process 5: Run multiqc on picard
 */

process multiqc_picard {

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_picard_out", mode:'copy'
   
    input:
        file("*") from picard_out_path.collect()
    
    output:
        file("multiqc_picard.html")

    script:
        """
        multiqc . -n multiqc_picard
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

process multiqc_star_alignments {

    conda '/home/ajan/.conda/envs/multiqc'

    publishDir "$multiqc_star_alignments_out", mode:'copy'
   
    input:
        file("*") from star_alignments_path.collect()
    
    output:
        file("multiqc_star_alignments.html")

    script:
        """
        multiqc .

        mv multiqc_report.html multiqc_star_alignments.html
	"""
}

workflow.onComplete {
log.info ( workflow.success ? "\n The workflow was complete!" : "Oops .. something went wrong" )
}
