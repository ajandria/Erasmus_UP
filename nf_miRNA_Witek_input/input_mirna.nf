
/*
 * Define output folders
 */ 

//astq_in = '/home/ajan/shared_mirna/subset'
fastq_in = '/home/ajan/shared_mirna/input'
fastqc_raw_out = '/home/ajan/shared_mirna/input_mirna/fastqc_raw'
fastp_out = '/home/ajan/shared_mirna/input_mirna/fastp'
fastqc_fastp_out = '/home/ajan/shared_mirna/input_mirna/fastqc_fastp_out'
star_index_in = '/home/ajan/shared_mirna/liam/index'
star_alignments_out = '/home/ajan/shared_mirna/input_mirna/alignments'
gtf_in = '/home/ajan/shared_mirna/chr_primary_assembly/mirna_Homo_sapiens.GRCh37.87.chr.gtf'
qualimap_out = '/home/ajan/shared_mirna/input_mirna/qualimap'
picard_out = '/home/ajan/shared_mirna/input_mirna/picard'
ref_flat_in = '/home/ajan/shared_mirna/chr_primary_assembly/refFlat_mirna.txt'
multiqc_fastqc_raw_out = '/home/ajan/shared_mirna/input_mirna/multiqc/fastqc_raw'
multiqc_fastqc_fastp_out = '/home/ajan/shared_mirna/input_mirna/multiqc/fastqc_fastp'
multiqc_fastp_out = '/home/ajan/shared_mirna/input_mirna/multiqc/fastp'
multiqc_star_alignments_out = '/home/ajan/shared_mirna/input_mirna/multiqc/star_alignments'
multiqc_qualimap_out = '/home/ajan/shared_mirna/input_mirna/multiqc/qualimap'
multiqc_picard_out = '/home/ajan/shared_mirna/input_mirna/multiqc/picard'
samtools_index_out = star_alignments_out
samtools_flagstat_out = '/home/ajan/shared_mirna/input_mirna/samtools_flagstat'
multiqc_samtools_flagstat_out = '/home/ajan/shared_mirna/input_mirna/multiqc/samtools_flagstat'
multiqc_featureCounts_out = '/home/ajan/shared_mirna/input_mirna/multiqc/featureCounts'
featureCounts_out = '/home/ajan/shared_mirna/input_mirna/featureCounts'

log.info """\


====================================================================================================================

 Transcriptomics 0.2a    

=====================================================================================================================
|                                                                                                                   
| .fastq files                           : $fastq_in                                                                
| STAR index                             : $star_index_in                                                           
| GTF                                    : $gtf_in
| refFlat                                : $ref_flat_in 
| Head output dir                        : /home/ajan/shared_mirna                             
|                                                                                                                     
=====================================================================================================================


"""

/*
 *  Parse the input parameters
 */
Channel
.fromPath("${fastq_in}/*.fastq.gz")
.into{ raw_fastq_fastqc; raw_fastq_fastp}

/*
 * Process 1: Fastqc on raw fastq
 */

process fastqc_raw {

    label 'fastqc_raw'

    conda '/home/ajan/.conda/envs/fastqc'

    publishDir "$fastqc_raw_out", mode:'copy'
   
    input:
        file(fastq_file) from raw_fastq_fastqc
    
    output:
        file ("*.zip") into fastqc_multiqc
        file "*.html"

    script:
        """
        fastqc -t ${task.cpus} ${fastq_file}
	"""
}

/*
 * Process 2: Filtering
 */
process fastp {

    label 'fastp'

    conda '/home/ajan/.conda/envs/fastp'
    
    publishDir "$fastp_out", mode:'copy', pattern: '*.fastq.gz'
    publishDir "$fastp_out", mode:'copy', pattern: '*fastp.json'
    publishDir "$fastp_out", mode:'copy', pattern: '*fastp.html'

    input:
        file(fastq_file) from raw_fastq_fastp
    
    output:
        file("${fastq_file}_fastp.fastq.gz") into (fastp_fastqc, fastp_star)
        file("${fastq_file}_fastp.json") into fastp_multiqc
        file "${fastq_file}_fastp.html"

    script:
        """
	fastp -i ${fastq_file} -o ${fastq_file}_fastp.fastq.gz -j ${fastq_file}_fastp.json -h ${fastq_file}_fastp.html -q 30 -e 25 -n 3 -l 19 -c -x 20 -w ${task.cpus}

	"""
}  

 /*
  * Process 3: Fastqc on filtered
  */

 process fastqc_fastp {

     label 'fastqc_fastp'

     conda '/home/ajan/.conda/envs/fastqc'

     publishDir "$fastqc_fastp_out", mode:'copy'
 
     input:
         file(fastq_file) from fastp_fastqc
  
     output:
         file("*.zip") into fastqc_fastp_multiqc
         file "*.html"

     script:
         """
         fastqc -t ${task.cpus} ${fastq_file}
 	"""
 }

 /*
  * Process 4: Run STAR 
  */

 process star_align {

     label 'star_align'

     conda '/home/ajan/.conda/envs/star'

     publishDir "$star_alignments_out", mode:'copy'
   
     input:
         file(fastq_file) from fastp_star
    
     output:
         file ("${fastq_file}*Aligned.out.bam") into (star_qualimap, star_samtools_index, star_samtools_flagstat, star_picard, star_subread, star_featureCounts)
         path ("${fastq_file}*") into star_alignments_path

     script:
         """
         STAR --runThreadN ${task.cpus} --runMode alignReads --genomeDir $star_index_in --readFilesIn ${fastq_file} --readFilesCommand zcat --readStrand Unstranded --outSAMtype BAM Unsorted --sjdbGTFfile $gtf_in --sjdbOverhang 1 --outFileNamePrefix ${fastq_file} --twopassMode Basic --outFilterMismatchNoverLmax 0.05 --alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMultimapScoreRange 0 --quantMode TranscriptomeSAM GeneCounts --outReadsUnmapped Fastx --outFilterMultimapNmax 10 --outSAMunmapped Within --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 16 --alignSJDBoverhangMin 1000 --alignIntronMax 1

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

      label 'qualimap'

     conda '/home/ajan/.conda/envs/qualimap'

      publishDir "$qualimap_out", mode:'copy'
   
      input:
          file(bam) from star_qualimap
    
      output:
          path ("*") into qualimap_out_path

      script:
          """
          export JAVA_OPTS="-Djava.io.tmpdir=/home/ajan/shared_mirna/input_mirna/trash"

          qualimap rnaseq -bam ${bam} -gtf ${gtf_in} -p non-strand-specific -outformat PDF:HTML --java-mem-size=16G
 	"""
  }

 process picard_matrix {

     label 'picard_matrix'

    conda '/home/ajan/.conda/envs/picard'

     publishDir "$picard_out", mode:'copy'
   
     input:
         file(bam) from star_picard
    
     output:
         path ("*") into picard_out_path

     script:
         """
         picard CollectRnaSeqMetrics -I ${bam} -O ${bam}.txt --REF_FLAT ${ref_flat_in} --STRAND NONE
	"""
 }

 process samtools_flagstat {

     label 'samtools_flagstat'

    conda '/home/ajan/.conda/envs/samtools'

     publishDir "$samtools_flagstat_out", mode:'copy'
   
     input:
         file(bam) from star_samtools_flagstat
    
     output:
         path ("*") into samtools_flagstat_out_path

     script:
         """
         samtools flagstat -@ ${task.cpus} ${bam} > ${bam}.txt
	"""
 }

 process featureCounts {

     label 'featureCounts'

    conda '/home/ajan/.conda/envs/subread'

     publishDir "$featureCounts_out", mode:'copy'
   
     input:
         file(bam) from star_featureCounts
    
     output:
         path ("*.summary") into featureCounts_multiqc
         file ("*")

     script:
         """
         featureCounts -t gene -s 0 -T ${task.cpus} --verbose -M -a ${gtf_in} -o ${bam}_featureCounts_matrix.txt ${bam}
	"""
 }

/*
 * Process 5: Run multiqc on qualimap
 */

process multiqc_featureCounts {

    label 'multiqc_featureCounts'

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

    label 'multiqc_samtools_flagstat'

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

     label 'multiqc_qualimap'

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

    label 'multiqc_picard'

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

    label 'multiqc_fastp'

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

    label 'multiqc_fastqc_raw'

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

    label 'multiqc_fastqc_fastp'

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

    label 'multiqc_star_alignments'

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
