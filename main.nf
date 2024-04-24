#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
      stageInMode 'symlink'
      stageOutMode 'move'

      script:
        """

        if [[ "${params.containers}" == "singularity" ]] ; 

          then

            cd ${params.image_folder}
            
            if [[ ! -f star-2.7.11b.sif ]] ;
              then
                singularity pull star-2.7.11b.sif docker://index.docker.io/mpgagebioinformatics/star:2.7.11b
            fi
        
            if [[ ! -f samtools-1.16.1.sif ]] ;
              then
                singularity pull samtools-1.16.1.sif docker://index.docker.io/mpgagebioinformatics/samtools:1.16.1
            fi

        fi


        if [[ "${params.containers}" == "docker" ]] ; 

          then

            docker pull mpgagebioinformatics/star:2.7.11b
            docker pull mpgagebioinformatics/samtools:1.16.1

        fi

        """
    }


process star_indexer {
    stageInMode 'symlink'
    stageOutMode 'move'

    when:
    ( ! file("${star_index}index.fa".exists() ) 
  
    script:
    
        """
        
        if [[ ! -e \$star_index ]] ; then mkdir -p \$star_index ; fi

        cd \$star_index
    
        STAR --runMode genomeGenerate \
             --genomeDir ${star_index} \ #output_folder
             --genomeFastaFiles ../genomes/${params.organism}/${params.release}/original.toplevel.fa \ #reference genome FASTA file
             --sjdbGTFfile ../genomes/${params.organism}/${params.release}/original.gtf \
             --runThreadN 12
        
        """
    }


process star_mapping {

    stageInMode 'symlink'
    stageOutMode 'move'

    input:
    tuple val(pair_id), path(fastq)

    output:
    val pair_id

    when:
    ( ! file("${params.project_folder}/star_output/${pair_id}.bam").exists() ) 
  
    script:
    
    def single = fastq instanceof Path
    
if ( single ) {

    """
      mkdir -p /$star_output
      cd /raw_data
      
      echo ${pair_id}
      
      
      STAR --readFilesCommand zcat --runThreadN 12 --genomeDir ${star_index} --outSAMtype BAM SortedByCoordinate \
           --limitBAMsortRAM 15000000000 --readFilesIn ${renamed}${sample} ${renamed}${sample%${read1_sufix}}${read2_sufix} \
           --outFileNamePrefix ${alignments}${sample%${read1_sufix}}. \
           --quantMode GeneCounts --genomeLoad NoSharedMemory --alignIntronMax ${max_intron} \
           --twopassMode Basic --outSAMattributes NH HI AS nM NM MD jM jI XS  --sjdbGTFfile ${gtf}
           
     """
}
    else {
    
    """
    
    mkdir -p /$star_output
    cd /raw_data
    
    echo ${pair_id}
      
    STAR --readFilesCommand zcat --runThreadN 12 --genomeDir ${star_index} --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM 15000000000 --readFilesIn ${renamed}${sample} \
         --outFileNamePrefix ${alignments}${sample%${read1_sufix}}. \
         --quantMode GeneCounts --genomeLoad NoSharedMemory --alignIntronMax ${max_intron} \
         --twopassMode Basic --outSAMattributes NH HI AS nM NM MD jM jI XS  --sjdbGTFfile ${gtf}
         
    """
  } 
  

}
