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
    ( ! file("${params.star_index}index.fa".exists() ) 
  
    script:
    
        """
        
        if [[ ! -e \$star_index ]] ; then mkdir -p \$star_index ; fi

        cd \$star_index
    
        STAR --runMode genomeGenerate \
             --genomeDir ${star_index} \ 
             --genomeFastaFiles ${params.genomes}${params.organism}/${params.release}/.toplevel.fa \ #reference genome FASTA file
             --sjdbGTFfile ../genomes/${params.organism}/${params.release}/original.gtf \
             --runThreadN 12
        
        """
    }


