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
            
            if [[ ! -f rnaseq.python-3.8-1.sif ]] ;
              then
                singularity pull rnaseq.python-3.8-1.sif docker://index.docker.io/mpgagebioinformatics/rnaseq.python:3.8-1
            fi

        fi


        if [[ "${params.containers}" == "docker" ]] ; 

          then

            docker pull mpgagebioinformatics/star:2.7.11b
            docker pull mpgagebioinformatics/samtools:1.16.1
            docker pull mpgagebioinformatics/rnaseq.python:3.8-1
        fi

        """
    }


process rename_sample {
    stageInMode 'symlink'
    stageOutMode 'move'

    script: 

      """
      #!/usr/local/bin/python3
      
      import os
      import pandas as pd
      import openpyxl
      
      samplesheet = pd.read_excel("/workdir/nf-star-test/sample_sheet.xlsx",engine="openpyxl")
      samplesheet['sample_name'] = ''
      samplesheet
      groups = dict()
      sample_name_colum = samplesheet.columns.get_loc("sample_name")

      for index, row in samplesheet.iterrows():
          if not row[1] in groups:
              groups[row[1]] = 1
          else: 
              groups[row[1]] += 1     
          print('ln -s %s%s %s%s_%s${params.read1_sufix}' %('${params.raw_data}', row[0], '${params.raw_renamed}', row[1], groups[row[1]]))
          os.system('ln -s %s%s %s%s_%s${params.read1_sufix}' %('${params.raw_data}', row[0], '${params.raw_renamed}', row[1], groups[row[1]]))
          row[sample_name_colum] = '%s_%s' %(row[1], groups[row[1]])
          if "${params.read2_sufix}" in row[0]:
              print('ln -s ${params.raw_data}%s ${params.raw_renamed}%s_%s${params.read2_sufix}' % (row[0].replace('${params.read1_sufix}', '${params.read2_sufix}'), row[1], groups[row[1]]))
              os.system('ln -s ${params.raw_data}%s ${params.raw_renamed}%s_%s${params.read2_sufix}' % (row[0].replace('${params.read1_sufix}', '${params.read2_sufix}'), row[1], groups[row[1]]))
      
      """
}


process star_indexer {
    stageInMode 'symlink'
    stageOutMode 'move'

    when:
    ( ! file("${params.star_index}index.fa").exists() ) 
  
    script:
    
        """
        
        if [[ ! -e ${params.star_index} ]] ; then mkdir -p ${params.star_index} ; fi

        cd ${params.star_index}
    
        STAR --runMode genomeGenerate \
             --genomeDir ${params.star_index} \
             --genomeFastaFiles "${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.fa" \
             --sjdbGTFfile "${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.gtf" \
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
    ( ! file("${params.project_folder}/star_output/${pair_id}.Aligned.sortedByCoord.out.bam").exists() ) 
   
    script:
    
    def single = fastq instanceof Path
    
    if ( single ) {

        """
        mkdir -p ${params.star_output}
        cd ${params.raw_renamed}

        echo ${pair_id}

        if [ ! -e ${params.star_output}${pair_id}.Aligned.sortedByCoord.out.bam ] ; then

            STAR --readFilesCommand zcat \
                --runThreadN 12 \
                --genomeDir ${params.star_index} \
                --outSAMtype BAM SortedByCoordinate \
                --limitBAMsortRAM 15000000000 \
                --readFilesIn ${params.raw_renamed}${pair_id}.READ_1.fastq.gz \
                --outFileNamePrefix ${params.star_output}${pair_id}. \
                --quantMode GeneCounts \
                --genomeLoad NoSharedMemory \
                --alignIntronMax ${params.max_intron} \
                --twopassMode Basic \
                --outSAMattributes NH HI AS nM NM MD jM jI XS \
                --sjdbGTFfile ${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.gtf
        fi
        """
    }
    else {
          
        """
        mkdir -p ${params.star_output}
        cd ${params.raw_renamed}

        echo ${pair_id}

        echo "pass paired"

        if [ ! -e ${params.star_output}${pair_id}.Aligned.sortedByCoord.out.bam ] ; then

            STAR --readFilesCommand zcat \
                --runThreadN 12 \
                --genomeDir ${params.star_index} \
                --outSAMtype BAM SortedByCoordinate \
                --limitBAMsortRAM 15000000000 \
                --readFilesIn ${params.raw_renamed}${pair_id}.READ_1.fastq.gz ${params.raw_renamed}${pair_id}.READ_2.fastq.gz \
                --outFileNamePrefix ${params.star_output}${pair_id}. \
                --quantMode GeneCounts \
                --genomeLoad NoSharedMemory \
                --alignIntronMax ${params.max_intron} \
                --twopassMode Basic \
                --outSAMattributes NH HI AS nM NM MD jM jI XS \
                --sjdbGTFfile ${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.gtf
        fi

        """
    }

} 

process samtools_index {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val pair_id
    tuple val(pair_id), path(fastq)

  when:
    ( ! file("${params.sajr_output}asplicing.merged.sorted.bam").exists() ) 


  script:
    """
    cd ${params.star_output}

    if [ ! -e ${params.star_output}${pair_id}.Aligned.sortedByCoord.out.bam ] ; then

    samtools index ${params.star_output}${pair_id}.Aligned.sortedByCoord.out.bam
    
    fi

    if [ ! -e ${params.sajr_output}asplicing.merged.sorted.bam ] ; then

    mkdir -p ${params.sajr_output}
  
    samtools merge -f -@ 10 -O BAM ${params.sajr_output}asplicing.merged.bam *.bam
    
    cd ${params.sajr_output}

    samtools sort -@ 10 -o asplicing.merged.sorted.bam asplicing.merged.bam
    samtools index asplicing.merged.sorted.bam

    fi
    """
}


workflow images {

  main:
  get_images()
}

workflow rename {
  rename_sample()
}


workflow index {
  main:
  star_indexer()
}


workflow map_reads {
  main:
    read_files=Channel.fromFilePairs( "${params.raw_renamed}/*.READ_{1,2}.fastq.gz", size: -1 )
    star_mapping( read_files )
    samtools_index( star_mapping.out.collect(), read_files )

}
