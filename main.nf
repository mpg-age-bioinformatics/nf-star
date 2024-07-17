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

            if [[ ! -f irfinder-1.3.1.sif ]] ;
              then
                singularity pull irfinder-1.3.1.sif docker://index.docker.io/mpgagebioinformatics/irfinder:1.3.1
            fi

        fi


        if [[ "${params.containers}" == "docker" ]] ; 

          then

            docker pull mpgagebioinformatics/star:2.7.11b
            docker pull mpgagebioinformatics/samtools:1.16.1
            docker pull mpgagebioinformatics/rnaseq.python:3.8-1
            docker pull mpgagebioinformatics/irfinder:1.3.1.sif
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
      
      samplesheet = pd.read_excel("${params.sample_sheet_xlsx}",engine="openpyxl")
      samplesheet['sample_name'] = ''
      samplesheet
      groups = dict()
      sample_name_colum = samplesheet.columns.get_loc("sample_name")

      for index, row in samplesheet.iterrows():
          if not row[1] in groups:
              groups[row[1]] = 1
          else: 
              groups[row[1]] += 1     
          print('ln -s %s%s %s%s_%s${params.read1_sufix}' %('${params.fastqc_raw_data}', row[0], '${params.raw_renamed}', row[1], groups[row[1]]))
          os.system('ln -s %s%s %s%s_%s${params.read1_sufix}' %('${params.fastqc_raw_data}', row[0], '${params.raw_renamed}', row[1], groups[row[1]]))
          row[sample_name_colum] = '%s_%s' %(row[1], groups[row[1]])
          if "${params.read2_sufix}" in row[0]:
              print('ln -s ${params.fastqc_raw_data}%s ${params.raw_renamed}%s_%s${params.read2_sufix}' % (row[0].replace('${params.read1_sufix}', '${params.read2_sufix}'), row[1], groups[row[1]]))
              os.system('ln -s ${params.fastqc_raw_data}%s ${params.raw_renamed}%s_%s${params.read2_sufix}' % (row[0].replace('${params.read1_sufix}', '${params.read2_sufix}'), row[1], groups[row[1]]))
      
      """
}


process star_indexer {
    stageInMode 'symlink'
    stageOutMode 'move'

    when:
    ( ! file("${params.star_index}index.fa").exists() ) 
  
    script:
    
        """
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
    ( ! file("${params.star_out}${pair_id}.Aligned.sortedByCoord.out.bam").exists() ) 
   
    script:
    
    def single = fastq instanceof Path
    
    if ( single ) {

        """
        cd ${params.raw_renamed}

        echo ${pair_id}

        if [ ! -e ${params.star_out}${pair_id}.Aligned.sortedByCoord.out.bam ] ; then

            STAR --readFilesCommand zcat \
                --runThreadN 12 \
                --genomeDir ${params.star_index} \
                --outSAMtype BAM SortedByCoordinate \
                --limitBAMsortRAM 15000000000 \
                --readFilesIn ${params.raw_renamed}${pair_id}.READ_1.fastq.gz \
                --outFileNamePrefix ${params.star_out}${pair_id}. \
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
        cd ${params.raw_renamed}

        echo ${pair_id}

        echo "pass paired"

        if [ ! -e ${params.star_out}${pair_id}.Aligned.sortedByCoord.out.bam ] ; then

            STAR --readFilesCommand zcat \
                --runThreadN 12 \
                --genomeDir ${params.star_index} \
                --outSAMtype BAM SortedByCoordinate \
                --limitBAMsortRAM 15000000000 \
                --readFilesIn ${params.raw_renamed}${pair_id}.READ_1.fastq.gz ${params.raw_renamed}${pair_id}.READ_2.fastq.gz \
                --outFileNamePrefix ${params.star_out}${pair_id}. \
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

  script:
    """
    cd ${params.star_out}

    if [ ! -e ${params.star_out}${pair_id}.Aligned.sortedByCoord.out.bam.bai ]; then
    samtools index ${params.star_out}${pair_id}.Aligned.sortedByCoord.out.bam
    fi

    if [ ! -e ${params.star_out}${pair_id}.Aligned.sortedByCoord.out.bed ] ; then
    samtools sort -@ ${task.cpus} ${params.star_out}${pair_id}.Aligned.sortedByCoord.out.bam | genomeCoverageBed -bga -split -ibam - > ${params.star_out}${pair_id}.Aligned.sortedByCoord.out.bed
    fi

    """
}

process samtools_merge {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.sajr_output}${params.series}.merged.sorted.bam.bai").exists() ) 


  script:
    """
    cd ${params.star_out}

    if [ ! -e ${params.sajr_output}${params.series}.merged.sorted.bam ] ; then
    samtools merge -f -@ ${task.cpus} -O BAM ${params.sajr_output}${params.series}.merged.bam *.bam
    fi
    
    cd ${params.sajr_output}

    if [ ! -e ${params.series}.merged.sorted.bam ] ; then 
    samtools sort -@ ${task.cpus} -o ${params.series}.merged.sorted.bam ${params.series}.merged.bam
    fi

    if [ ! -e ${params.series}.merged.sorted.bam.bai ]; then
    samtools index ${params.series}.merged.sorted.bam
    fi

    """
}


workflow images {

  main:
  get_images()
}

workflow rename {
  if ( ! file("${params.raw_renamed}").isDirectory() ) {
        file("${params.raw_renamed}").mkdirs()
      }
  rename_sample()
}

workflow index {
  main:
  if ( ! file("${params.star_index}").isDirectory() ) {
        file("${params.star_index}").mkdirs()
      }
  star_indexer()
}


workflow map_reads {

  if ( ! file("${params.star_out}").isDirectory() ) {
        file("${params.star_out}").mkdirs()
      }

  if ( ! file("${params.sajr_output}").isDirectory() ) {
        file("${params.sajr_output}").mkdirs()
      }

  read_files=Channel.fromFilePairs( "${params.raw_renamed}/*.READ_{1,2}.fastq.gz", size: -1 )
  star_mapping( read_files )

  samtools_index(star_mapping.out.collect(), read_files)
}

workflow merging {
  if ( ! file("${params.sajr_output}").isDirectory() ) {
        file("${params.sajr_output}").mkdirs()
      }

  samtools_merge()
}