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

    // input:
    //   val samplesheet

    script: 

      """
      #!/usr/bin/python3
      
      import os
      import pandas as pd
      import openpyxl
      
      samplesheet = pd.read_excel("/workdir/nf-star-test/sample_sheet.xlsx",engine="openpyxl")
      sample_sheet['sample_name'] = ''
      sample_sheet
      groups = dict()
      sample_name_colum = sample_sheet.columns.get_loc("sample_name")

      for index, row in sample_sheet.iterrows():
          if not row[1] in groups:
              groups[row[1]] = 1
          else: 
              groups[row[1]] += 1     
          print('ln -s %s/%s %s/%s_%s${params.read1_sufix}' %('${params.raw_data}', row[0], '${params.raw_renamed}', row[1], groups[row[1]]))
          os.system('ln -s %s/%s %s/%s_%s${params.read1_sufix}' %('${params.raw_data}', row[0], '${params.raw_renamed}', row[1], groups[row[1]]))
          row[sample_name_colum] = '%s_%s' %(row[1], groups[row[1]])
          if "${params.read2_sufix}" in row[0]:
              print('ln -s ${params.raw_data}/%s ${params.raw_renamed}/%s_%s${params.read2_sufix}' % (row[0].replace('${params.read1_sufix}', '${params.read2_sufix}'), row[1], groups[row[1]]))
              os.system('ln -s ${params.raw_data}/%s ${params.raw_renamed}/%s_%s${params.read2_sufix}' % (row[0].replace('${params.read1_sufix}', '${params.read2_sufix}'), row[1], groups[row[1]]))
      
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
             --genomeFastaFiles "${params.genomes}/${params.organism}/${params.release}/${params.organism}.${params.release}.fa" \
             --sjdbGTFfile "${params.genomes}/${params.organism}/${params.release}/${params.organism}.${params.release}.gtf" \
             --runThreadN 12
        
        """
    }






workflow images {

  main:
  get_images()
}

workflow rename {
  // if ( 'sample_sheet_xlsx' in params.keySet() ) {
  //   samplesheet="${params.sample_sheet_xlsx}"
  // } else {
  //   samplesheet=""
  // }
  // if ( ! file("${params.raw_renamed}").isDirectory() ) {
  //   file("${params.raw_renamed}").mkdirs()
  // }
  // rename_sample(samplesheet)
  rename_sample()
}


workflow index {
  main:
  star_indexer()
}

