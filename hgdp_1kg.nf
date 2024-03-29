/*
* Author: Peyton McClelland <peyton.mcclelland@mail.mcgill.ca>
* Version: 1.0
* Year: 2023
*/

params.nPCs = 20
params.min_prob = 0
params.seed = 11

process intersect {
  input:
  tuple val(chr), path(ref), path(ref_tbi), path(study), path(study_tbi)
  
  output:
  path("${chr}.ref.vcf.gz"), emit: reference_vcfs
  path("${chr}.study.vcf.gz"), emit: study_vcfs

  script:
  """
  bcftools isec -n=2 -p temp -Oz $ref $study

  mv temp/0000.vcf.gz ${chr}.ref.vcf.gz
  mv temp/0001.vcf.gz ${chr}.study.vcf.gz
  """
} 

process merge_ref {
  input:
  path "chr*.ref.vcf.gz" 

  output:
  path "reference.vcf.gz"

  script:
  """
  bcftools concat chr{1..22}.ref.vcf.gz -Oz -o reference.vcf.gz
  """
} 

process merge_study {
  input:
  path "chr*.study.vcf.gz"

  output:
  path "study.vcf.gz"

  script:
  """
  bcftools concat chr{1..22}.study.vcf.gz -Oz -o study.vcf.gz
  """
} 

process convert_geno {
  memory "64GB"
  time "1h"
  
  input:
  path "reference.vcf.gz"
  
  output:
  tuple path("reference.geno"), path("reference.site")

  publishDir "output/geno/", mode: "copy"

  script:
  """
  ${params.path_to_laser}/vcf2geno/vcf2geno --inVcf reference.vcf.gz --out reference
  """
}

process convert_geno2 {
  memory "64GB"
  time "2h"

  input:
  path "study.vcf.gz"
  
  output:
  tuple path("study.geno"), path("study.site")

  publishDir "output/geno/", mode: "copy"

  script:
  """
  ${params.path_to_laser}/vcf2geno/vcf2geno --inVcf study.vcf.gz --out study
  """

}

process reference_PCA {
  memory "32GB"
  time "1h"

  input:
  tuple path("reference.geno"), path("reference.site")

  output:
  path "reference.RefPC.coord"

  publishDir "output/", mode: "copy"

  script:
  """
  ${params.path_to_laser}/laser -g reference.geno -k ${params.nPCs} -pca 1 -o reference
  """

}

process project {
  memory "32GB"
  time "2d"

  input:
  path "reference.RefPC.coord"
  tuple path("reference.geno"), path("reference.site")
  tuple path("study.geno"), path("study.site")
  val job


  output:
  path("projection_job_${job}.ProPC.coord"), emit: parallel_study_proPC_coord

  """
  n_samples=\$(cat study.geno | wc -l )
  step=\$((\$n_samples/${params.n_jobs}))
  first=\$(( ((${job}-1)*\$step) + 1 ))
  last=\$((${job}*\$step))
  if (( ${job} == ${params.n_jobs} ))
  then
  	last=\$n_samples
  fi
  echo "STUDY_FILE study.geno" > trace.conf
  echo "GENO_FILE reference.geno" >> trace.conf
  echo "COORD_FILE reference.RefPC.coord" >> trace.conf
  echo "FIRST_IND \$first" >> trace.conf
  echo "LAST_IND \$last" >> trace.conf
  echo "OUT_PREFIX projection_job_${job}" >> trace.conf
  echo "DIM ${params.nPCs}" >> trace.conf

  ${params.path_to_laser}/trace -p trace.conf

  """
  
}

process project_merge {

 input:
 path "projection_job_*.ProPC.coord"

 output:
 path "trace.ProPC.coord"

 publishDir "output/", mode: "copy"

 script:
 """
 head -n 1 projection_job_1.ProPC.coord >> trace.ProPC.coord
 for i in projection_job_*.ProPC.coord
 do
 tail -n +2 \$i >> trace.ProPC.coord
 done
 """
}


process infer_ancestry {
  debug true 

  input:
  path "reference.RefPC.coord"
  path "trace.ProPC.coord"

  output:
  path "predicted_ancestry.txt"
  stdout emit: rf_log
  
  publishDir "output/", mode: "copy"

  script:
  """
  rf_model.py -r reference.RefPC.coord -s trace.ProPC.coord -k ${params.nPCs} -p ${params.reference_pop} --min_prob ${params.min_prob} --seed ${params.seed}
  """  

}

workflow {
ref_ch = Channel.fromPath(params.input_reference, checkIfExists:true) \
    | map { file -> [ file.name.toString().tokenize('.').get(0), file, file + ".tbi"] }

study_ch = Channel.fromPath(params.input_study, checkIfExists:true) \
    | map { file -> [ file.name.toString().tokenize('.').get(0), file, file + ".tbi"] }

input = ref_ch
   .join(study_ch)

vcfs = intersect(input)

reference_vcf = merge_ref(vcfs.reference_vcfs.collect())
study_vcf = merge_study(vcfs.study_vcfs.collect())

parallelization = Channel.from( 1..params.n_jobs )

reference_geno = convert_geno(reference_vcf)
study_geno = convert_geno2(study_vcf)

reference_PC_coord = reference_PCA(reference_geno)

parallel_study_proPC_coord = project(reference_PC_coord, reference_geno, study_geno, parallelization)

study_proPC_coord=project_merge(parallel_study_proPC_coord.collect())

infer_ancestry(reference_PC_coord, study_proPC_coord)

}
