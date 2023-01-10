// requires vcf.gz format ( with index ) 
params.input_ref = "/lustre06/project/rrg-vmooser/shared/gnomad.genomes.v3.1.2.hgdp_tgp/filtered_pruned/chr*.vcf.gz"
params.input_study = "../VCFs-BQC19/chr*.vcf.gz"
params.path_to_laser = "~/scratch/LASER-2.04"

process intersect {
  debug true
  
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
  debug true
  
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
  debug true

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
  debug true
  memory "64GB"
  
  input:
  path "reference.vcf.gz"

  script:
  """
  ${params.path_to_laser}/vcf2geno/vcf2geno --inVcf reference.vcf.gz --out reference
  """
}

process convert_geno2 {
  debug true
  memory "64GB"

  input:
  path "study.vcf.gz"

  script:
  """
  ${params.path_to_laser}/vcf2geno/vcf2geno --inVcf study.vcf.gz --out study
  """

}


workflow {
ref_ch = Channel.fromPath(params.input_ref, checkIfExists:true) \
    | map { file -> [ file.name.toString().tokenize('.').get(0), file, file + ".tbi"] }

study_ch = Channel.fromPath(params.input_study, checkIfExists:true) \
    | map { file -> [ file.name.toString().tokenize('.').get(0), file, file + ".tbi"] }

input = ref_ch
   .join(study_ch)

vcfs = intersect(input)

reference_vcf = merge_ref(vcfs.reference_vcfs.collect())
study_vcf = merge_study(vcfs.study_vcfs.collect())

ref_geno = convert_geno(reference_vcf)
study_geno = convert_geno2(study_vcf)

}
