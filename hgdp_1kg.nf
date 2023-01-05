// requires vcf.gz format ( with index ) 
params.input_ref = "/lustre06/project/rrg-vmooser/shared/gnomad.genomes.v3.1.2.hgdp_tgp/filtered_pruned/chr*.vcf.gz"
params.input_study = "../VCFs-BQC19/chr*.vcf.gz"

process intersect {
  debug true
  
  input:
  tuple val(chr), path(ref), path(ref_tbi), path(study), path(study_tbi)
  
  output:
  path("./${chr}-isec/000*.vcf.gz*"),

  script:
  """
  bcftools isec -n=2 -p $chr-isec -Oz $ref $study
  """
} 

workflow {
ref_ch = Channel.fromPath(params.input_ref, checkIfExists:true) \
    | map { file -> [ file.name.toString().tokenize('.').get(0), file, file + ".tbi"] }

study_ch = Channel.fromPath(params.input_study, checkIfExists:true) \
    | map { file -> [ file.name.toString().tokenize('.').get(0), file, file + ".tbi"] }

input = ref_ch
   .join(study_ch)

intersect(input)
}
