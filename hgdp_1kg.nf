// requires vcf.gz format
params.input_ref = "/lustre06/project/rrg-vmooser/shared/gnomad.genomes.v3.1.2.hgdp_tgp/filtered_pruned/chr*.vcf.gz"
params.input_study = "/lustre06/project/rrg-vmooser/shared/BQC19_HLA/VCFs/*.vcf.gz"

process intersect {
  debug true
  input:
  tuple val(chr), path(ref), path(study)

  script:
  """
  bcftools isec -n=2 -p isec -Oz -o $chr-isec $ref $study
  """
} 

workflow {
ref_ch = Channel.fromPath(params.input_ref, checkIfExists:true) \
    | map { file -> [ file.name.toString().tokenize('.').get(0), file] }

study_ch = Channel.fromPath(params.input_study, checkIfExists:true) \
    | map { file -> [ file.name.toString().tokenize('.').get(0), file] }

input = ref_ch
   .join(study_ch)

intersect(input)
}
