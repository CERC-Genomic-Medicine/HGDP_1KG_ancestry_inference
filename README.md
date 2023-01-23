# HGDP_1KG Ancestry Inference

## Requirements
- Nextflow
- bcftools
- LASER (download from <http://csg.sph.umich.edu//chaolong/LASER/> )
  *Extract the .tar.gz file with `tar -xzvf LASER-2.04.tar.gz`*


## Input
### Study and reference VCFs
- Input VCFs should be split by chromosome and named starting with "chr{#}."
- Input VCFs should be gzipped. 
- Inout VCFs should be indxed with corresponding '.tbi' files in the same location. 
- Study and reference files can be in different locations. 

### Parameters
Specify the following parameters in the nextflow.config file:
- `input_reference` - Path to reference VCFs (with index files)
- `input_study` - Path to study VCFs (with index files)
- `path_to_laser` - Path to unzipped LASER directory
