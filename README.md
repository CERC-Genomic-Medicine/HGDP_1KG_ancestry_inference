# HGDP_1KG Ancestry Inference

## Requirements
- Nextflow
- bcftools
- LASER (download from <http://csg.sph.umich.edu//chaolong/LASER/> )
  *Extract the .tar.gz file with `tar -xzvf LASER-2.04.tar.gz`*
- Python virtual enviroment with pandas, numpy, sklearn


## Input
### Study and reference VCFs
- Input VCFs should be split by chromosome and named starting with "chr{#}."
- Input VCFs should be gzipped. 
- Inout VCFs should be indxed with corresponding '.tbi' files in the same location. 
- Study and reference files can be in different locations. 
### Reference population labels
- Should be a CSV with an ID column and a population column ('genetic_region' for HGDP_TGP)


### Parameters
Specify the following parameters in the nextflow.config file:
- `input_reference` - Path to reference VCFs (with index files)
- `input_study` - Path to study VCFs (with index files)
- `path_to_laser` - Path to unzipped LASER directory
- `nPCs` - Number of PCs to compute/project, then use for ancestry inference (optional)
- `reference_pop` - Path to reference population file
- `min_prob` - Minimum probability threshold for assigning population label w RF model 
- `seed` - Random seed to use for random forest model (optional)
