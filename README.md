# HGDP_1KG Ancestry Inference

## Requirements
- Nextflow
- bcftools
- LASER (download from <http://csg.sph.umich.edu//chaolong/LASER/> )
  - Extract with `tar -xzvf LASER-2.04.tar.gz`
- Python virtual enviroment with pandas, numpy, sklearn


## To run:
Running the entire pipeline should take ~24h
### Input: Study and reference VCFs
- Input VCFs should be split by chromosome and named starting with "chr{#}."
- Input VCFs should be gzipped. 
- Input VCFs should be indexed with corresponding '.tbi' files in the same location. 
- Study and reference files can be in different locations. 
Input reference VCFs located at: `projects/rrg-vmooser/shared/HGDP_1KG/input_vcfs/`
### Input: Reference population labels
- Should be a CSV with an ID column and a population column ('genetic_region' for HGDP_TGP)
Reference population file located at: `projects/rrg-vmooser/shared/HGDP_1KG/hgdp_tgp_meta_genetic_region.csv`

### Parameters
Specify the following parameters in the nextflow.config file:
- `input_reference` - Path to reference VCFs (with index files)
- `input_study` - Path to study VCFs (with index files)
- `path_to_laser` - Path to unzipped LASER directory
- `nPCs` - Number of PCs to compute/project, then use for ancestry inference (optional)
- `reference_pop` - Path to reference population file
- `min_prob` - Minimum probability threshold for assigning population label w RF model 
- `seed` - Random seed to use for random forest model (optional)
- `n_jobs` - Number of jobs to use for parallelization of projection step


## Optional support 

#### {GRAPH}_Plot_PCA_Ancestry.py
options: 

| Short | Long        | Type    | Required | Default | Description                                        |
|-------|-------------|---------|----------|---------|----------------------------------------------------|
| `-P`  | `--Projected`  | `str`   | Yes      |         | Projected samples position                         |
| `-R`  | `--Reference` | `str`   | Yes      |         | Reference samples position                         |
| `-S`  | `--Study` | `str`   | Yes      |         | Infered Ethnicity projected samples               |
| `-A`  | `--Ancestry`  | `str`   | Yes      |         | Reference Samples ethnicity                        |
| `-c`  | `--selected` | `str`   | Yes      |         | Ethnicty to plot                                   |
| `-T`  | `--Threshold` | `float` | Yes      |         | Reference Samples ethnicity threshold              |
| `-n`  | `--PC       | `number`| `n`         | `int`   | Yes      |         | number of PC to plot                               |
| `-l`  | `--label` | `str`   | No       | `Study` | label of the study                                 |
|       | `--out`    | `str`   | No       | `output`| output                                             |
