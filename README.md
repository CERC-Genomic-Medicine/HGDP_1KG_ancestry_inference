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

## Preparing the HGDP + 1KG callset to use as input data
To use the HGDP + 1KG callset (gnomAD v3.1.2) as your reference data, download [the callset and sample metadata here](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg).
### Some pre-processing may be required:
- Filter per chromosome on preferred quality metrics.
  - E.g. `bcftools view -fPASS -q0.05 -Q0.95 -i 'F_MISSING < 0.001' <vcf> |  bcftools annotate -x INFO,^GT - -Oz -o <filtered_vcf>`
  - Additional QC filters from gnomAD sample metadata (also excludes related individuals) `bcftools view -S HGDP_1KG.PassedQC.id`
- Perform LD pruning per chromosome.
  - E.g. `plink --indep-pairwise 1000 100 0.9 --vcf <filtered_vcf> --double-id --id-delim , --out chr${i}.plink`
  - `plink --vcf <filtered_vcf> --extract chr${i}.plink.prune.in --recode vcf-iid --out <filtered_pruned_vcf>`
- _IMPORTANT:_ Rename gnomad files to begin with `chr#`. If needed, rename chromosomes to align with study data.
  - E.g. `bcftools annotate --rename-chrs rename_chrs.txt gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}.filtered_pruned.vcf.gz -Oz -o chr${i}.gnomad.genomes.v3.1.2.hgdp_tgp.filtered_pruned.vcf.gz`
  - An example `rename_chrs` file will look something like this:
    ```
    1 chr1
    2 chr2
    3 chr3
    ...
    ```
- Index all files.
  - `bcftools index -t chr${i}.gnomad.genomes.v3.1.2.hgdp_tgp.filtered_pruned.vcf.gz`
 ### Format the reference population data:
 - I.e. from the [sample metadata available here](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg), extract the population descriptors (we used the harmonized `genetic_region` label)
 - Format of the population label file should be a CSV (comma-separated), with two columns (Sample ID and Population label).
   - The pipeline currently assumes the file will look like this:
     ```
     ID,genetic_region
     HG00000,EUR
     ...
     ```
   - But can easily be updated to accomodate other ID and population column labels. 

## Optional support 

#### Script: `plot_PCA.py` to plot resulting TRACE study PCs against reference data.
This is a script to take the output of the ancestry inference pipeline and:
- Subset your study samples to the predicted ancestry group of interested (from `predicted_ancestry.txt` at your desired threshold of probability).
- Take the projected PC coordinates for these samples (from `trace.ProPC.coord`) and the reference PC coordinates against which they are projected (`reference.RefPC.coord`)
- Plot the top N PCs for your predicted ancestry data and reference data (with population labels).

**Parameters to specify:**
- `-P`,`--Projected` - Projected study data PCA coordinates (*output/trace.ProPC.coord*) 
- `-R`,`--Reference` - Reference data PCA coordinates (*output/reference.RefPC.coord*)
- `-S`,`--Study` - Predicted ancestry output for study samples (*output/predicted_ancestry.txt*)
- `-A`,`--Ancestry` - Population labels for reference data
- `-c`,`--selected` - Predicted ancestry group in study data to subset and plot
- `-T`,`--Threshold` - Threshold of probability to use when subsetting by predicted ancestry (e.g. 0.5)
- `-n`,`--PC` - Number of PCs to plot
- `-l`,`--label` - Study label to use when titling plot and legend
- `--out` - Prefix of output image file name (output.png)

**Example command:**
` python plot_PCA.py -P output/trace.ProPC.coord -R output/reference.RefPC.coord -S output/predicted_ancestry.txt -A genetic_region.csv -c CSA -T 0.80 -n4 -l STUDY --out figure `
