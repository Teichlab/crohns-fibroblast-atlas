## Activate conda environment before running the script

## Example usage
#conda activate create_cistarget_databases
#/software/team205/bin/jsub lsf -n cistarget_db_stim -q long -c 4 -m 100g "bash 002_create_cistarget_db.sh" | bsub -G teichlab

## Location of all peaks to be used as background for cistarget analyses
peaks_fasta='/lustre/scratch126/cellgen/team205/is10/fibroblasts/bulk_stim/cistarget/all_peaks.fa'
## Directory of cistarget database to be created
out_dir='/lustre/scratch126/cellgen/team205/is10/fibroblasts/bulk_stim/cistarget/cistarget_db'

## Set ${create_cistarget_databases_dir} variable to path where the repo was cloned to.
create_cistarget_databases_dir='/nfs/team205/is10/scripts/create_cisTarget_databases'
motif_collection='/nfs/team205/is10/resources/aerts_motifs/v10nr_clust_public'

## Create database
python ${create_cistarget_databases_dir}/create_cistarget_motif_databases.py \
-f $peaks_fasta \
-M ${motif_collection}/singletons/ \
-m ${motif_collection}/motifs.txt \
-o $out_dir \
-s 42 \
-t 14

## Create rankings
${create_cistarget_databases_dir}/convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py \
-i $out_dir.motifs_vs_regions.scores.feather \
-s 42