#!/usr/bin/env bash
#SBATCH --export ALL
#SBATCH --partition condo
#SBATCH --account ddp268
#SBATCH --qos condo
#SBATCH --job-name hprc_subset_minigraph_align
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 64
#SBATCH --mem 25G
#SBATCH --time 8:00:00
#SBATCH --error results/2024-09-12_subset_minigraph_align/slurm_logs/%j.err
#SBATCH --output results/2024-09-12_subset_minigraph_align/slurm_logs/%j.out
#SBATCH --array=1-95%16 # Run for indicies 1-95 inclusive and only let 16 jobs run at a time

set -e
function TIMESTAMP () { echo $(TZ=US/Pacific; date '+%Y-%m-%d_%H-%M-%S'); }

#################
### Constants ###
#################

# Inputs
RESULTS_DIR=results/2024-09-12_subset_minigraph_align
PANGENOME_GFA=data/hprc-v1.1-mc-grch38.gfa
PANGENOME_GFAB=data/hprc-v1.1-mc-grch38.gfab
ASSEMBLIES_AGC=data/HPRC-yr1.agc
REFERENCE_HG38=data/hg38.fa
RANGE=${CHROM}:${START}-${END}
RANGE_IN_PATH=${CHROM}_${START}-${END}

# Job Constants
JOB_INDEX=${SLURM_ARRAY_TASK_ID}
HAP_ID=$( agc listset ${ASSEMBLIES_AGC} | head -n ${JOB_INDEX} | tail -n 1 )

# Temporary File Location
TEMP_DIR=/var/tmp/wwford
TEMP_CLEAN_FASTA=${TEMP_DIR}/clean_${HAP_ID}.fa

ALIGN_DIR=${RESULTS_DIR}/hg38_alignments
ALIGN_BAM=${ALIGN_DIR}/hg38_align_${HAP_ID}.bam
SUBSET_RANGE_DIR=${RESULTS_DIR}/subsets/${RANGE_IN_PATH}
SUBSET_FASTA=${SUBSET_RANGE_DIR}/subset_${HAP_ID}_${JOB_INDEX}.fa
SUBSET_GFA=${SUBSET_RANGE_DIR}_mc.gfa

# Outputs
SLURM_LOG_FOLDER=${RESULTS_DIR}/slurm_logs
LOG_FOLDER=${RESULTS_DIR}/logs
OUTPUT_LOG=${LOG_FOLDER}/minigraph_${RANGE_IN_PATH}_${HAP_ID}_$(TIMESTAMP).log

OUTPUT_DIR=${RESULTS_DIR}/${RANGE_IN_PATH}
OUTPUT_GAF=${OUTPUT_DIR}/minigraph_${RANGE_IN_PATH}_${HAP_ID}.gaf

#################
### Functions ###
#################

# Exit if output exists 
output_exists () {
    if [ -e ${OUTPUT_GAF} ]; then
        echo "${OUTPUT_GAF} already exists. Exiting w/ code 1" >> ${OUTPUT_LOG}
        exit 1
    fi
}

# Make log, temp, and output dirs if they do not exist
make_dir_if () {
    if [ ! -d $1 ]; then
        mkdir -p $1 ;
    fi
}
make_dirs () {
    make_dir_if ${SLURM_LOG_FOLDER}
    make_dir_if ${LOG_FOLDER}
    make_dir_if ${TEMP_DIR}
    make_dir_if ${OUTPUT_DIR}
    make_dir_if ${SUBSET_RANGE_DIR}
    make_dir_if ${ALIGN_DIR}
}

# Delete likely contaminated contigs from HPRC assemblies
BAD_CONTIGS_TXT=data/y1_genbank_remaining_potential_contamination.txt
BAD_CONTIGS_REGEX=$(cut -f2 ${BAD_CONTIGS_TXT} | tail -n +2 \
    | sed ';s/$/|/' | tr -d '\n' | sed 's/|$//')
delete_contaminated_contigs () {
    awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } 
        else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } 
    }' \
    | grep -Ev ${BAD_CONTIGS_REGEX} \
    | tr "\t" "\n"
}

# Make clean fasta if it does not exist
clean_fasta () {
    if [ -e ${TEMP_CLEAN_FASTA} ]; then
        echo "${TEMP_CLEAN_FASTA} already exists" >> ${OUTPUT_LOG}
    else
        agc getset data/HPRC-yr1.agc ${HAP_ID} | delete_contaminated_contigs > ${TEMP_CLEAN_FASTA}
    fi
}

# align fasta to grch38
align_fasta () {
    if [ -e ${ALIGN_BAM} ]; then
        echo "${ALIGN_BAM} already exists" >> ${OUTPUT_LOG}
    else
        minimap2 -ax asm5 ${REFERENCE_HG38} ${TEMP_CLEAN_FASTA} -t ${SLURM_CPUS_PER_TASK} 2>> ${OUTPUT_LOG} | samtools sort -o ${ALIGN_BAM}
        samtools index ${ALIGN_BAM}
    fi
}

# subset fasta based on grch38 coords
subset_fasta () {
    samtools view ${ALIGN_BAM} ${RANGE} | awk -v id_col=1 -v seq_col=10 '{print ">"$id_col"\n"$seq_col}' > ${SUBSET_FASTA}
}

subset_graph () {
    if [ -e ${SUBSET_GFA} ]; then
        echo "${SUBSET_GFA} already exists" >> ${OUTPUT_LOG}
    else
        gfabase sub ${PANGENOME_GFAB} --range "GRCh38#${RANGE}" --view -o ${SUBSET_GFA}
    fi
}

# align subset fasta to subset gfa
align_to_graph () {
    minigraph \
        -x asm \
        ${SUBSET_GFA} \
        ${SUBSET_FASTA} \
        -t ${SLURM_CPUS_PER_TASK} > \
            ${OUTPUT_GAF} 2>> \
            ${OUTPUT_LOG}
}

#######################
### Main Entrypoint ###
#######################

make_dirs

echo "Job ID                ${SLURM_JOB_ID}" > ${OUTPUT_LOG}
echo "SLURM_ARRAY_JOB_ID    ${SLURM_ARRAY_JOB_ID}" >> ${OUTPUT_LOG}
echo "SLURM_ARRAY_TASK_ID   ${SLURM_ARRAY_TASK_ID}" >> ${OUTPUT_LOG}

# Check if output already exists
output_exists

echo "Saving clean temp fasta       $(TIMESTAMP)" >> ${OUTPUT_LOG}
clean_fasta
echo "Saved Fasta                   $(TIMESTAMP)" >> ${OUTPUT_LOG}

echo "Aligning fasta to grch38      $(TIMESTAMP)" >> ${OUTPUT_LOG}
align_fasta
echo "Aligned fasta to grch38       $(TIMESTAMP)" >> ${OUTPUT_LOG}

echo "Subsetting fasta              $(TIMESTAMP)" >> ${OUTPUT_LOG}
subset_fasta
echo "Subset fasta                  $(TIMESTAMP)" >> ${OUTPUT_LOG}

echo "Subsetting graph              $(TIMESTAMP)" >> ${OUTPUT_LOG}
subset_graph
echo "Subset graph                  $(TIMESTAMP)" >> ${OUTPUT_LOG}

echo "Start Align                   $(TIMESTAMP)" >> ${OUTPUT_LOG}
align_to_graph
echo "Done Align                    $(TIMESTAMP)" >> ${OUTPUT_LOG}
echo "Results at ${OUTPUT_GAF}"

# Delete temporary fasta
rm ${TEMP_CLEAN_FASTA}
echo "Deleted Temp Fasta            $(TIMESTAMP)" >> ${OUTPUT_LOG}
