#!/usr/bin/env bash
#SBATCH --export ALL
#SBATCH --partition condo
#SBATCH --account ddp268
#SBATCH --qos condo
#SBATCH --job-name hprc_minigraph_align
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 64
#SBATCH --mem 25G
#SBATCH --time 8:00:00
#SBATCH --output results/2024-08-30_HPRC_minigraph_alignments/slurm_logs/%j.out-%N
#SBATCH --output results/2024-08-30_HPRC_minigraph_alignments/slurm_logs/%j.err-%N
#SBATCH --array=1-2%4 # Run for indicies 1-95 inclusive and only let 4 jobs run at a time

# This slurm file iterates over all HPRC assemblies and aligns them to the pangenome.
# 1. Pull fasta from agc by index
# 2. Delete any contaminated contigs if they exist and save to a temporary fasta file
# 3. Align fasta to pangenome and save output

set -e

# Inputs
PANGENOME_GFA=data/hprc-v1.1-mc-grch38.gfa
ASSEMBLIES_AGC=data/HPRC-yr1.agc

BAD_CONTIGS_TXT=data/y1_genbank_remaining_potential_contamination.txt
BAD_CONTIGS_REGEX=$(cut -f2 ${BAD_CONTIGS_TXT} | tail -n +2 \
    | sed ';s/$/|/' | tr -d '\n' | sed 's/|$//')

# Job Constants
JOB_INDEX=${SLURM_ARRAY_TASK_ID}
HAP_ID=$( agc listset ${ASSEMBLIES_AGC} | head -n ${JOB_INDEX} | tail -n 1 )

# Temporary File Location
TEMP_DIR=/var/tmp/wwford
TEMP_FASTA=${TEMP_DIR}/minigraph_${HAP_ID}_${JOB_INDEX}.fa

# Outputs
SLURM_LOG_FOLDER=results/2024-08-30_HPRC_minigraph_alignments/slurm_logs
LOG_FOLDER=results/2024-08-30_HPRC_minigraph_alignments/logs/

function TIMESTAMP () { echo $(TZ=US/Pacific; date '+%Y-%m-%d_%H-%M-%S'); }
OUTPUT_GAF=results/2024-08-30_HPRC_minigraph_alignments/minigraph_${HAP_ID}.gaf
OUTPUT_LOG=${LOG_FOLDER}/minigraph_${HAP_ID}_$(TIMESTAMP).log


# - Reads fasta from stdin, 
# - deletes contaminated contigs and sends to stdout
# - not really a limiting step
delete_contaminated_contigs () {
    awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } 
        else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } 
    }' \
    | grep -Ev ${BAD_CONTIGS_REGEX} \
    | tr "\t" "\n"
}

### Main Entrypoint ###

echo "Job ID                ${SLURM_JOB_ID}" > ${OUTPUT_LOG}
echo "SLURM_ARRAY_JOB_ID    ${SLURM_ARRAY_JOB_ID}" >> ${OUTPUT_LOG}
echo "SLURM_ARRAY_TASK_ID   ${SLURM_ARRAY_TASK_ID}" >> ${OUTPUT_LOG}

# Make log dirs if they do not exist
if [ ! -d ${SLURM_LOG_FOLDER} ]; then
    mkdir -p ${SLURM_LOG_FOLDER} ;
fi

if [ ! -d ${LOG_FOLDER} ]; then
    mkdir -p ${LOG_FOLDER} ;
fi

# Check if output already exists
if [ -e ${OUTPUT_GAM} ]; then
    echo "${OUTPUT_GAM} already exists. Exiting w/o error" >> ${OUTPUT_LOG}
    exit 0
fi

# Make temp directory if it does not already exist
if [ ! -d ${TEMP_DIR} ]; then
  mkdir -p ${TEMP_DIR} ;
fi

# Generate temporary fasta from agc
echo "Saving Fasta to ${TEMP_FASTA} $(TIMESTAMP)" >> ${OUTPUT_LOG}
agc getset data/HPRC-yr1.agc ${HAP_ID} \
    | delete_contaminated_contigs > ${TEMP_FASTA}

echo "Saved Fasta                   $(TIMESTAMP)" >> ${OUTPUT_LOG}
echo "Start Align                   $(TIMESTAMP)" >> ${OUTPUT_LOG}

# Align fasta with minigraph
minigraph \
    -x asm \
    ${PANGENOME_GFA} \
    ${TEMP_FASTA} \
    -t ${SLURM_CPUS_PER_TASK} > \
        ${OUTPUT_GAF} 2>> \
        ${OUTPUT_LOG}

echo "Done Align to ${OUTPUT_GAF}   $(TIMESTAMP)" >> ${OUTPUT_LOG}

# Delete temporary fasta
rm ${TEMP_FASTA}
echo "Deleted Fasta                 $(TIMESTAMP)" >> ${OUTPUT_LOG}
