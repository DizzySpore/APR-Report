#!/bin/bash

# Function to get valid numeric input
get_numeric_input() {
    local prompt="$1"
    local min="$2"
    local max="$3"
    local value
    while true; do
        read -p "$prompt" value
        if [[ "$value" =~ ^[0-9]+$ ]] && [ "$value" -ge "$min" ] && [ "$value" -le "$max" ]; then
            echo "$value"
            return 0
        else
            echo "Please enter a valid number between $min and $max"
        fi
    done
}

# Function to get directory input with validation
get_directory_input() {
    local prompt="$1"
    local dir
    local return_dir
    while true; do
        read -p "$prompt" dir
        if [ -z "$dir" ]; then
            dir="$(pwd)"
            echo "Using current working directory: $dir" >&2
            return_dir="$dir"
            break
        elif [ -d "$dir" ]; then
            return_dir="$dir"
            break
        else
            read -p "Directory does not exist. Create it? (y/n): " create
            if [[ "$create" =~ ^[Yy]$ ]]; then
                mkdir -p "$dir"
                return_dir="$dir"
                break
            fi
        fi
    done
    echo "$return_dir"
    return 0
}

# Interactive inputs
echo "=== Metagenomics Workflow Configuration ==="

# Get working directory
DIRECTORY=$(get_directory_input "Enter working directory path (or press Enter for current directory): ")

# Get file extension
read -p "Enter file extension (e.g., .fastq.gz): " EXTENSION
if [ -z "$EXTENSION" ]; then
    EXTENSION=".fastq.gz"
    echo "Using default extension: $EXTENSION"
fi

# Get SLURM resource requirements
echo -e "\n=== SLURM Resource Configuration ==="
WALL_TIME=$(get_numeric_input "Enter wall time in hours (1-240): " 1 240)
# Pad wall time to three digits
WALL_TIME_PADDED=$(printf "%03d" "$WALL_TIME")
MEM_SIZE=$(get_numeric_input "Enter memory in GB (4-450): " 4 450)
NUM_CORES=$(get_numeric_input "Enter number of CPU cores (1-50): " 1 50)
CONCURRENT=$(get_numeric_input "Enter maximum concurrent jobs (1-50): " 1 50)

# Get email for notifications
read -p "Enter email for SLURM notifications: " EMAIL
if [ -z "$EMAIL" ]; then
    EMAIL="lxm697@student.bham.ac.uk"
    echo "Using default email: $EMAIL"
fi

# Confirm settings
echo -e "\n=== Configuration Summary ==="
echo "Working Directory: $DIRECTORY"
echo "File Extension: $EXTENSION"
echo "Wall Time: ${WALL_TIME_PADDED}:00:00"
echo "Memory: ${MEM_SIZE}G"
echo "CPU Cores: $NUM_CORES"
echo "Max Concurrent Jobs: $CONCURRENT"
echo "Email: $EMAIL"

read -p "Proceed with these settings? (y/n): " confirm
if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Configuration cancelled. Exiting..."
    exit 1
fi

# Set search directory
SEARCH_DIR="$DIRECTORY"

# Extract unique base names from files matching the extension using find
# This identifies sample names by removing common suffixes (_R1, _R2, etc.)
echo "$(date): Searching for files in $SEARCH_DIR with extension $EXTENSION..."
SAMPLE_LIST=($(find "$SEARCH_DIR" -maxdepth 1 -type f -name "*$EXTENSION" -exec basename {} "$EXTENSION" \; | sed 's/_[12]\(_sequence\)*$//; s/_R[12]\(_[0-9]\{3\}\)*$//' | grep -v '\.sh$' | sort -u))
if [ ${#SAMPLE_LIST[@]} -eq 0 ]; then
    echo "$(date): No samples found matching extension '$EXTENSION' in $SEARCH_DIR"
    exit 1
fi

# Set number of samples
# This determines how many unique samples were found
N=${#SAMPLE_LIST[@]}
if [ -z "$N" ] || [ "$N" -lt 1 ]; then
    echo "$(date): Error: Number of samples (N) is invalid: $N"
    exit 1
fi

# Log the number of samples and concurrency settings
echo "$(date): Found $N samples: ${SAMPLE_LIST[@]}"
echo "$(date): Queueing $N jobs and setting $CONCURRENT concurrent jobs"
echo "$(date): Working directory: $DIRECTORY"

# Define a variable for the base database path
DATABASES="/rds/projects/h/hallly-microbiome/software/databases"

# Generate the SLURM script
# This script will be submitted to the SLURM workload manager to process the samples
cat > metagenomics_workflow.sh << EOF
#!/bin/bash
#SBATCH --time=${WALL_TIME_PADDED}:00:00
#SBATCH --nodes=1
#SBATCH --mem=${MEM_SIZE}G
#SBATCH --ntasks=${NUM_CORES}
#SBATCH --mail-type=all
#SBATCH --mail-user=${EMAIL}
#SBATCH --array=1-${N}%${CONCURRENT}
#SBATCH -J metagenomics_workflow
#SBATCH --error=slurm-%A_%a.err
#SBATCH --output=slurm-%A_%a.out
#SBATCH --chdir=${DIRECTORY}

# Ensure the SLURM script is not empty
echo "\$(date): SLURM script generated successfully."

# Log the start of the workflow
echo "\$(date): Running combined bioinformatic workflow..."

# Activate Conda environment
CONDA_BASE=\$(conda info --base)
source \$CONDA_BASE/etc/profile.d/conda.sh

# Set trap for errors to mark sample as failed with line number
trap 'echo "Error occurred for \${SAMPLE} at line \$LINENO"; touch "\${STATUS_DIR}/\${SAMPLE}.failed"; exit 1' ERR

#set -x # Enable debugging output

# Retrieve SLURM resource allocations
NUM_THREADS=\$SLURM_NTASKS
MEMORY=\$SLURM_MEM_PER_NODE

# Log resource usage
echo "\$(date): Using \$NUM_THREADS threads and \$MEMORY MB total memory"

# Retrieve the sample name for the current array task
SAMPLE_LIST=(${SAMPLE_LIST[@]})
SAMPLE=\${SAMPLE_LIST[\$SLURM_ARRAY_TASK_ID - 1]}
echo "\$(date): Processing sample: \$SAMPLE"

# Verify SAMPLE is set
if [ -z "\$SAMPLE" ]; then
    echo "\$(date): Error: SAMPLE is not set" >&2
    exit 1
fi

# Identify input files for the sample
# Match paired-end files (_R1 and _R2) based on the sample name
PATTERN="${EXTENSION}"
SEARCH_DIR="${SEARCH_DIR}"
for file in \${SEARCH_DIR}/\${SAMPLE}*\${PATTERN}; do
    if [[ "\$file" =~ _R1(_[0-9]{3})?\${PATTERN}\$ || "\$file" =~ _1(_sequence)?\${PATTERN}\$ ]]; then
        READ1="\$file"
    elif [[ "\$file" =~ _R2(_[0-9]{3})?\${PATTERN}\$ || "\$file" =~ _2(_sequence)?\${PATTERN}\$ ]]; then
        READ2="\$file"
    fi
done

if [ -z "\$READ1" ] || [ -z "\$READ2" ]; then
    echo "Error: Could not find paired files for \$SAMPLE matching pattern '\${PATTERN}'" >&2
    exit 1
fi

# Log the identified input files
echo "\$(date): Input files: \$READ1 and \$READ2"

COMPLETENESS=70
CONTAMINATION=10

# Define directories for workflow outputs
# These directories will store intermediate and final results
QC_READS="${DIRECTORY}/qc_reads"
DEHUMAN_READS="${DIRECTORY}/dehuman_reads"
ASSEMBLY_DIR="${DIRECTORY}/megahit_assembly-test/\${SAMPLE}assembly"
BINNING_DIR="${DIRECTORY}/metawrap-out/metawrap_binref-checkm-\${SAMPLE}-out"
QUAST_ASSEMBLY_OUT="${DIRECTORY}/quast_assembly/\${SAMPLE}/\${SAMPLE}_quast_results"
QUAST_BINS_OUT="${DIRECTORY}/quast_bins/\${SAMPLE}/\${SAMPLE}_quast_results"
KRAKEN2_OUT="${DIRECTORY}/kraken2_bracken_out"
METAPHLAN_OUT="${DIRECTORY}/metaphlan_out"
KAIJU_OUT="${DIRECTORY}/kaiju_out"
GTDBTK_OUT="${DIRECTORY}/gtdbtk-out/gtdbtk-\${SAMPLE}-out"
MAPPING_DIR="${DIRECTORY}/mapping_output"
HUMANN_READS_DIR="${DIRECTORY}/humann_reads_output/\${SAMPLE}"
HUMANN_CONTIGS_DIR="${DIRECTORY}/humann_contigs_output/\${SAMPLE}"
STATUS_DIR="${DIRECTORY}/status"
KNEADDATA_OUT="\${DEHUMAN_READS}/kneaddata_output_\${SAMPLE}"
ASSEMBLY_FILE="\${ASSEMBLY_DIR}/final.contigs.fa"
BINNING_OUT="${DIRECTORY}/metawrap-out/metawrap-checkm-\${SAMPLE}-out"
BIN_DIR="\${BINNING_DIR}/metawrap_\${COMPLETENESS}_\${CONTAMINATION}_bins"
REASSEMBLED_BINS_DIR="\${BINNING_DIR}/reassembled_bins"
BAKTA_OUT="${DIRECTORY}/bakta_out/\${SAMPLE}"
EUKFINDER_OUT="${DIRECTORY}/eukfinder_out/\${SAMPLE}/"
BAKTA_QUAST_OUT="${DIRECTORY}/bakta_quast_out/\${SAMPLE}"
ANTIMASH_OUT="${DIRECTORY}/antismash_out/\${SAMPLE}"
EUK_READS="\${EUKFINDER_OUT}/gunzipped_reads/Intermediate_data/Classified_reads/"
STATS_FILE="${DIRECTORY}/metawrap-out/metawrap_binref-checkm-\${SAMPLE}-out/reassembled_bins/reassembled_bins.stats"
GOOD_BINS_DIR="${DIRECTORY}/good_bins/\${SAMPLE}"
METAQUAST_REFS_FILE="\${REASSEMBLED_BINS_DIR}/metaquast_references.txt"

# Create necessary directories for the workflow
mkdir -p \${QC_READS} \${DEHUMAN_READS} \${ASSEMBLY_DIR} \${QUAST_ASSEMBLY_OUT} \${BINNING_DIR} \\
         \${BINNING_DIR}/abundance \${QUAST_BINS_OUT} \${KRAKEN2_OUT} \${METAPHLAN_OUT} \\
         \${KAIJU_OUT} \${GTDBTK_OUT} \${MAPPING_DIR} \${HUMANN_READS_DIR} \${HUMANN_CONTIGS_DIR} \${STATUS_DIR} \\
         \${BAKTA_OUT} \${EUKFINDER_OUT} \${BAKTA_QUAST_OUT} \${ANTIMASH_OUT} \${GOOD_BINS_DIR}

# Remove any existing status files for the sample
rm -f "\${STATUS_DIR}/\${SAMPLE}.*"

# Set temporary directory for intermediate files
export TMPDIR=/rds/projects/h/hallly-microbiome/tmp/\$SAMPLE
mkdir -p \$TMPDIR

# Function to check the resume point for a sample
# Determines the last completed step in the workflow
check_resume_point() {
    local DIRECTORY="\$1"
    local SAMPLE="\$2"

    echo "\$(date): Checking resume point for \$SAMPLE in $DIRECTORY"

    if [ -z "\$SAMPLE" ] || [ -z "\$DIRECTORY" ]; then
        echo "\$(date): Error: SAMPLE or DIRECTORY is empty in check_resume_point" >&2
        exit 1
    fi

    # Ensure DIRECTORY does not end with a slash
    DIRECTORY="\${DIRECTORY%/}"

    local CHECKPOINTS=(
        "FASTQC:qc_reads/multiqc/multiqc_report.html"
        "FASTP:qc_reads/\${SAMPLE}_fastp.html"
        "KNEADDATA:dehuman_reads/kneaddata_output_\${SAMPLE}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq"
        "KNEADDATA_INSPECTION:qc_reads/multiqc_kneaddata/multiqc_report.html"
        "ASSEMBLY:megahit_assembly-test/\${SAMPLE}assembly/final.contigs.fa"
        "QUAST_ASSEMBLY:quast_assembly/\${SAMPLE}/\${SAMPLE}_quast_results/report.html"
        "BINNING:metawrap-out/metawrap_binref-checkm-\${SAMPLE}-out/reassembled_bins/reassembled_bins/*.fa"
        "TAXONOMY_GTDBTK:gtdbtk-out/gtdbtk-\${SAMPLE}-out/gtdbtk.bac120.summary.tsv"
        "QUAST_BINS:quast_bins/\${SAMPLE}/\${SAMPLE}_quast_results/report.html"
        "TAXONOMY_KRAKEN2:kraken2_bracken_out/\${SAMPLE}_euk-test_bracken_abundance.txt"
        "FOOD_cFMD_TAXONOMY_KRAKEN2:kraken2_bracken_out/\${SAMPLE}_cFMD_bracken_abundance.txt"
        "TAXONOMY_METAPHLAN:metaphlan_out/\${SAMPLE}_metaphlan_profile.txt"
        "TAXONOMY_KAIJU:kaiju_out/\${SAMPLE}.kaiju.out"
        "BIN_FILTER:good_bins/\${SAMPLE}/*.fa"
        "TAXONOMY_EUKFINDER:eukfinder_out/\${SAMPLE}/gunzipped_reads/Eukfinder_results/scf_\${SAMPLE}_eukfinder_short_seq_.Unk.fasta"
        "EUK_TAXONOMY_KRAKEN2:kraken2_bracken_out/\${SAMPLE}_euk-fin_bracken_abundance.txt"
        "ANNOTATION_BAKTA:bakta_out/\${SAMPLE}/\${SAMPLE}_Bakta_annotation/loop_completed.txt"
        "QUAST_BAKTA:bakta_quast_out/\${SAMPLE}/\${SAMPLE}_Bakta_annotation/"
        "ANTISMASH:antismash_out/\${SAMPLE}/antismash_loop_completed.txt"
        "HUMANN_READS:humann_reads_output/\${SAMPLE}/out.\${SAMPLE}_combined_pathabundance.tsv"
        "HUMANN_CONTIGS:humann_contigs_output/humann_loop_completed.txt"
        "COMBINING:humann_contigs/merged_pathabundance.tsv"
        "VISUALIZATION:status/abundance_visualization_complete.txt"
        "R_ANALYSES:maaslin3_contigs/all_results.tsv"
    )

    local RESUME_SET=false
    local LAST_FOUND=""
    
    # First, check all files to establish the last completed step
    for checkpoint in "\${CHECKPOINTS[@]}"; do
        local step="\${checkpoint%%:*}"
        local path="\${checkpoint#*:}"
        local full_path="\${DIRECTORY}/\${path}"

        echo "\$(date): Checking \$step: \$full_path"
        if [[ "\$path" =~ \*\. ]]; then
            if ls \$full_path 1>/dev/null 2>&1; then
                echo "\$(date): Found output for \$step: \$full_path"
                LAST_FOUND="\$step"
            else
                echo "\$(date): No output found for \$step"
            fi
        else
            if [ -f "\$full_path" ]; then
                echo "\$(date): Found output for \$step: \$full_path"
                LAST_FOUND="\$step"
            else
                echo "\$(date): No output found for \$step"
            fi
        fi
    done
    
    # After checking all files, determine where to resume based on the last found checkpoint
    if [ -n "\$LAST_FOUND" ]; then
        case "\$LAST_FOUND" in
            "R_ANALYSES")
                echo "\$(date): All steps completed for \$SAMPLE. No further action needed."
                touch "\${STATUS_DIR}/\${SAMPLE}.done"
                RESUME_FROM=""
                ;;
            "VISUALIZATION") RESUME_FROM="R_ANALYSES" ;;
            "COMBINING") RESUME_FROM="VISUALIZATION" ;;
            "HUMANN_CONTIGS") RESUME_FROM="COMBINING" ;;
            "HUMANN_READS") RESUME_FROM="HUMANN_CONTIGS" ;;
            "ANTISMASH") RESUME_FROM="HUMANN_READS" ;;
            "QUAST_BAKTA") RESUME_FROM="ANTISMASH" ;;
            "ANNOTATION_BAKTA") RESUME_FROM="QUAST_BAKTA" ;;
            "EUK_TAXONOMY_KRAKEN2") RESUME_FROM="ANNOTATION_BAKTA" ;;
            "TAXONOMY_EUKFINDER") RESUME_FROM="EUK_TAXONOMY_KRAKEN2" ;;
            "BIN_FILTER") RESUME_FROM="TAXONOMY_EUKFINDER" ;;
            "TAXONOMY_KAIJU") RESUME_FROM="BIN_FILTER" ;;
            "TAXONOMY_METAPHLAN") RESUME_FROM="TAXONOMY_KAIJU" ;;
            "FOOD_cFMD_TAXONOMY_KRAKEN2") RESUME_FROM="TAXONOMY_METAPHLAN" ;;
            "TAXONOMY_KRAKEN2") RESUME_FROM="FOOD_cFMD_TAXONOMY_KRAKEN2" ;;
            "QUAST_BINS") RESUME_FROM="TAXONOMY_KRAKEN2" ;;
            "TAXONOMY_GTDBTK") RESUME_FROM="QUAST_BINS" ;;
            "BINNING") RESUME_FROM="TAXONOMY_GTDBTK" ;;
            "QUAST_ASSEMBLY") RESUME_FROM="BINNING" ;;
            "ASSEMBLY") RESUME_FROM="QUAST_ASSEMBLY" ;;
            "KNEADDATA_INSPECTION") RESUME_FROM="ASSEMBLY" ;;
            "KNEADDATA") RESUME_FROM="KNEADDATA_INSPECTION" ;;
            "FASTP") RESUME_FROM="KNEADDATA" ;;
            "FASTQC") RESUME_FROM="FASTP" ;;
            *) RESUME_FROM="FASTQC" ;;  # Default if none matched
        esac
    else
        # If no checkpoint was found, start from the beginning
        RESUME_FROM="FASTQC"
    fi
    
    echo "Determined resume point: \$RESUME_FROM based on last found: \$LAST_FOUND"
    echo "export RESUME_FROM=\"\$RESUME_FROM\"" > /tmp/resume_\${SAMPLE}.sh
}

# Call the function and source the result to set RESUME_FROM
check_resume_point "${DIRECTORY}" "\$SAMPLE"
echo "\$(date): Sourcing resume point from /tmp/resume_\${SAMPLE}.sh..."
cat /tmp/resume_\${SAMPLE}.sh  # Debug: Show what's being sourced
source /tmp/resume_\${SAMPLE}.sh
rm -f /tmp/resume_\${SAMPLE}.sh

echo "\$(date): Workflow will resume from: \$RESUME_FROM"

# Execute workflow steps based on RESUME_FROM
case "\$RESUME_FROM" in
    "FASTQC")
        echo "\$(date): Starting FastQC and MultiQC..."
        conda activate fastqc-env
        mkdir -p \${QC_READS}/fastqc_raw \${QC_READS}/multiqc

        # Run FastQC on raw reads
        echo "\$(date): Running FastQC on raw reads for \${SAMPLE}..."
        fastqc -t \$NUM_THREADS -o \${QC_READS}/fastqc_raw "\$READ1" "\$READ2"

        # Run MultiQC to summarize FastQC results
        echo "\$(date): Running MultiQC to summarize FastQC results..."
        multiqc -o \${QC_READS}/multiqc \${QC_READS}/fastqc_raw

        conda deactivate
        ;&
    "FASTP")
        echo "\$(date): Starting fastp..."
        echo "\$(date): Running fastp on \$READ1 \$READ2"
        module purge; module load bluebear
        module load bear-apps/2023a
        module load fastp/0.24.0-GCC-12.3.0
        fastp -v
        if [ -z "\$NUM_THREADS" ] || [ -z "\$READ1" ] || [ -z "\$READ2" ]; then
            echo "\$(date): Error: NUM_THREADS, READ1, or READ2 is not set" >&2
            exit 1
        fi
        fastp --thread "\$NUM_THREADS" -q 20 \\
            -i "\$READ1" -I "\$READ2" \\
            -o "\${QC_READS}/out.\${SAMPLE}.R1.fq.gz" -O "\${QC_READS}/out.\${SAMPLE}.R2.fq.gz" \\
            --html "\${QC_READS}/\${SAMPLE}_fastp.html" --json "\${QC_READS}/\${SAMPLE}_fastp.json"
        ;&
    "KNEADDATA")
        if [ ! -f "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq" ] || [ ! -f "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq" ]; then
            echo "\$(date): Starting kneaddata..."
            kneaddata --version
            echo "\$(date): Running kneaddata on \${QC_READS}/out.\${SAMPLE}.R1.fq.gz \${QC_READS}/out.\${SAMPLE}.R2.fq.gz"
            kneaddata --threads \$NUM_THREADS --bypass-trim --reorder \\
                --input1 \${QC_READS}/out.\${SAMPLE}.R1.fq.gz --input2 \${QC_READS}/out.\${SAMPLE}.R2.fq.gz \\
                -db ${DATABASES}/human-GRCh38_no-alt/GRCh38_noalt_as/ \\
                --output \${KNEADDATA_OUT}
            if [ -f "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq" ] && [ -f "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq" ]; then
                rm -f "\${QC_READS}/out.\${SAMPLE}.R1.fq.gz" "\${QC_READS}/out.\${SAMPLE}.R2.fq.gz"
            fi
        else
            echo "\$(date): kneaddata output already exists, skipping..."
        fi

        echo "\$(date): Cleaning up Kneaddata outputs for \$SAMPLE"
        KNEADDATA_LOG="\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata.log"
        ;&
    "KNEADDATA_INSPECTION")
        echo "\$(date): Starting FastQC and MultiQC on Kneaddata outputs..."
        conda activate fastqc-env
        mkdir -p \${QC_READS}/fastqc_kneaddata \${QC_READS}/multiqc_kneaddata

        # Run FastQC on Kneaddata paired outputs
        echo "\$(date): Running FastQC on Kneaddata paired outputs for \${SAMPLE}..."
        fastqc -t \$NUM_THREADS -o \${QC_READS}/fastqc_kneaddata \\
            "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq" \\
            "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq"

        # Run MultiQC to summarize FastQC results
        echo "\$(date): Running MultiQC to summarize FastQC results for Kneaddata outputs..."
        multiqc -o \${QC_READS}/multiqc_kneaddata \${QC_READS}/fastqc_kneaddata

        conda deactivate
        ;&
    "ASSEMBLY")
        echo "\$(date): Starting megahit assembly..."
        megahit --version
        echo "\$(date): Running megahit assembly on kneaddata_output_\${SAMPLE}"
        megahit \\
            -1 \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq \\
            -2 \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq \\
            --num-cpu-threads \$NUM_THREADS \\
            -o \${ASSEMBLY_DIR} --force
        if [ -f "\${ASSEMBLY_FILE}" ]; then
            echo "\$(date): Cleaning up intermediate MEGAHIT files for \${SAMPLE}..."
            rm -rf "\${ASSEMBLY_DIR}/intermediate_contigs" \\
                   "\${ASSEMBLY_DIR}/k\*.gfa" \\
                   "\${ASSEMBLY_DIR}/k\*.fa" \\
                   "\${ASSEMBLY_DIR}/tmp.\*" \\
                   "\${ASSEMBLY_DIR}/\*.bin" 2>/dev/null || true
            echo "\$(date): Intermediate files removed, retaining final.contigs.fa and log/opts.txt"
        else
            echo "\$(date): Error: final.contigs.fa not found, retaining all MEGAHIT outputs" >&2
            exit 1
        fi
        
        echo "\$(date): Cleaning up Kneaddata outputs for \$SAMPLE"
        echo "\$(date): Removing kneaddata intermediate files for \${KNEADDATA_OUT}/out.\${SAMPLE}\* ..."
        rm -rf "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata.repeats.removed.1.fastq" \\
               "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata.repeats.removed.2.fastq" \\
               "\${KNEADDATA_OUT}/out.\${SAMPLE}_kneaddata_GRCh38_noalt_as_bowtie2_paired_contam_1.fastq" \\
               "\${KNEADDATA_OUT}/out.\${SAMPLE}_kneaddata_GRCh38_noalt_as_bowtie2_paired_contam_2.fastq" \\
               "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_contam_1.fastq \\
               "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_contam_2.fastq \\
               "\${KNEADDATA_OUT}/out.\${SAMPLE}_kneaddata_GRCh38_noalt_as_bowtie2_bowtie2_unmatched_1_contam.fastq" \\
               "\${KNEADDATA_OUT}/out.\${SAMPLE}_kneaddata_GRCh38_noalt_as_bowtie2_bowtie2_unmatched_1_contam.fastq" 2>/dev/null || true
        ;&
    "QUAST_ASSEMBLY")
        echo "\$(date): Starting QUAST on assembly..."
        conda activate quast-env
        quast --version
        echo "\$(date): Running QUAST on \${SAMPLE} assembly with coverage calculation"
        quast \${ASSEMBLY_FILE} \\
            -o \${QUAST_ASSEMBLY_OUT} \\
            --threads \$NUM_THREADS 
        conda deactivate
        ;&
    "BINNING")
        echo "\$(date): Starting metawrap binning..."
        conda activate metawrap-env
        metawrap --version
        echo "\$(date): Binning with metaWRAP for \${SAMPLE}assembly"
        metawrap binning --metabat2 --maxbin2 --concoct \\
            -a \${ASSEMBLY_FILE} \\
            -t \$NUM_THREADS -m \$MEMORY \\
            -o \${BINNING_OUT} \\
            \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq

        echo "\$(date): Bin refining with metaWRAP for \${SAMPLE} bins with \${COMPLETENESS}% completeness and \${CONTAMINATION}% contamination"
        metawrap bin_refinement \\
            -A \${BINNING_OUT}/concoct_bins \\
            -B \${BINNING_OUT}/maxbin2_bins \\
            -C \${BINNING_OUT}/metabat2_bins \\
            -t \$NUM_THREADS -m \$MEMORY \\
            -c \$COMPLETENESS -x \$CONTAMINATION \\
            -o \${BINNING_DIR}

        # Reassemble bins
        echo "\$(date): Starting bin reassembly with metaWRAP..."
        REASSEMBLED_BINS_DIR="\${BINNING_DIR}/reassembled_bins"
        metawrap reassemble_bins \\
            -o \${REASSEMBLED_BINS_DIR} \\
            -1 \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq \\
            -2 \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq \\
            -t \$NUM_THREADS -m \$MEMORY \\
            -c \$COMPLETENESS -x \$CONTAMINATION \\
            -b \${BINNING_DIR}/metawrap_\${COMPLETENESS}_\${CONTAMINATION}_bins

        # Update BIN_DIR to point to reassembled bins
        BIN_DIR="\${REASSEMBLED_BINS_DIR}"

        conda deactivate
        if [ -f "\${BIN_DIR}/reassembled_bins.stats" ]; then
            echo "\$(date): Cleaning up intermediate metaWRAP files for \${SAMPLE} within \${BIN_DIR}/..."
            rm -rf "\${BINNING_OUT}" \\
                   "\${BINNING_DIR}/work_files" \\
                   "\${BINNING_DIR}/metabat2_bins" \\
                   "\${BINNING_DIR}/maxbin2_bins" \\
                   "\${BIN_DIR}/work_files" \\
                   "\${BINNING_DIR}/concoct_bins" 2>/dev/null || true
            echo "\$(date): Intermediate bins removed, retaining reassembled bins and stats"
        else
            echo "\$(date): Error: \${BIN_DIR}/reassembled_bins/bin* not found, retaining all metaWRAP outputs" >&2
            exit 1
        fi
        ;&
    "TAXONOMY_GTDBTK")
        echo "\$(date): Starting GTDB-Tk..."
        export GTDBTK_DATA_PATH=${DATABASES}/gtdbtk/release220
        gtdbtk --version
        echo "\$(date): Classifying taxonomy with gtdbtk on \${SAMPLE} reassembled bins (\${COMPLETENESS}_\${CONTAMINATION} threshold)"
        gtdbtk classify_wf --cpus \$NUM_THREADS \
            --mash_db ${DATABASES}/gtdbtk/mash_db \
            --extension .fa \
            --genome_dir \${REASSEMBLED_BINS_DIR}/reassembled_bins/ \
            --out_dir \${GTDBTK_OUT}

        # Create a file for MetaQUAST references-list
        METAQUAST_REFS_FILE="\${REASSEMBLED_BINS_DIR}/metaquast_references.txt"
        echo "# References for MetaQUAST --references-list option" > \${METAQUAST_REFS_FILE}
        echo "# Generated from GTDB-Tk classifications on \$(date)" >> \${METAQUAST_REFS_FILE}

        # Rename bins to include species information
        echo "\$(date): Renaming bins to include taxonomic classification..."
        if [ -f "\${GTDBTK_OUT}/classify/gtdbtk.bac120.summary.tsv" ]; then
            while IFS=\$'\t' read -r user_genome classification other; do
                if [ "\$user_genome" != "user_genome" ]; then
                    # Extract species name (last part of the classification)
                    species=\$(echo \$classification | rev | cut -d';' -f1 | rev | sed 's/^s__//' | tr ' ' '_')
                    genus=\$(echo \$classification | rev | cut -d';' -f2 | rev | sed 's/^g__//' | tr ' ' '_')

                    # Check if genus and species are N/A
                    if [ "\$genus" = "N/A" ] && [ "\$species" = "N/A" ]; then
                        genus="Unknown"
                        species="organism"
                    fi

                    # If species is empty, use genus_sp format
                    if [ -z "\$species" ]; then
                        species="\${genus}_sp"
                    fi
                    
                    old_path="\${REASSEMBLED_BINS_DIR}/reassembled_bins/\${user_genome}.fa"
                    new_path="\${REASSEMBLED_BINS_DIR}/reassembled_bins/\${user_genome}_\${species}.fa"
                        
                    # Rename the file
                    mv "\$old_path" "\$new_path"
                    echo "Renamed: \${user_genome}.fa → \${user_genome}_\${species}.fa"
                    
                    # Add to MetaQUAST references file - replace underscores with spaces in genus and species for better NCBI search
                    if [[ \$genus != "Unknown" && ! \$species =~ "_sp" ]]; then
                        # Only include entries with proper species designations (not genus_sp)
                        metaquast_name=\$(echo "\$genus \$species" | tr '_' ' ')
                        echo "\$metaquast_name" >> \${METAQUAST_REFS_FILE}
                    fi
                fi
            done < "\${GTDBTK_OUT}/classify/gtdbtk.bac120.summary.tsv"
        fi

        if [ -f "\${GTDBTK_OUT}/classify/gtdbtk.ar53.summary.tsv" ]; then
            while IFS=\$'\t' read -r user_genome classification other; do
                if [ "\$user_genome" != "user_genome" ]; then
                    # Extract species name (last part of the classification)
                    species=\$(echo \$classification | rev | cut -d';' -f1 | rev | sed 's/^s__//' | tr ' ' '_')
                    genus=\$(echo \$classification | rev | cut -d';' -f2 | rev | sed 's/^g__//' | tr ' ' '_')

                    # If species is empty, use genus_sp format
                    if [ -z "\$species" ]; then
                        species="\${genus}_sp"
                    fi
                    
                    old_path="\${REASSEMBLED_BINS_DIR}/reassembled_bins/\${user_genome}.fa"
                    new_path="\${REASSEMBLED_BINS_DIR}/reassembled_bins/\${user_genome}_\${species}.fa"
                        
                    # Rename the file if it exists (might have already been renamed if both bac120 and ar53 results exist)
                    if [ -f "\$old_path" ]; then
                        mv "\$old_path" "\$new_path"
                        echo "Renamed: \${user_genome}.fa → \${user_genome}_\${species}.fa"
                        
                        # Add to MetaQUAST references file - replace underscores with spaces in genus and species for better NCBI search
                        if [[ \$genus != "Unknown" && ! \$species =~ "_sp" ]]; then
                            # Only include entries with proper species designations (not genus_sp)
                            metaquast_name=\$(echo "\$genus \$species" | tr '_' ' ')
                            echo "\$metaquast_name" >> \${METAQUAST_REFS_FILE}
                        fi
                    fi
                fi
            done < "\${GTDBTK_OUT}/classify/gtdbtk.ar53.summary.tsv"
        fi

        echo "\$(date): Created MetaQUAST references list at \${METAQUAST_REFS_FILE}"
        echo "Use with: metaquast.py --references-list \${METAQUAST_REFS_FILE} ..."
        ;&
    "QUAST_BINS")
        echo "\$(date): Starting metaQUAST on reassembled bins..."
        conda activate quast-env
        echo "\$(date): Running metaQUAST on \${SAMPLE} reassembled bins (\${COMPLETENESS}_\${CONTAMINATION} threshold) with coverage calculation"

        # Run metaQUAST on reassembled bins
        metaquast.py \${BINNING_DIR}/reassembled_bins/reassembled_bins/*.fa \\
            -o quast_bins/\${SAMPLE}/\${SAMPLE}_results/ \\
            --threads \$NUM_THREADS \\
            --references-list \${METAQUAST_REFS_FILE} \\
            --pe1 \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq \\
            --pe2 \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq 

        conda deactivate
        echo "\$(date): metaQUAST analysis completed for reassembled bins of \${SAMPLE}."
        
        # Removing the MEGAHIT assembly
        echo "\$(date): Cleaning up MEGAHIT assembly outputs for \$SAMPLE"
        rm -rf "\${ASSEMBLY_DIR}/final.contigs.fa" 2>/dev/null || true
        ;&
    "TAXONOMY_KRAKEN2")
        echo "\$(date): Starting Kraken2 and Bracken..."
        kraken2 --version
        bracken -v
        echo "\$(date): Running Kraken2 and Bracken on \${SAMPLE} reads (Bacteria and Archaea)"
        kraken2 \\
            --db ${DATABASES}/krakendb/standard/ \\
            --paired --threads \$NUM_THREADS \\
            --report \${KRAKEN2_OUT}/\${SAMPLE}.prok-kraken2_report.txt \\
            --output \${KRAKEN2_OUT}/\${SAMPLE}.prok-kraken2_output.txt \\
            \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq
        bracken \\
            -d ${DATABASES}/krakendb/standard/ \\
            -i \${KRAKEN2_OUT}/\${SAMPLE}.prok-kraken2_report.txt \\
            -o \${KRAKEN2_OUT}/\${SAMPLE}_prok-bracken_abundance.txt \\
            -r 100

        kraken2 \\
            --db ${DATABASES}/krakendb/Archaea/ \\
            --paired --threads \$NUM_THREADS \\
            --report \${KRAKEN2_OUT}/\${SAMPLE}.arc-test.kraken2_report.txt \\
            --output \${KRAKEN2_OUT}/\${SAMPLE}.arc-test.kraken2_output.txt \\
            \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq
        bracken \\
            -d ${DATABASES}/krakendb/Archaea/ \\
            -i \${KRAKEN2_OUT}/\${SAMPLE}.arc-test.kraken2_report.txt \\
            -o \${KRAKEN2_OUT}/\${SAMPLE}_arc-test_bracken_abundance.txt \\
            -r 100
        
        kraken2 \\
            --db ${DATABASES}/krakendb/Eukaryote/ \\
            --paired --threads \$NUM_THREADS \\
            --report \${KRAKEN2_OUT}/\${SAMPLE}.euk-test.kraken2_report.txt \\
            --output \${KRAKEN2_OUT}/\${SAMPLE}.euk-test.kraken2_output.txt \\
            \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq
        bracken \\
            -d ${DATABASES}/krakendb/Eukaryote/ \\
            -i \${KRAKEN2_OUT}/\${SAMPLE}.euk-test.kraken2_report.txt \\
            -o \${KRAKEN2_OUT}/\${SAMPLE}_euk-test_bracken_abundance.txt \\
            -r 100
        ;&
    "FOOD_cFMD_TAXONOMY_KRAKEN2")
        kraken2 \\
            --db ${DATABASES}/cFMD/kraken2-cFMD-db/cFMD/ \\
            --paired --threads \$NUM_THREADS \\
            --report \${KRAKEN2_OUT}/\${SAMPLE}.cFMD-kraken2_report.txt \\
            --output \${KRAKEN2_OUT}/\${SAMPLE}.cFMD-kraken2_output.txt \\
            \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq

        # Check if there are any classified reads
        UNCLASSIFIED_PCT=\$(head -n1 \${KRAKEN2_OUT}/\${SAMPLE}.cFMD-kraken2_report.txt | awk '{print \$1}')
        if (( \$(echo "\$UNCLASSIFIED_PCT < 100" | bc -l) )); then
            echo "\$(date): Classified reads found (\$UNCLASSIFIED_PCT% unclassified). Running Bracken..."
            bracken \\
                -d ${DATABASES}/cFMD/kraken2-cFMD-db/cFMD/ \\
                -i \${KRAKEN2_OUT}/\${SAMPLE}.cFMD-kraken2_report.txt \\
                -o \${KRAKEN2_OUT}/\${SAMPLE}_cFMD-bracken_abundance.txt \\
                -r 100
        else
            echo "\$(date): No classified reads found (100% unclassified). Skipping Bracken."
            touch \${KRAKEN2_OUT}/\${SAMPLE}_cFMD-bracken_abundance.txt
        fi
        ;&
    "TAXONOMY_METAPHLAN")
        echo "\$(date): Starting MetaPhlAn..."
        conda activate metaphlan-env
        metaphlan --version
        echo "\$(date): Running MetaPhlAn on \${SAMPLE} reads"
        metaphlan \\
            \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq,\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq \\
            --input_type fastq \\
            -o \${METAPHLAN_OUT}/\${SAMPLE}_metaphlan_profile.txt \\
            --nproc \$NUM_THREADS \\
            --db_dir ${DATABASES}/metaphlandb/ \\
            --mapout \${METAPHLAN_OUT}/\${SAMPLE}_bowtie2out.bowtie2.bz2
        conda deactivate
        ;&
    "TAXONOMY_KAIJU")
        echo "\$(date): Starting Kaiju..."
        conda activate kaiju-env
        echo "\$(date): Running Kaiju on \${SAMPLE} reads"
        kaiju \\
            -t ${DATABASES}/kaijudb/nodes.dmp \\
            -f ${DATABASES}/kaijudb/nr_euk/kaiju_db_nr_euk.fmi \\
            -i \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq \\
            -j \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq \\
            -o \${KAIJU_OUT}/\${SAMPLE}.kaiju.out \\
            -v -z \$NUM_THREADS
        echo "\$(date): Transforming Kaiju output to Krona format..."
        kaiju2krona \\
            -t ${DATABASES}/kaijudb/nodes.dmp \\
            -n ${DATABASES}/kaijudb/names.dmp \\
            -i \${KAIJU_OUT}/\${SAMPLE}.kaiju.out -o \${KAIJU_OUT}/\${SAMPLE}.kaiju.out.krona
        conda deactivate
        echo "\$(date): Generating Krona visualization..."
        conda activate krona-env
        ktImportText \\
            -o \${KAIJU_OUT}/\${SAMPLE}.kaiju.out.html \\
            \${KAIJU_OUT}/\${SAMPLE}.kaiju.out.krona
        echo "\$(date): Adding names to Kaiju output..."
        conda deactivate
        conda activate kaiju-env
        kaiju-addTaxonNames \\
            -t ${DATABASES}/kaijudb/nodes.dmp \\
            -n ${DATABASES}/kaijudb/names.dmp \\
            -i \${KAIJU_OUT}/\${SAMPLE}.kaiju.out \\
            -o \${KAIJU_OUT}/\${SAMPLE}.kaiju.out.names
        echo "\$(date): Summarizing Kaiju results..."
        kaiju2table \\
            -t ${DATABASES}/kaijudb/nodes.dmp \\
            -n ${DATABASES}/kaijudb/names.dmp \\
            -o \${KAIJU_OUT}/\${SAMPLE}.kaiju.report.txt \\
            -r species \\
            -m 0.5 -v \\
            \${KAIJU_OUT}/\${SAMPLE}.kaiju.out.names
        echo "\$(date): Combining Kaiju results with Kraken2 results..."
        # Sort Kraken2 and Kaiju results
        sort -k2,2 \${KAIJU_OUT}/\${SAMPLE}.kaiju.out >> \${KAIJU_OUT}/\${SAMPLE}.kaiju.out.sort
        sort -k2,2 \${KRAKEN2_OUT}/\${SAMPLE}.prok-kraken2_output.txt >> \${KRAKEN2_OUT}/\${SAMPLE}.prok-kraken2_output.sort
        # Combine Kraken2 and Kaiju results
        kaiju-mergeOutputs \\
            -i \${KAIJU_OUT}/\${SAMPLE}.kaiju.out.sort \\
            -j \${KRAKEN2_OUT}/\${SAMPLE}.prok-kraken2_output.sort \\
            -o \${KAIJU_OUT}/\${SAMPLE}.kaiju.kraken2.combined.txt \\
            -t ${DATABASES}/kaijudb/nodes.dmp \\
            -v -c lowest
        
        conda deactivate

        #bracken -v
        #echo "\$(date): Running Bracken on combined Kaiju and Kraken2 results..."
        #bracken \\
        #    -d ${DATABASES}/krakendb/standard/ \\
        #    -i \${KAIJU_OUT}/\${SAMPLE}.kaiju.kraken2.combined.txt \\
        #    -o \${KAIJU_OUT}/\${SAMPLE}_combined_bracken_abundance.txt \\
        #    -r 100
        ;&
    "BIN_FILTER")

        echo "\$(date): Starting bin filtering..."
        if [ ! -f "\${STATS_FILE}" ]; then
            echo "\$(date): Error: Stats file not found: \${STATS_FILE}" >&2
            exit 1
        fi

        # Ensure good bins directory exists
        mkdir -p "\${GOOD_BINS_DIR}"

        # Filter bins based on completeness >80 and contamination <10
        echo "\$(date): Filtering bins with >80% completeness and <10% contamination..."
        GOOD_BINS=\$(awk -F'\t' 'NR > 1 && \$2 > 80 && \$3 < 10 {print \$1}' "\${STATS_FILE}")
        echo "\$(date): Found \$(echo "\${GOOD_BINS}" | wc -l) good bins: \${GOOD_BINS}"
        echo "\${GOOD_BINS}" | while read -r bin; do
            # Find the renamed bin file (with taxonomy in the name)
            BIN_PATTERN="${DIRECTORY}/metawrap-out/metawrap_binref-checkm-\${SAMPLE}-out/reassembled_bins/reassembled_bins/\${bin}_*.fa"
            BIN_FILE=\$(ls \$BIN_PATTERN 2>/dev/null | head -n 1)
            if [ -n "\$BIN_FILE" ] && [ -f "\$BIN_FILE" ]; then
                cp "\$BIN_FILE" "\${GOOD_BINS_DIR}/"
                if [ \$? -eq 0 ]; then
                    echo "\$(date): Successfully copied \$BIN_FILE to \${GOOD_BINS_DIR}/"
                else
                    echo "\$(date): Error: Failed to copy \$BIN_FILE to \${GOOD_BINS_DIR}/" >&2
                fi
            else
                echo "\$(date): Warning: Bin file not found for pattern: \$BIN_PATTERN" >&2
            fi
        done
        echo "\$(date): Bin filtering completed for \${SAMPLE}. Filtered bins are in \${GOOD_BINS_DIR}."

        ;&
    "TAXONOMY_EUKFINDER")
        set +e
        echo "\$(date): Starting EukFinder..."
        conda activate eukfinder-env
        echo "\$(date): Running EukFinder on \${SAMPLE} reads"

        # Gunzip read1 and read2 into the temporary directory
        mkdir -p "\${EUKFINDER_OUT}/gunzipped_reads/"

        cd "\${EUKFINDER_OUT}/gunzipped_reads/"
        cp "\$READ1" .
        cp "\$READ2" .

        echo "\$(date): Gunzipping reads..."
        gunzip "\${SAMPLE}"*
        for file in \${SAMPLE}*.fastq; do
            if [[ "\$file" =~ _R1(_[0-9]{3})?\.fastq\$ || "\$file" =~ _1(_sequence)?\.fastq\$ ]]; then
                GUNZIP1="\$file"
            elif [[ "\$file" =~ _R2(_[0-9]{3})?\.fastq\$ || "\$file" =~ _2(_sequence)?\.fastq\$ ]]; then
                GUNZIP2="\$file"
            fi
        done

        eukfinder read_prep \\
            --r1 \$GUNZIP1 \\
            --r2 \$GUNZIP2 \\
            -o \${SAMPLE}_eukfinder_read_prep \\
            -n \$NUM_THREADS \\
            --hcrop 10 -l 15 -t 15 --wsize 40 --qscore 25 --mlen 20 --mhlen 20 \\
            -i ${DATABASES}/eukfinder-db/illumina_adapters/TrueSeq2_NexteraSE-PE.fa \\
            --hg ${DATABASES}/eukfinder-db/GRCh38.p13/GCF_000001405.39_GRCh38.p13_human_genome.fna \\
            --cdb ${DATABASES}/eukfinder-db/CentrifugeDB/Centrifuge_DB_2024
        if [ \$? -ne 0 ]; then
            echo "\$(date): eukfinder read_prep failed for \$SAMPLE" >> \${EUKFINDER_OUT}/error.log
        fi

        rm -f "\$GUNZIP1" "\$GUNZIP2"

        eukfinder short_seqs \\
            --r1 \${SAMPLE}_eukfinder_read_prep_p.1.fastq \\
            --r2 \${SAMPLE}_eukfinder_read_prep_p.2.fastq \\
            --un \${SAMPLE}_eukfinder_read_prep_un.fastq \\
            -o \${SAMPLE}_eukfinder_short_seq_ \\
            -t T --max_m 100 -e 0.01 --pid 60 --cov 30 --mhlen 50 -n \$NUM_THREADS -z 10\\
            --pclass \${SAMPLE}_eukfinder_read_prep_centrifuge_P \\
            --uclass \${SAMPLE}_eukfinder_read_prep_centrifuge_UP \\
            --plast-database ${DATABASES}/eukfinder-db/PlastDB/PlastDB.fasta \\
            --plast-id-map ${DATABASES}/eukfinder-db/PlastDB/PlastDB_map.txt \\
            --cdb ${DATABASES}/eukfinder-db/CentrifugeDB/Centrifuge_DB_2024
        conda deactivate
        if [ \$? -ne 0 ]; then
            echo "\$(date): eukfinder short_seqs failed for \$SAMPLE" >> \${EUKFINDER_OUT}/error.log
        fi
        set -e
        cd "${DIRECTORY}"
        ;&
    "EUK_TAXONOMY_KRAKEN2")
        echo "˘\$(date): Concatenating EukFinder output files..."
        cat \${EUK_READS}/\${SAMPLE}_eukfinder_short_seq_.{Euk,EUnk,Unk}.R1.fq > \${EUK_READS}/combined_target_R1.fq
        cat \${EUK_READS}/\${SAMPLE}_eukfinder_short_seq_.{Euk,EUnk,Unk}.R2.fq > \${EUK_READS}/combined_target_R2.fq
        cat \${EUK_READS}/\${SAMPLE}_eukfinder_short_seq_.{Euk,EUnk,Unk}.un.fq > \${EUK_READS}/combined_target_un.fq
        echo "\$(date): Combined EukFinder output files into combined_target_R1.fq, combined_target_R2.fq, and combined_target_un.fq"
        echo "\$(date): Starting Kraken2 and Bracken for Eukaryotes..."
        kraken2 --version
        bracken -v
        echo "\$(date): Running Kraken2 and Bracken on \${SAMPLE} reads (Eukaryotes)"
        kraken2 \\
            --db ${DATABASES}/krakendb/Eukaryote/ \\
            --paired \\
            --report \${KRAKEN2_OUT}/\${SAMPLE}.euk-fin.kraken2_report.txt \\
            --output \${KRAKEN2_OUT}/\${SAMPLE}.euk-fin.kraken2_output.txt \\
            \${EUK_READS}/combined_target_R1.fq \${EUK_READS}/combined_target_R2.fq 
        bracken \\
            -d ${DATABASES}/krakendb/EuPathDB46/ \\
            -i \${KRAKEN2_OUT}/\${SAMPLE}.euk-fin.kraken2_report.txt \\
            -o \${KRAKEN2_OUT}/\${SAMPLE}_euk-fin_bracken_abundance.txt \\
            -r 100
        ;&
    "ANNOTATION_BAKTA")
        echo "\$(date): Starting Bakta annotation..."
        conda activate bakta-env
        bakta --version

        # Check if GTDB-Tk taxonomy file exists
        if [ ! -f "\${GTDBTK_OUT}/gtdbtk.bac120.summary.tsv" ]; then
            echo "\$(date): Error: GTDB-Tk taxonomy file not found. Skipping Bakta annotation." >&2
            exit 1
        fi

        # Process each FASTA file in the specified directory
        for fasta_file in \${GOOD_BINS_DIR}/*.fa; do
            # Skip if not a regular file
            if [ ! -f "\$fasta_file" ]; then
                echo "Skipping \$fasta_file: not a regular file"
                continue
            fi
            
            # Extract basename from filename
            basename=$(basename "\$fasta_file")
            
            # Extract genus and species from filename
            if [[ "\$basename" =~ bin\.[0-9]+\.(strict|orig|permissive)_([^_]+)_?([^\.]*)?\.fa$ ]]; then
                bin_type=\${BASH_REMATCH[1]}
                genus=\${BASH_REMATCH[2]}
                species=\${BASH_REMATCH[3]}
                
                # Handle special case where species might not be present
                if [ -z "\$species" ]; then
                    organism="\$genus"
                    genus_arg="--genus \$genus"
                    species_arg=""
                else
                    # Handle cases where genus might have prefix like "NHYM01"
                    organism="\${genus}_\${species}"
                    genus_arg="--genus \$genus"
                    species_arg="--species \$species"
                fi
            else
                # Fallback for unrecognized naming pattern
                organism=$(echo "\$basename" | sed -E 's/^bin\.[0-9]+\.[^_]+_(.+)\.fa$/\1/')
                genus_arg=""
                species_arg="--species \$organism"
            fi
            
            echo "\$(date): Processing: \$fasta_file"
            echo "\$(date): Extracted genus: \$genus, species: \$species"
            echo "\$(date): Using organism name: \$organism"
            
            # Run bakta with the extracted genus and species
            echo "Running bakta with organism: \$organism"
            echo "\$(date): Annotating bin \${fasta_file}..."
            
            # Run bakta with the extracted species name
            echo "\$(date): Annotating bin \${fasta_file} with taxonomy \${organism}..."
                bakta "\${fasta_file}" --db ${DATABASES}/bakta-db/db \\
                      --output "bakta_out/\${SAMPLE}/\${SAMPLE}_Bakta_annotation" \\
                      --threads \$NUM_THREADS \\
                      --species "\${species}" \\
                      --genus "\${genus}"\\
                      --verbose --force
        done
        touch "bakta_out/\${SAMPLE}/\${SAMPLE}_Bakta_annotation/loop_completed.txt"

        conda deactivate
        echo "\$(date): Bakta annotation completed for \${SAMPLE}."
        ;&
    "QUAST_BAKTA")
        echo "\$(date): Starting QUAST on Bakta results..."
        # Run QUAST on Bakta annotated genomes
        for fna_file in \${BAKTA_OUT}/\${SAMPLE}_Bakta_annotation/*fna; do
            GFF3_file="\${fna_file%.fna}.gff3"
            conda activate quast-env
            echo "\$(date): Running QUAST on \$fna_file with GFF file \$GFF3_file..."
            quast \${fna_file} \\
                -o bakta_quast_out/\${SAMPLE}/\${SAMPLE}_Bakta_annotation/\${fna_file}/ \\
                --threads \$NUM_THREADS \\
                --features \${GFF3_file}
            conda deactivate
        done
        echo "\$(date): QUAST analysis completed for Bakta results of \${SAMPLE}."
        ;&
    "ANTISMASH")
        echo "\$(date): Starting antiSMASH on Bakta GBFF outputs..."
        conda activate antismash-env
        antismash --version

        # Run antiSMASH on Bakta GBFF files
        for gbk_file in \${BAKTA_OUT}/\${SAMPLE}_Bakta_annotation/*gbff; do
            if [ -f "\${gbk_file}" ]; then
                echo "\$(date): Running antiSMASH on \${gbk_file}..."
                
                # Extract just the base filename without path
                BASE_FILENAME=\$(basename "\${gbk_file}" .gbff)

                # Define a unique output directory for each gbk_file
                OUTPUT_DIR="\${ANTIMASH_OUT}/\${BASE_FILENAME}"
                mkdir -p "\${OUTPUT_DIR}"  # Create the directory if it doesn't exist
                
                GFF3_file="\${gbk_file%.gbff}.gff3"
                antismash \${gbk_file} \\
                        --output-dir "\${OUTPUT_DIR}" \\
                        --cpus \$NUM_THREADS \\
                        --genefinding-tool prodigal \\
                        --verbose --pfam2go
            else
                echo "\$(date): Warning: GBK file \${gbk_file} not found. Skipping."
            fi
        done

        touch "\${ANTIMASH_OUT}/antismash_loop_completed.txt"

        conda deactivate
        echo "\$(date): antiSMASH analysis completed for Bakta results of \${SAMPLE}."
        ;&
    "HUMANN_READS")
        echo "\$(date): Starting HUMAnN3 reads mode..."
        conda activate humann3-env
        humann --version
        echo "\$(date): Running HUMAnN3 in reads mode for \${SAMPLE}"
        cat \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq > \${KNEADDATA_OUT}/out.\${SAMPLE}_combined.fastq
        humann --input \${KNEADDATA_OUT}/out.\${SAMPLE}_combined.fastq \\
            --output \${HUMANN_READS_DIR} \\
            --output-basename out.\${SAMPLE}_combined \\
            --threads \$NUM_THREADS \\
            --memory-use minimum \\
            --nucleotide-database ${DATABASES}/chocophlan/chocophlan/ \\
            --protein-database ${DATABASES}/uniref90/uniref \\
            --verbose
        conda deactivate
        if [ -f "\${HUMANN_READS_DIR}/out.\${SAMPLE}_combined_humann_temp/out.\${SAMPLE}_combined_bowtie2_aligned.sam" ]; then
            echo "\$(date): Cleaning up intermediate HUMAnN3 files for \${SAMPLE}..."
            rm -rf "\${HUMANN_READS_DIR}/out.\${SAMPLE}_combined_humann_temp/out.\${SAMPLE}_combined_*" 2>/dev/null || true
            echo "\$(date): Intermediate files removed, retaining pathabundance tables"
        else
            echo "\$(date): Error: \${HUMANN_READS_DIR}/out.\${SAMPLE}_combined_humann_temp not found, retaining all HUMAnN3 outputs" >&2
            exit 1
        fi
        ;&
    "HUMANN_CONTIGS")
        echo "\$(date): Starting HUMAnN3 contigs mode..."
        conda activate humann3-env
        echo "\$(date): Running HUMAnN3 in contigs mode on each bin..."
        for fasta_file in \${GOOD_BINS_DIR}/*.fa; do
            if [ -f "\$fasta_file" ]; then
                echo "\$(date): Processing bin: \$fasta_file"
                humann --input "\$fasta_file" \\
                       --output \${HUMANN_CONTIGS_DIR} \\
                       --threads \$NUM_THREADS \\
                       --memory-use minimum \\
                       --nucleotide-database ${DATABASES}/chocophlan/chocophlan/ \\
                       --protein-database ${DATABASES}/uniref90/uniref \\
                       --verbose
            else
                echo "\$(date): Warning: FASTA file \${fasta_file} not found. Skipping."
            fi
        done
        touch "\${HUMANN_CONTIGS_DIR}/humann_loop_completed.txt"
        conda deactivate
        ;&
    "COMBINING")
        echo "\$(date): Combining HUMAnN3 results..."
        touch "\${STATUS_DIR}/\${SAMPLE}.done"

        # Check if this is the last job to finish
        NUM_DONE=\$(ls \${STATUS_DIR}/*.done 2>/dev/null | wc -l)
        NUM_FAILED=\$(ls \${STATUS_DIR}/*.failed 2>/dev/null | wc -l)
        NUM_PROCESSED=\$((\$NUM_DONE + \$NUM_FAILED))

        echo "\$(date): Number of samples processed: \$NUM_PROCESSED out of $N (Done: \$NUM_DONE, Failed: \$NUM_FAILED)"

        if [ "\$NUM_PROCESSED" -eq "$N" ]; then
            LOCK_FILE="${DIRECTORY}/combine_ready.lock"
            if [ ! -f "\$LOCK_FILE" ]; then
                touch "\$LOCK_FILE"
                echo "\$(date): This is the last job. Proceeding with combining results..."
                conda activate humann3-env

                # Combine reads mode
                if [ -d "\${HUMANN_READS_DIR}" ] && [ -n "\$(ls \${HUMANN_READS_DIR}/*_pathabundance.tsv 2>/dev/null)" ]; then
                    humann_join_tables --input "${DIRECTORY}/humann_reads_output/" --output "${DIRECTORY}/humann_reads_output/merged_pathabundance.tsv" --file_name pathabundance
                fi

                # Combine contigs mode
                if [ -d "\${HUMANN_CONTIGS_DIR}" ] && [ -n "\$(ls \${HUMANN_CONTIGS_DIR}/*/*_pathabundance.tsv 2>/dev/null)" ]; then
                    humann_join_tables --input "${DIRECTORY}/humann_contigs_output/" --output "${DIRECTORY}/humann_contigs/merged_pathabundance.tsv" --file_name pathabundance
                fi
                conda deactivate
                echo "\$(date): Combining step completed."
            else
                echo "\$(date): Another job is already combining results. Skipping."
            fi
        else
            echo "\$(date): Waiting for all jobs to finish before combining results."
        fi
        ;&
    "VISUALIZATION")
        echo "\$(date): Starting Kraken2/Bracken visualization..."
        NUM_DONE=\$(ls \${STATUS_DIR}/*.done 2>/dev/null | wc -l)
        NUM_FAILED=\$(ls \${STATUS_DIR}/*.failed 2>/dev/null | wc -l)
        NUM_PROCESSED=\$((\$NUM_DONE + \$NUM_FAILED))
        echo "\$(date): Number of samples processed: \$NUM_PROCESSED out of ${N} (Done: \$NUM_DONE, Failed: \$NUM_FAILED)"

        if [ "\$NUM_PROCESSED" -eq "${N}" ]; then
            VISUALIZATION_LOCK_FILE="${DIRECTORY}/visualization_ready.lock"
            if [ ! -f "\$VISUALIZATION_LOCK_FILE" ]; then
                touch "\$VISUALIZATION_LOCK_FILE"
                echo "\$(date): This is the last job. Proceeding with visualization..."
                
                # Visualization logic here
                # ...existing visualization code...
                # Set custom R library path (same as in R_ANALYSES)
                export R_LIBS="/rds/projects/h/hallly-microbiome/software/R_tools/R_LIB/lxm697_lib:\$R_LIBS"
                echo "\$(date): Using R library path: \$R_LIBS"
                
                # Load the same modules as in R_ANALYSES
                module purge
                module load bluebear
                module load bear-apps/2024a
                module load R-bundle-Bioconductor/3.21-foss-2024a-R-4.5.0
                
                # Create directory for visualizations
                mkdir -p ${DIRECTORY}/visualizations
                
# Create R script for visualization
cat > ${DIRECTORY}/abundance_visualization.R << 'RSCRIPT'
# Combined Bracken and MetaPhlAn Data Processing and Visualization Script
# This script processes both Bracken and MetaPhlAn files, creates relative abundance visualizations
# with original sample names ordered by age from Malawi_Samples_age.csv, and uses a shared color
# palette for species across both plots.
# Load required libraries explicitly
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(htmlwidgets)
library(patchwork)
# Set up for headless environment (no X11 display)
options(bitmapType = 'cairo')
if (!interactive()) {
  # Force cairo device for PNG output on headless systems
  if (capabilities("cairo")) {
    options(device = function(...) cairo_pdf(...))
  }
}
# ============================================================================
# DATA PROCESSING FUNCTIONS
# ============================================================================
# Function to process a single Bracken file
process_bracken_file <- function(file_path) {
  tryCatch({
    # Extract sample name from file path
    sample_name <- gsub("_prokaryote_bracken_abundance\\.txt$", "", basename(file_path))
    sample_name <- gsub("_prok-bracken_abundance\\.txt$", "", sample_name)
   
    # Read the data with more robust parameters
    data <- read.table(
      file_path,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      quote = "", # Disable quote parsing
      fill = TRUE, # Handle uneven columns
      check.names = FALSE # Preserve column names
    )
   
    # Verify expected columns are present
    required_cols <- c("name", "taxonomy_lvl", "fraction_total_reads")
    if (!all(required_cols %in% colnames(data))) {
      warning(sprintf("Missing required columns in file: %s", basename(file_path)))
      return(NULL)
    }
   
    # Filter for species-level data and remove NA values
    data <- data %>%
      filter(
        taxonomy_lvl == "S",
        !is.na(fraction_total_reads)
      )
   
    # Add sample column
    data$sample <- sample_name
   
    # Reorder columns to have sample first
    data <- data[, c("sample", names(data)[names(data) != "sample"])]
   
    return(data)
   
  }, error = function(e) {
    warning(sprintf("Error processing file %s: %s", basename(file_path), e$message))
    return(NULL)
  })
}
# Function to process a single MetaPhlAn profile file
process_metaphlan_file <- function(file_path) {
  tryCatch({
    # Extract sample name from file path
    sample_name <- gsub("_metaphlan_profile\\.txt$", "", basename(file_path))
   
    # Read the file lines to find where the data starts
    lines <- readLines(file_path)
    # Find the header line that contains column names
    header_line <- grep("^#clade_name", lines, value = TRUE)[1]
    # Find where the actual data starts (first taxonomic entry)
    data_start <- which(grepl("^k__", lines))[1]
   
    if (is.na(data_start)) {
      warning(paste("No taxonomic data found in:", basename(file_path)))
      return(NULL)
    }
   
    # Read only the data portion of the file
    data <- read.table(
      text = lines[data_start:length(lines)],
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE,
      quote = "", # Disable quote parsing
      fill = TRUE, # Handle uneven columns
      comment.char = "" # Don't treat # as comments
    )
   
    # Set column names based on the actual MetaPhlAn output format
    colnames(data) <- c("clade_name", "NCBI_tax_id", "relative_abundance", "additional_species")[1:ncol(data)]
   
    # Filter for species-level taxa and process the names
    species_data <- data %>%
      filter(grepl("\\|s__[^|]+$", clade_name)) %>% # Only get species-level entries
      mutate(
        sample = sample_name,
        # Extract species name from the full clade name
        name_clean = sub(".*\\|s__([^|]+)$", "\\1", clade_name),
        # Clean up the name
        name_clean = gsub("_", " ", name_clean),
        # Convert abundance to fraction (handling potential NA values)
        fraction_total_reads = as.numeric(relative_abundance) / 100
      ) %>%
      filter(!is.na(fraction_total_reads)) %>%
      select(sample, name_clean, fraction_total_reads)
   
    # Check if we have any valid species data
    if (nrow(species_data) == 0) {
      warning(paste("No valid species-level data found in:", basename(file_path)))
      return(NULL)
    }
   
    return(species_data)
   
  }, error = function(e) {
    warning(paste("Error processing file:", basename(file_path), "-", e$message))
    return(NULL)
  })
}
# ============================================================================
# COLOR PALETTE FUNCTION
# ============================================================================
# Function to create unified color palette with multiple options
create_unified_color_palette <- function(data, top_n = 20, palette_option = "qualitative") {
  # Get top species by total abundance
  top_species <- data %>%
    group_by(name_clean) %>%
    summarise(total_abundance = sum(fraction_total_reads)) %>%
    arrange(desc(total_abundance)) %>%
    head(top_n) %>%
    pull(name_clean)
  # Number of colors needed
  n_colors <- length(top_species)
  cat(sprintf("Selected %d top species for palette\n", n_colors))
  # Print the top species for debugging
  cat("Top species:\n")
  print(top_species)
  # Define different palette options
  palette <- switch(palette_option,
                    "default" = {
                      base_colors <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Accent"))
                      colorRampPalette(base_colors)(n_colors)
                    },
                    "vibrant" = {
                      base_colors <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Dark2"))
                      colorRampPalette(base_colors)(n_colors)
                    },
                    "pastel" = {
                      base_colors <- c(brewer.pal(8, "Pastel1"), brewer.pal(8, "Pastel2"))
                      colorRampPalette(base_colors)(n_colors)
                    },
                    "rainbow" = {
                      colorRampPalette(rainbow(10))(n_colors)
                    },
                    "qualitative" = {
                      base_colors <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Set3"))
                      colorRampPalette(base_colors)(n_colors)
                    },
                    "spectral" = {
                      colorRampPalette(rev(brewer.pal(11, "Spectral")))(n_colors)
                    },
                    {
                      base_colors <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Accent"))
                      colorRampPalette(base_colors)(n_colors)
                    }
  )
  # Add gray for "Other" category
  palette <- c(palette, "gray70")
  names(palette) <- c(top_species, "Other")
  return(list(palette = palette, top_species = top_species))
}
# Function to create normalized abundance plots
create_abundance_plot <- function(processed_data, top_species, color_palette, title_prefix, interactive = FALSE, age_data) {
  top_n <- min(20, length(top_species))
  top_categories <- top_species[1:top_n]
  all_samples <- unique(processed_data$sample)
  
  # Get actual abundances for top species (sum if multiple entries, but unlikely)
  actual_tops <- processed_data %>%
    filter(name_clean %in% top_categories) %>%
    group_by(sample, name_clean) %>%
    summarise(abundance = sum(fraction_total_reads), .groups = "drop") %>%
    rename(species_category = name_clean)
  
  # Create full grid for top categories, filling missing with 0
  top_grid <- expand_grid(sample = all_samples, species_category = top_categories) %>%
    left_join(actual_tops, by = c("sample", "species_category")) %>%
    mutate(abundance = ifelse(is.na(abundance), 0, abundance))
  
  # Calculate sum of top abundances per sample
  sum_tops <- top_grid %>%
    group_by(sample) %>%
    summarise(sum_tops = sum(abundance), .groups = "drop")
  
  # Create Other data
  other_data <- expand_grid(sample = all_samples) %>%
    left_join(sum_tops, by = "sample") %>%
    mutate(
      species_category = "Other",
      abundance = pmax(1.0 - sum_tops, 0)
    )
  
  # Combine all data
  plot_data <- bind_rows(top_grid, other_data) %>%
    left_join(age_data, by = "sample") %>%
    mutate(
      sample = factor(sample, levels = age_data$sample[order(age_data$age)]),
      species_category = factor(species_category, levels = c(top_categories, "Other"))
    )
  
  # Debug print
  cat(sprintf("For %s plot: %d unique species categories (including Other)\n", title_prefix, nlevels(plot_data$species_category)))
  cat(sprintf("  Top species count: %d\n", top_n))
  
  # Create the plot with normalized values
  p <- ggplot(plot_data, aes(x = sample, y = abundance, fill = species_category)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = color_palette, name = "Species", drop = FALSE) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
      title = paste(title_prefix, ": Relative Abundance of Top", top_n, "Species Across Samples"),
      subtitle = "Remaining species are grouped as 'Other', samples ordered by age (months)",
      x = "Sample",
      y = "Relative Abundance (%)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 1, drop = FALSE))
  if (interactive) {
    # Ensure plotly is loaded before calling ggplotly
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("plotly package is not installed. Skipping interactive plot.")
      return(p)
    }
    p <- tryCatch({
      ggplotly(p, tooltip = c("x", "y", "fill"))
    }, error = function(e) {
      warning("Error in ggplotly: ", e$message, ". Returning static plot.")
      return(p)
    })
  }
  return(p)
}
# ============================================================================
# MAIN PROCESSING PIPELINE
# ============================================================================
cat("Starting combined Bracken and MetaPhlAn data processing...\n")
# --- 1. PROCESS BRACKEN DATA ---
cat("\n=== 1. PROCESSING BRACKEN DATA ===\n")
bracken_files <- list.files(path = "kraken2_bracken_out", pattern = "*_prokaryote_bracken_abundance\\.txt$|*_prok-bracken_abundance\\.txt$", full.names = TRUE)
bracken_processed <- NULL
if (length(bracken_files) > 0) {
  cat("Found", length(bracken_files), "Bracken files.\n")
  all_bracken_data <- map_dfr(bracken_files, process_bracken_file)
  bracken_processed <- all_bracken_data %>%
    mutate(name_clean = gsub("_", " ", gsub("s__", "", sub(".*\\|", "", name)))) %>%
    filter(!is.na(fraction_total_reads) & name_clean != "") %>%
    group_by(sample) %>%
    mutate(fraction_total_reads = fraction_total_reads / sum(fraction_total_reads, na.rm = TRUE)) %>%
    ungroup()
  cat("- Bracken data processed and normalized to 100% per sample.\n")
  cat(sprintf("- Unique species in Bracken: %d\n", length(unique(bracken_processed$name_clean))))
} else {
  cat("No Bracken files found.\n")
}
# --- 2. PROCESS METAPHLAN DATA ---
cat("\n=== 2. PROCESSING METAPHLAN DATA ===\n")
metaphlan_files <- list.files(path = "metaphlan_out", pattern = "*_metaphlan_profile\\.txt$", full.names = TRUE)
metaphlan_processed <- NULL
if (length(metaphlan_files) > 0) {
  cat("Found", length(metaphlan_files), "MetaPhlAn files.\n")
  metaphlan_processed <- map_dfr(metaphlan_files, process_metaphlan_file) %>%
    filter(!is.na(fraction_total_reads) & name_clean != "") %>%
    group_by(sample) %>%
    mutate(fraction_total_reads = fraction_total_reads / sum(fraction_total_reads, na.rm = TRUE)) %>%
    ungroup()
  cat("- MetaPhlAn data processed and normalized to 100% per sample. Skipped files (if any) will be noted in warnings above.\n")
  cat(sprintf("- Unique species in MetaPhlAn: %d\n", length(unique(metaphlan_processed$name_clean))))
} else {
  cat("No MetaPhlAn files found.\n")
}
# --- 3. LOAD AGE DATA AND PREPARE COLOR MAPPING ---
cat("\n=== 3. LOADING AGE DATA AND CREATING COLOR MAPPING ===\n")
# Safely combine sample names and species data
all_sample_names <- unique(c(
  if (!is.null(bracken_processed)) bracken_processed$sample else NULL,
  if (!is.null(metaphlan_processed)) metaphlan_processed$sample else NULL
))
# Load age data
age_data <- read.csv("Malawi_Samples_age.csv", stringsAsFactors = FALSE)
colnames(age_data) <- c("sample", "age")
age_data <- age_data %>%
  filter(sample %in% all_sample_names) %>%
  arrange(age)
if (nrow(age_data) > 0) {
  cat("- Age data loaded for", nrow(age_data), "samples from 'Malawi_Samples_age.csv'.\n")
} else {
  cat("- No valid age data found.\n")
  stop("Age data is required for ordering samples.")
}
# Create a single, shared color palette for all species
data_to_combine <- list()
if (!is.null(bracken_processed)) {
  data_to_combine$bracken <- select(bracken_processed, name_clean, fraction_total_reads)
}
if (!is.null(metaphlan_processed)) {
  data_to_combine$metaphlan <- select(metaphlan_processed, name_clean, fraction_total_reads)
}
all_species_data <- bind_rows(data_to_combine)
if (nrow(all_species_data) > 0) {
  # Filter to only species with positive total abundance
  all_species_data <- all_species_data %>%
    group_by(name_clean) %>%
    filter(sum(fraction_total_reads) > 0) %>%
    ungroup()
  cat(sprintf("- Total unique species with >0 abundance across all data: %d\n", length(unique(all_species_data$name_clean))))
  color_info <- create_unified_color_palette(all_species_data, top_n = 20, palette_option = "qualitative")
  unified_color_palette <- color_info$palette
  top_species_list <- color_info$top_species
} else {
  cat("- No species data to generate a color palette.\n")
  stop("No species data available for visualization.")
}
# ============================================================================
# CREATE AND SAVE VISUALIZATIONS
# ============================================================================
cat("\n=== 4. CREATING VISUALIZATIONS ===\n")
if (!is.null(bracken_processed) && !is.null(metaphlan_processed) &&
    nrow(bracken_processed) > 0 && nrow(metaphlan_processed) > 0) {
  cat("Creating abundance plots...\n")
  # Create static plots
  bracken_plot <- create_abundance_plot(
    processed_data = bracken_processed,
    top_species = top_species_list,
    color_palette = unified_color_palette,
    title_prefix = "Bracken",
    age_data = age_data
  )
 
  metaphlan_plot <- create_abundance_plot(
    processed_data = metaphlan_processed,
    top_species = top_species_list,
    color_palette = unified_color_palette,
    title_prefix = "MetaPhlAn",
    age_data = age_data
  )
 
  # Save separate static PNGs
  png("bracken_relative_abundance_stacked_bar.png",
      width = 16, height = 9,
      units = "in", res = 300, type = "cairo")
  print(bracken_plot)
  dev.off()
 
  png("metaphlan_relative_abundance_stacked_bar.png",
      width = 16, height = 9,
      units = "in", res = 300, type = "cairo")
  print(metaphlan_plot)
  dev.off()
 
  # Create and save separate interactive plots
  bracken_interactive <- create_abundance_plot(
    processed_data = bracken_processed,
    top_species = top_species_list,
    color_palette = unified_color_palette,
    title_prefix = "Bracken",
    interactive = TRUE,
    age_data = age_data
  )
  saveWidget(bracken_interactive, "bracken_abundance_interactive.html",
             selfcontained = TRUE)
 
  metaphlan_interactive <- create_abundance_plot(
    processed_data = metaphlan_processed,
    top_species = top_species_list,
    color_palette = unified_color_palette,
    title_prefix = "MetaPhlAn",
    interactive = TRUE,
    age_data = age_data
  )
  saveWidget(metaphlan_interactive, "metaphlan_abundance_interactive.html",
             selfcontained = TRUE)
 
  cat("- Separate plots saved: .png and .html files for Bracken and MetaPhlAn\n")
} else {
  cat("- Insufficient data to create combined plots\n")
  # Create individual plots if available
  if (!is.null(bracken_processed) && nrow(bracken_processed) > 0) {
    bracken_plot <- create_abundance_plot(
      processed_data = bracken_processed,
      top_species = top_species_list,
      color_palette = unified_color_palette,
      title_prefix = "Bracken",
      age_data = age_data
    )
    png("bracken_relative_abundance_stacked_bar.png",
        width = 16, height = 9,
        units = "in", res = 300, type = "cairo")
    print(bracken_plot)
    dev.off()
    cat("- Bracken plot saved as 'bracken_relative_abundance_stacked_bar.png'\n")
  }
 
  if (!is.null(metaphlan_processed) && nrow(metaphlan_processed) > 0) {
    metaphlan_plot <- create_abundance_plot(
      processed_data = metaphlan_processed,
      top_species = top_species_list,
      color_palette = unified_color_palette,
      title_prefix = "MetaPhlAn",
      age_data = age_data
    )
    png("metaphlan_relative_abundance_stacked_bar.png",
        width = 16, height = 9,
        units = "in", res = 300, type = "cairo")
    print(metaphlan_plot)
    dev.off()
    cat("- MetaPhlAn plot saved as 'metaphlan_relative_abundance_stacked_bar.png'\n")
  }
}
cat("\nScript execution completed successfully!\n")
RSCRIPT

                # Replace _DIRECTORY_ placeholder with actual directory path
                sed -i "s|_DIRECTORY_|${DIRECTORY}|g" ${DIRECTORY}/abundance_visualization.R
                
                # Run the R script
                echo "\$(date): Running R visualization script..."
                Rscript ${DIRECTORY}/abundance_visualization.R
                
                # Create marker file to indicate visualization is complete
                touch "\${STATUS_DIR}/abundance_visualization_complete.txt"
                
                # Remove the lock file
                rm -f "\$VISUALIZATION_LOCK_FILE"
                echo "\$(date): Kraken2/Bracken visualization completed."
                echo "\$(date): Cleaning up temporary files..."
                rm -f "\${KRAKEN2_OUT}/\${SAMPLE}.kraken2_output.txt" 2>/dev/null || true
            fi
        else
            echo "\$(date): Waiting for all ${N} samples to be processed (currently \$NUM_PROCESSED done)."
        fi
        ;&
    "R_ANALYSES")
        echo "\$(date): Starting R analyses..."
        # Set custom R library path
        export R_LIBS="/rds/projects/h/hallly-microbiome/software/R_tools/R_LIB/lxm697_lib:$R_LIBS"
        echo "\$(date): Using R library path: \$R_LIBS"
        module purge
        module load bluebear
        module load bear-apps/2023a
        module load GCC/12.3.0
        module load R/4.4.1-gfbf-2023a
        mkdir -p ${DIRECTORY}/maaslin3_reads ${DIRECTORY}/maaslin3_contigs ${DIRECTORY}/melonnpan_reads ${DIRECTORY}/melonnpan_contigs

        # Set custom R library path
        export R_LIBS="/rds/projects/h/hallly-microbiome/software/R_tools/R_LIB/lxm697_lib:$R_LIBS"
        echo "\$(date): Using R library path: \$R_LIBS"

        # Run MelonnPan and MaAsLin 3 using the custom library
        Rscript -e "cat('Using library paths:', .libPaths(), '\n'); library(melonnpan); if (!'melonnpan' %in% rownames(installed.packages())) { stop('melonnpan package not found in library path') }; melonnpan.predict(metag='${DIRECTORY}/humann_reads_output/merged_genefamilies.tsv', output='${DIRECTORY}/melonnpan_reads')" || echo "MelonnPan for reads skipped due to error"
        
        Rscript -e "cat('Using library paths:', .libPaths(), '\n'); library(melonnpan); if (!'MelonnPan' %in% rownames(installed.packages())) { stop('MelonnPan package not found in library path') }; library(MelonnPan); melonnpan(input_data='${DIRECTORY}/humann_contigs/merged_pathabundance.tsv', output_dir='${DIRECTORY}/melonnpan_contigs', threads=$NUM_THREADS)" || echo "MelonnPan for contigs skipped due to error"

        Rscript -e "cat('Using library paths:', .libPaths(), '\n'); if (!'Maaslin3' %in% rownames(installed.packages())) { stop('Maaslin3 package not found in library path') }; library(Maaslin3); Maaslin3(input_data='${DIRECTORY}/humann_reads/merged_pathabundance.tsv', output='${DIRECTORY}/maaslin3_reads', metadata='${DIRECTORY}/metadata.tsv', analysis_method='LM')" || echo "MaAsLin 3 for reads skipped due to error"

        Rscript -e "cat('Using library paths:', .libPaths(), '\n'); if (!'Maaslin3' %in% rownames(installed.packages())) { stop('Maaslin3 package not found in library path') }; library(Maaslin3); Maaslin3(input_data='${DIRECTORY}/humann_contigs/merged_pathabundance.tsv', output='${DIRECTORY}/maaslin3_contigs', metadata='${DIRECTORY}/metadata.tsv', analysis_method='LM')" || echo "MaAsLin 3 for contigs skipped due to error"
        
        # Removing the dehumnanized reads
        #echo "\$(date): Cleaning up Kneaddata outputs for \$SAMPLE"
        #rm -rf "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_unmatched_1.fastq" \\
        #       "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_unmatched_2.fastq" 2>/dev/null || true

        ;;
    *)
        echo "\$(date): Error: Invalid RESUME_FROM value: '\$RESUME_FROM'" >&2
        exit 1
        ;;
esac

echo "\$(date): Workflow completed for \$SAMPLE."
echo "\$(date): Workflow summary - Completed: \$NUM_DONE, Failed: \$NUM_FAILED, Total: $N"
EOF

# Submit the SLURM job and capture the job ID
JOB_ID=$(sbatch metagenomics_workflow.sh | awk '/Submitted batch job/ {print $4}')  # Ensure job ID is captured correctly

if [[ -z "$JOB_ID" ]]; then
    echo "$(date): Error: Failed to submit SLURM job. Exiting."
    exit 1
fi

# Wait for the job to start
echo "$(date): Waiting for job $JOB_ID to start..."
while true; do
    # Check the job state using squeue
    JOB_STATE=$(squeue --job ${JOB_ID}_1 --noheader --format=%T)
    
    # If the job is running, break the loop
    if [[ "$JOB_STATE" == "RUNNING" ]]; then
        echo -e "\n$(date): Job ${JOB_ID}_1 has started."
        break
    elif [[ -z "$JOB_STATE" ]]; then
        echo -e "\n$(date): Job $JOB_ID is no longer in the queue. Exiting."
        exit 1
    fi

    # Use scontrol to get the scheduled start time
    SCHEDULED_START=$(scontrol show job $JOB_ID | grep -oP 'StartTime=\K[^ ]+')
    if [[ "$SCHEDULED_START" == "Unknown" ]]; then
        echo -ne "$(date): Job $JOB_ID is waiting to be scheduled...\r"
    else
        echo -ne "$(date): Job $JOB_ID is scheduled to start at $SCHEDULED_START.\r"
    fi

    sleep 5  # Check every 5 seconds
done

# Combine .out, .err, and .stats files into a single log file
LOG_FILE="${DIRECTORY}/first-job-run.log"
echo "$(date): Monitoring logs for job $JOB_ID..." | tee -a "$LOG_FILE"

# Monitor all relevant SLURM log files and append to the combined log file
(
    tail -f slurm-${JOB_ID}.stats slurm-${JOB_ID}_1.out slurm-${JOB_ID}_1.err 2>/dev/null | while read -r line; do
        echo "$line" | tee -a "$LOG_FILE"
    done
) || echo "$(date): Error: Log files not found for job $JOB_ID." | tee -a "$LOG_FILE"

# After the first job completes, switch to monitoring the SLURM queue
echo "$(date): Job $JOB_ID completed. Switching to monitoring the SLURM queue..."

# Monitor the SLURM queue until the final job in the array is submitted
while true; do
    # Get the list of jobs in the queue for the current user
    JOB_LIST=$(squeue --user=$USER --noheader --format="%.18i %.2t" | grep "${JOB_ID}_" | awk '{print $1, $2}')
    
    # Check if the final job in the array is in the queue
    FINAL_JOB="${JOB_ID}_${N}"  # Assuming $N is the total number of jobs in the array
    FINAL_JOB_STATE=$(echo "$JOB_LIST" | grep "$FINAL_JOB" | awk '{print $2}')

    if [[ "$FINAL_JOB_STATE" == "R" ]]; then
        echo "$(date): Final job $FINAL_JOB is running. Switching to monitor its logs..."
        break
    elif [[ -z "$FINAL_JOB_STATE" ]]; then
        echo "$(date): Final job $FINAL_JOB is not yet in the queue. Waiting..."
    else
        echo "$(date): Final job $FINAL_JOB is in state $FINAL_JOB_STATE. Waiting for it to start..."
    fi
    
    sleep 5  # Check every 5 seconds
done

# Monitor the logs of the final job
(
    tail -f slurm-${FINAL_JOB}.out slurm-${FINAL_JOB}.err 2>/dev/null | while read -r line; do
        echo "$line" | tee -a "$LOG_FILE"
    done
) || echo "$(date): Error: Log files not found for final job $FINAL_JOB." | tee -a "$LOG_FILE"