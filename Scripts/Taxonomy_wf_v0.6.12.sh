#!/bin/bash

# Check if at least a file pattern is provided as an argument
if [ -z "$1" ]; then
    echo "Usage: $0 <file_extension> [working_directory]"
    echo "Example: $0 '.fastq.gz' [/path/to/directory]"
    echo "If no directory is specified, uses current working directory"
    exit 1
fi

# Assign the provided extension (e.g., '.fastq.gz') to a variable
EXTENSION="$1"

# Set directory to either provided argument or current directory, removing trailing slash
if [ -n "$2" ]; then
    DIRECTORY="$2"
    DIRECTORY="${DIRECTORY%/}" # Remove trailing slash if present
else
    DIRECTORY="$(pwd)" # Use current directory if none is provided
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

# Set concurrency
# Limit the number of concurrent jobs to 10 or the total number of samples, whichever is smaller
if [ "$N" -le 15 ]; then 
    CONCURRENT=$N
else
    CONCURRENT=15
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
#SBATCH --time=55:00:00
#SBATCH --nodes=1
#SBATCH --mem=350G
#SBATCH --ntasks=50
#SBATCH --mail-type=all
#SBATCH --mail-user=lxm697@student.bham.ac.uk
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

# Check if the sample is already completed or failed
if [ -f "\${STATUS_DIR}/\${SAMPLE}.done" ]; then
    echo "\$(date): Sample \${SAMPLE} is already completed. Skipping..."
    exit 0
elif [ -f "\${STATUS_DIR}/\${SAMPLE}.failed" ]; then
    echo "\$(date): Sample \${SAMPLE} previously failed. Proceeding with the workflow..."
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

COMPLETENESS=50
CONTAMINATION=10

# Define directories for workflow outputs
# These directories will store intermediate and final results
QC_READS="${DIRECTORY}/qc_reads"
DEHUMAN_READS="${DIRECTORY}/dehuman_reads"
ASSEMBLY_DIR="${DIRECTORY}/megahit_assembly-test/\${SAMPLE}assembly"
BINNING_DIR="${DIRECTORY}/metawrap-out/metawrap_binref-checkm-\${SAMPLE}-out"
QUAST_ASSEMBLY_OUT="quast_assembly/\${SAMPLE}/\${SAMPLE}_results/"
QUAST_BINS_OUT="quast_bins/\${SAMPLE}/\${SAMPLE}_results/"
KRAKEN2_OUT="${DIRECTORY}/kraken2_bracken_out"
GTDBTK_OUT="${DIRECTORY}/gtdbtk-out/gtdbtk-\${SAMPLE}-out"
STATUS_DIR="${DIRECTORY}/status"
KNEADDATA_OUT="\${DEHUMAN_READS}/kneaddata_output_\${SAMPLE}"
ASSEMBLY_FILE="\${ASSEMBLY_DIR}/final.contigs.fa"
BINNING_OUT="${DIRECTORY}/metawrap-out/metawrap-checkm-\${SAMPLE}-out"
BIN_DIR="\${BINNING_DIR}/metawrap_\${COMPLETENESS}_\${CONTAMINATION}_bins"
REASSEMBLED_BINS_DIR="\${BINNING_DIR}/reassembled_bins"
METAQUAST_REFS_FILE="\${REASSEMBLED_BINS_DIR}/metaquast_references.txt"

# Create necessary directories for the workflow
mkdir -p \${DEHUMAN_READS} \${ASSEMBLY_DIR} \${QUAST_ASSEMBLY_OUT} \${BINNING_DIR} \\
         \${QUAST_BINS_OUT} \${KRAKEN2_OUT} \${GTDBTK_OUT} \${STATUS_DIR} \${QC_READS}

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
        "FASTP:qc_reads/\${SAMPLE}_fastp.html"
        "KNEADDATA:dehuman_reads/kneaddata_output_\${SAMPLE}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq"
        "ASSEMBLY:megahit_assembly-test/\${SAMPLE}assembly/final.contigs.fa"
        "QUAST_ASSEMBLY:megahit_assembly-test/\${SAMPLE}assembly/quast_results/report.html"
        "BINNING:metawrap-out/metawrap_binref-checkm-\${SAMPLE}-out/reassembled_bins/reassembled_bins/*.fa"
        "TAXONOMY_GTDBTK:gtdbtk-out/gtdbtk-\${SAMPLE}-out/gtdbtk.bac120.summary.tsv"
        "TAXONOMY_KRAKEN2:kraken2_bracken_out/\${SAMPLE}_bracken_abundance.txt"
        "QUAST_BINS:quast_bins/\${SAMPLE}/\${SAMPLE}_results/report.html"
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
            "QUAST_BINS")
                echo "\$(date): All steps completed for \$SAMPLE. No further action needed."
                touch "\${STATUS_DIR}/\${SAMPLE}.done"
                RESUME_FROM=""
                ;;
            "TAXONOMY_KRAKEN2") RESUME_FROM="QUAST_BINS" ;;
            "TAXONOMY_GTDBTK") RESUME_FROM="TAXONOMY_KRAKEN2" ;;
            "BINNING") RESUME_FROM="TAXONOMY_GTDBTK" ;;
            "QUAST_ASSEMBLY") RESUME_FROM="BINNING" ;;
            "ASSEMBLY") RESUME_FROM="QUAST_ASSEMBLY" ;;
            "KNEADDATA") RESUME_FROM="ASSEMBLY" ;;
            "FASTP") RESUME_FROM="KNEADDATA" ;;
            *) RESUME_FROM="FASTP" ;;  # Default if none matched
        esac
    else
        # If no checkpoint was found, start from the beginning
        RESUME_FROM="FASTP"
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
              # Cleanup for fastp
            if [ -f "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq" ] && [ -f "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq" ]; then
                rm -f "\${QC_READS}/out.\${SAMPLE}.R1.fq.gz" "\${QC_READS}/out.\${SAMPLE}.R2.fq.gz"
            fi
        else
            echo "\$(date): kneaddata output already exists, skipping..."
            if [ -f "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq" ] && [ -f "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq" ]; then
                rm -f "\${QC_READS}/out.\${SAMPLE}.R1.fq.gz" "\${QC_READS}/out.\${SAMPLE}.R2.fq.gz"
            fi
        fi

        echo "\$(date): Cleaning up Kneaddata outputs for \$SAMPLE"
        # Define the correct kneaddata log file path
        KNEADDATA_LOG="\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata.log"
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
                   "\${ASSEMBLY_DIR}/k*.gfa" \\
                   "\${ASSEMBLY_DIR}/k*.fa" \\
                   "\${ASSEMBLY_DIR}/tmp.*" \\
                   "\${ASSEMBLY_DIR}/*.bin" 2>/dev/null || true
            echo "\$(date): Intermediate files removed, retaining final.contigs.fa and log/opts.txt"
        else
            echo "\$(date): Error: final.contigs.fa not found, retaining all MEGAHIT outputs" >&2
            exit 1
        fi
        
        echo "\$(date): Cleaning up Kneaddata outputs for \$SAMPLE"
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
            -o quast_assembly/\${SAMPLE}/\${SAMPLE}_results/ \\
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
            -A \${BINNING_OUT}/concoct_bins/ \\
            -B \${BINNING_OUT}/maxbin2_bins/ \\
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
            echo "\$(date): Error: \${BIN_DIR}/reassembled_bins/bin\* not found, retaining all metaWRAP outputs" >&2
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
        echo "# Generated from GTDB-Tk classifications on $(date)" >> \${METAQUAST_REFS_FILE}

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
    "TAXONOMY_KRAKEN2")
        echo "\$(date): Starting Kraken2 and Bracken..."
        kraken2 --version
        bracken -v
        echo "\$(date): Running Kraken2 and Bracken on \${SAMPLE} reads (Bacteria and Archaea)"
        kraken2 \\
            --db ${DATABASES}/krakendb/standard/ \\
            --paired --threads \$NUM_THREADS \\
            --report \${KRAKEN2_OUT}/\${SAMPLE}.prokaryote.kraken2_report.txt \\
            --output \${KRAKEN2_OUT}/\${SAMPLE}.prokaryote.kraken2_output.txt \\
            \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq
        bracken \\
            -d ${DATABASES}/krakendb/standard/ \\
            -i \${KRAKEN2_OUT}/\${SAMPLE}.prokaryote.kraken2_report.txt \\
            -o \${KRAKEN2_OUT}/\${SAMPLE}_prokaryote_bracken_abundance.txt \\
            -r 100
        
        kraken2 \\
            --db ${DATABASES}/krakendb/Eukaryote/ \\
            --paired --threads \$NUM_THREADS \\
            --report \${KRAKEN2_OUT}/\${SAMPLE}.eukaryote.kraken2_report.txt \\
            --output \${KRAKEN2_OUT}/\${SAMPLE}.eukaryote.kraken2_output.txt \\
            \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq
        bracken \\
            -d ${DATABASES}/krakendb/Eukaryote/ \\
            -i \${KRAKEN2_OUT}/\${SAMPLE}.eukaryote.kraken2_report.txt \\
            -o \${KRAKEN2_OUT}/\${SAMPLE}_eukaryote_bracken_abundance.txt \\
            -r 100

        kraken2 \\
            --db ${DATABASES}/krakendb/Archaea/ \\
            --paired --threads \$NUM_THREADS \\
            --report \${KRAKEN2_OUT}/\${SAMPLE}.archaea.kraken2_report.txt \\
            --output \${KRAKEN2_OUT}/\${SAMPLE}.archaea.kraken2_output.txt \\
            \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq \${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq
        bracken \\
            -d ${DATABASES}/krakendb/Archaea/ \\
            -i \${KRAKEN2_OUT}/\${SAMPLE}.archaea.kraken2_report.txt \\
            -o \${KRAKEN2_OUT}/\${SAMPLE}_archaea_bracken_abundance.txt \\
            -r 100

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

        # Removing the dehumnanized reads
        echo "\$(date): Cleaning up Kneaddata outputs for \$SAMPLE"
        rm -rf "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_1.fastq" \\
               "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_paired_2.fastq" \\
               "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_unmatched_1.fastq" \\
               "\${KNEADDATA_OUT}/out.\${SAMPLE}.R1_kneaddata_unmatched_2.fastq" 2>/dev/null || true
        
        # Removing the MEGAHIT assembly
        echo "\$(date): Cleaning up MEGAHIT assembly outputs for \$SAMPLE"
        rm -rf "\${ASSEMBLY_DIR}/final.contigs.fa" 2>/dev/null || true
        ;;
esac

echo "\$(date): Workflow completed for \$SAMPLE."
touch "\${STATUS_DIR}/\${SAMPLE}.done"

NUM_DONE=\$(find \${STATUS_DIR} -name "*.done" 2>/dev/null | wc -l || echo 0)
NUM_FAILED=\$(find \${STATUS_DIR} -name "*.failed" 2>/dev/null | wc -l || echo 0)
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