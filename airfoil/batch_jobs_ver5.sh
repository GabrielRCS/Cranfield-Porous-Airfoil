#!/bin/bash
##
## Batch FLUENT submission script for PBS on DELTA
## ------------------------------------------------
## Runs multiple CFD jobs on a single node (2 cores each)
##
## STEP 1: Enter a job name after the -N on the line below:
##
#PBS -N palabos_batch_AoA_sweep
##
## STEP 2: Select resources (adjust based on number of jobs)
##
#PBS -l select=1:ncpus=12:mpiprocs=12
##
## STEP 3: Select the correct queue
##
#PBS -q half_day
##
## STEP 4: Replace with your Cranfield email address
##
#PBS -m abe
#PBS -M gabriel.rochette.306@cranfield.ac.uk
##
## ====================================
## DO NOT CHANGE THE LINES BETWEEN HERE
## ====================================
#PBS -l application=palabos
#PBS -j oe
#PBS -W sandbox=PRIVATE
#PBS -k n
ln -s $PWD $PBS_O_WORKDIR/$PBS_JOBID
## Change to working directory
cd $PBS_O_WORKDIR
## ========
## AND HERE
## ========

## STEP 5: Load modules
module use /gpfs/apps/modules/all/
module load ANSYS/2024R1

## Configuration - MODIFY THESE FOR YOUR PROJECT
PROJECT_DIR="porous_airfoil_project/Porous_NACA0012" # To change to my local path
CORES_PER_JOB=2

# List of subdirectories to process (one job per directory)
CASE_DIRS=("AoA0" "AoA10" "AoA20" "AoA28" "AoA4") # Different names for the cases

LOG_DIR="${PBS_O_WORKDIR}/fluent_logs_${PBS_JOBID}"
mkdir -p "$LOG_DIR"

echo "========================================" 
echo "Batch Palabos Job Started at $(date)"
echo "Job ID: $PBS_JOBID"
echo "Node: $(hostname)"
echo "Working directory: $PBS_O_WORKDIR"
echo "Project directory: $PROJECT_DIR"
echo "========================================"

# Function to run a single Fluent job
run_fluent_job() {
    case_dir=$1
    job_index=$2
    
    # Full path to job directory
    job_dir="${PBS_O_WORKDIR}/${PROJECT_DIR}/${case_dir}"
    
    # Check if job directory exists
    if [ ! -d "$job_dir" ]; then
        echo "ERROR: Directory $job_dir does not exist!" >> "${LOG_DIR}/job_${case_dir}.log"
        return 1
    fi
    
    echo "========================================" >> "${LOG_DIR}/job_${case_dir}.log"
    echo "Starting Palabos job for ${case_dir} at $(date)" >> "${LOG_DIR}/job_${case_dir}.log"
    echo "Working in directory: $job_dir" >> "${LOG_DIR}/job_${case_dir}.log"
    echo "========================================" >> "${LOG_DIR}/job_${case_dir}.log"
    
    # This is not needed for Palabos: we stay in the main directory.

    # # Change to job-specific directory
    # cd "$job_dir" || {
    #     echo "ERROR: Failed to cd to $job_dir" >> "${LOG_DIR}/job_${case_dir}.log"
    #     return 1
    # }
    
    ###### MODIFY BELOW FOR PALABOS JOB SETUP ######
    # Find the .cas.h5 file (assumes one case file per directory)
    case_file=$(ls -1 *.cas.h5 2>/dev/null | head -1)
    
    if [ -z "$case_file" ]; then
        echo "ERROR: No .cas.h5 file found in $job_dir" >> "${LOG_DIR}/job_${case_dir}.log"
        cd "$PBS_O_WORKDIR"
        return 1
    fi
    ################################################
    
    # Check if data file exists
    data_file="${case_file%.cas.h5}.dat.h5"
    
    echo "Case file: $case_file" >> "${LOG_DIR}/job_${case_dir}.log"
    if [ -f "$data_file" ]; then
        echo "Data file: $data_file (will continue from previous solution)" >> "${LOG_DIR}/job_${case_dir}.log"
        read_data_option="-i"
    else
        echo "No data file found - starting from initialization" >> "${LOG_DIR}/job_${case_dir}.log"
        read_data_option=""
    fi
    
    # Create a temporary node file for this job
    job_nodefile="${job_dir}/fluent_nodes_${case_dir}.$$"
    head -n $CORES_PER_JOB $PBS_NODEFILE | sort -u > "$job_nodefile"
    
    # Create a journal file for this run
    journal_file="run_${case_dir}_${PBS_JOBID}.jou"
    cat > "$journal_file" << 'EOF'
; Automated journal file for transient simulation
; Read case and data files
/file/read-case-data
EOF
    echo "${case_file}" >> "$journal_file"
    cat >> "$journal_file" << 'EOF'

; Continue transient calculation
/solve/dual-time-iterate

; Save results
/file/write-case-data
EOF
    echo "${case_file}" >> "$journal_file"
    cat >> "$journal_file" << 'EOF'

; Exit Fluent
/exit
yes
EOF
    
    echo "Created journal file: $journal_file" >> "${LOG_DIR}/job_${case_dir}.log"
    
    # Run Fluent
    # Change 2ddp/3ddp based on your simulation dimension
    # 2ddp = 2D double precision, 3ddp = 3D double precision
    fluent 2ddp -ssh -g -cflush -pib -pib.ofed -t${CORES_PER_JOB} \
        -cnf="$job_nodefile" -i "$journal_file" >> "${LOG_DIR}/job_${case_dir}.log" 2>&1
    
    exit_code=$?
    
    # Clean up temporary files
    rm -f "$job_nodefile"
    
    echo "========================================" >> "${LOG_DIR}/job_${case_dir}.log"
    echo "Finished Fluent job for ${case_dir} at $(date)" >> "${LOG_DIR}/job_${case_dir}.log"
    echo "Exit code: $exit_code" >> "${LOG_DIR}/job_${case_dir}.log"
    echo "========================================" >> "${LOG_DIR}/job_${case_dir}.log"
    
    # Return to original directory
    cd "$PBS_O_WORKDIR"
    
    return $exit_code
}

export -f run_fluent_job
export LOG_DIR
export CORES_PER_JOB
export PROJECT_DIR
export PBS_O_WORKDIR
export PBS_NODEFILE

# Verify all case directories exist
echo "Checking case directories..."
all_dirs_exist=true
for case_dir in "${CASE_DIRS[@]}"; do
    full_path="${PBS_O_WORKDIR}/${PROJECT_DIR}/${case_dir}"
    if [ ! -d "$full_path" ]; then
        echo "WARNING: Directory $full_path does not exist!"
        all_dirs_exist=false
    else
        echo "  Found: $case_dir"
        # Check if case file exists
        case_count=$(ls -1 "$full_path"/*.cas.h5 2>/dev/null | wc -l)
        if [ $case_count -eq 0 ]; then
            echo "  WARNING: No .cas.h5 file found in $case_dir"
            all_dirs_exist=false
        fi
    fi
done

if [ "$all_dirs_exist" = false ]; then
    echo "ERROR: Not all case directories or files exist. Please check above."
    rm $PBS_O_WORKDIR/$PBS_JOBID
    exit 1
fi

echo "All directories verified. Starting jobs..."
echo ""

# Calculate how many jobs to run in parallel based on available cores
jobs_per_batch=$((12 / CORES_PER_JOB))
echo "Running $jobs_per_batch jobs in parallel (${CORES_PER_JOB} cores each)"
echo ""

# Run jobs in batches
total_jobs=${#CASE_DIRS[@]}
job_index=0

while [ $job_index -lt $total_jobs ]; do
    batch_end=$((job_index + jobs_per_batch))
    if [ $batch_end -gt $total_jobs ]; then
        batch_end=$total_jobs
    fi
    
    echo "Starting batch: jobs $((job_index + 1)) to ${batch_end}..."
    
    # Launch jobs in this batch
    for i in $(seq $job_index $((batch_end - 1))); do
        run_fluent_job "${CASE_DIRS[$i]}" $((i + 1)) &
    done
    
    # Wait for this batch to complete
    wait
    
    echo "Batch complete."
    echo ""
    
    job_index=$batch_end
done

echo "========================================"
echo "All Fluent jobs completed at $(date)"
echo "========================================"
echo "Logs saved to: $LOG_DIR/"
echo "Results saved in: ${PROJECT_DIR}/*/

## Tidy up the log directory
## DO NOT CHANGE THE LINE BELOW
## ============================
rm $PBS_O_WORKDIR/$PBS_JOBID
#
