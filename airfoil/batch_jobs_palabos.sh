#!/bin/bash
##
## Batch submission script for PBS on DELTA
## Runs multiple MPI jobs for AoA sweep with Palabos/Airfoil
##
## STEP 1: Enter a job name after the -N on the line below:
#PBS -N airfoil_AoA_sweep
##
## STEP 2: Select resources (adjust based on number of jobs)
#PBS -l select=1:ncpus=64:mpiprocs=64
##
## STEP 3: Select the correct queue
#PBS -q half_day
##
## STEP 4: Replace with your Cranfield email address
#PBS -m abe
#PBS -M gabriel.rochette.306@cranfield.ac.uk
##
## ====================================
## DO NOT CHANGE THE LINES BETWEEN HERE
## ====================================
#PBS -j oe
#PBS -W sandbox=PRIVATE
#PBS -k n
ln -s $PWD $PBS_O_WORKDIR/$PBS_JOBID
## Change to working directory
cd $PBS_O_WORKDIR
## ========
## AND HERE
## ========

## Configuration - MODIFY THESE FOR YOUR PROJECT
EXECUTABLE="./airfoil"                # Your MPI executable
CORES_PER_JOB=12                      # Cores per job
# List of AoA subdirectories to process (one job per directory)
CASE_AOA=("0" "8" "10" "20" "28" )
LOG_DIR="${PBS_O_WORKDIR}/airfoil_logs_${PBS_JOBID}"
mkdir -p "$LOG_DIR"

echo "========================================"
echo "Batch Airfoil Job Started at $(date)"
echo "Job ID: $PBS_JOBID"
echo "Node: $(hostname)"
echo "Working directory: $PBS_O_WORKDIR"
echo "========================================"

# Function to run a single Airfoil job
run_airfoil_job() {
    case_aoa=$1
    job_index=$2

    # Full path to the parameter file
    param_file="$./xmlFiles/parameters_${case_aoa}.xml"

    # Check if parameter file exists
    if [ ! -f "$param_file" ]; then
        echo "ERROR: Parameter file $param_file does not exist!" >> "${LOG_DIR}/job_${case_aoa}.log"
        return 1
    fi

    echo "========================================" >> "${LOG_DIR}/job_${case_aoa}.log"
    echo "Starting Airfoil job for ${case_aoa} at $(date)" >> "${LOG_DIR}/job_${case_aoa}.log"
    echo "Parameter file: $param_file" >> "${LOG_DIR}/job_${case_aoa}.log"
    echo "========================================" >> "${LOG_DIR}/job_${case_aoa}.log"

    # Run your MPI executable from the root folder
    mpirun -np $CORES_PER_JOB "$EXECUTABLE" "$param_file" >> "${LOG_DIR}/job_${case_aoa}.log" 2>&1
    exit_code=$?

    echo "========================================" >> "${LOG_DIR}/job_${case_aoa}.log"
    echo "Finished Airfoil job for ${case_aoa} at $(date)" >> "${LOG_DIR}/job_${case_aoa}.log"
    echo "Exit code: $exit_code" >> "${LOG_DIR}/job_${case_aoa}.log"
    echo "========================================" >> "${LOG_DIR}/job_${case_dir}.log"

    return $exit_code
}

export -f run_airfoil_job
export LOG_DIR
export CORES_PER_JOB
export PBS_O_WORKDIR
export EXECUTABLE

# Verify all case directories and parameter files exist
echo "Checking case directories and parameter files..."
all_files_exist=true
for case_aoa in "${CASE_AOA[@]}"; do
    full_path="${PBS_O_WORKDIR}/xmlFiles/"
    param_file="${full_path}/parameters_${case_aoa}.xml"
    if [ ! -d "$full_path" ]; then
        echo "WARNING: Directory $full_path does not exist!"
        all_files_exist=false
    elif [ ! -f "$param_file" ]; then
        echo "WARNING: Parameter file $param_file does not exist!"
        all_files_exist=false
    else
        echo "  Found: $case_aoa"
    fi
done

if [ "$all_files_exist" = false ]; then
    echo "ERROR: Not all case directories or parameter files exist. Please check above."
    rm $PBS_O_WORKDIR/$PBS_JOBID
    exit 1
fi

echo "All directories and files verified. Starting jobs..."
echo ""

# Calculate how many jobs to run in parallel based on available cores
jobs_per_batch=$((64 / CORES_PER_JOB)) # To modify depending on user input above
echo "Running $jobs_per_batch jobs in parallel (${CORES_PER_JOB} cores each)"
echo ""

# Run jobs in batches
total_jobs=${#CASE_AOA[@]}
job_index=0
while [ $job_index -lt $total_jobs ]; do
    batch_end=$((job_index + jobs_per_batch))
    if [ $batch_end -gt $total_jobs ]; then
        batch_end=$total_jobs
    fi

    echo "Starting batch: jobs $((job_index + 1)) to ${batch_end}..."

    # Launch jobs in this batch
    for i in $(seq $job_index $((batch_end - 1))); do
        run_airfoil_job "${CASE_AOA[$i]}" $((i + 1)) &
    done

    # Wait for this batch to complete
    wait

    echo "Batch complete."
    echo ""

    job_index=$batch_end
done

echo "========================================"
echo "All Airfoil jobs completed at $(date)"
echo "========================================"
echo "Logs saved to: $LOG_DIR/"
echo "Results saved in: ${CASE_AOA[*]}"

## Tidy up the log directory
rm $PBS_O_WORKDIR/$PBS_JOBID
