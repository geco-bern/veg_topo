#! /usr/bin/bash -l
#SBATCH --job-name="twi450mClean"
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker
#SBATCH --nodes=1
#SBATCH --time=4:00:00     #
#SBATCH --ntasks=1
#SBATCH --partition=icpu-stocker
#SBATCH --cpus-per-task=1
#SBATCH --mem=400G  # > 400 G
#SBATCH --mail-user=ting.tan@students.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --chdir=/storage/homefs/tt22k003/veg_topo/analysis

# Set personal R library
export R_LIBS_USER=/storage/homefs/tt22k003/R/x86_64-pc-linux-gnu-library/4.4

# Load modules
module load foss/2024a
module load PROJ/9.4.1-GCCcore-13.3.0
module load GDAL/3.10.0-foss-2024a
module load R/4.4.2-gfbf-2024a

# Job information
echo "=================================================="
echo "Job started on: $(date --rfc-3339=seconds)"
echo "Job ID: $SLURM_JOB_ID"
echo "Job name: $SLURM_JOB_NAME"  # 输出job name
echo "Hostname: $(hostname)"
echo "Working directory: $PWD"
echo "R_LIBS_USER: $R_LIBS_USER"
echo "=================================================="

# Force Rscript to use the same library paths as RStudio Server
Rscript -e '.libPaths(c(
  "/storage/homefs/tt22k003/R/x86_64-pc-linux-gnu-library/4.4",
  "/storage/software/epyc2.9/software/R-bundle-CRAN/2024.11-foss-2024a",
  "/storage/software/epyc2.9/software/R/4.4.2-gfbf-2024a/lib64/R/library"
)); cat("Running script: 1_01_1_twi_450m_clean.R\n"); source("1_01_1_twi_450m_clean.R")'

# Capture the exit status
EXIT_STATUS=$?
echo "=================================================="
echo "Job finished on: $(date --rfc-3339=seconds)"
echo "Exit status: $EXIT_STATUS"
echo "Job name: $SLURM_JOB_NAME"
echo "=================================================="

# Exit with the same status as the R script
exit $EXIT_STATUS
