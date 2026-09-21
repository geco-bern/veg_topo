#!/bin/bash -l
#SBATCH --job-name="PaulH"
#SBATCH --time=96:00:00
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=icpu-stocker
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
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

# Set multi-threading environment variables
export OMP_NUM_THREADS=8
export GDAL_NUM_THREADS=8

echo "Job started on: $(date)"
BASE_URL="https://echosat.uni-muenster.de/"
OUTDIR="/storage/scratch/giub_geco/tting/data_raw/canopy_height_paul_2026"
mkdir -p "$OUTDIR"

# Error log file
ERROR_LOG="${OUTDIR}/errors.log"
> "$ERROR_LOG"

TEST_MODE=false
TEST_COUNT=1

echo "Fetching file list..."
FILES=$(wget -q -O - "$BASE_URL" | grep -oP 'href="[^"]+\.tif"' | cut -d'"' -f2)

if [ "$TEST_MODE" = true ]; then
    FILES=$(echo "$FILES" | head -n $TEST_COUNT)
    echo "TEST MODE: Processing $TEST_COUNT files"
fi

TOTAL=$(echo "$FILES" | wc -l)
COUNT=0

echo "$FILES" | while read -r file; do
    [ -z "$file" ] && continue
    COUNT=$((COUNT+1))
    url="${BASE_URL}${file}"
    fname=$(basename "$url")
    raw="$OUTDIR/$fname"
    out="$OUTDIR/${fname%.tif}_2020.tif"

    # If output file already exists, remove any residual raw file and skip
    if [ -f "$out" ]; then
        [ -f "$raw" ] && rm -f "$raw"
        continue
    fi

    # Download with resume support
    if ! wget -c -q -O "$raw" "$url" 2>/dev/null; then
        echo "[$COUNT/$TOTAL] Download failed: $fname" >> "$ERROR_LOG"
        continue
    fi

    # Extract Band 3
    # GDAL_CACHEMAX: 2 GB = 2147483648 bytes
    # GDAL_SWATH_SIZE: 100 MB = 104857600 bytes
    if ! gdal_translate -q -b 3 \
        -co COMPRESS=LZW \
        -co NUM_THREADS=8 \
        --config GDAL_CACHEMAX 2147483648 \
        --config GDAL_SWATH_SIZE 104857600 \
        "$raw" "$out" 2>/dev/null; then
        echo "[$COUNT/$TOTAL] GDAL failed: $fname" >> "$ERROR_LOG"
        continue
    fi

    # Remove raw file after successful processing
    rm -f "$raw"
done

echo "Job finished on: $(date)"
echo "Total files processed: $COUNT"
echo "Check errors.log for any failures"
