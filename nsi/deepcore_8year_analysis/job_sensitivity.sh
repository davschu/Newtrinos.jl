#!/bin/bash
#SBATCH --job-name=sensitivity_scan
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny        # single-node only; switch to cm4_std if/when multi-node distribution is wired up
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=66          # single hypothesis, single grid scan: threads ~= n-points budget, not the 5-hypothesis formula used by job_optimizer_tests.sh
#SBATCH --time=24:00:00
#SBATCH --array=0-5                 # one array task per hypothesis (eμ, eτ, μτ, ee_μμ, ττ_μμ, all) -- see HYPOTHESES below
# NOTE: the logs/ folder below must already exist when you run sbatch (Slurm opens the log files before the script starts; mkdir inside the script is too late)
#SBATCH --output=/dss/dsshome1/04/ge42niq2/Newtrinos.jl/nsi/deepcore_8year_analysis/results/sensitivity_results/logs/%x_%A_%a.out  # %x jobname, %A array job id, %a array task id
#SBATCH --error=/dss/dsshome1/04/ge42niq2/Newtrinos.jl/nsi/deepcore_8year_analysis/results/sensitivity_results/logs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=d.schultheiss@tum.de

module load slurm_setup

export JULIA_DEPOT_PATH="$HOME/.julia"
JULIA=/dss/lrzsys/sys/tools/julia_base/julia-1.12.6/bin/julia

# Canonical hypothesis order, indexed by $SLURM_ARRAY_TASK_ID (0-5) -- must match truth_grid_scan.jl's own all_cps Dict keys exactly.
HYPOTHESES=(eμ eτ μτ ee_μμ ττ_μμ all)
PARAMS=("ε_eμ_abs δ_eμ"  "ε_eτ_abs δ_eτ"  "ε_μτ_abs δ_μτ"  "Δ_eμ"  "Δ_τμ"  "ε_eμ_abs δ_eμ ε_eτ_abs δ_eτ ε_μτ_abs δ_μτ Δ_eμ Δ_τμ")
# optional restriction per task (empty = whole prior range); "name=lo,hi;name2=lo,hi"
RANGES=("ε_eμ_abs=0.01,0.3"  "ε_eτ_abs=0.01,0.35"  "ε_μτ_abs=0.005,0.07"  ""  ""  "")
NPOINTS=(64 64 64 20 20 100)

HYPOTHESIS="${HYPOTHESES[$SLURM_ARRAY_TASK_ID]}"
#for method scan
SCAN_PARAMS="${PARAMS[$SLURM_ARRAY_TASK_ID]}"
SCAN_RANGES="${RANGES[$SLURM_ARRAY_TASK_ID]}"
SCAN_NPOINTS="${NPOINTS[$SLURM_ARRAY_TASK_ID]}"

cd "$HOME/Newtrinos.jl/nsi"   # all relative paths below (setup, results, logs) are relative to nsi/ -- submit with sbatch from nsi/ too

# Don't comment sth out here, just delete the line or else it will cause a bug!
$JULIA ../src/analysis/sensitivity.jl \
    --setup setups/nsi_deepcore_8y.jl \
    --hypothesis "$HYPOTHESIS" \
    --strategy "mle optim-lbfgs" \
    --method scan \
    --params $SCAN_PARAMS \
    ${SCAN_RANGES:+--range "$SCAN_RANGES"} \
    --n-points  $SCAN_NPOINTS \
    --workers  64 \
    --results-dir deepcore_8year_analysis/results/sensitivity_results/ \
    --name "sensitivity_tests_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
