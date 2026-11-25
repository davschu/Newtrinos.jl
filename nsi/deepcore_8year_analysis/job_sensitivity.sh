#!/bin/bash
#SBATCH --job-name=sensitivity_fisher
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny        # single-node only; switch to cm4_std if/when multi-node distribution is wired up
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=27          # single hypothesis, single grid scan: threads ~= n-points budget, not the 5-hypothesis formula used by job_optimizer_tests.sh
#SBATCH --time=4:00:00
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
HYPOTHESIS="${HYPOTHESES[$SLURM_ARRAY_TASK_ID]}"

cd "$HOME/Newtrinos.jl/nsi"   # all relative paths below (setup, results, logs) are relative to nsi/ -- submit with sbatch from nsi/ too

# Don't comment sth out here, just delete the line or else it will cause a bug!
$JULIA ../src/analysis/sensitivity.jl \
    --setup setups/nsi_deepcore_8y.jl \
    --hypothesis "$HYPOTHESIS" \
    --method fisher \
    --truth grid \
    --n-truths 25 \
    --information expected \
    --drop-boundaries \
    --workers 25 \
    --results-dir deepcore_8year_analysis/results/sensitivity_results/ \
    --name "sensitivity_tests_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
