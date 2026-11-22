#!/bin/bash
#SBATCH --job-name=test_pso_5points_p800g125  # shared across all array tasks by design; squeue's JOBID column (e.g. 12345678_3) distinguishes tasks -- see %A_%a in --output/--error below
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny        # single-node only; switch to cm4_std if/when multi-node distribution is wired up
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=92          # single hypothesis, single grid scan: threads ~= n-points budget, not the 5-hypothesis formula used by job_optimizer_tests.sh
#SBATCH --time=24:00:00
#SBATCH --array=0-5                 # one array task per hypothesis (eμ, eτ, μτ, ee_μμ, ττ_μμ, all) -- see HYPOTHESES below
#SBATCH --output=pso_test_results/logs/%x_%A_%a.out  # %x jobname, %A array job id, %a array task id
#SBATCH --error=pso_test_results/logs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=d.schultheiss@tum.de

export JULIA_DEPOT_PATH="$HOME/.julia"
JULIA=/dss/lrzsys/sys/tools/julia_base/julia-1.12.6/bin/julia

# Canonical hypothesis order, indexed by $SLURM_ARRAY_TASK_ID (0-5) -- must match truth_grid_scan.jl's own all_cps Dict keys exactly.
HYPOTHESES=(eμ eτ μτ ee_μμ ττ_μμ all)
HYPOTHESIS="${HYPOTHESES[$SLURM_ARRAY_TASK_ID]}"
NAME="${SLURM_JOB_NAME}_${HYPOTHESIS}_$(date +%d%m%Y)_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

cd "$HOME/Newtrinos.jl/nsi"
mkdir -p pso_test_results/logs

# Don't comment sth out here, just delete the line or else it will cause a bug!
$JULIA --threads=$((SLURM_CPUS_PER_TASK - 2)) test_pso.jl \
    --hypothesis "$HYPOTHESIS" \
    --n-points 5 \
    --pso-subruns 5 \
    --fluctuate-nuisances \
    --vary-physics-truth \
    --global-iterations 125 \
    --pso-particles 800 \
    --suffix "p800g125"
