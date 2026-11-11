#!/bin/bash
#SBATCH --job-name=nsi_scan_posterior_1d_multinode
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_std
#SBATCH --qos=cm4_std
#SBATCH --mem=400G
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --array=0-1                 # 1D ("real") hypotheses only -- see HYPOTHESES below
#SBATCH --output=scan_results/logs/%x_%A_%a.out
#SBATCH --error=scan_results/logs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=d.schultheiss@tum.de

# Real 1D ("real") NSI profile scan, multi-node on cm4_std.
#
# Only ee_μμ (free param: Δ_eμ) and ττ_μμ (free param: Δ_τμ) are genuinely 1D --
# confirmed directly from nsi_scan.jl's cp_* conditioning Dicts: every other named
# hypothesis (eμ, eτ, μτ) leaves an ε_*_abs + δ_* pair free (2D -- see the 2D job
# script). eμ/eτ/μτ's "real" slice (ε_ab at δ_ab=0 or π) is extracted FROM the 2D
# scan's own grid after the fact -- it is not scanned separately here.
#
# Sizing: 4 nodes (cm4_std max) x 28 tasks/node x 4 cpus/task = 112 workers x 4
# threads = 448 cores (the cm4_std job-size ceiling). --nseeds 3 is handled internally
# via Threads.@threads across each worker's --threads=4, as already wired into
# nsi_scan.jl -- no separate seed-count flag needed here.
#
# --global-iterations 100 (not the earlier 200) per your own measured reduction
# (confirmed ~4.2h/point at --pso-particles 500 on cm4_tiny-equivalent per-core
# performance -- re-verify this holds on cm4_std's actual CPUs if timings look off).
#
# 51 grid points (odd, so quantile 0.5 of each prior lands exactly on a grid point --
# e.g. the midpoint of a symmetric/zero-centered NSI param) / 112 workers = at most 1
# point per worker (some idle) -- should comfortably finish well inside 24h on a
# single pass.

module load slurm_setup

export JULIA_DEPOT_PATH="$HOME/.julia"
JULIA=/dss/lrzsys/sys/tools/julia_base/julia-1.12.6/bin/julia

HYPOTHESES=(ee_μμ ττ_μμ)
HYPOTHESIS="${HYPOTHESES[$SLURM_ARRAY_TASK_ID]}"
NAME="${SLURM_JOB_NAME}_${HYPOTHESIS}_$(date +%d%m%Y)_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

cd "$HOME/Newtrinos.jl/nsi/deepcore_3year_analysis"
mkdir -p scan_results/logs

# Don't comment sth out here, just delete the line or else it will cause a bug!
$JULIA nsi_scan.jl \
    --name "$NAME" \
    --hypothesis "$HYPOTHESIS" \
    --experiments deepcore_3y \
    --task profile \
    --grid-points 51 \
    --objective posterior \
    --seed-strategy pso \
    --nseeds 3 \
    --pso-particles 500 \
    --global-iterations 100 \
    --workers 112 \
    --threads 4 \
    --multi-node
