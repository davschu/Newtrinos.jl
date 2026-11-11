#!/bin/bash
#SBATCH --job-name=nsi_scan_posterior_2d_multinode
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
#SBATCH --array=0-2                 # 2D ("complex") hypotheses only -- see HYPOTHESES below
#SBATCH --output=scan_results/logs/%x_%A_%a.out
#SBATCH --error=scan_results/logs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=d.schultheiss@tum.de

# Real 2D ("complex") NSI profile scan (21x21 = 441 points), multi-node on cm4_std.
#
# eμ, eτ, μτ are the genuinely 2D hypotheses -- confirmed directly from nsi_scan.jl's
# cp_* conditioning Dicts: cp_eμ leaves ε_eμ_abs+δ_eμ free, cp_eτ leaves ε_eτ_abs+δ_eτ
# free, cp_μτ leaves ε_μτ_abs+δ_μτ free (each a magnitude+phase pair). ee_μμ/ττ_μμ are
# each already 1D (free params Δ_eμ/Δ_τμ alone) -- see the 1D job script instead. The
# "real" ε_ab value for eμ/eτ/μτ is extracted from this scan's own δ_ab=0,π slices
# after the fact -- no separate 1D run is needed for them.
#
# Sizing: 4 nodes (cm4_std max) x 28 tasks/node x 4 cpus/task = 112 workers x 4
# threads = 448 cores (the cm4_std job-size ceiling). --nseeds 3 is handled internally
# via Threads.@threads across each worker's --threads=4, as already wired into
# nsi_scan.jl.
#
# --global-iterations 100, --pso-particles 500: same per-point budget as the 1D job
# and the validated 25-point timing run (~4.2h/point measured on cm4_tiny-equivalent
# per-core performance -- re-verify this holds on cm4_std's actual CPUs).
#
# 441 grid points / 112 workers = ceil(441/112) = 4 points per worker (sequential,
# same bucket as 400 -- no cost increase) -> ~4 x 4.2h = ~16.8h estimated wall-clock,
# leaving buffer under the 24h cap. This is the TIGHTEST case in the whole scaling
# investigation -- if real per-point cost on cm4_std exceeds ~6h, this will NOT fit
# and the deferred grid-slicing fallback (see the retained plan notes) needs to be
# revisited.
#
# --grid-points is the TOTAL points target, not per-dimension -- nsi_scan.jl computes
# n_per_dim = round(sqrt(grid_points)) for 2D hypotheses (nsi_scan.jl:355). 441 is used
# (not 400) so sqrt(441) = 21 exactly, an ODD per-dimension count -- this guarantees a
# grid point exactly at quantile 0.5 of each prior (e.g. the midpoint of the magnitude
# and the phase axes), so the SM point (NSI off) is covered without needing a separate
# dedicated run.

module load slurm_setup

export JULIA_DEPOT_PATH="$HOME/.julia"
JULIA=/dss/lrzsys/sys/tools/julia_base/julia-1.12.6/bin/julia

HYPOTHESES=(eμ eτ μτ)
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
    --grid-points 441 \
    --objective posterior \
    --seed-strategy pso \
    --nseeds 3 \
    --pso-particles 500 \
    --global-iterations 100 \
    --workers 112 \
    --threads 4 \
    --multi-node
