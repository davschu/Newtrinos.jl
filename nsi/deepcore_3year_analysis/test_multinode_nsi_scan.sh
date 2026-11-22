#!/bin/bash
# Minimal, fast-failing smoke test for --multi-node in nsi_scan.jl on LRZ cm4_std.
#
# nsi_scan.jl --task profile is the REAL deliverable (the final NSI profile-likelihood
# scan) -- roundtrips.jl is only a validation/truth-recovery tool. Its --multi-node path
# was implemented by mirroring roundtrips.jl but has NOT been smoke-tested at all yet
# (only parse-checked), so this must be validated independently before trusting it for
# the real 50/400-point jobs.
#
# Goal: confirm SlurmClusterManager.jl's SlurmManager() spans workers across BOTH
# allocated nodes for THIS script specifically, as cheaply as possible -- NOT to
# produce a usable physics result. Grid size and PSO/LBFGS budgets are deliberately
# tiny so this finishes (or fails) in minutes, not hours.
#
# cm4_std requires a minimum of 2 nodes AND a QOS floor of 112 physical cores per job
# (confirmed via QOSMinCpuNotSatisfied on a smaller first attempt -- cm4_std's minimum
# job size is NOT negotiable below 112 cores, unlike cm4_tiny). So core count can't be
# shrunk further than this; cost is kept low instead by keeping the actual workload
# (grid points, PSO budget, walltime) tiny -- this job will still finish in minutes.
# 4 tasks/node x 2 nodes x 14 threads/task = 112 cores total (the QOS floor, exactly).
#
#SBATCH --job-name=multinode_test_nsi_scan
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_std
#SBATCH --qos=cm4_std
#SBATCH --mem=50G
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=14
#SBATCH --time=00:15:00
#SBATCH --output=scan_results/logs/%x_%j.out
#SBATCH --error=scan_results/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=d.schultheiss@tum.de

module load slurm_setup

export JULIA_DEPOT_PATH="$HOME/.julia"
JULIA=/dss/lrzsys/sys/tools/julia_base/julia-1.12.6/bin/julia

cd "$HOME/Newtrinos.jl/nsi"
mkdir -p scan_results/logs

# Don't comment sth out here, just delete the line or else it will cause a bug!
$JULIA nsi_scan.jl \
    --name multinode_smoketest \
    --hypothesis μτ \
    --experiments deepcore \
    --task profile \
    --grid-points 4 \
    --objective posterior \
    --seed-strategy pso \
    --nseeds 1 \
    --pso-particles 5 \
    --global-iterations 2 \
    --iterations 5 \
    --workers 8 \
    --threads 14 \
    --multi-node
