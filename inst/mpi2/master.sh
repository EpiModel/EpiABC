
# Run batch queue
sbatch -p csde -A csde runsim.sh

# Run checkpoint queue (previously called backfill)
sbatch -p ckpt -A csde-ckpt runsim.sh

