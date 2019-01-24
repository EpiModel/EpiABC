
# Run batch queue
sbatch -p csde -A csde runsim.sh

# Run checkpoint queue
sbatch -p ckpt -A csde-ckpt runsim.sh

