#!/bin/bash

module load matlab/R2023b
# Loop over the range 2 to 32
for i in {1..32}; do
    # Construct the script name dynamically
    matlab_script="sAIF_BirthDeath_mu${i}"
    
    # Submit the job to SLURM using sbatch
    sbatch --ntasks=1 --cpus-per-task=16 --mem-per-cpu=2G --time=48:00:00 --wrap="matlab -nodisplay -r ${matlab_script}"
done

sbatch --ntasks=1 --cpus-per-task=16 --mem-per-cpu=2G --time=48:00:00 --wrap="matlab -nodisplay -r FP_BirthDeath"
sbatch --ntasks=1 --cpus-per-task=16 --mem-per-cpu=2G --time=48:00:00 --wrap="matlab -nodisplay -r FP_BirthDeath_FixedHightheta"
sbatch --ntasks=1 --cpus-per-task=16 --mem-per-cpu=2G --time=48:00:00 --wrap="matlab -nodisplay -r FP_BirthDeath_widedelta"