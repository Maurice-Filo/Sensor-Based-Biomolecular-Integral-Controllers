#!/bin/bash

module load matlab/R2023b

sbatch --ntasks=1 --cpus-per-task=16 --mem-per-cpu=2G --time=48:00:00 --wrap="matlab -nodisplay -r FP_BirthDeath_delta"
sbatch --ntasks=1 --cpus-per-task=16 --mem-per-cpu=2G --time=48:00:00 --wrap="matlab -nodisplay -r sAIF_BirthDeath_mu_eta"
sbatch --ntasks=1 --cpus-per-task=16 --mem-per-cpu=2G --time=48:00:00 --wrap="matlab -nodisplay -r rAIF_BirthDeath_mu_eta_k"

sbatch --ntasks=1 --cpus-per-task=16 --mem-per-cpu=2G --time=48:00:00 --wrap="matlab -nodisplay -r FP_GeneExp_delta"
sbatch --ntasks=1 --cpus-per-task=16 --mem-per-cpu=2G --time=48:00:00 --wrap="matlab -nodisplay -r sAIF_GeneExp_mu_eta"
sbatch --ntasks=1 --cpus-per-task=16 --mem-per-cpu=2G --time=48:00:00 --wrap="matlab -nodisplay -r rAIF_GeneExp_mu_eta_k"

sbatch --ntasks=1 --cpus-per-task=16 --mem-per-cpu=2G --time=48:00:00 --wrap="matlab -nodisplay -r FP_GeneExpMat_delta"
sbatch --ntasks=1 --cpus-per-task=16 --mem-per-cpu=2G --time=48:00:00 --wrap="matlab -nodisplay -r sAIF_GeneExpMat_mu_eta"
sbatch --ntasks=1 --cpus-per-task=16 --mem-per-cpu=2G --time=48:00:00 --wrap="matlab -nodisplay -r rAIF_GeneExpMat_mu_eta_k"