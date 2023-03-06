#!/bin/bash
#
#SBATCH --ntasks 1            # number of tasks
#SBATCH --mem 20G            # memory pool per process
#SBATCH -o LOCI_KC_B1_B2.out    # STDOUT
#SBATCH -t 03:00:00            # time (D-HH:MM)

./simulate_traits.sh LOCI KC B1 B2
./run_plink.sh LOCI KC_B1_B2
./run_gem.sh LOCI KC_B1_B2
./run_SH.sh LOCI KC_B1_B2
