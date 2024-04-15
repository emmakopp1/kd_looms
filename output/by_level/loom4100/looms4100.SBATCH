#!/bin/sh
# Fichier looms4100.SBATCH
#SBATCH --nodes=1
#SBATCH -c 30
#SBATCH --output=looms4100.out

~/beast/bin/beast -overwrite /home/users/kopp/work/loom4100/loom4100.xml