#!/bin/sh
# Fichier looms.SBATCH
#SBATCH --nodes=1
#SBATCH -c 30
#SBATCH --output=looms.out

~/beast/bin/beast -overwrite /home/users/kopp/work/loom1000/loom1000.xml