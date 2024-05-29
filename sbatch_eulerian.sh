#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --cpus-per-task 1
#SBATCH --time 3:00:00
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --mail-type=BEGIN,END,FAIL   # Send email on job start
#SBATCH --mail-user=jacopo.bilotto@epfl.ch

#source ./venv/bin/activate

python main_obstruction_shear.py -c 0.3 -a 3. -cw 0.6 -vw 10.0 -f 2 -vf 0.5
