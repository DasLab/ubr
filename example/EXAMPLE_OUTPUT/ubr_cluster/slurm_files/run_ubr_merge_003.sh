#!/bin/bash
#SBATCH --job-name=ubr_merge
#SBATCH --output=ubr_merge.o%j
#SBATCH --error=ubr_merge.e%j
#SBATCH --partition=biochem,owners
#SBATCH --time=2:00:00
#SBATCH -n 16
#SBATCH -N 1

ubr_merge.py UBR/*/ --merge_files RTB008_Twist_PK50_1M7.GA.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB010_CustomArray_PK50_1M7.GA.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB012_Twist_PK50_nomod.GA.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB014_CustomArray_PK50_nomod.GA.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB008_Twist_PK50_1M7.GC.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB010_CustomArray_PK50_1M7.GC.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB012_Twist_PK50_nomod.GC.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB014_CustomArray_PK50_nomod.GC.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB008_Twist_PK50_1M7.GT.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB010_CustomArray_PK50_1M7.GT.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB012_Twist_PK50_nomod.GT.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB014_CustomArray_PK50_nomod.GT.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB008_Twist_PK50_1M7.TA.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB010_CustomArray_PK50_1M7.TA.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB012_Twist_PK50_nomod.TA.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB014_CustomArray_PK50_nomod.TA.txt.gz &

wait
echo "DONE"
