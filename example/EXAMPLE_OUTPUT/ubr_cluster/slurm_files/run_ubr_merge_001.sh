#!/bin/bash
#SBATCH --job-name=ubr_merge
#SBATCH --output=ubr_merge.o%j
#SBATCH --error=ubr_merge.e%j
#SBATCH --partition=biochem,owners
#SBATCH --time=2:00:00
#SBATCH -n 16
#SBATCH -N 1

ubr_merge.py UBR/*/ --merge_files RTB008_Twist_PK50_1M7.coverage.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB010_CustomArray_PK50_1M7.coverage.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB012_Twist_PK50_nomod.coverage.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB014_CustomArray_PK50_nomod.coverage.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB008_Twist_PK50_1M7.muts.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB010_CustomArray_PK50_1M7.muts.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB012_Twist_PK50_nomod.muts.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB014_CustomArray_PK50_nomod.muts.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB008_Twist_PK50_1M7.AC.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB010_CustomArray_PK50_1M7.AC.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB012_Twist_PK50_nomod.AC.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB014_CustomArray_PK50_nomod.AC.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB008_Twist_PK50_1M7.AG.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB010_CustomArray_PK50_1M7.AG.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB012_Twist_PK50_nomod.AG.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB014_CustomArray_PK50_nomod.AG.txt.gz &

wait
echo "DONE"
