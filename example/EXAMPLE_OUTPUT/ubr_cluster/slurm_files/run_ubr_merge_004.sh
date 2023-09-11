#!/bin/bash
#SBATCH --job-name=ubr_merge
#SBATCH --output=ubr_merge.o%j
#SBATCH --error=ubr_merge.e%j
#SBATCH --partition=biochem,owners
#SBATCH --time=2:00:00
#SBATCH -n 16
#SBATCH -N 1

ubr_merge.py UBR/*/ --merge_files RTB008_Twist_PK50_1M7.TC.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB010_CustomArray_PK50_1M7.TC.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB012_Twist_PK50_nomod.TC.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB014_CustomArray_PK50_nomod.TC.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB008_Twist_PK50_1M7.TG.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB010_CustomArray_PK50_1M7.TG.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB012_Twist_PK50_nomod.TG.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB014_CustomArray_PK50_nomod.TG.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB008_Twist_PK50_1M7.ins.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB010_CustomArray_PK50_1M7.ins.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB012_Twist_PK50_nomod.ins.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB014_CustomArray_PK50_nomod.ins.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB008_Twist_PK50_1M7.del.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB010_CustomArray_PK50_1M7.del.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB012_Twist_PK50_nomod.del.txt.gz &
ubr_merge.py UBR/*/ --merge_files RTB014_CustomArray_PK50_nomod.del.txt.gz &

wait
echo "DONE"
