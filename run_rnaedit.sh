#!/bin/bash

#Arguments 

echo "Sample: $1"

sbatch \
--export=sample=$1 \
-J $1 \
-D /home/a1mark/logs \
-c 4 \
/expanse/lustre/projects/csd691/kfisch/RNA_editing_pipeline/scripts/rnaedit.sh
