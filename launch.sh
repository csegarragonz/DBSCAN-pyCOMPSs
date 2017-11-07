#/bin/bash 

scriptDir=$(pwd)/$(dirname $0)
LOCAL_CLASSPATH=${scriptDir}
EXEC_FILE=${scriptDir}/DBSCAN.py
WORK_DIR=${scriptDir}/results/outputs/

echo $EXEC_FILE

if [ ! -d $WORK_DIR ]; then
    mkdir $WORK_DIR
fi


enqueue_compss \
    --job_dependency=$1
    --num_nodes=$2
    --exec_time=$3
    --cpus_per_node=$4
    --tracing=$5
    --lang=python \
    --worker_working_dir=scratch \
    --master_working_dir=. \
    $EXEC_FILE $6 $7 $8
