#/bin/bash 

scriptDir=$(pwd)/$(dirname $0)
LOCAL_CLASSPATH=${scriptDir}
EXEC_FILE=${scriptDir}/DBSCAN.py
WORK_DIR=${scriptDir}/results/outputs/

echo $EXEC_FILE

if [ ! -d $WORK_DIR ]; then
    mkdir $WORK_DIR
fi


schedulers=(es.bsc.compss.scheduler.fifoScheduler.FIFOScheduler \
es.bsc.compss.scheduler.lifoScheduler.LIFOScheduler \
es.bsc.compss.scheduler.loadBalancingScheduler.LoadBalancingScheduler \
es.bsc.compss.scheduler.fifoDataScheduler.FIFODataScheduler
) 

enqueue_compss \
    --job_dependency=$1
    --num_nodes=$2
    --exec_time=$3
    --cpus_per_node=$4
    --tracing=$5
    --lang=python \
    --worker_working_dir=scratch \
    $EXEC_FILE $6 $7 $8
