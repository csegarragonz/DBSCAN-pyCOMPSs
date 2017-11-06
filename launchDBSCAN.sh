#Script to launch a batch of tests


bID() {
    sleep 6
        jobID=$(bjobs | tail -n 1 | cut -c -7)
        echo "jobID = $jobID"
}

wait_and_get_jobID_MT() {
    sleep 6
    jobID=$(squeue | grep $(whoami) | sort -n - | tail -n 1 | awk '{ print $1  }')
    echo
}


#Test Parameters:
    #Number of executions of the same test
    NUM_EXPERIMENTS=3
    #Diferent Datasets to be read from /gpfs/projects/bsc19/COMPSs_DATASETS/dbscan/DATASETS[i]/*.txt
    DATASETS=(1 2 3 4 5)
    DATASETS_LEN=${#DATASETS[@]}
    STRONG_SCALING=true
    WEAK_SCALING=false

#Arguments for the enqueue_compss
    jobID=None
    NUM_NODES=(2 3 4 5 6 7 8 9 10)
    NUM_NODES_LEN=${#NUM_NODES[@]}
    TIME=(100 100 100 100 200 250 300 150 100)
    CPUS_PER_NODE=(48)
    TRACING=false

#Arguments for the DBSCAN
    EPSILON=0.1
    MIN_POINTS=10

for (( j=0; j<$DATASETS_LEN; j++ ))
do
    for (( tp=0; tp<$NUM_NODES_LEN; tp++ ))
    do
        for (( i=0; i<$NUM_EXPERIMENTS; i++ ))
            do
                ./launch.sh $jobID ${NUM_NODES[tp]} ${TIME[tp]} ${CPUS_PER_NODE[i]} $TRACING $EPSILON $MIN_POINTS ${DATASETS[j]}
                wait_and_get_jobID_MT
            done
    done
done
