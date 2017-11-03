
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

jobID=None

NUM_EXPERIMENTS=1
POINTS=(6400000)
DIMENSIONS=(128)
ITERATIONS=(10)
FRAGMENTS=(512)
CENTERS=(100)
TASKS_PER_NODE=(56)
NUM_NODES_LEN=4
NUM_NODES=(4 2 1 6)
TIME=(200 250 300 150)
QUEUE=("xeon" "xeon" "xeon" "xeon")

for (( tp=0; tp<$NUM_NODES_LEN; tp++ ))
do
    for (( i=0; i<$NUM_EXPERIMENTS; i++ ))
        do
            ./launch.sh $jobID ${NUM_NODES[tp]} ${TIME[tp]} ${TASKS_PER_NODE[i]} true ${QUEUE[tp]} ${POINTS[i]} ${DIMENSIONS[i]} ${CENTERS[i]} ${FRAGMENTS[i]} ${ITERATIONS[i]}
            echo $jobID ${NUM_NODES[tp]} ${TIME[tp]} ${TASKS_PER_NODE[i]} true ${QUEUE[tp]} ${POINTS[i]} ${DIMENSIONS[i]} ${CENTERS[i]} ${FRAGMENTS[i]} ${ITERATIONS[i]}
            wait_and_get_jobID_MT
        done
done
