runcompss --lang=python -g ./launchDBSCAN.py ./data/moons.txt 3 0.1 10 1 2D
#python ./launchDBSCAN.py ./data/ds1s.txt 3 0.015 10 1 2D
#enqueue_compss --lang=python --num_nodes=4 --exec_time=10 -t  --log_level=debug --worker_working_dir=gpfs ./launchDBSCAN.py ./data/5k.txt 0.1 0.015 10
#enqueue_compss --lang=python --num_nodes=8 --exec_time=40 -t  ./launchDBSCAN.py ./data/100k.txt 0.1 0.015 10
