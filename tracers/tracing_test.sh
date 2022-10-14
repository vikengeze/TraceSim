#!/bin/bash

cd ..

INSTR_NUM=100

trace=/home/users/vgezekel/bfs-10.trace.gz

echo "Creating trace file for application $(basename $trace) and $INSTR_NUM million isntructions ..."
./run_champsim.sh trace_d-bimodal-no-no-no-lru-1core 1 $INSTR_NUM $trace test
echo "Done, results in statistics/gap/"$(basename $trace)"_"$INSTR_NUM"M_trace.txt"

cd tracers