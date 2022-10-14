#!/bin/bash

cd ..

INSTR_NUM=100

fp='/various/dstyliar/ML-DPC/ChampSimTraces/gap/b*'

for trace in $fp; do
	echo "Creating trace file for application $(basename "$trace") and $INSTR_NUM million isntructions ..."
	./run_champsim.sh trace_a-bimodal-no-no-no-lru-1core 1 $INSTR_NUM $trace gap
	echo "Done, results in statistics/gap/"$(basename "$trace")"_"$INSTR_NUM"M_trace.txt"	
done


sfp='/various/dstyliar/ML-DPC/ChampSimTraces/gap/sssp*'

for trace in $sfp; do
	echo "Creating trace file for application $(basename "$trace") and $INSTR_NUM million isntructions ..."
	./run_champsim.sh trace_a-bimodal-no-no-no-lru-1core 1 $INSTR_NUM $trace gap
	echo "Done, results in statistics/gap/"$(basename "$trace")"_"$INSTR_NUM"M_trace.txt"	
done

tfp='/various/dstyliar/ML-DPC/ChampSimTraces/gap/cc*'

for trace in $tfp; do
	echo "Creating trace file for application $(basename "$trace") and $INSTR_NUM million isntructions ..."
	./run_champsim.sh trace_a-bimodal-no-no-no-lru-1core 1 $INSTR_NUM $trace gap
	echo "Done, results in statistics/gap/"$(basename "$trace")"_"$INSTR_NUM"M_trace.txt"	
done

cd tracers