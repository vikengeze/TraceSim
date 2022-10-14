#!/bin/bash

cd ..

INSTR_NUM=100

fp='/various/dstyliar/ML-DPC/ChampSimTraces/spec06/*'

for trace in $fp; do
	if [ $(cut -d. -f1-1 <<<$(basename $trace)) -eq 429 ]; then
		echo "Creating trace file for application $(basename "$trace") and $INSTR_NUM million isntructions ..."
		./run_champsim.sh trace_c-bimodal-no-no-no-lru-1core 1 $INSTR_NUM $trace spec06
		echo "Done, results in statistics/spec06/"$(basename "$trace")"_"$INSTR_NUM"M_trace.txt"	
	fi
done

cd tracers