#!/bin/bash

cd ..

INSTR_NUM=100

fp='/various/dstyliar/ML-DPC/ChampSimTraces/spec17/*'

for trace in $fp; do
	if [ $(cut -d. -f1-1 <<<$(basename $trace)) -eq 605 ] || [ $(cut -d. -f1-1 <<<$(basename $trace)) -eq 607 ]; then
		echo "Creating trace file for application $(basename "$trace") and $INSTR_NUM million isntructions ..."
		./run_champsim.sh trace_f-bimodal-no-no-no-lru-1core 1 $INSTR_NUM $trace spec17
		echo "Done, results in statistics/spec17/"$(basename "$trace")"_"$INSTR_NUM"M_trace.txt"	
	fi
done

cd tracers