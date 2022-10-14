#!/bin/bash

cd ..

INSTR_NUM=100

fp='/various/dstyliar/ML-DPC/ChampSimTraces/spec17/*'

for trace in $fp; do
	if [ $(cut -d. -f1-1 <<<$(basename $trace)) -eq 602 ] || [ $(cut -d. -f1-1 <<<$(basename $trace)) -ge 649 ]; then
		echo "Creating trace file for application $(basename "$trace") and $INSTR_NUM million isntructions ..."
		./run_champsim.sh trace_g-bimodal-no-no-no-lru-1core 1 $INSTR_NUM $trace spec17
		echo "Done, results in statistics/spec17/"$(basename "$trace")"_"$INSTR_NUM"M_trace.txt"	
	fi
done

cd tracers