#!/bin/bash
TRACE_DIR='../traces_gap'
#TRACE_DIR='/various/dstyliar/ML-DPC/ChampSimTraces/'

binary=${1}
n_warm=${2}
n_sim=${3}
trace=${4}
option=${5}
extra=${6}

DESTINATION_FOLDER='statistics'
mkdir -p ${DESTINATION_FOLDER}/${option}

(./bin/${binary} -warmup_instructions ${n_warm}000000 -simulation_instructions ${n_sim}000000 ${extra} -traces ${TRACE_DIR}/${trace}) &> ${DESTINATION_FOLDER}/${option}/$(basename $trace)_${n_sim}M_trace.txt
