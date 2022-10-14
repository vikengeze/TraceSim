#!/bin/bash

echo "Building all tracing binaries... "

./build_champsim.sh bimodal no no no no lru 1 trace _a 
echo "Binary no. 1 is done"

./build_champsim.sh bimodal no no no no lru 1 trace _b > /dev/null
echo "Binary no. 2 is done"

./build_champsim.sh bimodal no no no no lru 1 trace _c > /dev/null
echo "Binary no. 3 is done"

./build_champsim.sh bimodal no no no no lru 1 trace _d > /dev/null
echo "Binary no. 4 is done"

./build_champsim.sh bimodal no no no no lru 1 trace _e > /dev/null
echo "Binary no. 5 is done"

./build_champsim.sh bimodal no no no no lru 1 trace _f > /dev/null
echo "Binary no. 6 is done"

./build_champsim.sh bimodal no no no no lru 1 trace _g > /dev/null
echo "Binary no. 7 is done"

./build_champsim.sh bimodal no no no no lru 1 trace _h > /dev/null
echo "Binary no. 8 is done"
echo "All 8 binaries have been built"

