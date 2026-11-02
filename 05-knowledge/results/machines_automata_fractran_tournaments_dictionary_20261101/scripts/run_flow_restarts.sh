#!/bin/sh
# usage: run_flow_restarts.sh N maxseeds -- rerun the flow solver with seeds 1..maxseeds until a path is found
N=$1; K=$2; s=1
while [ $s -le $K ]; do
  python3 -u scripts/prime_path_flow.py $N $s && exit 0
  s=$((s+1))
done
echo "no path found at N=$N after $K seeds"; exit 1
