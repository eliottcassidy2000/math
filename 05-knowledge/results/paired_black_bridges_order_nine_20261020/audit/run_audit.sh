#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python - <<'PY'
import json
with open('audit/expected_nodes.tsv','w') as f:
 for x in json.load(open('census_n9.json'))['nodes']:
  f.write('\t'.join(map(str,[x['node'],x['representative'],x['mass'],x['blue'],x['H'],x['aut'],int(x['SC']),int(x['strong']),x['c3'],x['degree']]))+'\n')
PY
g++ -O3 -std=c++17 -Wall -Wextra audit/verify_map.cpp -o audit/verify_map
./audit/verify_map > audit/map_audit.json 2> audit/map_audit.log
g++ -O3 -std=c++17 -Wall -Wextra audit/verify_graph.cpp -o audit/verify_graph
./audit/verify_graph > audit/graph_audit.out 2> audit/graph_audit.log
python audit/verify_certificates.py
python audit/verify_controls_and_squares.py
python audit/verify_score_construction.py
python audit/verify_producer_construction.py
python audit/verify_pendant_consequence.py
