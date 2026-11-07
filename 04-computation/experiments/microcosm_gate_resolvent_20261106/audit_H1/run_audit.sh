#!/bin/sh
# audit_H1: regenerate every audit output (independent code; needs python3 + sympy + numpy, cc, and nauty gentourng at /opt/homebrew/bin).
set -e
cd "$(dirname "$0")"
cc -O2 -o hp_count hp_count.c
python3 audit_marked_pair.py > audit_marked_pair.out 2>&1
python3 audit_towers_numeric.py > audit_towers_numeric.out 2>&1
python3 audit_tower_ocf.py ./hp_count > audit_tower_ocf.out 2>&1
python3 audit_hyp9490_tables.py > audit_hyp9490_tables.out 2>&1
python3 audit_data.py ./hp_count > audit_data.out 2>&1
python3 audit_mod_p_groups.py > audit_mod_p_groups.out 2>&1
python3 audit_rate_and_k1_shape.py > audit_rate_and_k1_shape.out 2>&1
python3 audit_consequence_kappa0.py > audit_consequence_kappa0.out 2>&1
sh audit_priority_hits.sh > audit_priority_hits.out 2>&1
rm -f hp_count
rm -rf __pycache__
grep -H "FAILURES" *.out
echo "audit_H1 done"
