"""Referee: analyse census_n7.json produced by census.py."""
import json
from collections import Counter

import sys
data = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "census_n7.json"))
tot = Counter()
for n in sorted(data, key=int):
    recs = data[n]
    full = [r for r in recs if "crit_all" in r]
    big = [r for r in recs if "crit_all" not in r]
    # linear iff (C) (all subgroups), cyclic version, intersection-closure version
    mism_all = [r for r in full if r["linear"] != r["crit_all"]]
    mism_cyc = [r for r in full if r["crit_cyclic"] != r["crit_all"]]
    mism_int = [r for r in recs if "crit_all" in r and r["crit_intersection"] != r["crit_all"]]
    big_ok = all(r["linear"] and r["crit_intersection"] and r["tree"] for r in big)
    nonlin = [r for r in recs if not r["linear"]]
    inv_def = sum(1 for r in nonlin if r["involution_defect"])
    nofix = [r for r in recs if not r["A_fixes_vertex"]]
    lin_not_cheap = [r for r in recs if r["linear"] and not r["linear_cheap"]]
    weak_only = [r for r in full if r["weak"] and not r["linear"]]
    weak_only_odd = [r for r in weak_only if r["odd_sub_defect"]]
    weak_only_odd_elem = [r for r in weak_only if r["odd_elem_defect"]]
    neither = [r for r in full if not r["weak"] and not r["linear"]]
    count_coinc = [r for r in full if r["count_matches_some_X"] and not r["weak"]]
    print(f"n={n}: carriers {len(recs)} (full {len(full)}, big {len(big)}, big all tree+linear+C: {big_ok})")
    print(f"   linear != crit_all: {len(mism_all)}; crit_cyclic != crit_all: {len(mism_cyc)}; "
          f"crit_intersection != crit_all: {len(mism_int)}")
    print(f"   non-linear {len(nonlin)}, with involution defect {inv_def}")
    print(f"   A fixes no vertex: {len(nofix)} -> {[r['arcs'] for r in nofix]}")
    print(f"   linear but X=V-r fails: {len(lin_not_cheap)}")
    print(f"   weak-only {len(weak_only)}; of which with an odd-order subgroup defect {len(weak_only_odd)}"
          f" (odd element defect {len(weak_only_odd_elem)})")
    for r in weak_only_odd:
        print(f"      weak-only with odd defect: |A|={r['order']} arcs={r['arcs']}")
    print(f"   neither {len(neither)}:")
    for r in neither:
        print(f"      |A|={r['order']} count={r['count']} count_matches_some_X={r['count_matches_some_X']}"
              f" odd_sub_defect={r['odd_sub_defect']} arcs={r['arcs']}")
    print(f"   count equals some |Tour(X)/A| but A-sets differ: {len(count_coinc)}")
    tot["carriers"] += len(recs)
    tot["linear"] += sum(1 for r in recs if r["linear"])
print(tot)
