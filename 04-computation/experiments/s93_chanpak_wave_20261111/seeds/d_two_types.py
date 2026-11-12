#!/usr/bin/env python3
"""S93 lane seeds: types of the unit sets over two-apex bases (d_two_<N>.npz): which apexes are in O, how many
spine vertices; number of pairs and of distinct slopes per type and length."""
import sys
from collections import defaultdict
import numpy as np

N = int(sys.argv[1])
d = np.load(f"d_two_{N}.npz")
for sgn in (1, -1):
    print("sign", sgn)
    tot = defaultdict(set)
    for n in range(1, N + 1):
        R = d[f"rec_{n}"]
        types = defaultdict(list)
        for w, O, sg, M, H in R.tolist():
            if sg != sgn:
                continue
            u = (O >> n) & 1
            v = (O >> (n + 1)) & 1
            k = bin(O & ((1 << n) - 1)).count("1")
            types[("u" if u else "") + ("v" if v else "") + f"+{k}"].append(M)
        line = f"  n = {n}: "
        for t in sorted(types):
            line += f"{t}: {len(types[t])} ({len(set(types[t]))})  "
            tot[t] |= set(types[t])
        print(line)
    print("  distinct slopes by type over all lengths:", {t: len(s) for t, s in sorted(tot.items())})
    allm = set().union(*tot.values())
    for t in sorted(tot):
        others = set().union(*[s for tt, s in tot.items() if tt != t]) if len(tot) > 1 else set()
        print(f"    type {t}: {len(tot[t])} slopes, {len(tot[t] - others)} not given by another type")
