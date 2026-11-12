#!/usr/bin/env python3
"""print the unit pairs of d_one_<N>.npz for strong words (y_1 = 1, y_n = 0), |O| >= 2."""
import sys
import numpy as np
N = int(sys.argv[1]); nshow = int(sys.argv[2]); sgn = int(sys.argv[3])
d = np.load(f"d_one_{N}.npz")
for n in range(2, nshow + 1):
    R = d[f"rec_{n}"]
    rows = []
    for w, O, sg, M, H in R.tolist():
        if sg != sgn: continue
        y = format(w, f"0{n}b")
        if sgn == 1 and not (y[0] == "1" and y[-1] == "0"): continue
        vs = [("a" if i == n else str(i + 1)) for i in range(n + 1) if (O >> i) & 1]
        if len(vs) < 2: continue
        rows.append((y, ",".join(vs), M, H))
    rows.sort()
    print(f"n = {n}: {len(rows)}")
    for r in rows: print("   ", r[0], "O={" + r[1] + "}", "M=%d" % r[2], "H=%d" % r[3])
