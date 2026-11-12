#!/usr/bin/env python3
"""Skeptic check 3b: S_w and N(K) for w = 17..20 by the eps enumeration (method 2 of s3_sw.py)."""
import sys, time
import numpy as np
import s3_sw
out = []
def P(*a):
    s = " ".join(str(t) for t in a)
    print(s, flush=True)
    out.append(s)
w0 = int(sys.argv[1]); w1 = int(sys.argv[2])
for w in range(w0, w1 + 1):
    t0 = time.time()
    cnt = s3_sw.eps_counts(w)
    s3_sw.report(w, cnt, P, "EPS")
    P("    time %.1fs" % (time.time() - t0))
    np.save("cnt_w%d.npy" % w, cnt)
with open(sys.argv[3], "w") as f:
    f.write("\n".join(out) + "\n")
