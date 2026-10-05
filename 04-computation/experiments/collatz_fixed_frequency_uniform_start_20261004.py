#!/usr/bin/env python3
"""
E5r (opus, 2026-10-05): the typical rate of the UNIFORM-START q-tower (Fourier note 4k/4m).

The uniform-start tower has the deepest level's digit string uniform to the full window depth and every shallower
level given by the exact q-th-root (Pascal) law; its second moment is exactly 3^-n (4k, PROVED).  Question: is its
TYPICAL rate the Collatz value 0.5700 (then the excess is the cube-root coupling with any start, and the Jensen gap
against the exactly known rms 0.5774 is 1.3%) or the i.i.d. value 0.5723 (then the excess needs the arithmetic,
integer start)?  Eight seeds each for q = 3 and q = 5, N = 1500, A = 40, rates over 300..1500 and 750..1500.
Usage: python collatz_fixed_frequency_uniform_start_20261004.py [N=1500] [A=40] [S=8]
"""
import sys, math, time
import numpy as np
sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from collatz_fixed_frequency_random_digits_20261004 import run

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    S = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"uniform-start q-towers, N={N}, A={A}, {S} seeds; reference: integer-start q=3 0.5700 +- 0.0009, i.i.d. 0.5723 +- 0.0008")
    t0 = time.time()
    for q in (3, 5):
        rates = []
        for s in range(S):
            v = run(N, A, "uniformstart", q=q, seed=4000 + s)
            rs = []
            for a_ in (300, 750):
                ns = np.arange(a_, N + 1); y = np.log(np.maximum(v[a_:N + 1], 1e-300)) - ns * math.log(3) / 2
                rs.append(math.exp(np.polyfit(ns, y, 1)[0]))
            rates.append(rs)
            P(f"  q={q} seed={s}: rate 300..{N} = {rs[0]:.4f}   750..{N} = {rs[1]:.4f}   [{time.time()-t0:.0f}s]")
        arr = np.array(rates)
        P(f"  q={q} uniform-start: mean rate 300..{N} = {arr[:,0].mean():.4f} +- {arr[:,0].std(ddof=1):.4f}; 750..{N} = {arr[:,1].mean():.4f} +- {arr[:,1].std(ddof=1):.4f}")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
