#!/usr/bin/env python3
"""
E5k (opus, 2026-10-04): is the q = 3 cold-rate anomaly (0.5700 against the universal 0.5725) seeded by the low
levels (where the phase string of 3^-n is periodic inside the window, n <= 10) or a steady-state property?

Hybrid towers: (A) real digits of -u 3^-n for n > K, i.i.d. digits for n <= K (K = 20, 100);
               (B) i.i.d. digits for n > K, real digits for n <= K (K = 20, 100);
               (C) real digits of 3^-n at every level but with a RANDOM extra unit per level (phases of -u_n 3^-n
                   with u_n fresh random odd 40-bit): levels unrelated by the cube-root law but each level a genuine
                   3^-n string.
Rates over 300..N for u = 1, 5, 7 (A, B) and three seeds (C).  Reference: q = 3 real 0.5700 +- 0.0009 (six units),
i.i.d. 0.5723 +- 0.0008 (twelve seeds).
Usage: python collatz_fixed_frequency_hybrid_20261004.py [N=1500] [A=40]
"""
import sys, math, time
import numpy as np
sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from collatz_fixed_frequency_random_digits_20261004 import phases_from_R

def run_hybrid(N, A, mode, K=20, u=1, seed=0):
    rng = np.random.default_rng(seed)
    w = math.sqrt(3.0) * 2.0 ** (-np.arange(1, A + 1))
    prev_lo = -N * A - A
    prev = np.ones(-prev_lo + 1, dtype=np.complex128)
    v0 = np.zeros(N + 1); v0[0] = 1.0
    for n in range(1, N + 1):
        lo = -(N - n) * A
        M = -(lo - A)
        real_here = (mode == "A" and n > K) or (mode == "B" and n <= K) or mode == "C"
        if real_here:
            uu = u if mode != "C" else (int(rng.integers(1, 1 << 40)) | 1)
            R = (-uu * pow(3, -n, 1 << (M + 1))) % (1 << (M + 1))
        else:
            R = int.from_bytes(rng.bytes((M + 1) // 8 + 1), "little") % (1 << (M + 1)) | 1
        ph = np.empty(M + 1, dtype=np.complex128); ph[:M] = phases_from_R(R, M); ph[M] = 1.0
        Pp = ph * prev[(lo - A) - prev_lo: (0 - prev_lo) + 1]
        cur = np.zeros(-lo + 1, dtype=np.complex128)
        for a in range(1, A + 1):
            cur += w[a - 1] * Pp[A - a: A - a + (-lo + 1)]
        prev, prev_lo = cur, lo
        v0[n] = abs(cur[-lo])
    return v0

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"hybrid towers for q = 3, N={N}, A={A}; rates over 300..{N} (reference: real 0.5700 +- 0.0009, i.i.d. 0.5723 +- 0.0008)")
    t0 = time.time()
    def rate(v):
        ns = np.arange(300, N + 1); y = np.log(np.maximum(v[300:N + 1], 1e-300)) - ns * math.log(3) / 2
        return math.exp(np.polyfit(ns, y, 1)[0])
    for mode, K in (("A", 20), ("A", 100), ("B", 20), ("B", 100)):
        rs = [rate(run_hybrid(N, A, mode, K=K, u=u, seed=7 + u)) for u in (1, 5, 7)]
        P(f"  mode {mode} K={K:3d} ({'real above K, iid below' if mode == 'A' else 'iid above K, real below'}): rates u=1,5,7: " + " ".join(f"{r:.4f}" for r in rs) + f"  mean {np.mean(rs):.4f}   [{time.time()-t0:.0f}s]")
    rs = [rate(run_hybrid(N, A, "C", seed=s)) for s in range(3)]
    P(f"  mode C (every level a genuine 3^-n string, fresh random unit per level): rates: " + " ".join(f"{r:.4f}" for r in rs) + f"  mean {np.mean(rs):.4f}   [{time.time()-t0:.0f}s]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
