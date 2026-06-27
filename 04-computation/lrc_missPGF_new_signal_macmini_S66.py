#!/usr/bin/env python3
"""NEW SIGNAL: the miss-count PGF  G_N(z) = E[z^N] = sum_t q_t z^t  and its ROOT structure (S66).

Research finding (INVESTIGATION-BACKLOG:1520, "Route 5"): LRC and tournament-H are two evaluations of ONE
fugacity partition function  G_x = sum_k x^k MISS_k, consec=argmax at x=-1 (LRC) but WORST at x=+2 (tourn).
Since MISS_k = S_k (factorial moments) and sum_k S_k x^k = E[(1+x)^N], this IS the PGF of the miss-count
N = #empty inner sectors, in the variable z = 1+x:
   G_N(z) = sum_{t=0}^{6} q_t z^t ,   q_t = P(N=t).
   z=0  -> q_0 = p0 (LRC coverage);   z=2 -> E[2^N];   z=3 -> E[3^N] (the tournament fugacity).

The project measures q_t and the moments S_r and the single value p0=G_N(0). It has NEVER measured the
WHOLE CURVE / the ROOT STRUCTURE of G_N(z). New signals proposed here:
  (1) the 6 roots of G_N(z) in C; real-rootedness; the root NEAREST z=0 (the 'coverage pole proximity').
  (2) log-concavity / Newton defect of q (q_t^2 - q_{t-1} q_{t+1}).
  (3) the RANK-FLIP: rank of each config by G_N(z) as z runs 0 -> 3 (does consec flip best->worst?).
Bold guess under test: consec's extremality at z=0 is governed by an EXTREMAL root structure of G_N.
"""
import sys, itertools, cmath
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def sector_of(p): return int((p % 1) * 7)
def missdist(E):
    E = sorted(set(E)); b = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): b.add(F(m, 7 * e))
    b = sorted(b); q = [F(0)] * 7
    for i in range(len(b) - 1):
        x0, x1 = b[i], b[i + 1]
        if x1 <= x0: continue
        t = 7 - len(set(sector_of(e * ((x0 + x1) / 2)) for e in E))
        if 0 <= t <= 6: q[t] += x1 - x0
    return q

def pgf_roots(q):
    import numpy as np
    coeffs = [float(q[t]) for t in range(6, -1, -1)]  # numpy.roots wants highest degree first
    # strip leading zeros
    while len(coeffs) > 1 and abs(coeffs[0]) < 1e-15: coeffs = coeffs[1:]
    r = np.roots(coeffs)
    return sorted(r, key=lambda z: abs(z))

def evalz(q, z): return sum(float(q[t]) * z**t for t in range(7))

# configs at the binding row k=8 (the hardest), and a few others
K = 8
CONFIGS = {
    "consec_8":        tuple(range(8)),
    "even-AP 2i":      tuple(2*i for i in range(8)),
    "break {1,5,7,8,9}+0": (0,1,5,7,8,9,11,13),     # middle-spread flavor, size 8
    "top-cluster":     (0,7,8,9,10,11,12,13),
    "random-spread":   (0,1,3,7,12,20,33,54),
    "covering 0,14..": (0,14,1,2,3,5,8,13),
}
print("=" * 96)
print(f" MISS-COUNT PGF  G_N(z)=sum q_t z^t   (k={K})   -- a NEW signal: roots, log-concavity, rank-flip")
print("=" * 96)
rows = {}
for name, E in CONFIGS.items():
    q = missdist(E)
    roots = pgf_roots(q)
    nearest = min((r for r in roots), key=lambda z: abs(z))
    n_real = sum(1 for r in roots if abs(r.imag) < 1e-9)
    # log-concavity / Newton defect: min_t (q_t^2 - q_{t-1} q_{t+1}); >=0 means log-concave
    lc = min(float(q[t])**2 - float(q[t-1])*float(q[t+1]) for t in range(1, 6))
    g0, g2, g3 = evalz(q, 0.0), evalz(q, 2.0), evalz(q, 3.0)
    rows[name] = dict(q=q, roots=roots, nearest=nearest, n_real=n_real, lc=lc, g0=g0, g2=g2, g3=g3)
    print(f"\n {name}:")
    print(f"   q = [{', '.join(f'{float(x):.4f}' for x in q)}]   (P(N=t), t=0..6)")
    print(f"   G_N(0)=p0={g0:.5f}   G_N(2)={g2:.4f}   G_N(3)={g3:.4f}")
    print(f"   roots |z| (6): [{', '.join(f'{abs(r):.3f}' for r in roots)}]   #real={n_real}/6   "
          f"nearest root |z*|={abs(nearest):.4f}")
    print(f"   log-concavity defect min(q_t^2 - q_(t-1)q_(t+1)) = {lc:+.5f}  {'(LOG-CONCAVE)' if lc>=0 else '(NOT log-concave)'}")

print("\n" + "=" * 96)
print(" RANK-FLIP: order the configs by G_N(z) at z=0 (LRC), z=2, z=3 (tournament fugacity)")
print("=" * 96)
for z, lbl in [(0.0, "z=0 (LRC p0)"), (2.0, "z=2"), (3.0, "z=3 (tourn)")]:
    order = sorted(rows, key=lambda n: -evalz(rows[n]["q"], z))
    pos = order.index("consec_8") + 1
    print(f"  {lbl:16s}: argmax={order[0]:18s}  consec_8 rank={pos}/{len(order)}   "
          f"order={[n[:8] for n in order]}")
print("\n  -> if consec flips from a top rank at z=0 to a low rank at z=3, the fugacity rank-flip (Route 5)")
print("     is reproduced; the NEW content is whether consec's ROOT structure (#real, nearest |z*|,")
print("     log-concavity) is extremal -- the analytic fingerprint behind 'consec is extremal at z=0'.")
