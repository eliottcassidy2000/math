#!/usr/bin/env python3
"""
LRC(14) S3 -- AIRTIGHT verification of the cluster-collapse single-gap lemma's claims,
and the multi-band structural theorem (monotone sub-band split).

CLAIM 1 (cluster window safety). For integers Vmin<=u<=Vmax and integer K, on the open
interval W_K = ((K+1/14)/Vmin, (K+13/14)/Vmax), every runner u in [Vmin,Vmax] satisfies
||u t|| >= 1/14.
  Proof: for t in W_K, u*t in (u(K+1/14)/Vmin, u(K+13/14)/Vmax). Since u>=Vmin:
  u(K+1/14)/Vmin >= K+1/14. Since u<=Vmax: u(K+13/14)/Vmax <= K+13/14. So u*t in
  [K+1/14, K+13/14], i.e. frac(u t) in [1/14,13/14] (mod 1, single integer K), giving
  ||u t|| >= 1/14. QED.
  We VERIFY by sampling exact rationals t in W_K and ALL u in [Vmin,Vmax].

CLAIM 2 (nonemptiness). W_K nonempty iff (K+13/14)/Vmax > (K+1/14)/Vmin
  iff (14K+13) Vmin > (14K+1) Vmax  iff 13 Vmin - Vmax > 14 K (Vmax - Vmin) = 14 K s.
  We VERIFY the iff arithmetically.

CLAIM 3 (multi-band structural theorem). For ANY t>0 and speeds u<u', floor(u t) <=
floor(u' t) -- so at a fixed witness the gap-index map u -> floor(u t) is monotone
non-decreasing. Hence the cluster partitions into CONTIGUOUS sub-bands by gap index.
  Trivial (floor is monotone in u for t>0). Already verified empirically 2445/2445.

We re-verify CLAIM 1 and CLAIM 2 exactly on many (Vmin,Vmax,K) including adversarial.
"""
from fractions import Fraction as F
from math import gcd, floor
import random

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def verify_window_safety(Vmin, Vmax, K, samples=40):
    lo = F(14*K+1, 14) / Vmin
    hi = F(14*K+13, 14) / Vmax
    if lo >= hi: return None  # empty
    # sample t exactly across (lo,hi), check ALL u in [Vmin,Vmax]
    bad = []
    for s in range(1, samples):
        t = lo + (hi - lo) * F(s, samples)
        for u in range(Vmin, Vmax + 1):
            if nrm(u * t) < F(1, 14):
                bad.append((u, t)); break
    # also check endpoints' interior limit: midpoint
    tm = (lo + hi) / 2
    for u in range(Vmin, Vmax + 1):
        if nrm(u * tm) < F(1, 14): bad.append((u, tm)); break
    return (lo, hi, len(bad), bad[:3])

if __name__ == "__main__":
    print("CLAIM 2 (nonemptiness iff 13 Vmin - Vmax > 14 K s) -- arithmetic identity:")
    rng = random.Random(1); viol2 = 0; tested2 = 0
    for _ in range(200000):
        Vmin = rng.randint(14, 2000); s = rng.randint(0, 120); Vmax = Vmin + s
        K = rng.randint(0, 50)
        lo = F(14*K+1, 14)/Vmin; hi = F(14*K+13, 14)/Vmax
        nonempty = lo < hi
        cond = (13*Vmin - Vmax > 14*K*s)
        tested2 += 1
        if nonempty != cond: viol2 += 1
    print(f"   tested {tested2}, mismatches: {viol2}  (0 = identity holds)")

    print("\nCLAIM 1 (cluster window safety) -- exact sampling, all u in [Vmin,Vmax]:")
    rng = random.Random(2); tested1 = 0; bad1 = 0; examples = 0
    for _ in range(4000):
        Vmin = rng.randint(14, 400); s = rng.randint(1, 60); Vmax = Vmin + s
        if Vmax - Vmin > 80: continue
        K = rng.randint(0, 30)
        r = verify_window_safety(Vmin, Vmax, K, samples=25)
        if r is None: continue
        lo, hi, nbad, bad = r
        tested1 += 1
        if nbad > 0:
            bad1 += 1
            if examples < 5:
                print(f"   *** UNSAFE: Vmin={Vmin} Vmax={Vmax} K={K} bad={bad}")
                examples += 1
    print(f"   tested {tested1} nonempty windows, windows with an unsafe sampled point: {bad1}")
    print(f"   (0 = cluster window safety verified)")

    print("\nCLAIM 3 (gap-index monotone in speed for t>0) -- trivial (floor monotone).")
    rng = random.Random(3); viol3 = 0
    for _ in range(100000):
        t = F(rng.randint(1, 999), 1000)
        u = rng.randint(1, 5000); up = u + rng.randint(1, 5000)
        if floor(u*t) > floor(up*t): viol3 += 1
    print(f"   floor(u t) <= floor(u' t) for u<u', t>0: violations {viol3}/100000")
