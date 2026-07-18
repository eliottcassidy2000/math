#!/usr/bin/env python3
"""
The fill-1 perturbation lemma -- base case  (boxeph-2026-07-17-S82)
==================================================================

CRUX (from THM-998 / the resonance-fill reflection): a covering family has no
EMPTY circle (b<=14), so loneliness is immediate unless the binding obstruction
is an UNDER-FILLED circle f_b=1 -- one runner v* stranded on the origin at the
resonance center a/b, which must be perturbed off without dropping another
runner below 1/N.

THE PERTURBATION.  Sit at t0 = a/b (gcd(a,b)=1).  b|v* => v* at 0; every other
runner (b does not divide it) is at distance >= 1/b.  Kick by the MINIMAL amount
that lifts v* to exactly 1/N:  t = a/b + 1/(N v*).  Then:
   * ||v* t|| = ||integer + 1/N|| = 1/N            (exactly on threshold)
   * body runner v:  ||v t|| >= ||v a/b|| - ||v/(N v*)|| >= 1/b - v/(N v*)
                     (reverse triangle inequality on R/Z)
   >= 1/N   iff   v <= (N-b)/b * v*.

So the WITNESS CERTIFIES iff  B := max{body v} satisfies  b*B <= (N-b)*v*.
This script (exact rationals) verifies the lemma, maps its reach, checks tightness
at the boundary, and confirms it CLOSES the single-killer far-element families
(incl. the deep-well extremal) while leaving the bounded/compact core as residual.
"""
from fractions import Fraction as F
from math import gcd

def norm(x):
    """distance from rational x to nearest integer, exact."""
    r = x - int(x)               # in (-1,1)
    r = r % 1                    # in [0,1)
    return min(r, 1 - r)

def min_gap(V, t):
    return min(norm(v * t) for v in V)

def fill(V, b):
    return sum(1 for v in V if v % b == 0)

def is_covering(V, N=14):
    return all(any(v % b == 0 for v in V) for b in range(2, N+1))

def witness(a, b, vstar, N=14):
    return F(a, b) + F(1, N * vstar)

def certifies(V, b, N=14, a=1):
    """Return (ok, t, mingap) for the fill-1 perturbation at circle b, center a/b."""
    mults = [v for v in V if v % b == 0]
    if len(mults) != 1:
        return (None, None, None)   # not a fill-1 circle
    vstar = mults[0]
    B = max(v for v in V if v % b != 0)
    cond = (b * B <= (N - b) * vstar)
    t = witness(a, b, vstar, N)
    mg = min_gap(V, t)
    return (cond, t, mg)

N = 14
print("="*76)
print("FILL-1 PERTURBATION LEMMA (base case)   threshold 1/%d" % N)
print("witness t = a/b + 1/(N v*);  certifies iff  b*B <= (N-b)*v*   (B = body max)")
print("="*76)

# ---- 1. The deep-well extremal and the single-killer family ----------------
print("\n[1] Single-killer covering families {1..12, w}, 182|w  (use circle b=13):")
for w in [182, 364, 546, 728, 1820]:
    V = list(range(1,13)) + [w]
    cond, t, mg = certifies(V, 13, N)
    print(f"   w={w:5d}: covering={is_covering(V)} cond(13*12<=1*{w})={13*12}<={w}:{cond}"
          f"  witness t={t}  min_gap={mg} (>=1/14: {mg>=F(1,14)})")

# ---- 2. Reach map: for body {1..12}, smallest killer per circle b ----------
print("\n[2] Body {1..12}: at circle b, minimal killer v* with b*12<=(14-b)v*:")
for b in range(2, 14):
    B = 12
    vmin = -(-b*B // (N-b)) if N-b>0 else None   # ceil(bB/(N-b))
    # round up to a multiple of b (v* must be divisible by b)
    if vmin is not None:
        vstar = ((vmin + b - 1)//b)*b
        vstar = max(vstar, b)  # at least b
        # ensure fill_b=1: vstar not creating another multiple; body 1..12 has mult of b if b<=12
        print(f"   b={b:2d}: (14-b)={N-b:2d}  min killer >= ceil({b}*12/{N-b})={vmin}"
              f"  -> smallest valid multiple-of-b v*={vstar}")
    else:
        print(f"   b={b:2d}: (14-b)=0  NO slack at threshold denominator (b=14 undodgeable this way)")

# ---- 3. Tightness: at the boundary the naive witness just fails ------------
print("\n[3] Tightness check at circle b=13, body {1..12} (need v* >= 156):")
for w in [143, 156, 169, 182]:  # 143=11*13,156=12*13,169=13*13,182=14*13 (multiples of 13)
    V = list(range(1,13)) + [w]
    cond, t, mg = certifies(V, 13, N)
    verdict = "CERT" if (mg is not None and mg>=F(1,14)) else "naive-fail"
    print(f"   w={w}: cond={cond}  min_gap(naive t)={mg}  -> {verdict}"
          f"   (w>=156: {w>=156})")

# ---- 4. Does the lemma reach REAL covering families beyond single-killer? --
print("\n[4] Assorted covering families -- which have a fill-1 circle the lemma dodges:")
fams = {
    "deep well {1..12,182}   ": list(range(1,13))+[182],
    "residue {1..11,13,84}   ": list(range(1,12))+[13,84],
    "{1..10,13,22,84} (k=10) ": list(range(1,11))+[13,22,84],
    "missing-2 {1,3..14}     ": [1]+list(range(3,15)),
    "{1..12,14} bounded core ": list(range(1,13))+[14],
}
for name, V in fams.items():
    cov = is_covering(V)
    hit = None
    for b in range(2, 14):
        cond, t, mg = certifies(V, b, N)
        if cond is True and mg is not None and mg >= F(1,14):
            hit = (b, t, mg); break
    print(f"   {name} covering={cov}: "
          + (f"DODGED at b={hit[0]} (t={hit[1]}, min={hit[2]})" if hit
             else "NO fill-1 dodge -> residual (bounded/compact core or multi-killer)"))

print("\n" + "="*76)
print("READING: the lemma CLOSES the isolated-far-element regime (single killer")
print("v* >= b*B/(N-b), incl. the deep-well and residue extremals) by an ELEMENTARY")
print("direct witness. Residual = bounded core (v* too small, e.g. {1..12,14}) and")
print("multi-killer (no single fill-1 circle with a dominant stranded runner).")
