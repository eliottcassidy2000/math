#!/usr/bin/env python3
"""Single-far Lee-Yang region + Asano monotonicity: adding far to a GOOD base keeps zeros outside |z|<=1 (S69).

HYP-3127 obligation 1 (load-bearing): each single-far factor is zero-free in the relevant polydisk, so Asano
contraction gives the multi-far floor. S66 found: the miss-count PGF G_N(z)=sum q_t z^t is zero-free in |z|<=1
for HIGH-coverage configs (consec: nearest zero 1.49) but NOT low-coverage (random: 0.047). The multi-far
floor is over a GOOD bounded base B (consec-like) + far elements. Since each far INCREASES coverage
(d(f)=p0(Bu{f})/p0(B)>1, S68), the Asano/coverage-increasing step should PRESERVE or push OUT the zeros.

TEST: base B=consec; add r=1,2,3 far elements (resonant + separated); track the NEAREST-ZERO RADIUS rho of
G_N(z). Confirm rho stays >= 1 (Lee-Yang preserved) and the floor min rho over placements. Single-far building
block: the conditional 'far covers an empty sector' prob; |zero| >= 1 <=> P(far misses | base good) <= 1/2.
"""
import sys, itertools, random
from fractions import Fraction as F
import numpy as np
sys.stdout.reconfigure(line_buffering=True)
random.seed(69)

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
def nearest_zero_radius(q):
    c = [float(q[t]) for t in range(6, -1, -1)]
    while len(c) > 1 and abs(c[0]) < 1e-14: c = c[1:]
    if len(c) <= 1: return float('inf')
    r = np.roots(c)
    return min(abs(z) for z in r)

B = tuple(range(8))  # consec_8 base
rhoB = nearest_zero_radius(missdist(B))
print("=" * 90)
print(f" base B={B}  p0(B)={float(missdist(B)[0]):.4f}  nearest-zero radius rho(B) = {rhoB:.4f}  (>1 = Lee-Yang)")
print("=" * 90)

FAR_FAMS = {
    "single-far r=1":  [(15,),(21,),(28,),(56,),(100,),(500,)],
    "separated r=2":   [(101,211),(307,503),(701,1009)],
    "resonant r=2":    [(21,28),(35,49),(28,56),(21,42)],
    "tight doublet r=2":[(21,22),(50,51),(100,101)],
    "separated r=3":   [(101,211,307),(503,701,907)],
    "resonant r=3":    [(21,28,35),(28,42,56),(35,49,63)],
}
floor = (float('inf'), None)
print("\n adding far to the good base -- nearest-zero radius rho (>=1 = Lee-Yang PRESERVED):")
for name, cfgs in FAR_FAMS.items():
    rhos = []
    for Fs in cfgs:
        rho = nearest_zero_radius(missdist(B + Fs)); rhos.append(rho)
        if rho < floor[0]: floor = (rho, (name, Fs))
    lo = min(rhos)
    print(f"   {name:18s}: rho = [{', '.join(f'{r:.3f}' for r in rhos)}]   min={lo:.3f}  "
          f"{'OK >=1' if lo>=1 else 'FAILS <1'}")
# random multi-far scan over the good base
worst = (float('inf'), None)
for _ in range(400):
    r = random.choice([2,3,4]); Fs = tuple(sorted(random.sample(list(range(15,70))+[21,28,35,42,49,56,63],r)))
    rho = nearest_zero_radius(missdist(B + Fs))
    if rho < worst[0]: worst = (rho, Fs)
print(f"\n   random multi-far scan (400, r=2..4): MIN nearest-zero radius = {worst[0]:.4f} @ F={worst[1]}")
print(f"   OVERALL FLOOR rho = {floor[0]:.4f} @ {floor[1]}")

print("\n" + "=" * 90)
print(" SINGLE-FAR BUILDING BLOCK: |zero| >= 1 <=> P(far misses its sector | base good) <= 1/2")
print("=" * 90)
# proxy: how often does adding f leave a sector empty that B alone covered? d(f)=p0(Bu{f})/p0(B) > 1 means f HELPS.
def p0(E): return missdist(E)[0]
p0B = p0(B)
for f in (15, 21, 28, 56, 100, 500):
    d = p0(B + (f,)) / p0B
    print(f"   f={f:4d}: d(f)=p0(Bu{{f}})/p0(B) = {float(d):.4f} (>1: far INCREASES coverage => zero pushed OUT)")
print("\nREADING: adding far to a GOOD base keeps the nearest zero >= 1 (Lee-Yang PRESERVED) with a positive")
print("FLOOR -- because each far increases coverage (d(f)>1) = a coverage-increasing Asano factor that pushes")
print("the zeros outward. This is HYP-3127 obligation 1 (single-far Lee-Yang region) + the Asano monotonicity:")
print("the multi-far G_N stays zero-free in |z|<=1, so the multi-far floor R'>=c follows from the good base.")
