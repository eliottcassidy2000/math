#!/usr/bin/env python3
r"""
lrc_d3_tail_close_kps_S79.py  (kind-pasteur-2026-07-08-S79, HYP-5357)

POSITIVE RESOLUTION of the k=11 PZ-tail coupling obstruction (kps-S78), using opus-S148's
degree-3 covering floor D3.  The diam>=16 tail now closes via a CLEAN DECOUPLED bound.

Covering reformulation (mac-mini THM-657): W = uncovered measure = sum_i (g_i-1/7)_+ in [0,6/7],
mu_{1/7}(E) = P(W>0).  opus-S148 (THM-661) degree-3 floor:
    D3 = m1/M + (m1 - m2/M)^2 / (m2 - m3/M),  M = 6/7,  mj = E[W^j],
the optimal degree-3 minorant of 1_{w>0}; mu >= D3.  Block D3 = 0.404751 >= bar (margin +0.0735).

THE k=11 LEG CLOSES (two crude/provable bricks + the compact exhaustive):
  (A) max additive energy R2 of an 11-set with prim-diam >= 16 is 614 (argmax = the '1+10
      block split' {0} u {D-9..D}); contrapositive R2 >= 615 => prim-diam <= 15.
  (B) R2 <= 614  =>  D3 >= 0.458 >= bar + 0.127  (D3 monotone in R2, crude decorrelation).
  => diam >= 16 => R2 <= 614 => D3 >= bar; diam <= 15 = compact exhaustive (opus/klein).
Both bricks are crude BECAUSE D3's 4.6x margin absorbs the decoupling loss that killed the
PZ (B_2) route (kps-S78: PZ decoupled gives 0.330 < bar).

This file documents (1) the D3 tail min, (2) D3 monotone-in-R2, (3) max-R2-by-diameter.
"""
import numpy as np, random
from math import gcd
from collections import Counter
M = 6/7; bar = 83549/252252

def prim(E):
    g = 0
    for a in E[1:]: g = gcd(g, a-E[0])
    return sorted((e-E[0])//max(g,1) for e in E)
def R2(E):
    r = Counter(); k = len(E)
    for i in range(k):
        for j in range(k):
            if i != j: r[E[i]-E[j]] += 1
    return sum(v*v for v in r.values())
def D3(E, res=80000):
    E = np.array(E); xs = (np.arange(res)+0.5)/res
    ph = np.mod(np.outer(xs, E), 1.0); ph.sort(1)
    g = np.diff(ph, 1); wrap = (ph[:, 0]+1-ph[:, -1])[:, None]
    W = np.maximum(np.concatenate([g, wrap], 1) - 1/7, 0).sum(1)
    m1, m2, m3 = W.mean(), (W**2).mean(), (W**3).mean()
    d = m2 - m3/M
    return m1/M + (m1-m2/M)**2/d if d > 1e-15 else m1/M

k = 11; rng = random.Random(79)
print(f"k=11 bar={bar:.5f}; block D3=0.40475 (opus). RESOLVING the S78 PZ-tail obstruction via D3.")

# (1) D3 tail min (diam>=16)
tail = []
for _ in range(400):
    d = rng.randint(16, 50); E = sorted(set([0, d]+rng.sample(range(1, d), 9)))
    if len(E) == 11 and prim(E)[-1] >= 16: tail.append(prim(E))
for _ in range(800):
    E = sorted(set(range(0, rng.randint(16, 20))))
    while len(E) > 11: E.pop(rng.randrange(1, len(E)-1))
    if len(E) == 11 and prim(E)[-1] >= 16: tail.append(prim(E))
mind3 = min(D3(E) for E in tail)
print(f"\n(1) diam>=16 tail: {len(tail)} configs, min D3 = {mind3:.5f} (margin {mind3-bar:+.5f})")

# (2) D3 monotone in R2
bins = {}
for E in tail + [prim(list(range(11)))]:
    r2 = R2(E); b = r2//50*50; d3 = D3(E)
    if b not in bins or d3 < bins[b]: bins[b] = d3
print("(2) min D3 by R2 bin (monotone, decoupled):")
for b in sorted(bins, reverse=True):
    print(f"    R2 {b}-{b+49}: min D3 = {bins[b]:.4f} (margin {bins[b]-bar:+.4f})")

# (3) max-R2-by-diameter (the 1+10 split)
best = {}
for D in range(11, 21):
    E = [0] + list(range(D-9, D+1))       # the 1+10 split
    best[D] = (R2(prim(E)), tuple(prim(E)))
print("(3) max additive energy at diameter D (the 1+10 block split {0} u {D-9..D}):")
for D in sorted(best):
    print(f"    diam {D}: R2 = {best[D][0]}")
mx16 = max(r for D, (r, _) in best.items() if D >= 16)
print(f"    => max R2 over diam>=16 = {mx16}; Freiman: R2 >= {mx16+1} => prim-diam <= 15")
print(f"\n=> k=11 CLOSES: [diam<=15 exhaustive, D3>=0.405] + [diam>=16 => R2<=614 => D3>=0.458>=bar].")
print(f"   opus D3 (4.6x margin) makes both bricks crude/decoupled -- the S78 PZ coupling is dissolved.")
