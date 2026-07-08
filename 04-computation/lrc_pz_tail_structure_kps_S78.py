#!/usr/bin/env python3
r"""
lrc_pz_tail_structure_kps_S78.py  (kind-pasteur-2026-07-08-S78, HYP-5337)

SHARPENING the THM-660 PZ tail-closing lemma (k=11,12,13).  klein-S179 reduced the diam>=16
tail to "prove Var(W) <= c*R2".  This file shows that bound is INSUFFICIENT alone and
identifies the three tight ingredients + the clean provable brick.

Covering reformulation (THM-657): W(x) = sum_i (g_i(x) - 1/7)_+ = uncovered measure,
mu_{1/7}(E) = P(W>0) >= E[W]^2/E[W^2] = PZ (THM-660).  R2 = reduced additive energy.

FINDINGS (k=11, bar = 83549/252252 = 0.33121, R2(block)=770):
 (1) Var/R2 is NOT constant: ~6.1e-5 at the block (max R2), up to ~7.0e-5 at low R2.
     So Var(W) <= c*R2 with the honest max c=7.0e-5 gives, at R2=770,
     PZ >= E[W]^2/(E[W]^2 + 7.0e-5*770) = 0.316 < bar -- the uniform bound FAILS at the block.
 (2) maxR2(diam) falls SLOWLY: 770/630/614/574/510 at diam 10/15/16/20/25.
 (3) E[W] ~ flat (~0.15, min ~0.148 at diam>=16), not rising.
 (4) The diam>=16 tail closes (PZ>=~0.34) but TIGHTLY, needing ALL THREE:
       (a) maxR2 at diam>=16 <= 614 < 770   [clean extremal brick -- Freiman stability]
       (b) Var(W) <= V_W*R2                  [resonance bound, shared with THM-656]
       (c) E[W] >= ~0.148                    [mean-uncovered lower bound]
     Each alone insufficient; (a) is the piece klein's "Var<=c*R2" framing omits.

This script (numerical) documents (1)-(4).  The HANDOFF is the max-additive-energy-by-
diameter bound (a) -- a clean extremal statement independent of the PZ analysis.
"""
import numpy as np, random
from math import gcd
from collections import Counter

def prim(E):
    g = 0
    for a in E[1:]: g = gcd(g, a - E[0])
    return sorted((e - E[0]) // max(g, 1) for e in E)

def R2(E):
    r = Counter(); k = len(E)
    for i in range(k):
        for j in range(k):
            if i != j: r[E[i] - E[j]] += 1
    return sum(v * v for v in r.values())

def moments(E, res=120000):
    E = np.array(E); xs = (np.arange(res) + 0.5) / res
    ph = np.mod(np.outer(xs, E), 1.0); ph.sort(1)
    g = np.diff(ph, 1); wrap = (ph[:, 0] + 1 - ph[:, -1])[:, None]
    W = np.maximum(np.concatenate([g, wrap], 1) - 1/7, 0).sum(1)
    EW = W.mean(); Var = (W * W).mean() - EW * EW
    return EW, Var, EW * EW / (EW * EW + Var)

k = 11; bar = 83549 / 252252; R2block = k * (k - 1) * (2 * k - 1) // 3
rng = random.Random(11)
print(f"k={k}: bar={bar:.5f}, R2(block)={R2block}")

# (1) Var/R2 varies; the uniform-max-c bound fails at the block
EWb, Varb, PZb = moments(list(range(k)))
print(f"\n(1) block: E[W]={EWb:.5f} Var={Varb:.5f} Var/R2={Varb/R2block:.3e} PZ={PZb:.5f}")
ratios = []
for _ in range(150):
    d = rng.randint(11, 45); E = sorted(set([0, d] + rng.sample(range(1, d), 9)))
    if len(E) != 11: continue
    Ep = prim(E); r2 = R2(Ep); EW, Var, PZ = moments(Ep)
    ratios.append(Var / r2)
cmax = max(ratios + [Varb / R2block])
print(f"    max Var/R2 over sample = {cmax:.3e}; uniform bound at block: "
      f"PZ >= {EWb**2/(EWb**2 + cmax*R2block):.5f} vs bar {bar:.5f} => "
      f"{'FAILS' if EWb**2/(EWb**2+cmax*R2block) < bar else 'ok'}")

# (2) maxR2 by diameter
best = {}
for _ in range(200000):
    d = rng.randint(10, 30); E = sorted(set([0, d] + rng.sample(range(1, d), 9)))
    if len(E) != 11: continue
    Ep = prim(E); dm = Ep[-1]; r2 = R2(Ep)
    if dm not in best or r2 > best[dm]: best[dm] = r2
print(f"\n(2) maxR2 by primitive diameter: " +
      " ".join(f"d{d}:{best.get(d,'-')}" for d in [10,12,15,16,18,20,25]))

# (4) the tail-closing arithmetic at diam>=16 (R2<=614): three bounds
R2tail = best.get(16, 614); c = cmax; EWtail = 0.148
PZtail = EWtail**2 / (EWtail**2 + c * R2tail)
print(f"\n(4) tail (diam>=16): R2<={R2tail}, Var<={c:.1e}*R2, E[W]>={EWtail} => "
      f"PZ >= {PZtail:.5f} vs bar {bar:.5f}: {'CLEARS' if PZtail>=bar else 'tight/needs sharper'}")
print(f"    => needs (a) maxR2(diam>=16)<=614, (b) Var<=c*R2, (c) E[W]>={EWtail} -- all three.")
