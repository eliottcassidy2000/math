#!/usr/bin/env python3
r"""
lrc_three_heresies_kps_S63.py   (kind-pasteur-2026-07-07-S63, HYP-4897)

THREE HERESIES -- creative probes that challenge standing assumptions (owner directive).

PROBE A (SCHUR / EXCHANGE-MONOTONICITY -- would restructure (A')):
  Assumption challenged: (A') "the AP minimizes mu" is treated as a global extremal
  statement.  Heresy: maybe mu_theta is monotone under POSITION-FIXED Robin-Hood step
  transfers toward equality: word w with delta_i >= delta_j + 2 vs w' = transfer one unit
  from i to j.  If mu(w') <= mu(w) ALWAYS, then iterating transfers reaches the AP and
  (A') becomes a LOCAL two-step exchange lemma (Schur-convexity shape).
  Test: EXHAUSTIVE and EXACT at diameters D = 13, 14, 15 (all C(D-1,11) compositions of D
  into 12 positive parts), exact order-cell mu engine (opus-S136's method, reimplemented
  compactly).  Count violations; if nonzero, characterize by position pair (Wiener weight?).

PROBE B (THE LATTICE-MINIMUM LAW -- the sparse-lane theorem shape):
  Assumption challenged: "decorrelation" is treated as analytic mush.  Heresy: the deficit
  iid - E[maxgap] should be controlled by the SHORTEST zero-sum weight->=3 relation of E
  (Poisson summation over the relation lattice; pair modes vanish by S59).  Measure
  deficit vs min-relation size over a zoo + random families; fit the envelope.

PROBE C (IS THE MEAN FLOOR A 2-POINT THEOREM? -- methodologically decisive either way):
  Assumption challenged: everyone hunts proofs using specific arithmetic structure.  If the
  mean floor E[maxgap] > 1/7 followed from the pair-uniformity alone (every integer family
  has EXACTLY uniform unlabeled pair-distance density), a Cohn-Elkies-style 2-point proof
  would exist.  Test: try to CONSTRUCT a synthetic cyclic 13-point process (not from
  integers!) with uniform aggregate pair-distance density and E[maxgap] < 1/7.
  Construction family: mixtures of jittered lattices + compensators, weights optimized by
  projected least squares (numpy).  If a valid construction dips below 1/7: NO 2-point
  proof exists; the weight->=3 usage is NECESSARY.  If it stalls above: conjecture the
  2-point LP floor.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import random

# ---------------------------------------------------------------- exact mu engine
def mu_exact(E, theta=F(1, 7)):
    """Exact mu_theta for a set of positive integers E (order-cell method)."""
    E = sorted(E)
    k = len(E)
    dmax = E[-1] - E[0] if E[0] != 0 else E[-1]
    dmax = max(dmax, max(E[i+1]-E[i] for i in range(k-1)) if k > 1 else 1, max(E))
    bps = {F(0), F(1)}
    for d in range(1, dmax + 1):
        for m in range(1, d):
            bps.add(F(m, d))
    bps = sorted(bps)
    tot = F(0)
    for a, b in zip(bps[:-1], bps[1:]):
        mid = (a + b) / 2
        fl = {e: (e * mid).__floor__() for e in E}
        order = sorted(E, key=lambda e: e * mid - fl[e])
        segs = []
        for s in range(k):
            if s < k - 1:
                e1, e2 = order[s], order[s + 1]
                c = e2 - e1; b0 = -(fl[e2] - fl[e1])
            else:
                c = order[0] - order[-1]; b0 = -(fl[order[0]] - fl[order[-1]]) + 1
            segs.append((c, F(b0)))
        # maxgap(x) = max_s (c x + b0); superlevel > theta within (a,b)
        # collect crossing points of each seg with theta and pairwise argmax switches
        pts = {a, b}
        for c, b0 in segs:
            if c != 0:
                xc = (theta - b0) / c
                if a < xc < b: pts.add(xc)
        for i in range(len(segs)):
            for j in range(i + 1, len(segs)):
                ci, bi = segs[i]; cj, bj = segs[j]
                if ci != cj:
                    xc = (bj - bi) / (ci - cj)
                    if a < xc < b: pts.add(xc)
        pts = sorted(pts)
        for u, v in zip(pts[:-1], pts[1:]):
            m2 = (u + v) / 2
            val = max(c * m2 + b0 for c, b0 in segs)
            if val > theta:
                tot += v - u
    return tot

def compositions(D, parts):
    """All compositions of D into `parts` positive integers."""
    if parts == 1:
        yield (D,); return
    for first in range(1, D - parts + 2):
        for rest in compositions(D - first, parts - 1):
            yield (first,) + rest

def fam(word, v1=1):
    v = [v1]
    for s in word: v.append(v[-1] + s)
    return v

# ---------------------------------------------------------------- PROBE A
print("=" * 96)
print("PROBE A -- exchange-monotonicity (Schur shape) of mu_1/7, EXACT, exhaustive small D")
print("=" * 96)
for D in (13, 14, 15):
    words = list(compositions(D, 12))
    vals = {}
    for w in words:
        vals[w] = mu_exact(fam(w))
    viol = []
    checked = 0
    for w in words:
        for i in range(12):
            for j in range(12):
                if i == j or w[i] < w[j] + 2: continue
                w2 = list(w); w2[i] -= 1; w2[j] += 1; w2 = tuple(w2)
                checked += 1
                if vals[w2] > vals[w]:
                    viol.append((w, w2, (i, j), vals[w], vals[w2]))
    print(f"  D={D}: {len(words)} words, {checked} transfers checked, VIOLATIONS: {len(viol)}"
          f"  ({100*len(viol)/max(checked,1):.1f}%)")
    if viol:
        # characterize: positions and sizes
        from collections import Counter
        cpos = Counter((min(v[2]), max(v[2])) for v in viol)
        worst = max(viol, key=lambda v: v[4] - v[3])
        print(f"    top position pairs: {cpos.most_common(4)}")
        print(f"    worst: {worst[0]} -> {worst[1]} at (i,j)={worst[2]}: mu {float(worst[3]):.4f} -> {float(worst[4]):.4f} (+{float(worst[4]-worst[3]):.4f})")
    ap = tuple([1]*11 + [D-11]) if D > 12 else tuple([1]*12)
    uni = min(vals, key=lambda w: vals[w])
    print(f"    global min word: {uni} (mu={float(vals[uni]):.4f}); AP-like (1^12 word exists only at D=12)")

# ---------------------------------------------------------------- PROBE B
print()
print("=" * 96)
print("PROBE B -- deficit vs shortest zero-sum relation (the Poisson/lattice-minimum law)")
print("=" * 96)
def Emg_numeric(E, res=30000):
    t = 0.0
    for r in range(res):
        x = (r + .5) / res
        ph = sorted((e * x) % 1.0 for e in E)
        g = max(ph[i+1] - ph[i] for i in range(len(ph)-1))
        t += max(g, ph[0] + 1.0 - ph[-1])
    return t / res

def shortest_relation(E, max_support=5, max_coeff=3):
    """min L1-norm of a zero-sum integer relation sum m_i e_i = 0, sum m_i = 0,
    support >= 3 (pairs impossible); brute force small support/coefficients."""
    E = sorted(E); k = len(E)
    best = None
    idx = range(k)
    from itertools import combinations as C_, product
    for sup in (3, 4, 5):
        if sup > max_support: break
        for S in C_(idx, sup):
            vals = [E[i] for i in S]
            for coeffs in product(range(-max_coeff, max_coeff + 1), repeat=sup):
                if 0 in coeffs: continue
                if sum(coeffs) != 0: continue
                if sum(c * v for c, v in zip(coeffs, vals)) != 0: continue
                l1 = sum(abs(c) for c in coeffs)
                if best is None or l1 < best: best = l1
        if best is not None and sup >= 4: break
    return best  # None = no relation within the search box

H13_13 = sum(1/j for j in range(1, 14)) / 13
zoo = {
    "AP": list(range(1, 14)),
    "record": [2,4,6,8,10,11,12,13,14,16,18,20,22],
    "primes": [2,3,5,7,11,13,17,19,23,29,31,37,41],
    "Sidon-greedy": [1,2,4,8,13,21,31,45,66,81,97,123,148],
    "geometric": [1,2,4,8,16,32,64,128,256,512,1024,2048,4096],
    "GW": list(range(1,12)) + [13,24],
    "3adic cascade": [30,36,45,50,54,60,63,70,72,81,90,108,120],
}
rng = random.Random(63)
for _ in range(8):
    E = sorted(rng.sample(range(1, 400), 13))
    zoo[f"rand{_}"] = E
print(f"  {'family':16s} {'E[mg]':>8} {'deficit':>9} {'shortest rel L1 (sup<=5,|c|<=3)':>32}")
pts = []
for nm, E in zoo.items():
    v = Emg_numeric(E, 20000)
    d = H13_13 - v
    r = shortest_relation(E)
    pts.append((r if r else 99, d, nm))
    print(f"  {nm:16s} {v:8.4f} {d:9.4f} {str(r) if r else '>box (none found)':>32}")
pts.sort()
print("\n  envelope check (sorted by relation size): deficit should DECAY as min-relation grows:")
for r, d, nm in pts:
    print(f"    minrel={r if r<99 else 'none':>4}: deficit={d:+.4f}  ({nm})")

# ---------------------------------------------------------------- PROBE C
print()
print("=" * 96)
print("PROBE C -- can uniform PAIR data alone go below 1/7?  (synthetic process construction)")
print("=" * 96)
import numpy as np
M = 273           # circle resolution (multiple of 13 and 7)
target = np.full(M, 156.0 / M)   # aggregate pair-distance density: 13*12=156 ordered pairs uniform
# component library: rotation-averaged processes; each = (pair-density vector, E[maxgap])
comps = []
labels = []
rng2 = random.Random(99)
def add_component(gaps, label):
    gaps = np.array(gaps, dtype=float)
    gaps = gaps / gaps.sum()
    pos = np.cumsum(gaps) - gaps[0]
    # rotation-averaged pair density: depends only on gap multiset arrangement
    dens = np.zeros(M)
    for i in range(13):
        for j in range(13):
            if i == j: continue
            dcirc = (pos[j] - pos[i]) % 1.0
            dens[int(dcirc * M) % M] += 1.0 / 1.0
    # rotation-average = same histogram (rotations shift all distances equally? no --
    # pair distances are rotation-INVARIANT, so no averaging needed)
    mg = 0.0
    srt = np.sort(pos)
    mg = max(np.max(np.diff(srt)), srt[0] + 1 - srt[-1])
    comps.append((dens, mg))
    labels.append(label)

# jittered lattices (many jitter draws), one-big-gap, two-cluster, iid draws
for w in (0.0, 0.1, 0.25, 0.5, 0.75, 1.0):
    for rep in range(30):
        base = np.ones(13)
        jit = np.array([rng2.random() for _ in range(13)])
        g = (1 - w) * base / 13 + w * jit / jit.sum()
        add_component(g, f"jitter{w}")
for rep in range(60):
    g = np.array([rng2.random() for _ in range(13)]); add_component(g, "iid")
for big in (0.16, 0.2, 0.25, 0.3):
    for rep in range(15):
        rest = np.array([rng2.random() for _ in range(12)])
        g = np.concatenate([[big], (1 - big) * rest / rest.sum()])
        rng2.shuffle(g)
        add_component(g, f"biggap{big}")
# anti-lattice: gaps alternating small/large in structured patterns
for pat in ((2,1),(3,1),(3,2),(5,1),(4,3)):
    a, b = pat
    for rep in range(10):
        g = np.array(([a, b] * 7)[:13], dtype=float) + 0.05*np.array([rng2.random() for _ in range(13)])
        add_component(g, f"alt{pat}")

A = np.stack([c[0] for c in comps]).T          # M x K pair densities
m = np.array([c[1] for c in comps])            # K maxgaps
K = A.shape[1]
# minimize m.w  s.t.  A w = target, w >= 0, sum w = 1  -- projected gradient on the simplex
# with heavy penalty for the density constraint
wgt = np.full(K, 1.0 / K)
lam = 400.0
for it in range(4000):
    resid = A @ wgt - target
    grad = m + lam * (A.T @ resid) * 2
    wgt -= 0.002 * grad
    wgt = np.maximum(wgt, 0)
    s = wgt.sum()
    if s > 0: wgt /= s
resid = A @ wgt - target
Emax = float(m @ wgt)
print(f"  components: {K}; achieved E[maxgap] = {Emax:.4f}  (1/7 = {1/7:.4f})")
print(f"  pair-density constraint residual: L2 = {float(np.linalg.norm(resid)):.4f}, Linf = {float(np.abs(resid).max()):.4f} (target per-bin {156/M:.4f})")
top = np.argsort(-wgt)[:6]
print("  top components:", [(labels[t], round(float(wgt[t]), 3), round(float(m[t]), 3)) for t in top])
if Emax < 1/7 and np.abs(resid).max() < 0.1 * 156 / M:
    print("  ==> CONSTRUCTION SUCCEEDED: uniform-pair data does NOT force E[maxgap] > 1/7;")
    print("      any proof of the mean floor MUST use weight->=3 information.  (Verify harder.)")
else:
    print("  ==> construction stalls above 1/7 (or constraint unmet) -- evidence toward a")
    print("      2-point (Cohn-Elkies-style) LP floor; needs a true LP to decide.")
print()
print("DONE.")
