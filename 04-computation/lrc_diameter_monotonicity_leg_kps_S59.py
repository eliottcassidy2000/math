#!/usr/bin/env python3
r"""
lrc_diameter_monotonicity_leg_kps_S59.py   (kind-pasteur-2026-07-07-S59, HYP-4797)

THE DIAMETER-MONOTONICITY LEG of the density-floor mean target inf_E E[maxgap] > 1/7.

BRICK 1 (pointwise, one line, rigorous):
  For integer sets E \subseteq F, the point sets nest: {frac(ex)} \subseteq {frac(fx)},
  so the widest F-gap (an arc with no F-points, hence no E-points) sits inside an
  E-gap:   maxgap(E, x) >= maxgap(F, x)   FOR EVERY x.

BRICK 2 (translation/dilation invariance, already canon):
  maxgap(E+c, x) = maxgap(E, x) (all points rotate together by cx);
  E[maxgap(gE)] = E[maxgap(E)] (x -> gx preserves Lebesgue on the circle).

CONSEQUENCE (the small-diameter leg): translate E to min(E)=0, then E ⊆ {0..D},
  D = diam(E) = max(E)-min(E), and {0..D} is a translate of AP_{D+1} = {1..D+1}:
      E[maxgap(E)]  >=  A(D+1),   A(n) := E[maxgap(AP_n)] = Farey-roof sum
                                          sum_{p/q < p'/q' consecutive in F_n} 1/(q q'^2)
  (opus-S134 Farey roof, verified against all canon constants).
  So EVERY 13-element E with A(diam(E)+1) > 1/7 satisfies the mean target,
  RIGOROUSLY, with no structure analysis.  The leg covers diam(E) <= D0 where
  D0+1 = the largest n with A(n) > 1/7.

  For the QUANTITATIVE bar T* = 1/7 + (6/7) m_P (monad-explorer HYP-4787 / MISTAKE-118),
  the leg covers diam(E) <= D0* where D0*+1 = largest n with A(n) > T*.

This script:
  PART 1: verify Brick 1 pointwise on random/structured (E,F,x) (sanity; proof is 1 line).
  PART 2: exact A(n) for n = 13..34 via the Farey-roof sum, independently cross-checked
          against the death-star exact engine at n = 13, 14, 21 and numerics; find the
          crossings vs 1/7 and vs T*.
  PART 3: place every known minimizer (AP, GW, opus-stretch, death-star prim-sat,
          monad record) against its diameter bound.
  PART 4: per-diameter structured descent (dilated-AP+splitters, GAP interlacings,
          parity families) at fixed primitive diameter D = 18..32: does the
          per-diameter empirical inf dip under 1/7 anywhere the leg doesn't cover?
  PART 5: the deficit frame ledger: deficit(E) = H_13/13 - E[maxgap(E)] vs short-relation
          counts (3APs / Schur triples / additive energy), for the zoo — the forgotten
          pair-uniformity factoid says PAIR statistics are universal, so the deficit is
          carried entirely by >=3-term relations; measure how much each contributes.
"""
from fractions import Fraction as F
import math, random, itertools

# ---------------------------------------------------------------- exact engines
def farey(n):
    fr = set()
    for q in range(1, n + 1):
        for p in range(0, q + 1):
            fr.add(F(p, q))
    return sorted(fr)

def A_farey_roof(n):
    """E[maxgap(AP_n)] = sum over consecutive Farey(n) pairs of 1/(q q'^2)."""
    Fs = farey(n)
    S = F(0)
    for a, b in zip(Fs[:-1], Fs[1:]):
        S += F(1, a.denominator * b.denominator ** 2)
    return S

def Emaxgap_exact(E, kdenom=None):
    """death-star S1 corrected exact engine: E[maxgap of {frac(j x): j in E}].
    Breakpoints at all Farey m/d, d <= kdenom >= max|e_i-e_j|, max(|e|)."""
    E = list(E)
    if kdenom is None:
        kdenom = max(max(abs(j) for j in E),
                     max(abs(a - b) for a in E for b in E))
    bps = set([F(0), F(1)])
    for d in range(1, kdenom + 1):
        for m in range(1, d):
            bps.add(F(m, d))
    bps = sorted(bps)
    total = F(0)
    for a, b in zip(bps, bps[1:]):
        mid = (a + b) / 2
        fl = {j: (j * mid).__floor__() for j in E}
        order = sorted(E, key=lambda j: (j * mid - fl[j]))
        gaps = []
        for s in range(len(order)):
            if s < len(order) - 1:
                j1, j2 = order[s], order[s + 1]
                c = F(j2 - j1); b0 = F(-(fl[j2] - fl[j1]))
            else:
                c = F(order[0] - order[-1]); b0 = F(-(fl[order[0]] - fl[order[-1]]) + 1)
            gaps.append((c, b0))
        subbp = set([a, b])
        for i in range(len(gaps)):
            for jx in range(i + 1, len(gaps)):
                ci, bi = gaps[i]; cj, bj = gaps[jx]
                if ci != cj:
                    xc = (bj - bi) / (ci - cj)
                    if a < xc < b:
                        subbp.add(xc)
        subbp = sorted(subbp)
        for u, v in zip(subbp, subbp[1:]):
            m2 = (u + v) / 2
            cb, bb = max(gaps, key=lambda cb: cb[0] * m2 + cb[1])
            total += ((cb * u + bb) + (cb * v + bb)) / 2 * (v - u)
    return total

def maxgap_float(E, x):
    ph = sorted((e * x) % 1.0 for e in E)
    g = max(ph[i + 1] - ph[i] for i in range(len(ph) - 1)) if len(ph) > 1 else 0.0
    return max(g, ph[0] + 1.0 - ph[-1])

def Emaxgap_numeric(E, res=60000):
    return sum(maxgap_float(E, (r + .5) / res) for r in range(res)) / res

# ---------------------------------------------------------------- PART 1
print("=" * 88)
print("PART 1 -- BRICK 1 pointwise monotonicity  maxgap(E,x) >= maxgap(F,x) for E subset F")
print("=" * 88)
rng = random.Random(59)
viol = 0; trials = 0
for _ in range(400):
    D = rng.randrange(14, 60)
    Fset = sorted(rng.sample(range(0, D + 1), rng.randrange(14, min(D + 1, 30))))
    ksub = rng.randrange(3, len(Fset))
    Eset = sorted(rng.sample(Fset, ksub))
    for _ in range(40):
        x = rng.random()
        trials += 1
        if maxgap_float(Eset, x) < maxgap_float(Fset, x) - 1e-12:
            viol += 1
print(f"  random nested pairs: {viol} violations / {trials} pointwise checks")
# structured: every 13-subset bound vs its {0..D} superset, numeric mean
E_rec = [2, 4, 6, 8, 10, 11, 12, 13, 14, 16, 18, 20, 22]          # monad record
E0 = [e - 2 for e in E_rec]                                        # translate to 0
D = max(E0)
lhs = Emaxgap_numeric(E_rec, 30000); rhs = Emaxgap_numeric(list(range(0, D + 1)), 30000)
print(f"  record family: numeric E[mg(E)]={lhs:.6f} >= numeric E[mg({{0..{D}}})]={rhs:.6f}: {lhs>=rhs-1e-4}")

# ---------------------------------------------------------------- PART 2
print()
print("=" * 88)
print("PART 2 -- exact A(n)=E[maxgap(AP_n)] (Farey roof), crossings vs 1/7 and vs T*")
print("=" * 88)
mP = F(14249, 252252)
Tstar = F(1, 7) + F(6, 7) * mP
print(f"  1/7 = {float(F(1,7)):.6f}    T* = 1/7+(6/7)m_P = {Tstar} = {float(Tstar):.6f}")
# cross-check the roof against the independent death-star engine
for n in (13, 14, 21):
    ap = list(range(1, n + 1))
    e1 = A_farey_roof(n); e2 = Emaxgap_exact(ap)
    print(f"  cross-check n={n}: roof={e1}  engine={e2}  match={e1 == e2}")
Avals = {}
cross17 = crossT = None
for n in range(13, 35):
    Avals[n] = A_farey_roof(n)
    a = Avals[n]
    m17 = "> 1/7" if a > F(1, 7) else "<= 1/7  ***"
    mT = "> T*" if a > Tstar else "<= T*"
    if cross17 is None and a <= F(1, 7): cross17 = n
    if crossT is None and a <= Tstar: crossT = n
    print(f"  n={n:2d} (diam D={n-1:2d}):  A(n) = {str(a):>24} = {float(a):.6f}   {m17:>12}  {mT}")
print(f"\n  ==> A(n) stays > 1/7 through n = {cross17-1} (diam <= {cross17-2}); first n <= 1/7 at n = {cross17}")
print(f"  ==> A(n) stays > T*  through n = {crossT-1} (diam <= {crossT-2}); first n <= T*  at n = {crossT}")
print(f"  LEG (positivity bar): every 13-set E with diam(E) <= {cross17-2} has E[maxgap] > 1/7.")
print(f"  LEG (m_P bar T*):     every 13-set E with diam(E) <= {crossT-2} has E[maxgap] > T*.")

# ---------------------------------------------------------------- PART 3
print()
print("=" * 88)
print("PART 3 -- known minimizers vs their diameter bound")
print("=" * 88)
zoo = {
    "AP {1..13}":               list(range(1, 14)),
    "GW {1..11,13,24}":         list(range(1, 12)) + [13, 24],
    "opus stretch {0,2..12,17,28}": [0] + list(range(2, 13)) + [17, 28],
    "death-star 2*{1..12}u{13}": [2,4,6,8,10,12,13,14,16,18,20,22,24],
    "monad record 2*{1..11}u{11,13}": E_rec,
    "kps-S57 adversarial":      [2,6,8,10,11,12,14,16,18,20,22,26,42],
}
for nm, Eset in zoo.items():
    d = max(Eset) - min(Eset)
    g = math.gcd(*[abs(a - b) for a in Eset for b in Eset if a != b])
    dprim = d // g
    bound = Avals.get(dprim + 1)
    ex = Emaxgap_exact(Eset)
    btxt = f"A({dprim+1})={float(bound):.5f}" if bound else f"A({dprim+1}) beyond table"
    cov = "covered by leg(1/7)" if bound and bound > F(1, 7) else "NOT covered"
    print(f"  {nm:34s} diam={d:3d} gcd={g} primdiam={dprim:3d}  E[mg]={float(ex):.6f}  {btxt}  {cov}")

# ---------------------------------------------------------------- PART 4
print()
print("=" * 88)
print("PART 4 -- per-primitive-diameter structured descent (does inf dip near/below 1/7?)")
print("=" * 88)
def neighbors(Eset, D):
    """structured moves at fixed diameter: keep 0 and D pinned, move interior points."""
    E = sorted(Eset); out = []
    interior = E[1:-1]
    allowed = [v for v in range(1, D) ]
    for i, e in enumerate(interior):
        for step in (-3,-2,-1,1,2,3):
            v = e + step
            if 0 < v < D and v not in E:
                cand = sorted([0, D] + interior[:i] + [v] + interior[i+1:])
                out.append(cand)
    return out

def descend(seed, D, iters=250, rs=1):
    rng2 = random.Random(rs)
    cur = sorted(seed); best = cur; bval = Emaxgap_numeric(cur, 4000)
    for _ in range(iters):
        cands = neighbors(cur, D)
        if not cands: break
        vals = [(Emaxgap_numeric(c, 2500), c) for c in rng2.sample(cands, min(18, len(cands)))]
        v, c = min(vals)
        if v < bval - 1e-6:
            bval, best, cur = v, c, c
        else:
            cur = rng2.choice(cands)  # sideways kick
    return bval, best

print(f"  {'D':>3} {'best numeric E[mg]':>20}  {'exact (refined)':>16}  family")
per_diam = {}
for D in range(18, 33, 2):
    seeds = []
    # parity interlacing generalizations (monad mechanism): evens + odd splitters
    ev = [2 * i for i in range(1, D // 2 + 1)]
    if len(ev) >= 11:
        for spl in ([D//2-1, D//2+1], [D//2, D//2+2], [D//2-3, D//2-1]):
            spl = [s if s % 2 == 1 else s + 1 for s in spl]
            cand = sorted(set(ev[:11] + spl))
            cand = [c - min(cand) for c in cand]
            if len(cand) == 13 and max(cand) <= D:
                base = cand + [D] if D not in cand else cand
                base = sorted(set(base))[:13]
                if len(base) == 13: seeds.append(sorted(base))
    # dilated AP + boundary defects
    for d in (2, 3):
        L = D // d
        if L >= 11:
            core = [d * i for i in range(0, 12)]
            cand = sorted(set(core + [D]))
            if len(cand) == 13 and max(cand) == D:
                seeds.append(cand)
    # AP with stretched tail (opus stretch shape scaled)
    seeds.append(sorted(set([0] + list(range(2, 13)) + [D])))
    best = (9.9, None)
    for s in seeds:
        s = sorted(set(s))
        if len(s) != 13 or max(s) != D or min(s) != 0:
            continue
        v, fam = descend(s, D, rs=D)
        if v < best[0]: best = (v, fam)
    if best[1] is None:
        continue
    ex = Emaxgap_exact(best[1]) if D <= 30 else None
    per_diam[D] = (best[0], best[1], ex)
    extxt = f"{float(ex):.6f}" if ex is not None else "--"
    print(f"  {D:3d} {best[0]:>20.6f}  {extxt:>16}  {best[1]}")
print(f"  (1/7 = {1/7:.6f}; T* = {float(Tstar):.6f})")

# ---------------------------------------------------------------- PART 5
print()
print("=" * 88)
print("PART 5 -- deficit frame: deficit = H_13/13 - E[mg] vs short-relation counts")
print("=" * 88)
H13 = sum(F(1, j) for j in range(1, 14)) / 13
print(f"  iid mean H_13/13 = {H13} = {float(H13):.6f};  budget to 1/7 = {float(H13 - F(1,7)):.6f}")
def relation_counts(Eset):
    E = sorted(Eset); S = set(E)
    ap3 = sum(1 for a, b in itertools.combinations(E, 2) if (a + b) % 2 == 0 and (a + b) // 2 in S and (a+b)//2 not in (a,b))
    schur = sum(1 for a, b in itertools.combinations(E, 2) if a + b in S)
    # additive energy E_+ = #{(a,b,c,d): a+b=c+d} via difference multiset
    diffs = {}
    for a in E:
        for b in E:
            diffs[a - b] = diffs.get(a - b, 0) + 1
    energy = sum(v * v for v in diffs.values())
    return ap3, schur, energy
zoo2 = dict(zoo)
zoo2["primes 13"] = [2,3,5,7,11,13,17,19,23,29,31,37,41]
zoo2["Sidon-ish (greedy B2)"] = [0,1,3,7,12,20,30,44,65,80,96,122,147]
zoo2["geometric 2^i"] = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096]
print(f"  {'family':34s} {'E[mg]':>9} {'deficit':>9} {'#3AP':>5} {'#Schur':>7} {'energy':>7}")
for nm, Eset in zoo2.items():
    if max(Eset) - min(Eset) <= 60:
        v = float(Emaxgap_exact(Eset))
    else:
        v = Emaxgap_numeric(Eset, 40000)
    ap3, sch, en = relation_counts(Eset)
    print(f"  {nm:34s} {v:9.6f} {float(H13)-v:9.6f} {ap3:5d} {sch:7d} {en:7d}")
print()
print("  Pair-uniformity factoid: for ANY distinct integers, each ||(e_i-e_j)x|| is EXACTLY")
print("  uniform on [0,1/2] in x -- pair statistics are UNIVERSAL; the deficit above is")
print("  carried entirely by >=3-term relations (loop-vs-torus Fourier: E_x[mg] = H13/13 +")
print("  sum_{k in L(E), k!=0} mg-hat(k), and pure pair frequencies m(delta_i-delta_j)")
print("  never lie in L(E) since m(e_i-e_j) != 0).")
print("DONE.")
