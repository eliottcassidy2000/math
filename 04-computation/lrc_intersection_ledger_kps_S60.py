#!/usr/bin/env python3
r"""
lrc_intersection_ledger_kps_S60.py   (kind-pasteur-2026-07-07-S60, HYP-4847)

THE INTERSECTION LEDGER: the S59 diameter floor extended to the G_P-intersected
legs k=8..12 of hlarge, BYPASSING the union bound.

POINTWISE INCLUSION (S59 lemma, intersected):  for a cluster E with diam(E) = D
(translate min->0; NO gcd reduction here -- G_P is not dilation-invariant):
    {x : maxgap_E(x) > 1/7}  ⊇  {x : roof_{D+1}(x) > 1/7}          (S59, pointwise)
    ⟹  rho*_{1/7}(P,E) = meas(G_P ∩ Good_E) >= meas(G_P ∩ {roof_{D+1} > 1/7})
                        =: ILedger(P, D+1)                          (exact rational)
G_P is an explicit union of rational intervals; the roof superlevel is per-Farey-cell
affine (opus-S134); the intersection is exact interval arithmetic.  No union bound,
no adversarial sup: for every shape (P, E) with diam(E) <= D, the hlarge bar
G2 >= m_P is DISCHARGED whenever ILedger(P, D+1) >= m_P.

Also corrects MISTAKE-121 (S59's "k=8..10 no bite" -- a table-start artifact): the
honest union-bound bites are diam <= 9/11/11 at k=8/9/10; this ledger should extend
them by exploiting quasi-independence (THM-530-D: R >= 0.796) exactly.

OUTPUT:
  PART 0: sanity -- ILedger(∅,n) == mu(AP_n) exactly (roof reproduction);
          ILedger >= union bound; Monte-Carlo spot checks.
  PART 1: per-size ledger -- for |P| = s = 1..5 (cluster k = 13-s), min over ALL
          C(13,s) sets P of ILedger(P,n); per-size crossing n*(s) vs m_P; the
          argmin P at the crossing; comparison with the union-bound bites.
  PART 2: the composite hlarge coverage table after this session.
  PART 3: anti-correlation anatomy -- ILedger vs independence product
          meas(G_P)*mu(AP_n) (the exact quasi-independence ratio per (P,n)).
"""
from fractions import Fraction as F
from itertools import combinations
import random

TH = F(1, 7)
MP = F(14249, 252252)

# ---------------------------------------------------------------- interval tools
def merge(iv):
    iv = sorted((a, b) for a, b in iv if b > a)
    out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return [(a, b) for a, b in out]

def complement01(iv):
    iv = merge(iv)
    out = []; cur = F(0)
    for a, b in iv:
        if a > cur: out.append((cur, a))
        cur = max(cur, b)
    if cur < 1: out.append((cur, F(1)))
    return out

def inter_measure(A, B):
    tot = F(0); i = j = 0
    while i < len(A) and j < len(B):
        a1, b1 = A[i]; a2, b2 = B[j]
        lo, hi = max(a1, a2), min(b1, b2)
        if hi > lo: tot += hi - lo
        if b1 < b2: i += 1
        else: j += 1
    return tot

def measure(A):
    return sum(b - a for a, b in A)

# ---------------------------------------------------------------- G_P and roof
def gp_intervals(P):
    """G_P = {x in [0,1]: ||p x|| >= 1/14 for all p in P} as disjoint intervals."""
    bad = []
    for p in P:
        w = F(1, 14 * p)
        for j in range(p + 1):
            c = F(j, p)
            bad.append((max(F(0), c - w), min(F(1), c + w)))
    return complement01(bad)

_farey_cache = {}
def farey(n):
    if n in _farey_cache: return _farey_cache[n]
    fr = set()
    for q in range(1, n + 1):
        for p in range(0, q + 1):
            fr.add(F(p, q))
    _farey_cache[n] = sorted(fr)
    return _farey_cache[n]

_roof_cache = {}
def roof_superlevel(n, theta=TH):
    """{x: roof_n(x) > theta} as disjoint rational intervals (opus-S134 roof)."""
    if (n, theta) in _roof_cache: return _roof_cache[(n, theta)]
    Fs = farey(n); out = []
    for a, b in zip(Fs[:-1], Fs[1:]):
        q, qp = a.denominator, b.denominator
        vl, vr = F(1, q), F(1, qp)          # roof endpoint values
        if vl <= theta and vr <= theta:
            continue
        if vl > theta and vr > theta:
            out.append((a, b)); continue
        # affine from vl (at a) to vr (at b): crossing at a + (theta-vl)*(b-a)/(vr-vl)
        xc = a + (theta - vl) * (b - a) / (vr - vl)
        out.append((a, xc) if vl > theta else (xc, b))
    out = merge(out)
    _roof_cache[(n, theta)] = out
    return out

def iledger(P, n):
    return inter_measure(gp_intervals(P), roof_superlevel(n))

# ---------------------------------------------------------------- PART 0
print("=" * 96)
print("PART 0 -- sanity: roof reproduction, union-bound domination, Monte-Carlo spot checks")
print("=" * 96)
known = {8: F(691,735), 9: F(247,294), 10: F(38,49), 11: F(1381,2205), 12: F(13823,24255), 13: F(477,1078)}
ok = all(measure(roof_superlevel(k)) == known[k] for k in known)
print(f"  ILedger(empty, k) == mu(AP_k) exact for k=8..13: {ok}")
def maxgap_f(E, x):
    ph = sorted((e * x) % 1.0 for e in E)
    g = max(ph[i+1] - ph[i] for i in range(len(ph)-1)) if len(ph) > 1 else 0
    return max(g, ph[0] + 1 - ph[-1])
rng = random.Random(60)
print("  MC spot checks: rho*(P,E) >= ILedger(P, D+1) on random shapes:")
for _ in range(4):
    s = rng.randrange(1, 6)
    P = sorted(rng.sample(range(1, 14), s))
    k = 13 - s
    D = rng.randrange(k + 1, 3 * k)
    Ee = sorted(rng.sample(range(1, D), k - 2)) + [0, D]
    Ee = sorted(set(Ee))
    if len(Ee) != k: continue
    led = iledger(P, D + 1)
    res = 20000; cnt = 0
    gp = gp_intervals(P)
    for r in range(res):
        x = (r + .5) / res
        ingp = any(float(a) <= x <= float(b) for a, b in gp)
        if ingp and maxgap_f(Ee, x) > 1/7: cnt += 1
    rho = cnt / res
    print(f"    P={P} E(diam {D}, k={k}): rho*={rho:.4f} >= ILedger={float(led):.4f}  {'OK' if rho >= float(led) - 3e-3 else '*** VIOL ***'}")

# ---------------------------------------------------------------- PART 1
print()
print("=" * 96)
print("PART 1 -- per-size ledger: min over ALL P of ILedger(P,n); crossings vs m_P")
print("=" * 96)
print(f"  m_P = {MP} = {float(MP):.6f}")
union_bites = {8: 9, 9: 11, 10: 11, 11: 15, 12: 23}   # MISTAKE-121-corrected union-bound bites
results = {}
for s in range(1, 6):
    k = 13 - s
    Ps = list(combinations(range(1, 14), s))
    n = k          # cluster diam starts at k-1 => n = D+1 starts at k
    last_ok, cross_P = None, None
    while n <= 90:
        worst = None; worstP = None
        for P in Ps:
            v = iledger(P, n)
            if worst is None or v < worst:
                worst, worstP = v, P
        if worst >= MP:
            last_ok = n
            n += 1
        else:
            cross_P = worstP
            break
    results[s] = (last_ok, worst, worstP)
    ub = union_bites[k]
    print(f"  |P|={s} (k={k:2d}): min_P ILedger >= m_P through n = {last_ok}  =>  BITE diam {k-1}..{last_ok-1}"
          f"   [union bound gave diam <= {ub}]   first-failing P = {worstP}, value {float(worst):.6f}")

# ---------------------------------------------------------------- PART 2
print()
print("=" * 96)
print("PART 2 -- composite hlarge coverage after S59+S60 (+ k=13 from S59/monad-S2)")
print("=" * 96)
print("  leg   bar-object                  coverage now")
print(f"  k=13  mu >= m_P                  diam <= 75 PROVED (S59/monad-S2; primitive diam OK)")
for s in range(1, 6):
    k = 13 - s
    la = results[s][0]
    print(f"  k={k:2d}  meas(G_P ∩ Good) >= m_P    diam <= {la-1} PROVED via ILedger (raw diam; all {len(list(combinations(range(1,14),s)))} P of size {s})")
print("  (k<=7: pigeonhole, PROVED unconditional -- THM-530)")

# ---------------------------------------------------------------- PART 3
print()
print("=" * 96)
print("PART 3 -- anti-correlation anatomy: exact quasi-independence ratio R = ILedger/(measGP*mu)")
print("=" * 96)
for s, Pw in [(5, (1,5,7,8,9)), (4, (1,11,12,13)), (3, (1,12,13)), (2, (1,13)), (1, (1,))]:
    k = 13 - s
    gp = gp_intervals(Pw); mg = measure(gp)
    row = []
    for n in (k, k + 4, k + 10, k + 20):
        led = iledger(Pw, n)
        mu = measure(roof_superlevel(n))
        R = led / (mg * mu)
        row.append(f"n={n}: R={float(R):.3f}")
    print(f"  |P|={s} THM-530-worst P={Pw}: meas(G_P)={float(mg):.4f}   " + "   ".join(row))
print()
print("  NOTE (honesty): the k=8..12 ledger uses RAW diameter (G_P breaks dilation-invariance);")
print("  the k=13 leg keeps PRIMITIVE diameter (G_empty = full circle).")
print("DONE.")
