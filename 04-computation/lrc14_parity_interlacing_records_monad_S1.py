#!/usr/bin/env python3
"""
monad-explorer-2026-07-07-S1 (part 2) -- parity-interlacing adversaries vs T*, and the
weakened per-k tail targets.

Q1: does 2-level (bisection) / 3-level (trisection) / 6-level structure drive the
    crux-class E[maxgap] below T* = 56291/294294 ~ 0.191275?  (If yes, the mean sidecar
    dies even at k=13; the tail route is unaffected.)
Q2: does the SAME adversary class dent mu_{1/7} at k=11,12,13 below the WEAKENED
    per-k requirements (union-bound + m_P)?  The DAG needs only:
      k=8: mu>=0.67502   k=9:  mu>=0.56223   k=10: mu>=0.45209
      k=11: mu>=0.33121  k=12: mu>=0.19934   k=13: mu>=0.05649
    (observed minima 0.940/0.840/0.776/0.626/0.570/0.4425 -- slack 1.39x..7.8x).
    Exact-AP-minimality is NOT needed; these six floors are.
Q3: M values of the record families (are mean-minimizers near the loneliness moat?).

Tournament Analysis declaration:
  vertices: adversary family classes (2-level, 3-level, 6-level, free);
  pairwise observable: which class produces the lower exact E[maxgap] under the
            crux constraints; switch/gauge: descend within the winning class;
  tie Hamiltonian path: bisection -> trisection -> mixed -> free descent.
"""
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path
import random

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_spec = spec_from_file_location(
    "ds", ROOT / "04-computation" / "lrc_maxgap_ap_minimality_check_deathstar_S1.py")
_ds = module_from_spec(_spec); _spec.loader.exec_module(_ds)
Emaxgap_exact = _ds.Emaxgap_exact

TSTAR = F(1, 7) + F(6, 7) * F(14249, 252252)

def stats(E, res=30000):
    xs = (np.arange(res) + 0.5) / res
    ph = np.mod(np.outer(xs, np.array(E, dtype=np.float64)), 1.0)
    ph.sort(axis=1)
    gaps = np.empty_like(ph)
    gaps[:, :-1] = np.diff(ph, axis=1)
    gaps[:, -1] = ph[:, 0] + 1.0 - ph[:, -1]
    mg = gaps.max(axis=1)
    return mg.mean(), (mg > 1.0 / 7).mean()

def saturated(v): return all(any(x % q == 0 for x in v) for q in range(2, 15))
def primitive(v):
    g = 0
    for x in v: g = gcd(g, x)
    return g == 1
def single_scale(v): return max(v) <= 13 * min(v)
def in_crux(v): return saturated(v) and primitive(v) and single_scale(v)

def M_exact_farey(v, qmax=1200):
    """M(v) = max over Farey t=a/q of min_i ||v_i t||, exact lower profile.
    Adequate for census-scale sanity (true M attained at bounded q for these heights)."""
    best = F(0)
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            t = F(a, q)
            m = min(abs(vi * t - round(vi * t)) for vi in v)
            if m > best: best = m
    return best

print("=" * 78)
print("Q1 -- structured descent on crux-class E[maxgap] vs T* = %s ~ %.6f"
      % (TSTAR, float(TSTAR)))
print("=" * 78)

random.seed(4787)
cands = {}

# 2-level: evens 2*{1..m} + odd bisectors
for nodd in (1, 2, 3, 4):
    m = 13 - nodd
    evens = [2 * i for i in range(1, m + 1)]
    for odds in combinations(range(1, 2 * m + 2, 2), nodd):
        v = sorted(set(evens) | set(odds))
        if len(v) == 13 and in_crux(v):
            cands[tuple(v)] = f"2-level m={m} odds={odds}"

# 3-level: 3*{1..m} + fillers != 0 mod 3
for nf in (2, 3, 4, 5, 6):
    m = 13 - nf
    base = [3 * i for i in range(1, m + 1)]
    pool = [x for x in range(1, 3 * m + 3) if x % 3 != 0]
    for _ in range(400):
        fill = random.sample(pool, nf)
        v = sorted(set(base) | set(fill))
        if len(v) == 13 and in_crux(v):
            cands[tuple(v)] = f"3-level m={m} fill={sorted(fill)}"

# 6-level combined: 6*{1..m} + 3*odd + 2*odd + units
for _ in range(600):
    m6 = random.choice([3, 4, 5])
    base = [6 * i for i in range(1, m6 + 1)]
    pool = [x for x in range(2, 6 * m6 + 4) if x not in base]
    fill = random.sample(pool, 13 - m6)
    v = sorted(set(base) | set(fill))
    if len(v) == 13 and in_crux(v):
        cands[tuple(v)] = f"6-level m={m6}"

print(f"  candidate pool: {len(cands)} crux-class families; numeric prefilter...")
scored = sorted(((stats(v, 12000)[0], v, tag) for v, tag in cands.items()))[:12]
print("  top-12 by numeric E[mg]; exact-verifying...")
best_ex, best_v = F(1), None
for mgnum, v, tag in scored:
    ex = Emaxgap_exact(list(v))
    mark = " *** BELOW T* ***" if ex < TSTAR else ""
    print(f"    {str(list(v)):>62s} {tag:24s} exact={ex} ~{float(ex):.6f}{mark}")
    if ex < best_ex:
        best_ex, best_v = ex, list(v)

# free descent seeded from the structured winner
v = best_v[:]
cur = stats(v, 16000)[0]
for step in range(1200):
    i = random.randrange(13)
    cand = random.randrange(max(1, min(v) // 2), min(3 * max(v), 400))
    w = sorted(set(v[:i] + v[i + 1:] + [cand]))
    if len(w) != 13 or not in_crux(w):
        continue
    c = stats(w, 8000)[0]
    if c < cur - 1e-5:
        v, cur = w, c
ex = Emaxgap_exact(v)
print(f"\n  free-descent from winner: {v}")
print(f"    exact E[maxgap] = {ex} ~ {float(ex):.6f}   margin over T*: {float(ex - TSTAR):+.6f}")
if ex < best_ex:
    best_ex, best_v = ex, v

print(f"\n  SESSION RECORD (crux class): {best_v}")
print(f"    E[maxgap] = {best_ex} ~ {float(best_ex):.6f}  vs T* {float(TSTAR):.6f} "
      f"margin {float(best_ex - TSTAR):+.6f}  vs 1/7 margin {float(best_ex - F(1,7)):+.6f}")

print()
print("=" * 78)
print("Q2 -- parity adversary vs the WEAKENED per-k tail floors (mu_{1/7})")
print("=" * 78)
need = {8: 0.67502, 9: 0.56223, 10: 0.45209, 11: 0.33121, 12: 0.19934, 13: 0.05649}
obs_ap = {8: 691/735, 9: 247/294, 10: 38/49, 11: 1381/2205, 12: 13823/24255, 13: 477/1078}
for k in (11, 12, 13):
    best_mu, best_E = 2.0, None
    # parity shapes at size k: 2*{1..m}+odds and 3-level, plus free descent
    pool = []
    for nodd in range(1, min(5, k - 3)):
        m = k - nodd
        evens = [2 * i for i in range(1, m + 1)]
        for odds in combinations(range(1, 2 * m + 2, 2), nodd):
            v = sorted(set(evens) | set(odds))
            if len(v) == k:
                pool.append(v)
    random.shuffle(pool)
    pool = pool[:250]
    for v in pool:
        mu = stats(v, 8000)[1]
        if mu < best_mu:
            best_mu, best_E = mu, v
    # free descent on mu from the parity winner
    v = best_E[:]
    for step in range(500):
        i = random.randrange(k)
        cand = random.randrange(1, 80)
        w = sorted(set(v[:i] + v[i + 1:] + [cand]))
        if len(w) != k:
            continue
        c = stats(w, 8000)[1]
        if c < best_mu - 1e-4:
            v, best_mu = w, c
    best_mu = stats(v, 40000)[1]
    ok = "OK (floor holds)" if best_mu > need[k] else "*** BELOW WEAKENED FLOOR ***"
    print(f"  k={k}: adversarial min mu ~ {best_mu:.4f} at {v}")
    print(f"        needed {need[k]:.4f} | AP value {obs_ap[k]:.4f} | {ok}")

print()
print("=" * 78)
print("Q3 -- M values (loneliness) of the record mean-minimizers")
print("=" * 78)
for v in [[2, 4, 6, 8, 10, 11, 12, 13, 14, 16, 18, 20, 22],
          [2, 4, 6, 8, 10, 12, 13, 14, 16, 18, 20, 22, 24],
          best_v]:
    Mv = M_exact_farey(v, qmax=400)
    print(f"  M({v}) >= {Mv} ~ {float(Mv):.5f}  (1/14={1/14:.5f}, 1/13={1/13:.5f}, "
          f"in crux={in_crux(v)})")
print("\nDONE.")
