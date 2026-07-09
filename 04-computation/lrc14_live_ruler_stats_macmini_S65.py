#!/usr/bin/env python3
"""
lrc14_live_ruler_stats_macmini_S65.py -- the SCHUR-BUDGET / live-ruler statistics (THM-668 follow-up)

Classify every pair-sum ruler q = v_i+v_j (i<=j) of a 13-set:
  DEAD        : q | v_l for some runner l  (Schur-kill; only possible for q <= Vmax)
  BAND-EMPTY  : no multiplier p in [1,q-1] puts all residues v_l*p mod q in [q/14, 13q/14]
  LIVE        : some p works  (hosts a witness with M >= 1/14 at t = p/q)

FREE FACTS (proved in-session, verified here):
  (a) every ruler with q > Vmax is UNDEAD (q | v_l impossible: 0 < v_l <= Vmax < q).
  (b) at p = 1, a large ruler q is banded iff 14*Vmax/13 <= q <= 14*Vmin  -- i.e. p=1 on pair-sum
      rulers IS spread13_lonely (r <= 13).  The open regime r > 13 lives on multipliers p >= 2.

QUESTIONS: how many live rulers do covering sets have?  which q hosts the witness (q/Vmax)?
does the live count ever get dangerously low on covering sets (the AP has live supply squeezed
to the single tangent ruler q=14)?
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce
import random

random.seed(65)

def covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def primitive(S):
    return reduce(gcd, S) == 1

def band_ok(S, q, p):
    """all residues v*p mod q in the CLOSED band [q/14, 13q/14] i.e. 14*r >= q and 14*(q-r) >= q"""
    for v in S:
        r = v * p % q
        if 14 * r < q or 14 * (q - r) < q:
            return False
    return True

def classify_rulers(S):
    """return dict q -> ('DEAD'|'EMPTY'|'LIVE', best_p or None), one entry per distinct pair sum."""
    out = {}
    for i in range(13):
        for j in range(i, 13):
            q = S[i] + S[j]
            if q in out:
                continue
            if any(v % q == 0 for v in S):
                out[q] = ('DEAD', None)
                continue
            live_p = next((p for p in range(1, q) if band_ok(S, q, p)), None)
            out[q] = (('LIVE', live_p) if live_p is not None else ('EMPTY', None))
    return out

def exact_M_and_q(S):
    bestn, bestd, bestt = 0, 1, None
    qseen = set()
    for i in range(13):
        for j in range(i, 13):
            q = S[i] + S[j]
            if q in qseen: continue
            qseen.add(q)
            for p in range(1, q):
                m = min(min(v * p % q, q - v * p % q) for v in S)
                if m * bestd > bestn * q:
                    bestn, bestd, bestt = m, q, F(p, q)
    return F(bestn, bestd), bestt

# ---------------------------------------------------------------- census over [1,18] covering sets
print("=" * 100)
print("LIVE-RULER SUPPLY over the 966 primitive covering 13-subsets of [1,18]")
print("=" * 100)
live_counts, host_ratio, min_live = [], [], None
verify_a = True
for S in combinations(range(1, 19), 13):
    if not (covering(S) and primitive(S)):
        continue
    S = list(S)
    cls = classify_rulers(S)
    # verify free fact (a): no dead ruler above Vmax
    if any(st == 'DEAD' for q, (st, _) in cls.items() if q > S[-1]):
        verify_a = False
    nlive = sum(1 for st, _ in cls.values() if st == 'LIVE')
    live_counts.append(nlive)
    if min_live is None or nlive < min_live[0]:
        min_live = (nlive, S)
    M, t = exact_M_and_q(S)
    host_ratio.append(F(t.denominator) / S[-1] if t else None)
n = len(live_counts)
print(f"{n} sets.  free fact (a) [q > Vmax never dead] verified: {verify_a}")
print(f"live rulers per set: min = {min(live_counts)}, mean = {sum(live_counts)/n:.1f}, "
      f"max = {max(live_counts)}  (of ~78-91 distinct pair sums)")
print(f"minimum-supply set: {min_live[1]} with {min_live[0]} live rulers")
hr = sorted(host_ratio)
print(f"witness modulus / Vmax: min = {float(hr[0]):.3f}, median = {float(hr[n//2]):.3f}, "
      f"max = {float(hr[-1]):.3f}")
inwin = sum(1 for x in host_ratio if 1 < x <= 2)
print(f"witness modulus in (Vmax, 2Vmax] (the undead large rulers): {inwin}/{n} = {inwin/n:.1%}")

# ---------------------------------------------------------------- the AP and the 91-cluster
print()
print("=" * 100)
print("CONTRAST: the tight AP vs the covering 7-structured 91-cluster")
print("=" * 100)
for name, S in [("tight AP {1..13} (NON-covering)", list(range(1, 14))),
                ("worst7Struct @91 (covering)", sorted(91 - e for e in [0,7,14,21,26,29,37,44,51,58,67,75,82]))]:
    cls = classify_rulers(S)
    live = [(q, p) for q, (st, p) in sorted(cls.items()) if st == 'LIVE']
    dead = [q for q, (st, _) in sorted(cls.items()) if st == 'DEAD']
    empty = [q for q, (st, _) in sorted(cls.items()) if st == 'EMPTY']
    print(f"{name}: {len(live)} LIVE {live[:6]}{'...' if len(live)>6 else ''}")
    print(f"    {len(dead)} DEAD (Schur-killed) {dead[:14]}{'...' if len(dead)>14 else ''}")
    print(f"    {len(empty)} BAND-EMPTY {empty[:10]}{'...' if len(empty)>10 else ''}")

# ---------------------------------------------------------------- adversarial: minimize live count
print()
print("=" * 100)
print("ADVERSARIAL: hill-climb MINIMIZING the live-ruler count over covering sets (cap 60)")
print("=" * 100)
def nlive_of(S):
    return sum(1 for st, _ in classify_rulers(S).values() if st == 'LIVE')

best = None
for restart in range(6):
    while True:
        S = sorted(random.sample(range(1, 61), 13))
        if covering(S) and primitive(S):
            break
    cur = nlive_of(S)
    for step in range(150):
        T = list(S)
        T[random.randrange(13)] = random.randrange(1, 61)
        T = sorted(set(T))
        if len(T) != 13 or not covering(T) or not primitive(T):
            continue
        nl = nlive_of(T)
        if nl <= cur:
            S, cur = T, nl
    if best is None or cur < best[0]:
        best = (cur, S)
    print(f"  restart {restart}: min live = {cur}  S = {S}")
print(f"ADVERSARIAL MIN live-ruler count on covering sets (cap 60): {best[0]}  at {best[1]}")
M, t = exact_M_and_q(best[1])
print(f"  its exact M = {M} at t* = {t}  (>= 1/14: {M >= F(1,14)})")
print()
print("If the adversarial min stays >> 1, the Schur-budget route has room: covering sets cannot")
print("squeeze the live supply to the AP's single tangent ruler.")
print("Done.")
