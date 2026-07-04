#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
PIN THE GLOBAL R'-MINIMIZER over covering families by greedy descent + many restarts (opus-S59).

R'(S) = meas(lonely S)/(m_R m_Q), R={s:14∤s}, Q={s/14:14|s}. S59 pass-1 showed drop-7 (0.315) is NOT
global; a random family hit 0.285. Since R' -> 1 as the far part grows (decorrelation, Lemma C proved),
the min lives on SMALL magnitude. Greedy 1-swap descent from many seeds. Descent uses FLOAT arc measure
(fast); the final global minimizer is RE-VERIFIED in EXACT rational arithmetic.
"""
import sys, random
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

MAXV = 42; BANDf = 1.0/14.0
def primitive(S): return reduce(gcd, S) == 1   # gcd=1: else scale-invariance reduces it (kps HYP-4060)
# ---- float arc measure (fast, for descent decisions) ----
def safe_arcs_f(v):
    return [((k + BANDf)/v, (k + 1 - BANDf)/v) for k in range(v)]
def inter_f(A, B):
    r = []; i = j = 0
    while i < len(A) and j < len(B):
        lo = A[i][0] if A[i][0] > B[j][0] else B[j][0]
        hi = A[i][1] if A[i][1] < B[j][1] else B[j][1]
        if lo < hi: r.append((lo, hi))
        if A[i][1] < B[j][1]: i += 1
        else: j += 1
    return r
_mc = {}
def meas_f(S):
    key = tuple(sorted(set(S)))
    if key in _mc: return _mc[key]
    a = safe_arcs_f(key[0])
    for v in key[1:]:
        a = inter_f(a, safe_arcs_f(v))
        if not a: _mc[key] = 0.0; return 0.0
    m = sum(h - l for l, h in a); _mc[key] = m; return m
def covers(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def Rprime_f(S):
    S = sorted(set(S))
    if len(S) != 13 or not covers(S) or not primitive(S): return None
    R = [s for s in S if s % 14 != 0]; Q = [s // 14 for s in S if s % 14 == 0]
    if not Q: return None
    mR, mQ, mS = meas_f(R), meas_f(Q), meas_f(S)
    if mR == 0 or mQ == 0: return None
    return mS / (mR * mQ)
# ---- exact arc measure (verify the winner) ----
def safe_arcs_e(v):
    B = Fr(1, 14); return [((Fr(k) + B)/v, (Fr(k+1) - B)/v) for k in range(v)]
def inter_e(A, B):
    r = []; i = j = 0
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0]); hi = min(A[i][1], B[j][1])
        if lo < hi: r.append((lo, hi))
        if A[i][1] < B[j][1]: i += 1
        else: j += 1
    return r
def meas_e(S):
    S = sorted(set(S));
    if not S: return Fr(1)
    a = safe_arcs_e(S[0])
    for v in S[1:]:
        a = inter_e(a, safe_arcs_e(v))
        if not a: return Fr(0)
    return sum(h - l for l, h in a)
def Rprime_e(S):
    R = [s for s in S if s % 14 != 0]; Q = [s // 14 for s in S if s % 14 == 0]
    return meas_e(S) / (meas_e(R) * meas_e(Q))

POOL = list(range(1, MAXV + 1))
def descend(S):
    S = sorted(set(S)); cur = Rprime_f(S)
    if cur is None: return None, None
    improved = True
    while improved:
        improved = False
        for i in range(13):
            for w in POOL:
                if w in S: continue
                T = S[:i] + S[i+1:] + [w]
                v = Rprime_f(T)
                if v is not None and v < cur - 1e-12:
                    S = sorted(set(T)); cur = v; improved = True; break
            if improved: break
    return cur, S

seeds = [sorted([1,2,3,4,5,6,8,9,10,11,12,13,14]),
         sorted([1,3,4,5,8,9,10,11,12,13,14,17,23])]
random.seed(7); tries = 0
while len(seeds) < 250 and tries < 500000:
    tries += 1
    S = random.sample([x for x in POOL if x <= 32], 13)
    if any(s % 14 == 0 for s in S) and covers(S) and primitive(S): seeds.append(sorted(S))

best = (9.9, None); locmins = {}
for s in seeds:
    v, S = descend(s)
    if v is None: continue
    locmins[tuple(S)] = v
    if v < best[0]: best = (v, S)

Sb = best[1]
print("="*100); print(" GLOBAL R'-MIN (greedy 1-swap descent, speeds <= %d, %d seeds)" % (MAXV, len(seeds))); print("="*100)
print(f"  GLOBAL MIN (float): R'={best[0]:.5f}  S={Sb}")
print(f"  EXACT re-verify:    R'={float(Rprime_e(Sb)):.5f}")
R = [s for s in Sb if s%14]; Q=[s//14 for s in Sb if s%14==0]
print(f"    R={R}   dropped from 1..13: {[x for x in range(1,14) if x not in Sb]}   Q(far/14)={Q}")
print(f"    m_R={float(meas_e(R)):.4f}  m_Q={float(meas_e(Q)):.4f}  meas(S)={float(meas_e(Sb)):.5f}")
print(f"    drop-7 reference R'={float(Rprime_e([1,2,3,4,5,6,8,9,10,11,12,13,14])):.5f}")
print("\n  distinct local minima (top 15 lowest, EXACT R'):")
for S, v in sorted(locmins.items(), key=lambda kv: kv[1])[:15]:
    dr = [x for x in range(1,14) if x not in S]
    print(f"    R'={float(Rprime_e(list(S))):.5f}  S={list(S)}  drop{dr}")
ng = sum(1 for v in locmins.values() if abs(v-best[0])<1e-5)
print(f"\n  #distinct local minima: {len(locmins)};  seeds reaching global min: {ng}/{len(locmins)}")
print("DONE.")
