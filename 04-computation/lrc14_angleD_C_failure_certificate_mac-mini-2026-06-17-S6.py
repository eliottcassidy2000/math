#!/usr/bin/env python3
"""
lrc14_angleD_C_failure_certificate — mac-mini-2026-06-17-S6  (ANGLE D RESULT)

RESULT: criterion C is NOT universal. It is SUFFICIENT but NOT NECESSARY for M>=1/14.

CERTIFIED C-FAILURE (covering 13-set with NO v satisfying W(S\\{v}) > 1/(7v)):
    S* = {1,2,3,5,7,8,9,10,11,12,13,38,42}
Every runner's removal leaves a level-1/14 safe set whose widest arc W(S*\\{v}) is
SMALLER than the danger-tooth width 1/(7v); the tightest is v=42 (margin -1/1430-ish).
Yet M(S*) = 2/23 = 0.08696 >= 1/14 (witness tau = 4/23), so LRC(14) is NOT violated.

=> The generalized arc-width criterion C does NOT prove LRC(14) by itself; it must be
   STRENGTHENED (e.g. allow PAIRS of runners / wider safe-arc certificates, or a
   two-tooth / SDR argument), because there exist hard covering 13-sets that the single-
   runner arc-width test cannot certify even though they are loose.

This script is the standalone reproducible certificate: it recomputes the criterion rows,
M by the candidate-set method, and M by an INDEPENDENT dense rational grid, and confirms:
  (1) C fails (all single-runner margins negative),
  (2) M(S*) = 2/23 >= 1/14 by two independent methods (so NOT an LRC break).
"""
from fractions import Fraction as F
from math import gcd

C = F(1, 14)

def Wsafe(A):
    """widest level-1/14 safe arc of A, exact, via integers on D=14*lcm(A)."""
    A = sorted(set(A)); L = 1
    for u in A: L = L * u // gcd(L, u)
    D = 14 * L; iv = []
    for u in A:
        cu = D // u; hw = D // (14 * u)
        for k in range(u):
            c = k * cu; iv.append((c - hw, c + hw))
    if not iv: return F(1)
    norm = []
    for lo, hi in iv:
        ln = hi - lo; a = lo % D; b = a + ln
        if b <= D: norm.append((a, b))
        else: norm.append((a, D)); norm.append((0, b - D))
    norm.sort(); mg = []; cl, ch = norm[0]
    for lo, hi in norm[1:]:
        if lo <= ch:
            if hi > ch: ch = hi
        else: mg.append((cl, ch)); cl, ch = lo, hi
    mg.append((cl, ch)); best = 0; n = len(mg)
    for i in range(n):
        hi = mg[i][1]; lo = mg[(i + 1) % n][0] + (D if i == n - 1 else 0)
        g = lo - hi
        if g > best: best = g
    return F(best, D) if best > 0 else F(0)

def covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); Cc = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1,2): Cc.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1,2): Cc.add(F(k, d)); k += 1
    Cc.add(F(1,2)); return Cc
def M_candidates(S):
    best = F(0); bt = None
    for t in cand(S):
        val = g(S, t)
        if val > best: best = val; bt = t
    return best, bt
def M_grid(S, maxden=3000):
    best = F(0); bt = None; seen = set()
    for q in range(1, maxden + 1):
        for p in range(0, q + 1):
            t = F(p, q)
            if t in seen: continue
            seen.add(t)
            val = g(S, t)
            if val > best: best = val; bt = t
    return best, bt

S = (1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 38, 42)
print("="*78)
print("CERTIFIED C-FAILURE (criterion C sufficient-not-necessary)")
print("="*78)
print(f"S* = {S}")
print(f"covering 2..14: {covering(S)};  |S|={len(set(S))}")
print()
print("Criterion C rows  (v, W(S*\\v), 1/(7v), margin = W - 1/(7v)):")
anyok = False
tightest = (F(-99), None)
for v in sorted(set(S)):
    A = [u for u in S if u != v]; W = Wsafe(A); thr = F(1, 7*v); m = W - thr
    if m > 0: anyok = True
    if m > tightest[0]: tightest = (m, v)
    print(f"   v={v:3d}:  W={str(W):>12s} = {float(W):.7f}   1/(7v)={float(thr):.7f}   margin={float(m):+.7f}  {'OK' if m>0 else ''}")
print(f"\nC holds (some v with W>1/(7v)): {anyok}   "
      f"{'==> C FAILS (all margins <= 0)' if not anyok else ''}")
print(f"least-negative (tightest) margin: {tightest[0]} at v={tightest[1]}")
print()
Mc, tc = M_candidates(S)
Mg, tg = M_grid(S, 3000)
print(f"M(S*) candidate-set method : {Mc} = {float(Mc):.7f}  (tau={tc})")
print(f"M(S*) independent grid<=3000: {Mg} = {float(Mg):.7f}  (tau={tg})")
print(f"methods agree: {Mc == Mg}")
print(f"1/14 = {F(1,14)} = {float(F(1,14)):.7f}")
print(f"M(S*) >= 1/14: {Mc >= F(1,14)}")
print()
print("CONCLUSION:")
if not anyok and Mc >= F(1,14) and Mc == Mg:
    print("  C FAILS on S* yet M(S*) = 2/23 >= 1/14.  => criterion C is SUFFICIENT but")
    print("  NOT NECESSARY.  Proving C universal CANNOT prove LRC(14); C must be strengthened.")
    print("  LRC(14) itself is NOT threatened by S* (it is loose, M = 2/23 > 1/14).")
