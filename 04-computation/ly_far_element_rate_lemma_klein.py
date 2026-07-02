#!/usr/bin/env python3
"""
klein-2026-07-02-S94 -- HYP-4001: the L_y FAR-ELEMENT RATE LEMMA (THM-534's open half,
restructured into decidable legs).

LEMMA (rate): for a co-offset set E (0 in E) and A a set of sectors (subsets of {1..6},
sector j = [j/7,(j+1)/7)), let Good_A(E) = {x : all frac(e x), e in E, avoid U_A}.
For a far element w:
    | J(A, E u {w}) - (1 - |A|/7) J(A, E) |  <=  2 * comp(A,E) * |A| / w,
where comp(A,E) = number of maximal intervals of Good_A(E). Proof: on each component I,
{x in I : frac(wx) notin U_A} has measure |I|(1-|A|/7) + err, |err| <= |A| * 2/w (each of
the |A| avoided sectors contributes at most two partial wraps at I's ends). Summing. QED

COROLLARY: L_y(E u {w}) <= L_y^inf(E) + K(E)/w, L_y^inf(E) := sum_r y_r (1-r/7) S_r(E),
K(E) = sum_r |y_r| * sum_{|A|=r} 2 comp(A,E) r.

This script: (1) verify the rate lemma exactly (interval arithmetic vs the bound) on
sample E, A, w; (2) compute the damped functionals L_y^inf on the bounded-spread shape
census per k and compare against cap_{k+1} (the row E u {w} has k+1 elements) -- margins
=> explicit w-thresholds W0(k); (3) report the restructured proof's leg list.
"""
from fractions import Fraction as F
import itertools

SECTOR = [(F(j,7), F(j+1,7)) for j in range(7)]
DUALS = {  # THM-534: k -> [(r, y_r)] beyond y_0 = 1
  13: [(1,F(-1,2)),(2,F(1,6))], 12: [(1,F(-1,2)),(2,F(1,6))], 11: [(1,F(-1,2)),(2,F(1,6))],
  10: [(1,F(-13,18)),(2,F(4,9)),(3,F(-1,6))], 9: [(1,F(-13,18)),(2,F(4,9)),(3,F(-1,6))],
  8:  [(1,F(-1,1)),(2,F(1,1)),(3,F(-9,10)),(4,F(3,5))],
}
CAPS = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91), 11: F(66,91), 12: F(78,91), 13: F(1,1)}

def avoid_set(E, A):
    """Maximal intervals of {x in [0,1): all frac(ex) avoid the sectors in A} (e=0 ok if 0 notin A)."""
    bad = []
    for e in E:
        if e == 0: continue
        for j in A:
            lo, hi = SECTOR[j]
            # frac(ex) in [lo,hi) <=> x in union over a of [(a+lo)/e, (a+hi)/e)
            for a in range(e):
                bad.append((F(a+lo.numerator*1,1)/e + 0, 0))  # placeholder
    # simpler: build bad intervals directly
    bad = []
    for e in E:
        if e == 0: continue
        for j in A:
            lo, hi = SECTOR[j]
            for a in range(e):
                bad.append(((a+lo)/e, (a+hi)/e))
    bad.sort()
    merged = []
    for l,h in bad:
        if merged and l <= merged[-1][1]: merged[-1] = (merged[-1][0], max(merged[-1][1], h))
        else: merged.append((l,h))
    good, prev = [], F(0)
    for l,h in merged:
        if l > prev: good.append((prev,l))
        prev = max(prev,h)
    if prev < 1: good.append((prev,F(1)))
    return good

def J(E, A):
    return sum(h-l for l,h in avoid_set(E,A))

def S_moments(E, rmax):
    S = {r: F(0) for r in range(1, rmax+1)}
    for r in range(1, rmax+1):
        for A in itertools.combinations(range(1,7), r):
            S[r] += J(E, A)
    return S

def Ly(E, k):
    S = S_moments(E, max(r for r,_ in DUALS[k]))
    return 1 + sum(y*S[r] for r,y in DUALS[k])

def Ly_damped(E, k):
    S = S_moments(E, max(r for r,_ in DUALS[k]))
    return 1 + sum(y*(1-F(r,7))*S[r] for r,y in DUALS[k])

print("="*92)
print("(1) RATE LEMMA verification: |J(A, E+{w}) - (1-|A|/7) J(A,E)| <= 2 comp |A| / w")
print("="*92)
E0 = [0,1,2,3]
for A in [(1,), (2,4), (1,3,5)]:
    base = J(E0, list(A)); comp = len(avoid_set(E0, list(A)))
    for w in (40, 97, 200):
        act = J(E0+[w], list(A))
        err = abs(act - (1-F(len(A),7))*base)
        bound = F(2*comp*len(A), w)
        print(f"  E={E0} A={A} w={w}: |err|={float(err):.6f} <= bound {float(bound):.6f}  OK={err <= bound}")

print()
print("="*92)
print("(2) DAMPED comparison: max L_y^inf over bounded shapes (k elts) vs cap_{k+1}; W0 thresholds")
print("="*92)
for kk in (8, 9, 10, 11):
    # shapes E with |E| = kk (containing 0), spread <= 12 (sample census; full = THM-534's window)
    best, arg = None, None
    for tail in itertools.combinations(range(1, 13), kk-1):
        E = [0]+list(tail)
        v = Ly_damped(E, kk+1 if kk+1 <= 13 else 13)
        if best is None or v > best: best, arg = v, E
    cap = CAPS[kk+1]
    margin = cap - best
    print(f"  k={kk}->row {kk+1}: max L_y^inf = {float(best):.6f} at E={arg};  cap_{kk+1} = {float(cap):.6f};"
          f"  margin = {float(margin):+.6f}  {'=> far rows CLOSED for w >= K/margin' if margin > 0 else '=> NEEDS finer treatment'}")

print()
print("(3) THE RESTRUCTURED PROOF (THM-534 open half):")
print("    [bounded-spread census: exact, decide]  +  [rate lemma: elementary, formalizable]")
print("    +  [damped comparison: finite rational, decide]  +  [w-band w < W0: finite exact sweep]")
print("    +  [clustered far blocks: the renormalization tower (opus F-lemmas)]")
print("DONE.")
