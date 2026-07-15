#!/usr/bin/env python3
"""lrc14_lenses_moebius_rigidity_macmini_S111.py -- mac-mini-2026-07-15-S111.
The S110 lenses pointed at LRC(14). Two parts:

(I) THM-874 -- THE MOEBIUS-LOG^2 GRAMMAR of the Farey profile (roots-of-unity lens).
    The depth-layer weights of THM-826's profile ladder assemble into
        F(x) = sum_{d>=1} (mu(d)/d^2) * log^2(1/(1-x^d)),
    i.e. [x^s] F = (2/s) * H*(s),  H*(s) = sum_{i<s, gcd(i,s)=1} 1/i
    (the corridor constants of THM-819/853(II) for EVERY depth s at once;
    prime depth = pure log^2 term, composite = Moebius-corrected).
    Also: H*(m) = sum_{d|m} (mu(d)/d) H_{m/d - 1} (Moebius form), and the
    LRC(15) corridor instance  m({1..13}; lam) = (H*(14)/7)(1 - 14 lam) on
    [1/15, 1/14], margin values at rational points -- all exact in Q.
    THM-868 contrast: the tournament/figurate ladders are GEOMETRIC in one
    substituted variable; the LRC ladder's grammar is MOEBIUS-LOG^2 -- multi-
    scale by construction (the d-sum is the scale tower). That absence of a
    single substitution variable is the structural face of LRC's difficulty.

(II) THE LOW-M RIGIDITY ASSEMBLY (shelf lens on the peel walk; klein-S308 handoff).
    Claim (assembled from THM-726 + THM-724): every covering-type 13-set with
    M < 1/13 has at most ONE element > 14, hence >= 12 elements in {1..14},
    hence is near-AP (>= 10 in {1..14}).  Probe: generate covering-type sets
    (m(S; 1/13) = 0) by deep-well mutations + random small-element sets;
    verify near-AP-ness and the outlier count for every hit; verify
    multi-outlier samples have m(S; 1/13) > 0 (THM-726's face).
"""
import sys, random
from fractions import Fraction as F
from math import comb, gcd, log
sys.stdout.reconfigure(line_buffering=True)

# ---------- exact good measure (interval sweep, kps-style)
def good_measure(speeds, lam):
    pieces = []
    for w in speeds:
        r = F(lam) / w
        for k in range(w):
            c = F(k, w); lo, hi = c - r, c + r
            if lo < 0: pieces.append((F(0), hi)); pieces.append((lo + 1, F(1)))
            elif hi > 1: pieces.append((lo, F(1))); pieces.append((F(0), hi - 1))
            else: pieces.append((lo, hi))
    pieces.sort(); tot = F(0); cur = F(0)
    for lo, hi in pieces:
        if lo > cur: tot += lo - cur
        cur = max(cur, hi)
    if cur < 1: tot += 1 - cur
    return tot

def mobius(n):
    m, p, cnt = 1, 2, 0
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0: return 0
            m = -m
        p += 1
    if n > 1: m = -m
    return m

H = [F(0)]
for i in range(1, 60): H.append(H[-1] + F(1, i))
def Hstar(m): return sum(F(1, i) for i in range(1, m) if gcd(i, m) == 1)

print("(I) THE MOEBIUS-LOG^2 GRAMMAR (all exact in Q)")
# (a) H* Moebius form
okM = all(Hstar(m) == sum(F(mobius(d), d) * H[m // d - 1] for d in range(1, m + 1) if m % d == 0)
          for m in range(2, 40))
print("   H*(m) = sum_{d|m} (mu(d)/d) H_{m/d-1}, m<40:", okM)
# (b) [x^s] sum_d (mu(d)/d^2) log^2(1/(1-x^d)) = (2/s) H*(s)
#     log^2(1/(1-y)) = sum_{n>=2} (2 H_{n-1}/n) y^n  (classic); check the assembled claim.
NS = 40
coeff = [F(0)] * (NS + 1)
for d in range(1, NS + 1):
    mu = mobius(d)
    if mu == 0: continue
    for n in range(2, NS // d + 1):
        coeff[n * d] += F(mu, d * d) * 2 * H[n - 1] / n
okG = all(coeff[s] == 2 * Hstar(s) / s for s in range(2, NS + 1))
print("   [x^s] sum_d (mu(d)/d^2) log^2(1/(1-x^d)) = (2/s) H*(s), s<=40:", okG)
# (c) THM-826 cross-check: corridor segment constants A_{k+1} = (2/(k+1)) H*(k+1)
def farey_pairs(k):
    """consecutive denominator pairs (i,j) of F_k on the circle, ordered."""
    fr = sorted(set(F(a, q) for q in range(1, k + 1) for a in range(q)))
    fr.append(F(1))
    return [(fr[t].denominator, fr[t + 1].denominator if fr[t + 1] != 1 else 1)
            for t in range(len(fr) - 1)]
okA = True
for k in range(3, 14):
    A = sum(F(1, i * j) for (i, j) in farey_pairs(k) if i + j == k + 1)
    okA &= A == 2 * Hstar(k + 1) / (k + 1)
print("   Farey corridor constants A_{k+1} = (2/(k+1)) H*(k+1), k=3..13:", okA)
# (d) the LRC(15) corridor (first composite instance, k=13, depth 14):
k = 13
c14 = 2 * Hstar(14) / 14
pts = [F(1, 15), F(15, 211), F(1, 14)]
okC = True
for lam in pts:
    pred = c14 * (1 - 14 * lam)
    got = good_measure(list(range(1, 14)), lam)
    okC &= got == pred
print(f"   LRC(15) corridor m({{1..13}};lam) = (H*(14)/7)(1-14lam) on [1/15,1/14]: {okC}")
print(f"      H*(14) = {Hstar(14)} = 1+1/3+1/5+1/9+1/11+1/13; constant = {c14}")
print(f"      [prime contrast k=12: constant 2H_12/13 = {2*H[12]/13} (THM-853 II)]")

print()
print("(II) LOW-M RIGIDITY ASSEMBLY PROBE (shelf lens; klein-S308 handoff)")
rng = random.Random(20260715)
LAM = F(1, 13)
def covering_type(S): return good_measure(S, LAM) == 0
hits = []
tested = 0
# family A: deep-well mutations {1..12}+w and core swaps
cands = []
for w in [170, 176, 180, 182, 183, 190, 196, 205, 238, 250]:
    cands.append(list(range(1, 13)) + [w])
for _ in range(60):
    core = list(range(1, 13))
    i = rng.randrange(12)
    core[i] = rng.randrange(13, 30)
    w = rng.choice([168, 182, 196, 210, 240])
    cands.append(sorted(set(core + [w]))[:13] if len(set(core + [w])) == 13 else core + [w])
# family B: random small-element 13-sets
for _ in range(120):
    S = sorted(rng.sample(range(1, 26), 13))
    cands.append(S)
# family C: two-outlier (multi-killer face) samples
multi = []
for _ in range(25):
    core = sorted(rng.sample(range(1, 15), 11))
    S = core + sorted(rng.sample(range(20, 300), 2))
    multi.append(S)
for S in cands:
    if len(set(S)) != 13: continue
    g = gcd(*S) if False else 0
    from functools import reduce
    if reduce(gcd, S) != 1: continue
    tested += 1
    if covering_type(S):
        hits.append(S)
near_ap = lambda S: sum(1 for v in S if v <= 14) >= 10
out_ct = lambda S: sum(1 for v in S if v > 14)
bad = [S for S in hits if not near_ap(S)]
multi_bad = []
for S in multi:
    if covering_type(S): multi_bad.append(S)
print(f"   candidates tested: {tested}; covering-type (m(1/13)=0) found: {len(hits)}")
print(f"   ALL covering-type hits near-AP (>=10 in {{1..14}}): {not bad}"
      + (f"  VIOLATIONS: {bad[:3]}" if bad else ""))
print(f"   outlier(>14) counts among hits: { {out_ct(S) for S in hits} } (assembly predicts <= 1)")
print(f"   two-outlier samples that are covering-type: {len(multi_bad)} of {len(multi)}"
      f" (THM-726 face predicts 0)")
if hits:
    ms = {}
    for S in hits[:6]:
        # bracket M(S) in [14/183-eps, 1/13): check m at covering-min level
        m183 = good_measure(S, F(14, 183))
        ms[tuple(S[-2:])] = m183 > 0
    print(f"   sample hits still good at 14/183 (M > covering-min): {ms}")
print("\nDONE")
