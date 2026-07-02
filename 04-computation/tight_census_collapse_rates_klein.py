#!/usr/bin/env python3
"""
klein-2026-07-01-S87 -- HYP-3835 follow-up: THE TIGHT CENSUS WITH COLLAPSE RATES (small n)

Question: is the collapse rate c(S) = lim_{r->1/n^-} Lambda_S(r)/(1-nr) UNIVERSAL over the
tight locus (as it is over {AP, GW} at n=14, both = (2/n)Hx(n)), or does it VARY?

Key structural dichotomy (mini-theorem, this session):
  If M(S)=1/n and no v in S is divisible by n, then EVERY k/n (gcd(k,n)=1) is a maximizer,
  S must contain a speed = +-k^{-1} mod n for each k, and
      c(S) = (1/n) * sum_k [ 1/max(S cap class(+k^{-1})) + 1/max(S cap class(-k^{-1})) ].
  If S CONTAINS a multiple of n, all k/n are dead (m=0) and witnesses live at deeper
  denominators -> c(S) escapes the unit-pair formula entirely.

Census: exhaustively find ALL tight sets (M = 1/n exactly) with k = n-1... wait, k runners,
threshold 1/(k+1); n := k+1. Enumerate primitive k-subsets of {1..Vmax}, filter by float
Lambda, confirm exactly, compute exact collapse rates.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools

# ---------- fast float Lambda ----------
def lam_float(S, r):
    ivs = []
    for v in S:
        rv = r / v
        step = 1.0 / v
        for a in range(v + 1):
            c = a * step
            lo, hi = c - rv, c + rv
            if hi <= 0 or lo >= 1:
                continue
            ivs.append((max(lo, 0.0), min(hi, 1.0)))
    ivs.sort()
    tot, cl, ch = 0.0, None, None
    for lo, hi in ivs:
        if ch is None:
            cl, ch = lo, hi
        elif lo <= ch:
            ch = max(ch, hi)
        else:
            tot += ch - cl
            cl, ch = lo, hi
    if ch is not None:
        tot += ch - cl
    return 1.0 - tot

# ---------- exact machinery (as in lonely_profile_farey_slope_klein.py) ----------
def Lambda(S, r):
    ivs = []
    for v in S:
        rv = r / v
        for a in range(v + 1):
            c = F(a, v)
            lo, hi = c - rv, c + rv
            if hi <= 0 or lo >= 1:
                continue
            ivs.append((max(lo, F(0)), min(hi, F(1))))
    ivs.sort()
    tot = F(0)
    cl = ch = None
    for lo, hi in ivs:
        if ch is None:
            cl, ch = lo, hi
        elif lo <= ch:
            ch = max(ch, hi)
        else:
            tot += ch - cl
            cl, ch = lo, hi
    if ch is not None:
        tot += ch - cl
    return 1 - tot

def m_at(S, t):
    best = None
    for v in S:
        x = (v * t) % 1
        d = min(x, 1 - x)
        if best is None or d < best:
            best = d
    return best

def M_exact(S):
    qs = set()
    for u, v in itertools.combinations(S, 2):
        qs.add(u + v); qs.add(abs(u - v))
    for v in S:
        qs.add(2 * v)
    qs.discard(0)
    best, wit = F(0), []
    for q in qs:
        for c in range(q + 1):
            t = F(c, q)
            m = m_at(S, t)
            if m > best:
                best, wit = m, [t]
            elif m == best:
                wit.append(t)
    return best, sorted(set(wit))

def collapse_rate(S, thresh):
    s1 = Lambda(S, thresh * (1 - F(1, 10**4))) / F(1, 10**4)
    s2 = Lambda(S, thresh * (1 - F(1, 2 * 10**4))) / F(1, 2 * 10**4)
    return s1, (s1 == s2)

def Hx(n):
    return sum(F(1, u) for u in range(1, n) if gcd(u, n) == 1)

# ---------- census ----------
def census(n, vmax):
    k = n - 1
    thr = 1.0 / n
    thrF = F(1, n)
    cands = []
    for c in itertools.combinations(range(1, vmax + 1), k):
        if reduce(gcd, c) != 1:
            continue
        # necessary: dead-or-witness test at 1/n (fast float): Lambda just above 1/n must be ~0,
        # just below must be small
        hi = lam_float(c, thr * (1 + 1e-9))
        if hi > 1e-12:      # M > 1/n strictly (loose) -- keep only M <= 1/n
            continue
        lo = lam_float(c, thr * (1 - 1e-7))
        # M < 1/n would be a counterexample (should not happen); M = 1/n gives tiny positive lo
        cands.append((c, lo))
    tight = []
    for c, lo in cands:
        Mv, wits = M_exact(list(c))
        if Mv == thrF:
            tight.append((c, wits))
        elif Mv < thrF:
            print(f"  !!! COUNTEREXAMPLE?? {c} M={Mv} < 1/{n}")
    return tight

print("=" * 90)
print("TIGHT CENSUS WITH COLLAPSE RATES -- k=n-1 runners, threshold 1/n")
print("=" * 90)
for n, vmax in [(4, 40), (5, 40), (6, 30), (7, 24)]:
    print(f"\n--- n={n} (k={n-1} runners), speeds <= {vmax}, primitive ---")
    uni = F(2, n) * Hx(n)
    print(f"    universal candidate (2/n)Hx(n) = {uni} = {float(uni):.6f}")
    tset = census(n, vmax)
    print(f"    tight sets found: {len(tset)}")
    seen_rates = {}
    for S, wits in tset:
        s, ok = collapse_rate(list(S), F(1, n))
        multn = any(v % n == 0 for v in S)
        # dilation-primitive label: is S a dilate of a smaller tight set? (c*S' with c>1)
        is_dilate = False
        for d in range(2, max(S) + 1):
            if all(v % d == 0 for v in S):
                is_dilate = True
                break
        tag = []
        if multn: tag.append(f"has-mult-of-{n}")
        if s == uni: tag.append("UNIVERSAL")
        wl = f"#wit={len(wits)}"
        if len(wits) <= 6:
            wl += " at " + ",".join(str(t) for t in wits)
        print(f"    {str(S):32s} c = {str(s):>12} = {float(s):.6f}  {wl:40s} {' '.join(tag)}")
        seen_rates.setdefault(s, []).append(S)
    print(f"    DISTINCT collapse rates at n={n}: {len(seen_rates)}")
    for s, Ss in sorted(seen_rates.items()):
        print(f"      c={str(s):>12} = {float(s):.6f}  x{len(Ss)} sets")

print("\nDONE.")

# ---------- extension: n=8 and deeper n=5 (refined law test) ----------
print("\n" + "=" * 90)
print("EXTENSION: n=8 (k=7, speeds <= 26) and n=5 deeper (speeds <= 75)")
print("=" * 90)
for n, vmax in [(5, 75), (8, 26)]:
    print(f"\n--- n={n}, speeds <= {vmax} ---")
    uni = F(2, n) * Hx(n)
    tset = census(n, vmax)
    print(f"    tight sets found: {len(tset)}   (2/n)Hx = {uni}")
    for S, wits in tset:
        s, ok = collapse_rate(list(S), F(1, n))
        resid = sorted(set(v % n for v in S))
        # class-max formula check
        pred = F(0)
        formula_ok = True
        for k in range(1, n):
            if gcd(k, n) != 1:
                continue
            kin = pow(k, -1, n)
            plus = [v for v in S if v % n == kin]
            minus = [v for v in S if v % n == (n - kin) % n]
            if plus and minus:
                pred += F(1, n) * (F(1, max(plus)) + F(1, max(minus)))
            else:
                formula_ok = False
        tag = "UNIVERSAL" if s == uni else f"c<uni ({float(s):.4f}<{float(uni):.4f})" if s < uni else "c>uni ??"
        fm = f"class-max-formula={'MATCH' if (formula_ok and pred == s) else ('pred '+str(pred)) }"
        print(f"    {str(S):36s} c={str(s):>10}  {tag:28s} {fm}  residues mod {n}: {resid}")
