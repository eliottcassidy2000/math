#!/usr/bin/env python3
"""claudebox-2026-06-11-S7 part C — do the deep joint failures break HYP-2197?
HYP-2197 (S622/THM-412): twisted-shell(m<=27) UNION Criterion B' has zero residual on
multiple-of-14 configs. Here: find the joint-failure configs whose minimal witness
modulus exceeds 27, then test (a) plain twisted shells m<=27, (b) Criterion B'
(THM-398 section-4 width form: some component of the level-1/14 safe set of S minus v
is wider than the danger-arc diameter 2/(14 v), for v a multiple of 14 — also tested
for every v in S as the generous form).
Exact arithmetic (Fraction) throughout."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

TARGET = F(1, 14)

def witness_at_q(S, q):
    B = q // 14
    prof = [v % q for v in S]
    for a in range(1, q):
        if gcd(a, q) != 1: continue
        ok = True
        for p in prof:
            r = (p * a) % q
            if min(r, q - r) <= B: ok = False; break
        if ok: return a
    return None

def shells_le_27(S):
    for q in range(2, 28):
        a = witness_at_q(S, q)
        if a: return q, a
    return None, None

def safe_components(U):
    """Components of {t in [0,1): ||ut|| > 1/14 for all u in U}, exact.
    Danger arcs (closed): [(14k-1)/(14u), (14k+1)/(14u)]."""
    arcs = []
    for u in U:
        for k in range(0, u + 1):
            a, b = F(14 * k - 1, 14 * u), F(14 * k + 1, 14 * u)
            arcs.append((max(a, F(0)), min(b, F(1))))
            if a < 0: arcs.append((F(0), b))          # wrap handled by clamping; arc at 0
    arcs = sorted(arcs)
    merged = []
    for a, b in arcs:
        if merged and a <= merged[-1][1]: merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    gaps = []
    for i in range(len(merged) - 1):
        gaps.append((merged[i][1], merged[i + 1][0]))
    # note: t=0 and t=1 are danger (k=0,k=u arcs), so no wraparound gap survives
    return [g for g in gaps if g[1] > g[0]]

def criterion_Bprime(S, only_mult14=True):
    """Width-form B': exists v (mult of 14, or any v if only_mult14=False) and a
    component I of safe(S minus v) with width(I) > 2/(14 v)."""
    cands = [v for v in S if v % 14 == 0] if only_mult14 else list(S)
    for v in cands:
        comps = safe_components([u for u in S if u != v])
        rho2 = F(2, 14 * v)
        if any(b - a > rho2 for a, b in comps): return True, v
    return False, None

def min_witness_modulus(S, qmax=400):
    for q in range(2, qmax + 1):
        if witness_at_q(S, q): return q
    return None

# regenerate joint failures with deep minimal moduli (reuse the family generators,
# but here just rescan and keep configs with min modulus > 27)
def head_B_fails(S):
    if all(v % 3 for v in S): return False
    if all(v % 9 for v in S): return False
    if any(v % 27 == 0 for v in S): return True
    units = [v % 27 for v in S if v % 3]
    for a in range(1, 27):
        if gcd(a, 27) != 1: continue
        if all((a * u) % 27 not in (1, 26) for u in units): return False
    return True

def safe_window(S_other, t0):
    L, R = t0 - 1, t0 + 1
    for v in S_other:
        x = v * t0; kf = int(x)
        cl, cr = [], []
        for k in (kf - 1, kf, kf + 1):
            lo, hi = k - TARGET, k + TARGET
            if hi <= x: cl.append(hi)
            if lo >= x: cr.append(lo)
            if lo < x < hi: return None
        L = max(L, max(cl) / v if cl else t0 - 1)
        R = min(R, min(cr) / v if cr else t0 + 1)
        if L >= R: return None
    return (L, R)

def core_dodge_in(core, L, R):
    arcs = []
    for v in core:
        for k in range(int(v * L) - 1, int(v * R) + 2):
            a, b = F(k - TARGET) / v, F(k + TARGET) / v
            if b > L and a < R: arcs.append((a, b))
    arcs.sort(); cur = L
    for a, b in arcs:
        if a > cur: return (cur + min(a, R)) / 2
        cur = max(cur, b)
        if cur >= R: return None
    return (cur + R) / 2 if cur < R else None

def head_A_fails(S):
    for d in (14, 7, 2):
        core = [v for v in S if v % d == 0]
        rest = [v for v in S if v % d]
        if not core: return False
        for b in range(1, d):
            if gcd(b, d) != 1: continue
            win = safe_window(rest, F(b, d))
            if win and core_dodge_in(core, *win) is not None: return False
    return True

deep = []
def consider(S):
    S = sorted(S)
    if len(set(S)) != 13 or reduce(gcd, S) != 1: return
    if not any(v % 14 == 0 for v in S): return
    if head_A_fails(S) and head_B_fails(S):
        mq = min_witness_modulus(S)
        if mq and mq > 27: deep.append((S, mq))

for r in range(1, 1200):
    if r % 7 == 0: continue
    consider([7 * k for k in range(1, 13)] + [r])
for r in range(1, 1200):
    if r % 7 == 0 or r == 54: continue
    consider([14 * k for k in range(1, 7)] + [7, 21, 35, 49, 77, 54] + [r])
random.seed(14)
for trial in range(4000):
    c = random.randint(6, 12)
    ks = random.sample(range(1, 19), c)
    if not any((7 * k) % 14 == 0 for k in ks): ks[0] = 2 * random.randint(1, 9)
    consider([7 * k for k in ks] +
             random.sample([v for v in range(1, 120) if v % 7], 13 - c))

print(f"deep joint failures (min witness modulus > 27): {len(deep)}")
for S, mq in deep:
    q27, a27 = shells_le_27(S)
    bp14, v14 = criterion_Bprime(S, only_mult14=True)
    bpany, vany = criterion_Bprime(S, only_mult14=False)
    hyp2197 = (q27 is not None) or bp14
    print(f"  S = {S}")
    print(f"    min witness modulus = {mq}; shells<=27: {q27}; B'(mult-14): {bp14} (v={v14}); "
          f"B'(any v): {bpany} (v={vany})")
    print(f"    HYP-2197 [shells<=27 u B'(mult14)] covers it: {hyp2197}")

# also push the same families further (r up to 4000) hunting more deep cases
extra = 0
for r in range(1200, 4000):
    if r % 7 == 0: continue
    S = sorted([7 * k for k in range(1, 13)] + [r])
    if reduce(gcd, S) != 1 or len(set(S)) != 13: continue
    if head_A_fails(S) and head_B_fails(S):
        mq = min_witness_modulus(S, 500)
        if mq and mq > 27:
            q27, _ = shells_le_27(S)
            bp14, v14 = criterion_Bprime(S, only_mult14=True)
            if q27 is None and not bp14:
                extra += 1
                print(f"  EXTRA candidate vs HYP-2197: {S} minq={mq} B'14={bp14}")
print(f"extended fam1 scan r<4000: HYP-2197-violating candidates = {extra}")
