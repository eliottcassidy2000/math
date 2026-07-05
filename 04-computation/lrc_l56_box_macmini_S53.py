#!/usr/bin/env python3
"""
mac-mini-2026-07-05-S53 -- HYP-4109 part 3: l=5,6 BOXED sweeps + structured corners.

Full chain-domain closure at l=5,6 is out of compute reach (volumes ~10^14+;
the pyramid kill rates DROP as free coords add orbit-window constraints).
Honest deliverables here:
  (1) BOXED exhaustive k <= 4 per coord (extends kps-S1's k <= 2 by two notches):
      l=5: C(12,5)*4^5 = 811,008;  l=6: C(12,6)*4^6 = 3,784,704 -- scan-first.
  (2) STRUCTURED CORNERS at l = 3..6 (where every extremal so far has lived):
      +13 blocks (all k=1 on a coordinate set), 14r towers (k_r = r), mixed
      block+deep, and the S52 floor family's extensions.  Exact M each.
  (3) The LEDGER anchors for what remains open (chain-domain minus box):
      A_5 = 1600/3 ~ 533, A_6 = 2291; ratio chain R_4..R_1 = 200/7, 300/13,
      1100/51, 24.  Plus the UNLIFTED-BLOCK RECIPROCAL sharpening: when C
      pushes min(base) = a up (1 in C), the anchor fee uses margin a/(a+12)
      via kps band_margin_reciprocal instead of the citation 1/(13-l).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import sys, time

sys.path.insert(0, '04-computation')
from lonely_profile import profile

T0 = time.time()
def log(m=""):
    print(m, flush=True)

BETA = F(2, 25)
AP = list(range(1, 13))

def M_exact(S):
    for cap in (11, 8, 6, 4, 3, 2):
        p = profile(sorted(S), F(1, cap))
        m = p.M()
        if m is not None:
            return m
    return None

def dq(x, q):
    x %= q
    return min(x, q - x)

def ok(x, q):
    return dq(x, q) * 25 >= 2 * q

def sieve_ok(W):
    return all(any(v % m == 0 for v in W) for m in range(2, 13))

QLIB = [25, 50] + [13 * u for u in range(2, 22)] + [q for q in range(8, 42) if q % 13]

def scan(W):
    for q in QLIB:
        if any(v % q == 0 for v in W):
            continue
        for a in range(1, q // 2 + 1):
            if gcd(a, q) != 1:
                continue
            if all(ok(a * v, q) for v in W):
                return (q, a)
    return None

# ---------------------------------------------------------------------------
log("=" * 78)
log("PART 1 -- boxed k <= 4 sweeps, l = 5 and l = 6 (beyond kps-S1's k <= 2)")
log("=" * 78)
for l in (5, 6):
    t0 = time.time()
    stats = dict(total=0, nonprim=0, sieved=0, scanned=0, exact=0, sub=0)
    worst = None
    K = 4 if l == 5 else 3
    kvecs = []
    def gen(pref):
        if len(pref) == l:
            kvecs.append(tuple(pref))
            return
        for k in range(1, K + 1):
            gen(pref + [k])
    gen([])
    for C in combinations(range(1, 13), l):
        base = [v for v in AP if v not in C]
        for kv in kvecs:
            W = sorted(base + [c + 13 * k for c, k in zip(C, kv)])
            stats['total'] += 1
            if not sieve_ok(W):
                stats['sieved'] += 1
                continue
            if reduce(gcd, W) != 1:
                stats['nonprim'] += 1
                continue
            hit = scan(W)
            if hit:
                stats['scanned'] += 1
                continue
            stats['exact'] += 1
            M = M_exact(W)
            if M < BETA:
                stats['sub'] += 1
                log(f"  << SUB-2/25 (l={l}): M={M} W={list(W)} C={C} k={kv}")
            if worst is None or M < worst[0]:
                worst = (M, tuple(W))
    log(f"l={l} boxed k<={K}: {stats}  t={time.time()-t0:.0f}s")
    if worst:
        log(f"   worst exact case: M = {worst[0]} at {list(worst[1])}")
    log(f"   VERDICT: {'>= 2/25 everywhere in the box' if stats['sub']==0 else 'SUB-2/25 FOUND'}")

# ---------------------------------------------------------------------------
log("\n" + "=" * 78)
log("PART 2 -- structured corners, l = 3..6, exact M")
log("=" * 78)
corners = []
# +13 blocks: all k = 1
for l in (3, 4, 5, 6):
    for C in combinations(range(1, 13), l):
        corners.append((C, tuple([1] * l), f"block+13 l={l}"))
# 14r towers: k_r = r for each lifted coord (the deep-well move per coord)
for l in (2, 3, 4):
    for C in combinations(range(7, 13), l):
        corners.append((C, tuple(C), f"14r-tower l={l}"))
# mixed: the S52 floor block {4,6} extended by one deep coordinate
for r in range(7, 13):
    corners.append(((4, 6, r), (1, 1, r), "block{4,6}+14r"))
# double blocks at height 2 (k=2): the 2nd-harmonic at doubled scale
for C in combinations(range(1, 13), 2):
    corners.append((C, (2, 2), "block+26 l=2"))

seen = set()
rows = []
t0 = time.time()
for C, kv, tag in corners:
    key = (C, kv)
    if key in seen:
        continue
    seen.add(key)
    base = [v for v in AP if v not in C]
    W = sorted(base + [c + 13 * k for c, k in zip(C, kv)])
    if reduce(gcd, W) != 1 or not sieve_ok(W):
        continue
    M = M_exact(W)
    rows.append((M, tuple(W), C, kv, tag))
rows.sort()
log(f"structured corners computed: {len(rows)}  [{time.time()-t0:.0f}s]")
log(f"{'M':>10}  {'tag':>16}  C / k / W (bottom 25)")
for M, W, C, kv, tag in rows[:25]:
    log(f"{str(M):>10}  {tag:>16}  C={list(C)} k={list(kv)}")
sub = [r for r in rows if r[0] < BETA]
log(f"\ncorners below 2/25: {len(sub)}")
for M, W, C, kv, tag in sub:
    log(f"   M = {M}  {tag}  C={list(C)} k={list(kv)}  W={list(W)}")

# ---------------------------------------------------------------------------
log("\n" + "=" * 78)
log("PART 3 -- the open-ledger anchors for l=5,6 (chain domain minus box)")
log("=" * 78)
for l in (5, 6):
    A = l * BETA * 12 / ((F(1, 13 - l) - BETA) * (1 - 2 * l * BETA))
    log(f"l={l}: anchor w_1 <= {A} (~{float(A):.0f}); chain R_4..R_1 = 200/7, 300/13, 1100/51, 24")
    # reciprocal sharpening per pattern: if min(base) = a >= 2, margin a/(a+12)
    best = []
    for C in combinations(range(1, 13), l):
        base = [v for v in AP if v not in C]
        a, b = min(base), max(base)
        m0 = F(a, a + b)
        if m0 > BETA:
            A2 = l * BETA * b / ((m0 - BETA) * (1 - 2 * l * BETA))
            best.append((A2, C, m0))
    best.sort()
    if best:
        log(f"   reciprocal-sharpened anchors: best {float(best[0][0]):.0f} at C={list(best[0][1])} "
            f"(base margin {best[0][2]}); {sum(1 for x in best if x[0] < A)} of {len(best)} patterns improve")
log(f"\nDONE [t = {time.time()-T0:.0f}s]")
