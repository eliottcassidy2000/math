#!/usr/bin/env python3
"""
Signed-LRC pairwise gap: the n-ASYMPTOTIC of  inf_S Gstar(S).
monad-compute-2026-06-06.

Open handoff (THM-429 / HYP-2296, monad-explorer-S3 SESSION-LOG):

    n * inf_S Gstar(S)  ->  ?     decay to 1/n^2  (n*inf -> 0)
                                  or true Theta(1/n) floor (n*inf -> c > 0)?

Exhaustive search (signed_lrc_inf_highB) is capped at tiny n (<=8) and modest B,
and at each n the inf is still DROPPING with B, so the exhaustive table only
gives a slowly-tightening UPPER bound on inf at a few n.  To see the n-TREND we
need the best-known inf upper bound across a WIDER n range.  This script does a
randomized + greedy + structured-family search for the SMALLEST Gstar(S) at each
n = 6..11, over gcd-1 speed sets with entries <= B, and reports:

    n,  best Gstar = k/q (exact),  n*Gstar,  the minimizer S, and the binding
    pair {a,b} with a+b = q (the cut-exposed shell-partner, HYP-2296).

Gstar(S) = max over cuts (A,B) of M(W),   W = {|v_i-v_j| same-side} u {v_i+v_j across}.
THM-426.  Exact rational arithmetic throughout; no floating point in the optimum.

This is an UPPER bound on inf (we may miss the true minimizer), but a DECREASING
upper bound across n that stays >> 1/n^2 is evidence for a Theta(1/n) floor,
while one that plunges toward 1/n^2 supports decay.  We also include the explicit
"small-cluster" families the S3 handoff flagged as candidate drivers of q->inf.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd
from functools import reduce
import sys, random

random.seed(20260606)   # reproducible (Math.random/time are fine in plain python)

def norm(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def candidate_times(W):
    W = sorted(set(abs(w) for w in W if w != 0))
    T = set()
    for w in W:
        for a in range(0, w):
            T.add(F(2 * a + 1, 2 * w))
    for x, y in combinations(W, 2):
        for d in (x + y, x - y):
            if d == 0: continue
            d = abs(d)
            for a in range(1, d):
                T.add(F(a, d))
    return T

def maximin_at_least(Ws, thresh):
    """max_t min_w ||w t||, abandoning early once it exceeds thresh (strict)."""
    best = F(0)
    for t in candidate_times(Ws):
        ok = True; m = F(1, 2)
        for w in Ws:
            nv = norm(w * t)
            if nv <= best:
                ok = False; break
            if nv < m: m = nv
        if ok and m > best:
            best = m
            if best > thresh:
                return best, True
    return best, False

def cut_Ws(V, side):
    r = len(V); W = []
    for i, j in combinations(range(r), 2):
        W.append(abs(V[i] - V[j]) if side[i] == side[j] else V[i] + V[j])
    return sorted(set(w for w in W if w != 0))

def gstar_capped(V, cap):
    """Gstar(V) = max over cuts of M(W); abandon set once some cut > cap."""
    r = len(V); best = F(0)
    for tail in product([0, 1], repeat=r - 1):
        side = (0,) + tail
        Ws = cut_Ws(V, side)
        m, reached = maximin_at_least(Ws, cap)
        if m > best: best = m
        if reached:
            return best, False     # exceeds running best; not a candidate
    return best, True              # fully evaluated

def binding_pair(V, g):
    """Find a pair {a,b} in the relative-speed system with a+b = denom(g)
    (the cut-exposed shell-partner, HYP-2296).  Best-effort report."""
    q = g.denominator
    for i, j in combinations(range(len(V)), 2):
        if V[i] + V[j] == q:
            return (V[i], V[j])
    # else maybe a same-side difference equals q (rare) -- report any |vi-vj|=q
    return None

def random_sets(r, B, count):
    seen = set()
    out = []
    tries = 0
    while len(out) < count and tries < count * 40:
        tries += 1
        s = tuple(sorted(random.sample(range(1, B + 1), r)))
        if s in seen: continue
        seen.add(s)
        if reduce(gcd, s) == 1:
            out.append(s)
    return out

def structured_families(r, B):
    """Candidate floor-breakers: small consecutive cluster + spread tail.
    The S3 minimizers (2,3,4,8,11),(4,5,8,10,15),(3,8,9,11,18) all carry a
    small near-consecutive low block.  Generate such shapes deterministically."""
    fams = []
    # consecutive low block of length L (2..r) then arithmetic-ish tail
    for L in range(2, r + 1):
        low = list(range(2, 2 + L))          # {2,3,...} avoid 1 (keeps gcd via tail)
        need = r - L
        if need == 0:
            s = tuple(sorted(set(low)))
            if len(s) == r: fams.append(s)
            continue
        # tail: spread values up to B
        step = max(1, (B - (low[-1] + 1)) // max(1, need))
        tail = [min(B, low[-1] + 1 + k * step) for k in range(1, need + 1)]
        s = sorted(set(low + tail))
        if len(s) == r:
            fams.append(tuple(s))
    # plain AP {1..r} (known Gstar = 1/r, the non-breaker baseline)
    fams.append(tuple(range(1, r + 1)))
    # dedup, gcd-1 only
    out = []
    seen = set()
    for s in fams:
        if s in seen: continue
        seen.add(s)
        if reduce(gcd, s) == 1:
            out.append(s)
    return out

def search_n(n, B, n_random):
    r = n - 1
    pool = structured_families(r, B) + random_sets(r, B, n_random)
    inf_g = F(1, 2); inf_sets = []
    for V in pool:
        g, is_cand = gstar_capped(V, inf_g)
        if not is_cand:
            continue
        if g < inf_g:
            inf_g = g; inf_sets = [V]
        elif g == inf_g:
            inf_sets.append(V)
    return inf_g, inf_sets, len(pool)

def main():
    print("=" * 78)
    print("SIGNED-LRC  n-ASYMPTOTIC of inf_S Gstar(S)   monad-compute")
    print("randomized+structured search (UPPER bound on the true inf)")
    print("=" * 78)
    # (n, B, n_random)  -- bigger n gets fewer random sets (cut cost 2^{n-2})
    plan = [(6, 22, 6000), (7, 20, 4000), (8, 18, 2500),
            (9, 16, 1200), (10, 15, 500), (11, 14, 250)]
    if len(sys.argv) > 1:
        plan = [(int(sys.argv[1]), int(sys.argv[2]),
                 int(sys.argv[3]) if len(sys.argv) > 3 else 2000)]
    print(f"\n{'n':>3} {'best Gstar':>12} {'float':>9} {'n*Gstar':>8} {'1/n':>8} "
          f"{'brk?':>5}  minimizer / binding pair")
    print("-" * 78)
    rows = []
    for n, B, nr in plan:
        g, sets, npool = search_n(n, B, nr)
        brk = "YES" if g < F(1, n) else "no"
        V = sets[0] if sets else None
        bp = binding_pair(V, g) if V else None
        rows.append((n, g))
        print(f"{n:>3} {str(g):>12} {float(g):>9.5f} {float(n*g):>8.4f} "
              f"{float(F(1,n)):>8.5f} {brk:>5}  {V}  pair={bp}  (B={B},pool={npool})")
        sys.stdout.flush()
    print("-" * 78)
    print("\nn*Gstar series (decreasing-with-B upper bounds; watch the n-trend):")
    for n, g in rows:
        print(f"   n={n:2d}: n*inf <= {float(n*g):.4f}   inf <= {str(g):>10} = {float(g):.5f}"
              f"   (1/n^2 = {1.0/n**2:.5f})")
    print("\nReading: if n*inf stays bounded away from 0 (and inf >> 1/n^2),")
    print("the Theta(1/n) floor is supported; if n*inf plunges toward 0 with")
    print("inf approaching 1/n^2, decay is supported.")

if __name__ == "__main__":
    main()
