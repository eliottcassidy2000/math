#!/usr/bin/env python3
"""
The DEFLATION/DESCENT approach to the compact-core residual  (boxeph-2026-07-17-S83)
====================================================================================

THM-1003 (fill-1 perturbation) is the f=1 case of a general DESCENT: perturbing
around a resonance center a/b rescales the STRANDED runners S_b={v: b|v} by 1/b
into a sub-family W = S_b/b.  Loneliness of W (a smaller LRC instance, SETTLED by
LRC(|W|+1) when |W|<=12) transfers back to loneliness of V IF the kick can be
kept small enough to protect the body (non-multiples of b).

Witness:  t = a/b + u*/b,  u* a lonely time of W.  Then
   stranded v=b*w:  ||v t|| = ||w u*|| >= 1/(|W|+1) >= 1/13 >= 1/14     (exact)
   body v (b∤v):    ||v t|| >= 1/b - v|u*|/b >= 1/14  iff  v|u*| <= (14-b)/14
So the base-case condition is: W has a lonely (>=1/14) time u* with
   |u*| <= (14-b)/(14*B),   B = body max = max{v: b∤v}.

This script: (1) computes M exactly via the pair-sum bound (THM-999: q|v_i+v_j,
q<=2Vmax); (2) tests whether the deflation witness certifies M>=1/14 on compact
covering families; (3) measures the required vs available window to find the
clean base-case condition; (4) maps where deflation reaches vs the residual.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

N = 14

def norm(x):                      # distance to nearest integer, exact
    r = x % 1
    return min(r, 1 - r)

def min_gap(V, t):
    return min(norm(F(v) * t) for v in V)

def is_covering(V, n=N):
    return all(any(v % b == 0 for v in V) for b in range(2, n + 1))

# --- exact M via the pair-sum denominator bound (THM-999 / klein THM-1002) ---
def exact_M(V):
    """M(V) = max over t=a/q (q | v_i+v_j, q<=2max) of min_v ||v a/q||, exact."""
    Vmax = max(V)
    best = F(0)
    qs = set()
    for i in range(len(V)):
        for j in range(i, len(V)):
            s = V[i] + V[j]
            for d in range(1, s + 1):
                if s % d == 0:
                    qs.add(d)
    qs.add(N)  # always include the threshold denominator
    for q in qs:
        for a in range(1, q):
            if gcd(a, q) == 1:
                m = min(min((v * a) % q, q - (v * a) % q) for v in V)
                cand = F(m, q)
                if cand > best:
                    best = cand
    return best

# --- W's lonely time closest to 0 (rational search, W bounded) -------------
def W_lonely_near_zero(W, thr=F(1, N)):
    """Smallest |u|>0, u=a/q rational (q|w_i+w_j, q<=2max), with min_W ||w u||>=thr."""
    if not W:
        return F(0)  # empty: any u works; return 0
    Wmax = max(W)
    cands = []
    qs = set([N])
    for i in range(len(W)):
        for j in range(i, len(W)):
            s = W[i] + W[j]
            for d in range(1, s + 1):
                if s % d == 0:
                    qs.add(d)
    for q in qs:
        for a in range(1, q):
            if gcd(a, q) == 1:
                u = F(a, q)
                if min(min((w * a) % q, q - (w * a) % q) for w in W) * 1 >= thr * q / q:
                    # exact check
                    if min(norm(F(w) * u) for w in W) >= thr:
                        uu = u if u <= F(1, 2) else 1 - u
                        cands.append(uu)
    return min(cands) if cands else None

def deflation_certificate(V):
    """Try every circle b; return (b, u*, t, mingap, window_ok) for the first that certifies."""
    Vmax = max(V)
    for b in range(2, N):
        S = [v for v in V if v % b == 0]
        if not S or len(S) > 12:
            continue
        body = [v for v in V if v % b != 0]
        B = max(body) if body else 0
        W = [v // b for v in S]
        # window the base-case condition allows for |u*|
        window = F(N - b, N * B) if B > 0 else F(1, 2)
        u0 = W_lonely_near_zero(W)          # closest lonely-for-W time to 0
        if u0 is None:
            continue
        within = (u0 <= window)
        # build witness t = 1/b + u0/b  (a=1)
        t = F(1, b) + u0 / b
        mg = min_gap(V, t)
        if mg >= F(1, N):
            return (b, u0, window, within, t, mg, len(W))
    return None

# ---------------------------------------------------------------------------
# Compact covering families: bounded, no isolated far element.
# Build a representative pool (13 speeds, covering, Vmax bounded).
# ---------------------------------------------------------------------------
def compact_covering_pool(vmax_cap=30, limit=60):
    """Greedy/structured pool of covering 13-sets with all speeds <= vmax_cap."""
    pool = []
    # start from {2..14} and its bounded perturbations
    base = list(range(2, 15))  # {2..14}, covering, Vmax=14
    pool.append(base)
    # families {1..k} + a few bounded extras to reach 13 speeds, staying covering
    import itertools
    smalls = list(range(1, 15))
    tried = 0
    for extra_cnt in [0, 1, 2]:
        for drop in itertools.combinations(range(1, 15), extra_cnt):
            fam = [x for x in range(1, 15) if x not in drop]
            # need exactly 13; {1..14} has 14, drop 'extra_cnt' -> 14-extra_cnt
            if len(fam) != 13:
                continue
            if is_covering(fam):
                pool.append(sorted(fam))
                tried += 1
                if tried > 25:
                    break
    # add some Vmax in (14, 30] compact families: replace one small by a mid value
    for rep in range(15, vmax_cap + 1):
        fam = list(range(1, 14))  # {1..13}
        fam = fam[:-1] + [rep]    # {1..12, rep}
        if is_covering(fam):
            pool.append(sorted(fam))
    # dedup
    uniq = []
    seen = set()
    for f in pool:
        key = tuple(f)
        if key not in seen and len(f) == 13:
            seen.add(key); uniq.append(f)
    return uniq[:limit]

if __name__ == "__main__":
    print("=" * 78)
    print("DEFLATION on compact covering families (Vmax bounded). LRC(14), thr 1/14.")
    print("witness t=1/b+u*/b, u*=lonely-for-(stranded/b) time near 0; certifies if")
    print("min_V||v t||>=1/14. Base-case cond: |u*| <= (14-b)/(14 B).")
    print("=" * 78)
    pool = compact_covering_pool(vmax_cap=30, limit=50)
    n_cert = 0; n_within = 0; n_total = 0; residual = []
    for V in pool:
        M = exact_M(V)
        lonely = (M >= F(1, N))
        cert = deflation_certificate(V)
        n_total += 1
        tag = ""
        if cert:
            b, u0, window, within, t, mg, wlen = cert
            n_cert += 1
            if within: n_within += 1
            tag = (f"DEFLATE b={b} |W|={wlen} u*={u0} win={window} "
                   f"within={within} min={mg}")
        else:
            residual.append(V)
            tag = "NO deflation witness"
        print(f"  V(max={max(V):2d}) M={str(M):8s} lonely={lonely}  {tag}")
    print("-" * 78)
    print(f"deflation certified {n_cert}/{n_total}; base-case window held {n_within}/{n_cert}")
    print(f"residual (no deflation witness found): {len(residual)}")
    for V in residual[:8]:
        print(f"    {V}  M={exact_M(V)}")
