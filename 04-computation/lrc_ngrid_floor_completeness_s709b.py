"""
monad-explorer-2026-06-06-S709b — Is the n-GRID complete for the LRC floor 1/n?

Finding from s709: M(tight) is attained at t*=m/n (denominator exactly n), NOT on the
signed-LRC shell grid G_{2n-1}. THM-401: M attained at reduced pair-sum s=(v_a+v_b)/gcd;
for tight configs s=n. So the n-shell (pair-sum ≡0 mod n) controls the floor 1/n, while
2n-1 is the NEXT Farey pinch 2/(2n-1) (gap ABOVE the floor).

Reformulation tested:
  GRID-WITNESS  G(S):  exists m in {1..n-1} with  m*v_i ≢ 0 (mod n)  for all i.
    (Sufficient for M(S) >= 1/n via witness t=m/n.)
  Q1: Is G(S) NECESSARY too?  M(S) >= 1/n  <=>  G(S)?
  Q2: tight configs: active maximizing pair = mod-n shell-partner (a+b≡0 mod n)?
  Q3: mod-n mirror redundancy.

Fast integer M (no Fraction objects in the hot loop): compares r/d via cross-mult.
"""
import sys
from itertools import combinations
from math import gcd
from fractions import Fraction as F


def pinch_denoms(V):
    ds = set(V)
    L = len(V)
    for i in range(L):
        for j in range(i + 1, L):
            ds.add(V[i] + V[j])
            ds.add(abs(V[i] - V[j]))
    ds.discard(0)
    return ds


def M_fast(V):
    """Return (num, den) reduced of M = max_{t=m/d} min_i ||v_i t||, and a list of
       (d, m) achieving it. Pure-int hot loop; compare a/b vs c/d by cross mult."""
    best_n, best_d = -1, 1
    args = []
    for d in pinch_denoms(V):
        for m in range(1, d):
            # min over i of min(r, d-r), r = v_i*m mod d
            mn = d
            for v in V:
                r = (v * m) % d
                rr = d - r
                if r < rr:
                    if r < mn:
                        mn = r
                else:
                    if rr < mn:
                        mn = rr
                if mn == 0:
                    break
            # value mn/d  vs best_n/best_d
            if mn * best_d > best_n * d:
                best_n, best_d = mn, d
                args = [(d, m)]
            elif mn * best_d == best_n * d:
                args.append((d, m))
    g = gcd(best_n, best_d)
    return (best_n // g, best_d // g), args


def grid_witness(V, n):
    for m in range(1, n):
        if all((m * v) % n != 0 for v in V):
            return m
    return None


def gcd_list(V):
    g = 0
    for v in V:
        g = gcd(g, v)
    return g


def nrm_frac(x):
    r = x % 1
    return min(r, 1 - r)


pr = print
def flush():
    sys.stdout.flush()


# -------------------------------------------------------------------------
pr("=" * 78)
pr("Q1: M(S) >= 1/n   <=>   exists m in 1..n-1 : all m*v_i ≢ 0 (mod n) ?")
pr("    (exhaustive over (n-1)-subsets of [1..W], gcd=1)")
pr("=" * 78)
flush()
for n in range(3, 9):
    W = 2 * n + 2
    both = only_M = only_G = neither = 0
    ex = []
    for V in combinations(range(1, W + 1), n - 1):
        if gcd_list(V) != 1:
            continue
        (Mn, Md), _ = M_fast(V)
        hasM = (Mn * n >= Md)        # M >= 1/n  <=> Mn/Md >= 1/n
        hasG = grid_witness(V, n) is not None
        if hasM and hasG:
            both += 1
        elif hasM and not hasG:
            only_M += 1
            if len(ex) < 6:
                ex.append((V, (Mn, Md)))
        elif hasG and not hasM:
            only_G += 1
        else:
            neither += 1
    pr(f" n={n} (W={W}): both={both}  M-only(M>=1/n,NO grid)={only_M}  "
       f"grid-only(should be 0)={only_G}  neither={neither}")
    if ex:
        pr("    NECESSITY-BREAKING (M>=1/n but no grid witness):")
        for V, (a, b) in ex:
            pr(f"       {V}  M={a}/{b}")
    flush()

# -------------------------------------------------------------------------
pr("")
pr("=" * 78)
pr("Q2: tight configs: active maximizing pair a mod-n shell-partner (a+b≡0 mod n)?")
pr("    vs signed-LRC (2n-1) shell-partner")
pr("=" * 78)
flush()
for n in range(4, 9):
    C = 2 * n - 1
    W = 2 * n
    thr = F(1, n)
    for V in combinations(range(1, W + 1), n - 1):
        if gcd_list(V) != 1:
            continue
        (Mn, Md), args = M_fast(V)
        if F(Mn, Md) != thr:
            continue
        d, m = args[0]
        t = F(m, d)
        active = [v for v in V if nrm_frac(v * t) == thr]
        modn_sp = [(a, b) for a in active for b in active if a < b and (a + b) % n == 0]
        modC_sp = [(V[i], V[j]) for i in range(len(V)) for j in range(i + 1, len(V))
                   if (V[i] + V[j]) % C == 0]
        pr(f" n={n} V={V}: t*={t} active={active}")
        pr(f"      mod-{n} shell(active)={modn_sp if modn_sp else 'NONE'};  "
           f"mod-{C} shell(all)={modC_sp if modC_sp else 'NONE'}")
    flush()

# -------------------------------------------------------------------------
pr("")
pr("=" * 78)
pr("Q3: mod-n mirror redundancy (a+b≡0 mod n => ||a k/n||=||b k/n|| all k)")
pr("=" * 78)
viol = checked = 0
for n in range(3, 40):
    for a in range(1, 4 * n):
        for b in range(1, 4 * n):
            if (a + b) % n != 0:
                continue
            for k in range(n):
                ra = (a * k) % n
                rb = (b * k) % n
                if min(ra, n - ra) != min(rb, n - rb):
                    viol += 1
                checked += 1
pr(f"  {checked} checks, {viol} violations => {'HOLDS' if viol == 0 else 'FAILS'}")
flush()
