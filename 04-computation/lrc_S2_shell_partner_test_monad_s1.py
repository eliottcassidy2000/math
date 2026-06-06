#!/usr/bin/env python3
"""
lrc_S2_shell_partner_test_monad_s1.py   (monad-explorer-2026-06-06-S1)

TEST of HYP-2290: is the signed second clock (q=2n-1) the S_2 (j=2, additive-energy)
face of the covering-depth distribution (THM-406)?  Concretely:
  (b) does the n=8 worry-set SPLIT (AP_8 shell-free vs the two shell-partner configs,
      all M=1/8) show as a MEASURABLE S_2 EXCESS on the shell-partner configs?
  + per-pair overlaps: is the shell-partner PAIR itself the locus of the excess?

Exact rational arithmetic. Danger arc D_i = {t: ||v_i t|| < δ}, δ=1/n, is a union of v_i
open intervals with rational endpoints (k±δ)/v_i = (kn±1)/(n v_i). S_j = Σ_{|I|=j} μ(∩ D_i).
"""
from fractions import Fraction as Fr
from itertools import combinations

def danger_intervals(v, n):
    """D_v = {t in [0,1): ||v t|| < 1/n} as a sorted list of (lo,hi) Fraction intervals."""
    d = Fr(1, n)
    ivs = []
    for k in range(0, v + 1):
        lo = (Fr(k) - d) / v
        hi = (Fr(k) + d) / v
        # clip to [0,1)
        lo = max(lo, Fr(0)); hi = min(hi, Fr(1))
        if lo < hi:
            ivs.append((lo, hi))
    # merge touching/overlapping
    ivs.sort()
    merged = []
    for lo, hi in ivs:
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged

def measure(ivs):
    return sum(hi - lo for lo, hi in ivs)

def intersect(A, B):
    """measure of intersection of two interval-lists, exact."""
    res = []
    i = j = 0
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0]); hi = min(A[i][1], B[j][1])
        if lo < hi:
            res.append((lo, hi))
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return res

def depth_profile(speeds, n):
    """Compute p_k = μ(N=k) exactly by overlaying all arc endpoints, and S_j."""
    Ds = [danger_intervals(v, n) for v in speeds]
    # build event points
    pts = set([Fr(0), Fr(1)])
    for D in Ds:
        for lo, hi in D:
            pts.add(lo); pts.add(hi)
    pts = sorted(pts)
    p = {}
    for a, b in zip(pts, pts[1:]):
        mid = (a + b) / 2
        depth = sum(1 for D in Ds for lo, hi in D if lo < mid < hi)
        p[depth] = p.get(depth, Fr(0)) + (b - a)
    return p, Ds

def Sj(Ds, j):
    """overlap sum S_j = Σ_{|I|=j} μ(∩_{i in I} D_i)."""
    tot = Fr(0)
    for I in combinations(range(len(Ds)), j):
        cur = Ds[I[0]]
        for idx in I[1:]:
            cur = intersect(cur, Ds[idx])
            if not cur:
                break
        tot += measure(cur)
    return tot

# ----- pinch M (to confirm worry-set membership) -----
def nint(x):
    f = x - (x.numerator // x.denominator); return min(f, 1 - f)
def M_pinch(speeds):
    best = Fr(-1)
    cand = set()
    for a, b in combinations(speeds, 2):
        s = a + b
        for m in range(1, s):
            cand.add(Fr(m, s))
    for t in cand:
        mn = min(nint((v * t) % 1) for v in speeds)
        if mn > best:
            best = mn
    return best

# ===========================================================================
n = 8; C = 2*n - 1; delta = Fr(1, n)
configs = {
    "AP_8        ": [1,2,3,4,5,6,7],          # shell-free
    "(1..5,7,12) ": [1,2,3,4,5,7,12],          # shell-partner (3,12): 3+12=15=C
    "(1,4,5,6,7,11,13)": [1,4,5,6,7,11,13],    # shell-partner (4,11): 4+11=15=C
}
print(f"n={n}, gap δ=1/{n}, C=2n-1={C}.  S_1 (config-free) = (n-1)*2δ = {(n-1)*2*delta}")
poisson_S2 = ( (n-1)*2*delta )**2 / 2
indep_S2   = Fr((n-1)*(n-2), 2) * (2*delta)**2     # Σ μ(D_i)μ(D_j) if independent
print(f"reference S_2: Poisson λ²/2={poisson_S2}={float(poisson_S2):.5f}; "
      f"independent C(7,2)(2δ)²={indep_S2}={float(indep_S2):.5f}")
print("="*86)
rows = []
for name, S in configs.items():
    M = M_pinch(S)
    p, Ds = depth_profile(S, n)
    p0 = p.get(0, Fr(0))
    S1 = Sj(Ds, 1); S2 = Sj(Ds, 2)
    sp = [(a,b) for a,b in combinations(S,2) if (a+b) % C == 0]
    # per-pair overlap, find max and the shell-partner pair's overlap
    pair_ov = {}
    for ia, ib in combinations(range(len(S)), 2):
        pair_ov[(S[ia],S[ib])] = measure(intersect(Ds[ia], Ds[ib]))
    maxpair = max(pair_ov.items(), key=lambda kv: kv[1])
    excess = S2 - indep_S2
    print(f"\n{name.strip()}: speeds={S}")
    print(f"   M={M} (=1/{n}? {M==delta});  shell-partners mod {C}: {sp if sp else 'NONE'}")
    print(f"   p_0(lonely measure) = {p0} (=0 worry-set? {p0==0})")
    print(f"   S_1={S1} (=7/4? {S1==Fr(7,4)});  S_2={S2}={float(S2):.6f}")
    print(f"   S_2 EXCESS over independent = {excess} = {float(excess):.6f}")
    print(f"   max pairwise overlap: pair {maxpair[0]} -> {maxpair[1]}={float(maxpair[1]):.5f}")
    if sp:
        for (a,b) in sp:
            print(f"   SHELL-PARTNER pair {(a,b)} overlap = {pair_ov[(a,b)]}={float(pair_ov[(a,b)]):.5f}"
                  f"  (vs avg pair {float(S2/Fr(len(list(combinations(S,2))))):.5f})")
    rows.append((name.strip(), float(S2), float(excess), bool(sp)))

print("\n" + "="*86)
print("SUMMARY (HYP-2290 test): does shell-partner => higher S_2?")
for nm, s2, ex, hassp in sorted(rows, key=lambda r: r[1]):
    print(f"   {nm:22s} S_2={s2:.6f}  excess={ex:+.6f}  shell-partner={hassp}")
print()
print("VERDICT: compare AP_8 (shell-free) S_2 vs the two shell-partner configs.")
print("  If shell-partner configs have strictly larger S_2 -> HYP-2290(b) supported.")
print("  If the shell-partner PAIR itself is the max-overlap pair -> HYP-2290(a) supported.")
