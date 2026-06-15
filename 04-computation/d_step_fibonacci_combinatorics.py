#!/usr/bin/env python3
"""
d-step Fibonacci family a_d(n) = a_d(n-1) + a_d(n-d).
Verify the combinatorial readings:
 (1) strip tilings by 1-tiles and d-tiles
 (2) Pascal-slope-d binomial: sum_k C(n-1-(d-1)k, k) counts tilings with exactly k d-tiles
 (3) independent sets in the (d-1)-th power of a path
 (4) lattice-path / no-(d-1)-adjacent-1s binary strings
Plus the project-flavored question: does the STAIRCASE tiling model
(fixed base path, m=C(n-1,2) flippable tiles) have a 1-and-d structure?
"""
from itertools import combinations, product
from math import comb

# ---------- (0) the recurrence ----------
def a_d(n, d, cache=None):
    # number of tilings of a 1 x n strip by tiles of length 1 and length d
    # convention: a_d(0)=1 (empty), a_d(n)=0 for n<0
    if cache is None: cache = {}
    if n < 0: return 0
    if n == 0: return 1
    if (n,d) in cache: return cache[(n,d)]
    v = a_d(n-1,d,cache) + a_d(n-d,d,cache)
    cache[(n,d)] = v
    return v

# ---------- (1) brute-force strip tilings ----------
def strip_tilings_bruteforce(L, d):
    """Count tilings of length-L strip with tiles of length 1 and length d.
    Return total and the distribution over (# d-tiles)."""
    count = 0
    bydtiles = {}
    def rec(remaining, kd):
        nonlocal count
        if remaining == 0:
            count += 1
            bydtiles[kd] = bydtiles.get(kd,0)+1
            return
        # place a 1-tile
        rec(remaining-1, kd)
        # place a d-tile
        if remaining >= d:
            rec(remaining-d, kd+1)
    rec(L, 0)
    return count, bydtiles

# ---------- (3) independent sets in (d-1)-th power of a path ----------
def path_power_independent_sets(num_vertices, h):
    """P_n^h: vertices 0..n-1, edges between i,j with |i-j|<=h.
    Count independent sets (incl empty). Independent = no two chosen within distance h.
    => chosen indices pairwise differ by > h, i.e. gaps >= h+1."""
    verts = list(range(num_vertices))
    cnt = 0
    for r in range(num_vertices+1):
        for S in combinations(verts, r):
            ok = all(S[i+1]-S[i] > h for i in range(len(S)-1))
            if ok: cnt += 1
    return cnt

# ---------- (4) binary strings no two 1s within distance < d (gap-d) ----------
def binary_strings_gap(L, d):
    """Count binary strings of length L where any two 1s are at least d positions apart
    (i.e., consecutive 1s differ by >= d). For d=2 this is the 'no 11' Fibonacci rule."""
    cnt = 0
    for bits in product([0,1], repeat=L):
        ones = [i for i,b in enumerate(bits) if b==1]
        ok = all(ones[i+1]-ones[i] >= d for i in range(len(ones)-1))
        if ok: cnt += 1
    return cnt

print("="*70)
print("d-STEP FIBONACCI a_d(n)=a_d(n-1)+a_d(n-d): COMBINATORIAL READINGS")
print("="*70)

for d in [2,3,4]:
    print(f"\n--- d = {d} ---")
    print(f"{'L':>3} {'recur a_d(L)':>13} {'strip(BF)':>10} {'binom sum':>10} "
          f"{'indep P_(L)^(d-1)':>18} {'binstr gap-d, len L-? ':>10}")
    for L in range(0, 11):
        rec = a_d(L, d)
        bf, dist = strip_tilings_bruteforce(L, d)
        # (2) binomial-slope-d: sum_k C(L - (d-1)k, k)
        binom = sum(comb(L-(d-1)*k, k) for k in range(0, L//d+1) if L-(d-1)*k >= k)
        # check per-k against distribution
        per_k_ok = all(dist.get(k,0) == comb(L-(d-1)*k, k) for k in range(0, L//d+1))
        # (3) indep sets in (d-1)-th power of path on L+? vertices.
        # a_d(L) = strip tiling of length L. independent sets of P_v^{d-1} = a_d(v + (d-1))? test:
        # Known: indep sets of P_n^{h} satisfy p_n = p_{n-1}+p_{n-h-1}, p_n=F_{n+2} when h=1.
        indep = path_power_independent_sets(L, d-1) if L<=10 else None
        print(f"{L:>3} {rec:>13} {bf:>10} {binom:>10} {indep:>18} "
              f"{'per-k OK' if per_k_ok else 'PER-K FAIL':>10}")

# Establish the exact index shifts between the four objects.
print("\n" + "="*70)
print("EXACT INDEX ALIGNMENT (find the shift that makes them coincide)")
print("="*70)
d = 3
print(f"d={d}")
print("strip length L : tilings = a_d(L)")
print("indep sets of P_v^{d-1} as function of v:")
for v in range(0, 9):
    ind = path_power_independent_sets(v, d-1)
    # match to a_d(L): find L with a_d(L)==ind
    Ls = [L for L in range(0,30) if a_d(L,d)==ind]
    print(f"  v={v}: indep={ind}  matches a_d(L) for L in {Ls[:3]}")

print("\nbinary strings gap>=d of length L (consecutive 1s >= d apart):")
for L in range(0,9):
    b = binary_strings_gap(L, d)
    Ls = [LL for LL in range(0,30) if a_d(LL,d)==b]
    print(f"  L={L}: count={b}  matches a_d(LL) for LL in {Ls[:3]}")

# ---------- The d=2 sanity: classic Fibonacci ----------
print("\n" + "="*70)
print("d=2 SANITY: classic Fibonacci F_{n} = 1,1,2,3,5,8,...")
print("="*70)
print("a_2(L) for L=0..9:", [a_2 := a_d(L,2) for L in range(10)])
print("These are Fibonacci F_{L+1} (with F_1=F_2=1).")

# ---------- (Project) STAIRCASE tiling model has a 1-and-d structure? ----------
print("\n" + "="*70)
print("PROJECT TILING MODEL CHECK")
print("="*70)
print("Staircase delta_{n-2}: m = C(n-1,2) flippable tiles, 2^m tilings total.")
for n in range(3,9):
    m = comb(n-1,2)
    print(f"  n={n}: m=C({n-1},2)={m}, total tilings 2^{m}={2**m}")
print("\nThese are 2^m, NOT a d-step Fibonacci number in general.")
print("So the FULL staircase tile-flip hypercube is NOT a_d(n).")
print("The Fibonacci/d-step appears in CONSTRAINED sub-counts (circular tournaments,")
print("hard-core/independent-set readings), not the unconstrained 2^m hypercube.")
