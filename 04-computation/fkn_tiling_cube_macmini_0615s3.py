#!/usr/bin/env python3
"""
fkn_tiling_cube_macmini_0615s3.py  (mac-mini-2026-06-15-S3, T823/THM-512)

FKN ON THE TILING CUBE.
Tiling space = {0,1}^m, m=C(n-1,2). Base path n->n-1->...->1 (arc k->k-1 fixed).
Tile (x,y), x-y>=2: bit 1 = "aligned" (x->y, transitive direction); bit 0 = reversed.
The ALL-ALIGNED tiling = the TRANSITIVE tournament = the ground state. Backward-arc
count (# reversed tiles) = Hamming distance from the ground state = energy.

Computes:
 (1) the iso-class / tournament map of the 2^m tilings (n=4,5,6), and the distribution
     of iso classes by BACKWARD-ARC Hamming weight (energy);
 (2) Walsh-Fourier spectra (degree-by-level weight) of key invariants on the cube:
     c3 (signed triangle count), H (Ham-path count), and the "is-transitive" indicator;
 (3) single-tile-flip ("wiggly") perturbation: how often a flip changes the iso class,
     by tile and by energy level (the FKN influence / dictator structure);
 (4) the MÖBIUS inclusion-exclusion check: for an additive invariant, codim-c sub-
     tournaments with count C(n-1,c) and sign (-1)^{c+1} (the user's A+B+C-D-E-F+G).
"""
import sys, itertools
from math import comb
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def tiles_of(n):
    # explorer order: for y=1..n-2, for x=n down to y+2: tile (x,y)
    return [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1) if x-y >= 2]

def tiling_to_adj(n, bits, TILES):
    # vertices 1..n. base path arc k->k-1. tile aligned(bit1): x->y; reversed(bit0): y->x
    A = [[0]*(n+1) for _ in range(n+1)]
    for k in range(n, 1, -1):
        A[k][k-1] = 1
    for i, (x, y) in enumerate(TILES):
        if (bits >> i) & 1:
            A[x][y] = 1
        else:
            A[y][x] = 1
    return A

def canon(n, A):
    # canonical form over all relabelings (n<=6 ok)
    best = None
    verts = list(range(1, n+1))
    for perm in itertools.permutations(verts):
        # build bitmask of arcs under this relabeling
        key = 0
        idx = {v: p for p, v in enumerate(perm)}
        b = 0
        for i in range(1, n+1):
            for j in range(1, n+1):
                if A[i][j]:
                    b |= 1 << (idx[i]*n + idx[j])
        if best is None or b < best:
            best = b
    return best

def c3_count(n, A):
    c = 0
    for a, b, d in itertools.combinations(range(1, n+1), 3):
        if (A[a][b] and A[b][d] and A[d][a]) or (A[b][a] and A[d][b] and A[a][d]):
            c += 1
    return c

def H_count(n, A):
    out = [0]*(n+1)
    for i in range(1, n+1):
        r = 0
        for j in range(1, n+1):
            if A[i][j]: r |= 1 << (j-1)
        out[i] = r
    full = (1 << n) - 1
    dp = [[0]*(n+1) for _ in range(1 << n)]
    for v in range(1, n+1): dp[1 << (v-1)][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last in range(1, n+1):
            cc = row[last]
            if cc:
                nx = out[last] & ~mask
                while nx:
                    bb = nx & (-nx); j = bb.bit_length(); dp[mask|bb][j] += cc; nx ^= bb
    return sum(dp[full][1:])

def walsh_levels(vals, m):
    """real Walsh transform on +-1 input cube; return weight per Hamming-level of S."""
    N = 1 << m
    f = [float(v) for v in vals]
    step = 1
    while step < N:
        for base in range(0, N, step*2):
            for i in range(base, base+step):
                a, b = f[i], f[i+step]; f[i], f[i+step] = a+b, a-b
        step *= 2
    f = [c/N for c in f]
    lev = [0.0]*(m+1)
    for s in range(N):
        lev[bin(s).count('1')] += f[s]**2
    return lev

print("="*72)
print("FKN ON THE TILING CUBE: ground state, energy, Walsh spectra, perturbations")
print("="*72)

for n in range(4, 7):
    TILES = tiles_of(n); m = len(TILES); N = 1 << m
    Hv = [0]*N; c3v = [0]*N; cls = [0]*N; bw = [0]*N
    for bits in range(N):
        A = tiling_to_adj(n, bits, TILES)
        Hv[bits] = H_count(n, A)
        c3v[bits] = c3_count(n, A)
        bw[bits] = m - bin(bits).count('1')   # backward arcs = reversed tiles
        cls[bits] = canon(n, A)
    nclasses = len(set(cls))
    print(f"\n--- n={n}: m={m} tiles, {N} tilings, {nclasses} iso classes (= A000568({n})) ---")
    # (1) iso classes by backward-arc energy: for each energy level, # distinct iso classes
    from collections import defaultdict
    cls_by_energy = defaultdict(set); tilings_by_energy = defaultdict(int)
    for bits in range(N):
        cls_by_energy[bw[bits]].add(cls[bits]); tilings_by_energy[bw[bits]] += 1
    print("  energy(backward arcs) -> (#tilings, #distinct iso classes):")
    for e in range(m+1):
        print(f"    e={e}: tilings={tilings_by_energy[e]:4d}, distinct iso={len(cls_by_energy[e])}")
    print(f"  ground state e=0: transitive (c3={c3v[N-1]}, H={Hv[N-1]}); "
          f"c3 range {min(c3v)}..{max(c3v)}, H range {min(Hv)}..{max(Hv)}")
    # (2) Walsh spectra of c3, H, is-transitive indicator (centered to +-1 inputs: tiles +-1)
    levH = walsh_levels(Hv, m)
    levc3 = walsh_levels(c3v, m)
    istrans = [1.0 if c3v[bits]==0 and bits==N-1 else (1.0 if c3v[bits]==0 else 0.0) for bits in range(N)]
    levT = walsh_levels(istrans, m)
    def topdeg(lev, tol=1e-9):
        return max(k for k in range(len(lev)) if lev[k] > tol)
    print(f"  Walsh degree (top level w/ mass): c3={topdeg(levc3)}, H={topdeg(levH)} "
          f"(THM-163 band D=2*floor((n-1)/2)={2*((n-1)//2)}), is-transitive={topdeg(levT)}")
    print(f"  c3 Walsh weight by level: {[round(x,3) for x in levc3]}")
    print(f"  H  Walsh weight by level: {[round(x,2) for x in levH[:8]]}{'...' if m>=8 else ''}")
    # (3) single-flip influence: how often flipping tile i changes the iso class
    inf = [0]*m
    for bits in range(N):
        for i in range(m):
            if cls[bits] != cls[bits ^ (1<<i)]:
                inf[i] += 1
    inf = [x/N for x in inf]
    print(f"  single-flip iso-class-change influence per tile (FKN dictator strength): "
          f"{[round(x,3) for x in inf]}")
    print(f"    total influence (sum) = {sum(inf):.3f}; apex tile (n,1)=index0 influence={inf[0]:.3f}")

print("\n" + "="*72)
print("MÖBIUS INCLUSION-EXCLUSION (the A+B+C-D-E-F+G structure)")
print("="*72)
print("codim-c sub-tournament count C(n-1,c), sign (-1)^(c+1):")
for n in range(4, 8):
    terms = [(c, comb(n-1, c), (-1)**(c+1)) for c in range(1, n)]
    desc = " ".join(f"{'+' if s>0 else '-'}{cnt}*sz{n-c}" for c, cnt, s in terms)
    alt = sum(s*cnt for c, cnt, s in terms)
    print(f"  n={n}: {desc}   (alt-sum of counts = {alt}; sum |counts| = {sum(cnt for _,cnt,_ in terms)} = 2^(n-1)-1)")
print("  Note: Σ_c C(n-1,c) = 2^(n-1)-1 (all nonempty subsets of the n-1 non-anchor vertices);")
print("  Σ_c (-1)^(c+1) C(n-1,c) = 1 - (1-1)^(n-1) = 1 (the Möbius/inclusion-exclusion sieve).")
print("\nDONE.")
