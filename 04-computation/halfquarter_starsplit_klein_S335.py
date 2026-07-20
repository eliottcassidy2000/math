#!/usr/bin/env python3
"""
klein-2026-07-20-S335 -- HOW THE HALF AND QUARTER FOLDS SIT INSIDE THE STAR (CUT/CYCLE) SPLIT.

Prior art this builds on (does NOT duplicate):
  THM-549/550  half tiling = fold by the grid reflection sigma(x,y)=(n+1-y,n+1-x);
               Fix(sigma) = {x+y=n+1}, |Fix| = floor((n-1)/2); h = floor((n-1)^2/4).
               sigma(tiling) = phi_relabel(T^op): the grid reflection IS the complement.
  kps S15/S16  quarter tiling = Q_m / <sigma, phi>, Klein 4-group, phi = antipodal
               (flip ALL tiles).  MISTAKE-033: phi != T^op.
  THM-1382 (klein-S333) / THM-1405 (mac-mini-S128): star flips span cut(H),
               invariants = cycle(H), H = K_n \ P; m-n+1 "holonomy" bits.

NEW QUESTION ANSWERED HERE: sigma and phi are the two generators of the quarter fold.
Where do they sit relative to cut(H) (+) cycle(H)?

Claims tested, n = 3..11:
 (A) sigma is the tile permutation induced by the VERTEX involution rho(v)=n+1-v,
     which is a graph automorphism of H.  Hence sigma(star(v)) = star(rho(v)):
     sigma NORMALISES the star group, so it preserves cut(H) and cycle(H) setwise.
     => the HALF fold is INTERNAL to the star algebra.
 (B) phi = all-ones vector 1.  1 in cut(H)  <=>  H bipartite  <=>  n <= 4.
     For n >= 5 the triangle {1,3,5} sits inside H and obstructs.
     => the QUARTER fold is TRANSVERSE to the star algebra.
 (C) For n >= 5, phi acts on the star-orbit space (= the holonomy word) as
     TRANSLATION by the nonzero functional eps(C) = |C| mod 2 (odd-cycle indicator).
     Translation by a nonzero vector is FREE => phi pairs star orbits with no fixed
     orbit, so the quarter fold costs EXACTLY ONE holonomy bit.
     phi-invariant conserved quantities = EVEN-cycle parities, dim = m-(n-1)-1.
 (D) The triangle's three sides in star language: the two LEGS are star_H(1) and
     star_H(n) (both are star-flip generators, i.e. cut side); the HYPOTENUSE is
     Fix(sigma), the rho-mirror matching; and sigma SWAPS the two legs.
 (E) Quantitative defect of phi from the star group: m - maxcut(H).
"""
import itertools, sys
from math import comb

# ---------- GF(2) linear algebra on ints (bitmasks over m tile coordinates) ----------
def rref(vecs):
    basis = []
    for v in vecs:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v); basis.sort(reverse=True)
    return basis

def rank(vecs):
    return len(rref(vecs))

def in_span(v, basis):
    for b in basis:
        v = min(v, v ^ b)
    return v == 0

# ---------- the staircase model ----------
def tiles_of(n):
    """Explorer enumeration order: for y=1..n-2: for x=n down to y+2: tile (x,y)."""
    T = []
    for y in range(1, n - 1):
        for x in range(n, y + 1, -1):
            T.append((x, y))
    return T

def build(n):
    T = tiles_of(n)
    m = len(T)
    assert m == comb(n - 1, 2)
    idx = {t: i for i, t in enumerate(T)}
    # H = K_n minus the base path; its edges ARE the tiles
    star = {}
    for v in range(1, n + 1):
        s = 0
        for (x, y) in T:
            if x == v or y == v:
                s |= 1 << idx[(x, y)]
        star[v] = s
    return T, m, idx, star

def sigma_perm(n, T, idx):
    """grid map (x,y) -> (n+1-y, n+1-x); as a permutation of tile indices."""
    p = [0] * len(T)
    for (x, y) in T:
        u, w = n + 1 - y, n + 1 - x       # u > w since x > y
        p[idx[(x, y)]] = idx[(u, w)]
    return p

def apply_perm(vec, p):
    out = 0
    for i in range(len(p)):
        if vec >> i & 1:
            out |= 1 << p[i]
    return out

def cycle_basis(n, T, idx):
    """fundamental cycles of H w.r.t. a spanning tree, as tile-index bitmasks + lengths."""
    parent, order, seen = {}, [], set()
    adj = {v: [] for v in range(1, n + 1)}
    for (x, y) in T:
        adj[x].append(y); adj[y].append(x)
    roots = []
    for s in range(1, n + 1):
        if s in seen: continue
        roots.append(s); seen.add(s); parent[s] = None; stack = [s]
        while stack:
            v = stack.pop(); order.append(v)
            for w in adj[v]:
                if w not in seen:
                    seen.add(w); parent[w] = v; stack.append(w)
    tree_edges = set()
    for v in parent:
        if parent[v] is not None:
            tree_edges.add(frozenset((v, parent[v])))
    def path_to_root(v):
        out = []
        while parent[v] is not None:
            out.append(frozenset((v, parent[v]))); v = parent[v]
        return out, v
    def emask(e):
        a, b = sorted(e, reverse=True)
        return 1 << idx[(a, b)]
    cyc = []
    for (x, y) in T:
        e = frozenset((x, y))
        if e in tree_edges: continue
        px, rx = path_to_root(x); py, ry = path_to_root(y)
        if rx != ry: continue
        sym = set(px) ^ set(py)
        mask = emask(e)
        for f in sym:
            mask ^= emask(f)
        cyc.append(mask)
    return cyc, len(roots)

def maxcut(n, T):
    best = 0
    for S in range(1 << (n - 1)):            # fix vertex n on one side
        bits = S
        side = [False] * (n + 1)
        for v in range(1, n):
            side[v] = bool(bits >> (v - 1) & 1)
        side[n] = False
        c = sum(1 for (x, y) in T if side[x] != side[y])
        best = max(best, c)
    return best

def is_bipartite(n, T):
    col = {}
    adj = {v: [] for v in range(1, n + 1)}
    for (x, y) in T:
        adj[x].append(y); adj[y].append(x)
    for s in range(1, n + 1):
        if s in col: continue
        col[s] = 0; st = [s]
        while st:
            v = st.pop()
            for w in adj[v]:
                if w not in col:
                    col[w] = col[v] ^ 1; st.append(w)
                elif col[w] == col[v]:
                    return False
    return True

# ---------------------------------- run ----------------------------------
print("=" * 78)
print("klein-S335: the half fold (sigma) and quarter fold (sigma,phi) vs the star split")
print("=" * 78)

rows = []
for n in range(3, 12):
    T, m, idx, star = build(n)
    p = sigma_perm(n, T, idx)
    cutb = rref([star[v] for v in range(1, n + 1)])
    dcut = len(cutb)
    cyc, ncomp = cycle_basis(n, T, idx)
    dcyc = m - dcut
    ones = (1 << m) - 1

    # (A) sigma = permutation induced by rho(v) = n+1-v, and it permutes the stars
    rho_ok = all(apply_perm(star[v], p) == star[n + 1 - v] for v in range(1, n + 1))
    # sigma preserves cut(H) setwise
    cut_stable = all(in_span(apply_perm(b, p), cutb) for b in cutb)

    # (B) is phi = all-ones inside cut(H)?
    phi_in_cut = in_span(ones, cutb)
    bip = is_bipartite(n, T)

    # (C) eps(C) = |C| mod 2 on the cycle basis; phi acts by translation by eps
    eps_bits = [bin(c).count("1") & 1 for c in cyc]
    eps_nonzero = any(eps_bits)
    # verify the translation law on random tilings: hol(t+1) = hol(t) + eps
    import random
    random.seed(1000 + n)
    trans_ok = True
    for _ in range(2000):
        t = random.getrandbits(m)
        h0 = [bin(t & c).count("1") & 1 for c in cyc]
        h1 = [bin((t ^ ones) & c).count("1") & 1 for c in cyc]
        if [(a ^ b) for a, b in zip(h0, h1)] != eps_bits:
            trans_ok = False; break
    # freeness of phi on star orbits <=> eps != 0
    free = eps_nonzero

    # (D) legs and hypotenuse
    legA = star[1] & ~0; legB = star[n]
    leg1 = sum(1 for (x, y) in T if y == 1)
    legn = sum(1 for (x, y) in T if x == n)
    legs_swapped = (apply_perm(star[1], p) == star[n])
    fixtiles = [i for i in range(m) if p[i] == i]
    hyp = [T[i] for i in fixtiles]
    hyp_is_antidiag = all(x + y == n + 1 for (x, y) in hyp) and len(hyp) == (n - 1) // 2
    # is the hypotenuse a matching?
    seenv = set(); matching = True
    for (x, y) in hyp:
        if x in seenv or y in seenv: matching = False
        seenv.add(x); seenv.add(y)

    # (E) defect
    mc = maxcut(n, T) if n <= 11 else None
    defect = m - mc

    # invariant dimensions
    dim_half_note = dcyc            # sigma permutes; does not drop dimension
    dim_quarter = dcyc - (1 if eps_nonzero else 0)

    rows.append((n, m, dcut, dcyc, phi_in_cut, bip, eps_nonzero, dim_quarter, defect))
    print(f"\n--- n={n}  m={m} ---")
    print(f" (A) sigma = rho-induced tile perm, sigma(star v)=star(n+1-v): {rho_ok};  cut(H) sigma-stable: {cut_stable}")
    print(f" (B) phi=all-ones in cut(H): {phi_in_cut}   H bipartite: {bip}   [must agree]")
    print(f" (C) eps != 0: {eps_nonzero}   translation law hol(t+1)=hol(t)+eps: {trans_ok}   phi free on orbits: {free}")
    print(f" (D) leg(y=1)={leg1} tiles = star_H(1); leg(x=n)={legn} tiles = star_H(n); sigma swaps legs: {legs_swapped}")
    print(f"     Fix(sigma) = {hyp}  (antidiagonal x+y=n+1: {hyp_is_antidiag}, a matching: {matching})")
    print(f" (E) dim cut={dcut} (=n-1: {dcut==n-1})  dim cycle={dcyc} (=m-n+1: {dcyc==m-n+1})")
    print(f"     star orbits = 2^{dcyc};  QUARTER-invariant dim = {dim_quarter} = m-(n-1)-1")
    print(f"     maxcut(H)={mc}, defect of phi from cut(H) = m-maxcut = {defect}")

print("\n" + "=" * 78)
print("SUMMARY TABLE")
print("=" * 78)
print(f"{'n':>3} {'m':>4} {'dim cut':>8} {'dim cyc':>8} {'phi in cut':>11} {'H bip':>7} {'eps!=0':>7} {'quarter dim':>12} {'defect':>7}")
for (n, m, dcut, dcyc, pic, bip, eps, dq, df) in rows:
    print(f"{n:>3} {m:>4} {dcut:>8} {dcyc:>8} {str(pic):>11} {str(bip):>7} {str(eps):>7} {dq:>12} {df:>7}")

print("\nREADING:")
print(" half fold  (sigma): a PERMUTATION of the holonomy word -- internal to the star algebra.")
print(" quarter fold (phi): a TRANSLATION of the holonomy word by eps -- transverse, free for n>=5,")
print("                     costing exactly ONE holonomy bit (the odd-cycle-size functional).")
