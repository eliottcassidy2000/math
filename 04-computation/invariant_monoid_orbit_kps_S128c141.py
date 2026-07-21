#!/usr/bin/env python3
"""
invariant / monoid / orbit — verify the checkable creative statements.
kind-pasteur-2026-07-21-S128c141.  Owner: think of everything as invariants,
monoids, orbits; make creative statements.  (bit-packed, fast to n=6.)

  (A) switching-class count = cut-space orbits = 2^{C(n-1,2)}; unlabeled = A002854 = V(E_n)
      (repo even-graph metagraph = classical Seidel two-graphs).
  (B) iso classes = S_n-orbits (A000568).
  (C) score-vector classes; Ryser: directed-3-cycle-reversal orbits = out-degree-vector classes.
  (D) invariant lattice: spectrum vs score refinement; first COspectral-non-iso n.
  (E) permanent vs determinant: finite vs positive-dimensional stabilizer (GCT hardness split).
"""
import itertools, sys
from fractions import Fraction
from collections import Counter, defaultdict
import numpy as np

# ---------- bit-packed tournaments ----------
def setup(n):
    pairs = list(itertools.combinations(range(n), 2))
    pos = {pr: k for k, pr in enumerate(pairs)}
    # for each perm, precompute source-bit -> (target-bit, invert)
    perms = list(itertools.permutations(range(n)))
    maps = []
    for p in perms:
        m = []
        for (i, j) in pairs:
            a, b = p[i], p[j]
            if a < b:
                m.append((pos[(a, b)], 0))
            else:
                m.append((pos[(b, a)], 1))
        maps.append(m)
    return pairs, pos, maps

def canon_bits(bits, maps, npairs):
    best = None
    for m in maps:
        v = 0
        for k in range(npairs):
            bit = (bits >> k) & 1
            tb, inv = m[k]
            if inv:
                bit ^= 1
            v |= bit << tb
        if best is None or v < best:
            best = v
    return best

def adj(bits, pairs, n):
    A = np.zeros((n, n), dtype=int)
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1:
            A[i, j] = 1
        else:
            A[j, i] = 1
    return A

def charpoly_int(A):
    n = A.shape[0]
    c = [1]; Mk = np.eye(n, dtype=object)
    for k in range(1, n + 1):
        Mk = A.dot(Mk)
        ck = Fraction(-1, k) * int(np.trace(Mk))
        c.append(int(ck)); Mk = Mk + ck * np.eye(n, dtype=object)
    return tuple(c)

def score_vec(A):
    return tuple(int(A[i].sum()) for i in range(A.shape[0]))

print("=" * 72)
print("(A) switching classes = cut-space orbits = 2^{C(n-1,2)}  vs  A002854 = V(E_n)")
print("=" * 72)
A002854 = {3: 2, 4: 3, 5: 7, 6: 16, 7: 54, 8: 243, 9: 2038}
for n in range(3, 9):
    Cn2 = n * (n - 1) // 2; Cn12 = (n - 1) * (n - 2) // 2
    cut_orbits = 2 ** (Cn2 - (n - 1))
    print(f"n={n}: 2^C(n,2)={2**Cn2:>8} / 2^(n-1) = cut-orbits {cut_orbits:>7} "
          f"= 2^C(n-1,2) {2**Cn12:>7}  match={cut_orbits==2**Cn12}   "
          f"unlabeled two-graphs (A002854=V(E_n)) = {A002854.get(n,'?')}")

print()
print("=" * 72)
print("(B)-(D) exhaustive census (bit-packed iso-canon)")
print("=" * 72)
A000568 = {3: 2, 4: 4, 5: 12, 6: 56}
first_cospectral = None
for n in range(3, 7):
    pairs, pos, maps = setup(n); npairs = len(pairs)
    iso = {}                                   # canon-bits -> a representative bits
    for bits in range(2 ** npairs):
        c = canon_bits(bits, maps, npairs)
        if c not in iso:
            iso[c] = bits
    scores, specs = {}, {}
    for c, bits in iso.items():
        A = adj(bits, pairs, n)
        scores[c] = tuple(sorted(score_vec(A)))
        specs[c] = charpoly_int(A)
    n_iso = len(iso); n_score = len(set(scores.values())); n_spec = len(set(specs.values()))
    spec_ct = Counter(specs.values())
    cospectral = sum(v - 1 for v in spec_ct.values() if v > 1)
    if cospectral and first_cospectral is None:
        first_cospectral = n
    spec_to_scores = defaultdict(set); score_to_specs = defaultdict(set)
    for c in iso:
        spec_to_scores[specs[c]].add(scores[c]); score_to_specs[scores[c]].add(specs[c])
    spec_refines_score = all(len(s) == 1 for s in spec_to_scores.values())
    score_refines_spec = all(len(s) == 1 for s in score_to_specs.values())
    tag = 'ok' if n_iso == A000568[n] else 'BAD'
    print(f"n={n}: iso(S_n-orbits)={n_iso}(A000568={A000568[n]} {tag})  "
          f"score-classes={n_score}  spectral-classes={n_spec}  "
          f"cospectral-non-iso-pairs={cospectral}")
    print(f"      spectrum refines score? {spec_refines_score}   "
          f"score refines spectrum? {score_refines_spec}")
print(f"  => first n with COspectral non-isomorphic tournaments: n = {first_cospectral} "
      f"(spectral map = STRICT quotient of the iso-class quotient)")

print()
print("=" * 72)
print("(C) Ryser: directed-3-cycle-reversal orbits == out-degree-VECTOR classes")
print("=" * 72)
def directed_3cycles(A):
    n = A.shape[0]; out = []
    for i, j, k in itertools.combinations(range(n), 3):
        for (a, b, c) in [(i, j, k), (i, k, j)]:
            if A[a, b] and A[b, c] and A[c, a]:
                out.append((a, b, c))
    return out
for n in range(4, 6):
    pairs, pos, maps = setup(n); npairs = len(pairs)
    N = 2 ** npairs
    Amats = [adj(b, pairs, n) for b in range(N)]
    def encode(A):
        v = 0
        for k, (i, j) in enumerate(pairs):
            if A[i, j]: v |= 1 << k
        return v
    parent = list(range(N))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for i in range(N):
        A = Amats[i]
        for cyc in directed_3cycles(A):
            B = A.copy(); a, b, c = cyc
            for (x, y) in [(a, b), (b, c), (c, a)]:
                B[x, y] = 0; B[y, x] = 1
            union(i, encode(B))
    comps = defaultdict(list)
    for i in range(N): comps[find(i)].append(i)
    svclass = defaultdict(list)
    for i in range(N): svclass[score_vec(Amats[i])].append(i)
    comp_of = {i: find(i) for i in range(N)}
    ok_const = all(len({score_vec(Amats[i]) for i in mem}) == 1 for mem in comps.values())
    ok_single = all(len({comp_of[i] for i in mem}) == 1 for mem in svclass.values())
    print(f"n={n}: #T={N}  #reversal-orbits={len(comps)}  #score-vector-classes={len(svclass)}  "
          f"orbits==score-classes? {ok_const and ok_single and len(comps)==len(svclass)}")

print()
print("=" * 72)
print("(E) permanent vs determinant — symmetry-group hardness split")
print("=" * 72)
rng = np.random.default_rng(12345); n = 5
A = rng.integers(0, 3, (n, n))
def permanent(M):
    tot = 0
    for p in itertools.permutations(range(M.shape[0])):
        pr = 1
        for i in range(M.shape[0]): pr *= M[i, p[i]]
        tot += pr
    return tot
def pm(p):
    P = np.zeros((len(p), len(p)), dtype=int)
    for i, pi in enumerate(p): P[i, pi] = 1
    return P
p = list(rng.permutation(n)); q = list(rng.permutation(n))
per0, per1 = permanent(A), permanent(pm(p).dot(A).dot(pm(q)))
g = rng.integers(-2, 3, (n, n)).astype(float)
while abs(np.linalg.det(g)) < 1e-6:
    g = rng.integers(-2, 3, (n, n)).astype(float)
det0 = np.linalg.det(A.astype(float))
det1 = np.linalg.det(g.dot(A.astype(float)).dot(np.linalg.inv(g)))
print(f"per(A)={per0}  per(PAQ)={per1}  equal={per0==per1}  (per stabilizer S_n x S_n, FINITE)")
print(f"det(A)={det0:.3f}  det(gAg^-1)={det1:.3f}  equal={abs(det0-det1)<1e-6}  "
      f"(det stabilizer contains GL-conjugation, POSITIVE-DIMENSIONAL)")
print("=> tractability tracks stabilizer DIMENSION: continuous (det/trace, P) vs finite (per/H, #P);")
print("   Mulmuley-Sohoni GCT split; THM-1780 'H leaves the ladder at n=6' is its tournament instance.")
