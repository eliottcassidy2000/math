#!/usr/bin/env python3
"""
THE DEFINITIVE INVARIANT LATTICE of tournaments, in the invariant/monoid/orbit frame.
kind-pasteur-2026-07-21-S128c142.  Owner: work reframings & conjectures, chase DEFINITIVE results.

Each invariant f is an Sn-orbit function (constant on iso classes). f REFINES g iff
(same f-value => same g-value) across all iso classes: f's orbits are >= g's. This is a partial
order; we compute the FULL Hasse diagram + first-separation n for every pair, exhaustively n<=6.

Invariants (all rigorous, computed from scratch on iso-class reps):
  score  = sorted out-degree sequence            (cut-space / hierarchy)
  specA  = char poly of A (adjacency spectrum)    (GL-conjugation invariant, THM-1945 (5))
  specS  = char poly of S=A-A^T (skew spectrum)   (Seidel-side; carries var(lambda^2))
  cyc    = (c3,...,cn) simple directed cycle counts
  H      = # Hamiltonian paths  (Redei; #P; finite-stabilizer per THM-1945 (5))
  R      = signed Ham-path count sum sgn(pi)  (mac-mini THM-1936)
  disc   = |det(I+K)|/2^{n-1}, K=A-A^T  (klein THM-1950; poly-time skew-determinant)
  arb    = # spanning arborescences rooted at 0 (Matrix-Tree; poly-time)
  aut    = |Aut(T)|  (stabilizer order)

Definitive questions answered:
  (Q1) the complete refinement Hasse diagram + first-separation n.
  (Q2) MEET = iso?  which minimal invariant-sets jointly determine the iso class.
  (Q3) mac-mini's open Q: does R distinguish cospectral / same-H tournaments? (is R finer there?)
  (Q4) where do the poly-time invariants (disc, arb, specA) sit vs the #P ones (H, R)?
"""
import itertools
from fractions import Fraction
from collections import defaultdict
import numpy as np

# ---------- bit-packed iso census ----------
def setup(n):
    pairs = list(itertools.combinations(range(n), 2)); pos = {pr: k for k, pr in enumerate(pairs)}
    maps = []
    for p in itertools.permutations(range(n)):
        m = []
        for (i, j) in pairs:
            a, b = p[i], p[j]
            m.append((pos[(a, b)], 0) if a < b else (pos[(b, a)], 1))
        maps.append(m)
    return pairs, maps

def canon_and_orbit(bits, maps, npairs):
    seen = set()
    for m in maps:
        v = 0
        for k in range(npairs):
            bit = (bits >> k) & 1
            tb, inv = m[k]
            if inv: bit ^= 1
            v |= bit << tb
        seen.add(v)
    return min(seen), len(seen)   # canonical rep, orbit size

def adj(bits, pairs, n):
    A = np.zeros((n, n), dtype=int)
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1: A[i, j] = 1
        else: A[j, i] = 1
    return A

# ---------- invariants ----------
def charpoly_int(M):
    n = M.shape[0]; c = [1]; Mk = np.eye(n, dtype=object)
    for k in range(1, n + 1):
        Mk = M.dot(Mk); ck = Fraction(-1, k) * int(np.trace(Mk))
        c.append(int(ck)); Mk = Mk + ck * np.eye(n, dtype=object)
    return tuple(c)

def score_seq(A): return tuple(sorted(int(A[i].sum()) for i in range(A.shape[0])))

def cycle_vector(A):
    n = A.shape[0]; cyc = []
    for k in range(3, n + 1):
        cnt = 0
        for S in itertools.combinations(range(n), k):
            s0 = S[0]; rest = S[1:]
            for perm in itertools.permutations(rest):
                seq = (s0,) + perm
                if all(A[seq[i], seq[(i + 1) % k]] for i in range(k)):
                    cnt += 1
        cyc.append(cnt)
    return tuple(cyc)

def ham_paths_count(A):
    n = A.shape[0]
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            if not (mask >> v) & 1 or dp[mask][v] == 0: continue
            for w in range(n):
                if (mask >> w) & 1 or not A[v, w]: continue
                dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def ham_paths_signed(A):
    """R = sum over Ham paths of sgn(pi); incremental inversion parity (mac-mini)."""
    n = A.shape[0]
    dp = [defaultdict(int) for _ in range(1 << n)]   # dp[mask][last] = signed count
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1 << n):
        if not dp[mask]: continue
        for v, val in list(dp[mask].items()):
            if val == 0: continue
            for w in range(n):
                if (mask >> w) & 1 or not A[v, w]: continue
                # appending w: inversions added = #{used vertices > w}
                inv = bin(mask >> (w + 1)).count("1")
                sign = -1 if inv & 1 else 1
                dp[mask | (1 << w)][w] += val * sign
    full = (1 << n) - 1
    return sum(dp[full].values())

def disc(A):
    n = A.shape[0]; K = A - A.T
    d = round(np.linalg.det(np.eye(n) + K))
    return Fraction(int(d), 2 ** (n - 1))

def arborescences(A):
    """# spanning arborescences rooted at 0 (edges away from root): Matrix-Tree on L=Dout - A^T."""
    n = A.shape[0]
    Aout = A.astype(float)
    L = np.diag(Aout.sum(axis=0)) - Aout.T     # in-degree Laplacian for out-trees
    M = np.delete(np.delete(L, 0, 0), 0, 1)
    return round(np.linalg.det(M))

# ---------- census + invariants per iso class ----------
def census(n):
    pairs, maps = setup(n); npairs = len(pairs)
    reps = {}
    for bits in range(2 ** npairs):
        c, osz = canon_and_orbit(bits, maps, npairs)
        if c not in reps: reps[c] = (bits, osz)
    data = {}
    import math
    for c, (bits, osz) in reps.items():
        A = adj(bits, pairs, n)
        data[c] = dict(
            score=score_seq(A), specA=charpoly_int(A), specS=charpoly_int(A - A.T),
            cyc=cycle_vector(A), H=ham_paths_count(A), R=abs(ham_paths_signed(A)),
            disc=disc(A), arb=arborescences(A), aut=math.factorial(n) // osz)
    return data

INVS = ["score", "specA", "specS", "cyc", "H", "R", "disc", "arb", "aut"]

def refines(data, f, g):
    """f refines g: same f-value => same g-value. Returns True/False."""
    buckets = defaultdict(set)
    for d in data.values():
        buckets[d[f]].add(d[g])
    return all(len(s) == 1 for s in buckets.values())

def n_classes(data, f):
    return len(set(d[f] for d in data.values()))

# ================= run =================
print("=" * 78)
print("DEFINITIVE INVARIANT LATTICE — exhaustive iso census n=3..6")
print("=" * 78)
CEN = {}
for n in range(3, 7):
    CEN[n] = census(n)
    print(f"n={n}: iso classes = {len(CEN[n])}   "
          + "  ".join(f"|{f}|={n_classes(CEN[n], f)}" for f in INVS))

print()
print("=" * 78)
print("(Q1) REFINEMENT MATRIX at n=6  (row f REFINES col g?  '<' = yes)")
print("=" * 78)
data6 = CEN[6]
hdr = "        " + "".join(f"{g[:5]:>6}" for g in INVS)
print(hdr)
for f in INVS:
    row = f"{f:>7} "
    for g in INVS:
        if f == g: row += "     ."
        else: row += "     <" if refines(data6, f, g) else "     -"
    print(row)

print()
print("=" * 78)
print("(Q1b) FIRST-SEPARATION n for selected pairs (smallest n where f does NOT refine g)")
print("=" * 78)
pairs_test = [("score","specA"),("specA","score"),("score","cyc"),("cyc","score"),
              ("specA","H"),("H","specA"),("score","H"),("H","score"),
              ("specA","specS"),("specS","specA"),("H","R"),("R","H"),
              ("disc","H"),("H","disc"),("disc","specS"),("specS","disc"),
              ("arb","specA"),("specA","arb"),("R","specA"),("specA","R"),
              ("cyc","H"),("H","cyc"),("specA","aut"),("cyc","specS"),("specS","cyc")]
for f, g in pairs_test:
    first = None
    for n in range(3, 7):
        if not refines(CEN[n], f, g): first = n; break
    verdict = f"n={first}" if first else "refines through n=6"
    print(f"  {f:>6} refines {g:<6}? first FAILS at {verdict}")

print()
print("=" * 78)
print("(Q2) DOES A SET OF INVARIANTS JOINTLY DETERMINE THE ISO CLASS?  (meet = iso?)")
print("=" * 78)
def joint_determines_iso(data, fs):
    buckets = defaultdict(int)
    for c, d in data.items():
        buckets[tuple(d[f] for f in fs)] += 1
    # determines iso iff every joint-value bucket has exactly ONE iso class
    return max(buckets.values()) == 1, sum(v - 1 for v in buckets.values() if v > 1)
SETS = [("score",), ("specA",), ("cyc",), ("H",),
        ("score","specA"), ("score","cyc"), ("specA","specS"),
        ("score","specA","cyc"), ("score","specA","specS","cyc"),
        ("score","specA","specS","cyc","H","R","disc","arb","aut")]
for n in range(4, 7):
    print(f" n={n}:")
    for fs in SETS:
        det, extra = joint_determines_iso(CEN[n], fs)
        tag = "DETERMINES iso" if det else f"MISSES ({extra} collisions)"
        print(f"    {{{', '.join(fs)}}} -> {tag}")

print()
print("=" * 78)
print("(Q3) mac-mini's open Q: does R distinguish COSPECTRAL / SAME-H tournaments?")
print("=" * 78)
for n in range(4, 7):
    d = CEN[n]
    # among classes sharing specA, does R vary?  among classes sharing H, does R vary?
    byspec = defaultdict(list); byH = defaultdict(list); byR = defaultdict(list)
    for c, dd in d.items():
        byspec[dd["specA"]].append(dd); byH[dd["H"]].append(dd); byR[dd["R"]].append(dd)
    R_splits_cospectral = any(len(set(x["R"] for x in g)) > 1 for g in byspec.values())
    R_splits_sameH = any(len(set(x["R"] for x in g)) > 1 for g in byH.values())
    H_splits_sameR = any(len(set(x["H"] for x in g)) > 1 for g in byR.values())
    specA_splits_sameR = any(len(set(x["specA"] for x in g)) > 1 for g in byR.values())
    print(f" n={n}: R distinguishes some cospectral pair? {R_splits_cospectral}   "
          f"R distinguishes some same-H pair? {R_splits_sameH}")
    print(f"        H distinguishes some same-R pair? {H_splits_sameR}   "
          f"specA distinguishes some same-R pair? {specA_splits_sameR}  "
          f"(=> R and H/specA are INCOMPARABLE)")

print()
print("=" * 78)
print("(Q4) poly-time (disc, arb, specA) vs #P (H, R): does any poly invariant determine H?")
print("=" * 78)
for n in range(4, 7):
    d = CEN[n]
    for f in ["disc", "arb", "specA", "specS", "score", "cyc"]:
        det, _ = joint_determines_iso(d, (f,)) if False else (None, None)
        r = refines(d, f, "H")
        print(f" n={n}: {f:>6} refines H (determines Redei count)? {r}", end="   ")
        print(f"| H refines {f}? {refines(d,'H',f)}")
    print()

print("DONE.")
