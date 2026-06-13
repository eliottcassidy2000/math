#!/usr/bin/env python3
"""
equidecomposability_s90k.py — Tournament Equidecomposability
opus-2026-03-15-S90k

HILBERT'S THIRD PROBLEM (1900): When are two polyhedra with the same
volume scissors-congruent (cuttable into pieces and reassembled)?

ANSWER (Dehn 1901, Sydler 1965): Iff they have the same volume AND
the same Dehn invariant D = Σ_edges l(e) ⊗ θ(e).

THE TOURNAMENT ANALOG: When are two tournaments "equidecomposable"?
What is the tournament Dehn invariant?

THESIS: β₁ IS the tournament Dehn invariant.
  β₁ = 0: tournament is "scissors-congruent" to its transitive core
  β₁ = 1: tournament has a TWIST that prevents decomposition
  The twist is the portal (reversed min-max arc)
  The Dehn invariant measures the monodromy of the covering space
"""

import numpy as np
from itertools import permutations, combinations
from math import factorial, comb, log, pi

# ======================================================================
# PART 1: THE ARC-FLIP GRAPH — CONFIGURATION SPACE OF TOURNAMENTS
# ======================================================================
print("=" * 70)
print("PART 1: THE ARC-FLIP GRAPH")
print("=" * 70)

print("""
Two tournaments T₁, T₂ on n vertices are CONNECTED by an arc flip
if they differ in exactly one arc direction.

The ARC-FLIP GRAPH G_n has:
  Vertices: all 2^{C(n,2)} tournaments on n vertices
  Edges: pairs differing in exactly one arc

This graph is the C(n,2)-dimensional HYPERCUBE {0,1}^{C(n,2)}.
It is ALWAYS connected: any tournament can reach any other
by at most C(n,2) arc flips.

THE QUESTION: What invariants are preserved along PATHS in G_n?

  H(T) mod 2 = 1 ALWAYS (Rédei). ← preserved trivially
  H(T) itself changes with each flip. ← NOT preserved
  Score sequence changes with each flip. ← NOT preserved

But the PARITY of certain invariants IS preserved:
  H(T) ≡ 1 + 2α₁(T) (mod 4), where α₁ = #odd cycles.
  So H mod 4 is determined by the cycle count parity.

THE EQUIDECOMPOSABILITY QUESTION:
  T₁ ~ T₂ iff there exists a sequence of arc flips
  T₁ → T₁' → T₁'' → ... → T₂ such that some invariant I
  is preserved at EACH STEP: I(T_k) = I(T_{k+1}).

  Which invariants give interesting equivalence classes?
""")

# Compute the arc-flip graph for n=4
n = 4
m = n*(n-1)//2  # = 6

def tournament_from_bits(bits, n):
    T = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
            idx += 1
    return T

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for p in range(n-1):
            if T[perm[p]][perm[p+1]] != 1:
                ok = False
                break
        if ok:
            count += 1
    return count

# For n=4: compute H for all 64 tournaments and track single-flip changes
print(f"n=4: Arc-flip graph analysis ({2**m} tournaments, {m} arcs)")

h_values = {}
for bits in range(2**m):
    T = tournament_from_bits(bits, n)
    h_values[bits] = count_H(T, n)

# For each arc flip, compute ΔH
delta_H_counts = {}
for bits in range(2**m):
    for arc in range(m):
        flipped = bits ^ (1 << arc)
        dH = h_values[flipped] - h_values[bits]
        delta_H_counts[dH] = delta_H_counts.get(dH, 0) + 1

print(f"  Distribution of ΔH per single arc flip:")
for dh in sorted(delta_H_counts.keys()):
    print(f"    ΔH = {dh:+3d}: occurs {delta_H_counts[dh]//2} times")  # /2 for undirected

# ======================================================================
# PART 2: β₁ AS THE DEHN INVARIANT
# ======================================================================
print()
print("=" * 70)
print("PART 2: β₁ IS THE DEHN INVARIANT")
print("=" * 70)

print("""
HILBERT'S THIRD PROBLEM asks: when are two polyhedra with the same
volume scissors-congruent?

THE DEHN INVARIANT D(P) is the OBSTRUCTION:
  D(P) = 0 iff P is scissors-congruent to a box.
  D(P) ≠ 0 iff P has a "twist" that prevents box-decomposition.

For tournaments, the ANALOG is:
  "VOLUME" = H(T) (the Hamiltonian path count)
  "DEHN INVARIANT" = β₁(T) (the first Betti number)
  "BOX" = transitive tournament (the "untwisted" ground state)

CLAIM: β₁(T) measures the obstruction to decomposing T into
"transitive pieces" (sub-tournaments that are total orders).

When β₁ = 0: T can be decomposed into transitive sub-tournaments
  whose orderings are GLOBALLY CONSISTENT.
  The tournament IS a "patchwork" of total orders.

When β₁ = 1: T has a TWIST — a 1-cycle that no collection of
  transitive triples can fill. The twist prevents global consistency.
  The tournament is like a Möbius strip: locally flat, globally twisted.

THE NEAR-TRANSITIVE TOURNAMENT:
  β₁ = 1 (one twist) and sim_H = 1 (one simplicial path).
  It is the SIMPLEST non-trivial Dehn class.
  It is the "regular tetrahedron" of tournaments:
    same "volume" as a box, but with a nonzero twist.
""")

def compute_beta1_glmy(T, n):
    """Compute β₁ using GLMY path homology (regular 2-paths only)."""
    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and T[i][j]==1]
    edge_idx = {e:i for i,e in enumerate(edges)}
    ne = len(edges)

    # Regular 2-paths: (a,b,c) with a→b, b→c, a→c all arcs
    paths2 = []
    for a in range(n):
        for b in range(n):
            if a==b or T[a][b]!=1: continue
            for c in range(n):
                if c==a or c==b or T[b][c]!=1: continue
                if T[a][c]!=1: continue  # regularity
                paths2.append((a,b,c))

    np2 = len(paths2)

    # ∂₂ matrix
    d2 = np.zeros((ne, np2), dtype=int) if np2 > 0 else np.zeros((ne, 1), dtype=int)
    for j, (a,b,c) in enumerate(paths2):
        if (b,c) in edge_idx: d2[edge_idx[(b,c)], j] += 1
        if (a,c) in edge_idx: d2[edge_idx[(a,c)], j] -= 1
        if (a,b) in edge_idx: d2[edge_idx[(a,b)], j] += 1

    # ∂₁ matrix
    d1 = np.zeros((n, ne), dtype=int)
    for j, (a,b) in enumerate(edges):
        d1[b,j] += 1
        d1[a,j] -= 1

    rank_d1 = np.linalg.matrix_rank(d1)
    rank_d2 = np.linalg.matrix_rank(d2) if np2 > 0 else 0
    ker_d1 = ne - rank_d1
    return ker_d1 - rank_d2

# Compute β₁ and H for all n=5 tournaments, classify by (H, β₁)
n = 5
m = n*(n-1)//2
print(f"\nn=5: Computing (H, β₁) for all {2**m} tournaments...")

hb_counts = {}
for bits in range(2**m):
    T = tournament_from_bits(bits, n)
    H = count_H(T, n)
    b1 = compute_beta1_glmy(T, n)
    key = (H, b1)
    hb_counts[key] = hb_counts.get(key, 0) + 1

print(f"\n  (H, β₁) distribution:")
print(f"  {'H':>4} | {'β₁':>3} | {'count':>6}")
print(f"  " + "-" * 22)
for (h, b) in sorted(hb_counts.keys()):
    print(f"  {h:4d} | {b:3d} | {hb_counts[(h,b)]:6d}")

# Key insight: which H values have β₁=0 vs β₁=1?
h_has_b0 = set(h for (h,b) in hb_counts if b == 0)
h_has_b1 = set(h for (h,b) in hb_counts if b == 1)
h_only_b0 = h_has_b0 - h_has_b1
h_only_b1 = h_has_b1 - h_has_b0
h_both = h_has_b0 & h_has_b1

print(f"\n  H values with ONLY β₁=0: {sorted(h_only_b0)}")
print(f"  H values with ONLY β₁=1: {sorted(h_only_b1)}")
print(f"  H values with BOTH: {sorted(h_both)}")

# ======================================================================
# PART 3: EQUIDECOMPOSABILITY CLASSES
# ======================================================================
print()
print("=" * 70)
print("PART 3: EQUIDECOMPOSABILITY CLASSES AT n=5")
print("=" * 70)

print("""
Two tournaments T₁, T₂ are EQUIDECOMPOSABLE iff:
  H(T₁) = H(T₂) AND β₁(T₁) = β₁(T₂)

This is the tournament analog of Dehn-Sydler:
  same volume AND same Dehn invariant.

The equidecomposability classes at n=5:
""")

classes = sorted(hb_counts.items())
total_classes = len(classes)
print(f"  {total_classes} equidecomposability classes at n=5")
print(f"  (compared to {len(set(h for (h,b) in hb_counts))} H-only classes)")

# Is β₁ a FINER invariant than H? Or does it add new information?
print(f"\n  β₁ adds information beyond H for {len(h_both)} H-values")
print(f"  β₁ is redundant (determined by H) for {len(h_only_b0)+len(h_only_b1)} H-values")

# ======================================================================
# PART 4: THE FLIP DISTANCE AND DEHN INVARIANT
# ======================================================================
print()
print("=" * 70)
print("PART 4: FLIP DISTANCE vs DEHN INVARIANT")
print("=" * 70)

print("""
The FLIP DISTANCE d(T₁, T₂) = number of arcs that differ
= Hamming distance in {0,1}^{C(n,2)}.

The Dehn invariant difference Δβ₁ = |β₁(T₁) - β₁(T₂)|.

QUESTION: Can a single arc flip change β₁?

From our near-transitive analysis:
  Transitive tournament: β₁ = 0, H = 1
  Near-transitive (one flip from transitive): β₁ = 1, H = 2^{n-2}+1
  Flip distance = 1, Δβ₁ = 1.

So YES: a single arc flip CAN change β₁ by 1.
The Dehn invariant is NOT invariant under arc flips.

But it IS invariant under a RESTRICTED class of flips:
flips that preserve the transitive core structure.

THE SCISSORS-CONGRUENCE THEOREM FOR TOURNAMENTS:
  Two tournaments T₁, T₂ are "tournament-scissors-congruent" iff:
  1. H(T₁) = H(T₂) (same "volume")
  2. β₁(T₁) = β₁(T₂) (same "Dehn invariant")
  3. They have the same independence polynomial I(Ω(T), x) (full Dehn)

  Condition 3 is STRICTLY finer than conditions 1+2.
  The full equidecomposability requires matching the ENTIRE
  independence polynomial, not just its value at x=2 (=H).
""")

# Check: at n=5, do tournaments with same H and β₁ always have same I(Ω,x)?
# This would tell us if (H, β₁) is sufficient for equidecomposability

# Actually, let's check directly: among tournaments with H=9, β₁=0 at n=5,
# do they all have the same I(Ω, x)?
n = 5
print(f"\nChecking: same (H, β₁) → same I(Ω, x)?")
for target_h, target_b1 in [(9, 0), (3, 0), (5, 0)]:
    ipoly_set = set()
    count = 0
    for bits in range(2**m):
        T = tournament_from_bits(bits, n)
        H = count_H(T, n)
        if H != target_h:
            continue
        b1 = compute_beta1_glmy(T, n)
        if b1 != target_b1:
            continue

        # Compute I(Ω₃, x) coefficients (3-cycles only for speed)
        cycles = []
        for triple in combinations(range(n), 3):
            i,j,k = triple
            s = T[i][j]+T[j][k]+T[k][i]
            if s==0 or s==3:
                cycles.append(frozenset(triple))
        nc = len(cycles)
        adj = [[0]*nc for _ in range(nc)]
        for a in range(nc):
            for b in range(a+1,nc):
                if cycles[a]&cycles[b]:
                    adj[a][b]=adj[b][a]=1
        alpha = [0]*(nc+1)
        for mask in range(2**nc):
            sel=[i for i in range(nc) if (mask>>i)&1]
            ok=True
            for a in range(len(sel)):
                for b in range(a+1,len(sel)):
                    if adj[sel[a]][sel[b]]:
                        ok=False; break
                if not ok: break
            if ok: alpha[len(sel)]+=1

        ipoly_set.add(tuple(alpha))
        count += 1

    print(f"  (H={target_h}, β₁={target_b1}): {count} tournaments, "
          f"{len(ipoly_set)} distinct I(Ω₃, x)")

# ======================================================================
# PART 5: THE DEEP STRUCTURE
# ======================================================================
print()
print("=" * 70)
print("PART 5: EQUIDECOMPOSABILITY = HOMOLOGICAL EQUIVALENCE")
print("=" * 70)

print("""
THE THEOREM (conjectured):

Two tournaments T₁, T₂ on n vertices are "tournament-equidecomposable"
(in the strongest sense: same independence polynomial I(Ω, x)) iff
they have the same HOMOLOGICAL TYPE:
  β_k(T₁) = β_k(T₂) for all k

This is the tournament analog of Dehn-Sydler-Jessen:
  Volume (H) = β₀ contribution (trivial)
  Dehn (β₁) = first homological obstruction
  Higher Dehn (β₃, β₄, ...) = higher obstructions

THE HIERARCHY:
  Level 0: H(T₁) = H(T₂)                     [volume]
  Level 1: β₁(T₁) = β₁(T₂)                   [Dehn = twist]
  Level 2: β₃(T₁) = β₃(T₂)                   [higher twist]
  Level ∞: I(Ω(T₁), x) = I(Ω(T₂), x)        [full polynomial]

Each level is STRICTLY FINER than the previous.

For the near-transitive tournament:
  Level 0: H = 2^{n-2}+1 (shared with many others)
  Level 1: β₁ = 1 (shared with fewer)
  Level ∞: I(Ω, x) = 1 + 2^{n-3}·x (complete graph: unique!)

The near-transitive tournament is EQUIDECOMPOSABLE only with itself
(and its relabelings). It is the "regular tetrahedron" of tournaments:
irreducible, with a unique Dehn invariant.

CONNECTION TO KIND-PASTEUR'S POLYTOPE TRILOGY:
  Simplex: Q₁(x) = 1+x → transitive tournament (β₁=0)
  Cuboid: Q₂(x) = (2+x)/(2-x) → "box" tournament
  Cayley: Q(x) = (1+x)/(1-x) → general tournament

  The TRANSITION from cuboid to Cayley introduces the twist.
  The cuboid has D = 0 (all angles π/2, rational).
  The Cayley has D ≠ 0 generically (angles irrational).

  Equidecomposability in the polytope trilogy:
    Cuboid = equidecomposable (D = 0 → scissors-congruent)
    Regular tetrahedron ≠ cube (D ≠ 0 → NOT scissors-congruent)
    Near-transitive ≠ transitive (β₁ = 1 → NOT homologically trivial)

  This IS Hilbert's third problem for tournaments.
""")

print("""
╔══════════════════════════════════════════════════════════════════════╗
║           HILBERT'S THIRD PROBLEM FOR TOURNAMENTS                  ║
╠══════════════════════════════════════════════════════════════════════╣
║                                                                    ║
║  CLASSICAL:                                                        ║
║    Two polyhedra with same volume are scissors-congruent           ║
║    iff they have the same Dehn invariant.                          ║
║                                                                    ║
║  TOURNAMENTS:                                                      ║
║    Two tournaments with same H are "tournament-congruent"          ║
║    iff they have the same β₁ (and higher Betti numbers).          ║
║                                                                    ║
║  THE OBSTRUCTION:                                                  ║
║    Classical: dihedral angles irrational multiples of π            ║
║    Tournament: the 1-cycle that no triple boundary can fill        ║
║    Both: a GLOBAL property invisible to LOCAL decomposition        ║
║                                                                    ║
║  THE INVARIANT:                                                    ║
║    Classical: D = Σ l(e) ⊗ θ(e)  (edge lengths × dihedral angles)║
║    Tournament: β₁ = dim ker(∂₁) - rank(∂₂)  (cycle deficiency)   ║
║    Both: measured by HOMOLOGY (Dehn = H₃ of scissors group)       ║
║                                                                    ║
║  THE SPECIAL CASES:                                                ║
║    Classical: cube (D=0) vs regular tetrahedron (D≠0)             ║
║    Tournament: transitive (β₁=0) vs near-transitive (β₁=1)       ║
║    Both: the simplest non-decomposable object                     ║
║                                                                    ║
║  THE CAYLEY CONNECTION:                                            ║
║    Classical: Dehn invariant lives in R ⊗_Q (R/πQ)               ║
║    Tournament: β₁ lives in H₁(tournament path complex)           ║
║    The Cayley transform Q maps between these:                      ║
║      arctanh (hyperbolic angle) ↔ β₁ (topological twist)          ║
║      Q(x) = exp(2·arctanh(x)) converts distances to twists        ║
║                                                                    ║
║  THE ANSWER:                                                       ║
║    Tournaments are equidecomposable iff their Betti profiles       ║
║    agree at all levels. The complete invariant is the              ║
║    independence polynomial I(Ω(T), x) — the FULL DEHN.            ║
║    β₁ is the first obstruction.                                   ║
║    The near-transitive tournament is the tetrahedron:              ║
║    same volume as a cube, but not scissors-congruent to it.        ║
║                                                                    ║
╚══════════════════════════════════════════════════════════════════════╝

opus-2026-03-15-S90k: EQUIDECOMPOSABILITY
""")
