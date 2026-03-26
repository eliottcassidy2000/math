#!/usr/bin/env python3
"""
CREATIVE SEESAW PROOF ATTEMPT

Instead of algebraic dimension counting, think GEOMETRICALLY.

IDEA 1: DELETION INDUCTION
  If β₁(T) = 1, then for every vertex v:
    β₁(T\v) ∈ {0, 1} (β₁ ≤ 1 is known)
  The LES gives: β₁(T) ≤ β₁(T\v) + dim H₁(T, T\v)
  And: β₃(T) ≤ β₃(T\v) + dim H₃(T, T\v)

  Key: if β₁(T) = 1, can we find a "good" vertex v where:
    - β₁(T\v) = 0 (the 1-hole goes through v)
    - dim H₃(T, T\v) = 0 (no new 3-holes from deleting v)
  Then β₃(T) ≤ β₃(T\v) ≤ [inductive bound for n-1]

IDEA 2: THE RELATIVE HOMOLOGY H_k(T, T\v)
  What IS H_k(T, T\v)? It measures the directed paths THROUGH v.
  A path in H_k(T, T\v) goes through v and its boundary misses v.

  For tournaments: every vertex v has d⁺ out-neighbors and d⁻ in-neighbors.
  The paths through v = combinations of in-arcs and out-arcs.

  dim H₁(T, T\v) counts the "independent directed triangles through v"
  that don't bound in the relative complex.

IDEA 3: THE β₁=1 CONDITION
  β₁ = 1 iff the tournament has a "harmonic 1-cycle" — a directed cycle
  that can't be expressed as a boundary.

  From THM-223: β₁ = 1 iff rank(TT) = C(n,2) - n (one deficiency).
  The transitive triple matrix TT has one fewer independent row.

  This means there's ONE "missing constraint" — one pair of vertices
  where the transitive triple structure doesn't force a relationship.

  IF this missing constraint FORCES β₃ = 0, we have the seesaw.

IDEA 4: COMPUTE H_k(T, T\v) EXPLICITLY
  For each vertex v in a tournament with β₁ = 1:
    Compute the relative homology H₃(T, T\v).
    If it's always 0 → the 3-hole "needs" v to exist, but v is
    the one supporting the 1-hole → both can't exist.

Let me test Idea 4 computationally.

kind-pasteur-2026-03-25-S27
"""
import sys, time, random
import numpy as np
from itertools import permutations
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)
random.seed(2026327)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def delete_vertex(A, n, v):
    """Delete vertex v from tournament."""
    B = []
    for i in range(n):
        if i == v: continue
        row = []
        for j in range(n):
            if j == v: continue
            row.append(A[i][j])
        B.append(row)
    return B

def enum_allowed(A, n, p):
    paths = []
    for perm in permutations(range(n), p+1):
        ok = True
        for i in range(p):
            if A[perm[i]][perm[i+1]] != 1:
                ok = False; break
        if ok: paths.append(perm)
    return paths

def compute_omega_and_rank(A, n, p):
    """Compute dim(Ω_p) and rank(∂_p)."""
    ap = enum_allowed(A, n, p)
    if not ap: return 0, 0, np.zeros((0,0))

    if p == 0: return len(ap), 0, np.eye(len(ap))

    ap_m1 = enum_allowed(A, n, p-1)
    ap_m1_set = set(ap_m1)

    non_allowed = {}; na_idx = 0
    for path in ap:
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face not in ap_m1_set and face not in non_allowed:
                non_allowed[face] = na_idx; na_idx += 1

    if not non_allowed:
        omega_dim = len(ap); omega_basis = np.eye(len(ap))
    else:
        P = np.zeros((len(non_allowed), len(ap)))
        for j, path in enumerate(ap):
            for i in range(p+1):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in non_allowed:
                    P[non_allowed[face], j] += sign
        _, S, Vt = np.linalg.svd(P)
        rk = sum(s > 1e-8 for s in S)
        omega_dim = len(ap) - rk
        omega_basis = Vt[rk:].T if omega_dim > 0 else np.zeros((len(ap), 0))

    # Boundary rank
    if ap_m1:
        idx = {path: i for i, path in enumerate(ap_m1)}
        D = np.zeros((len(ap_m1), len(ap)))
        for j, path in enumerate(ap):
            for i in range(p+1):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in idx: D[idx[face], j] += sign
        D_omega = D @ omega_basis
        _, S2, _ = np.linalg.svd(D_omega)
        rank_d = sum(s > 1e-8 for s in S2)
    else:
        rank_d = 0

    return omega_dim, rank_d, omega_basis

def compute_betti_fast(A, n, max_p=4):
    """Quick Betti computation up to dimension max_p."""
    results = {}
    for p in range(max_p + 1):
        omega_dim, rank_d, basis = compute_omega_and_rank(A, n, p)
        results[p] = {'omega': omega_dim, 'rank': rank_d, 'basis': basis}

    betti = []
    for p in range(max_p + 1):
        ker = results[p]['omega'] - results[p]['rank']
        # image from above
        if p + 1 <= max_p:
            im = results[p+1]['rank']  # rank(∂_{p+1}) = dimension of image
            # But we need image INTO Ω_p, not just rank of ∂_{p+1}
            # For now, approximate: the image rank is bounded by rank(∂_{p+1})
            # The exact computation needs projection into Ω_p
            pass
        im = results.get(p+1, {}).get('rank', 0) if p+1 <= max_p else 0
        betti.append(max(0, ker - im))
    return betti

# ============================================================
# MAIN: Search for tournaments with β₁=1 and check β₃
# ============================================================

print("=" * 80)
print("  CREATIVE SEESAW PROOF: Deletion analysis of β₁=1 tournaments")
print("  kind-pasteur-2026-03-25-S27")
print("=" * 80)

n = 7
print(f"\n  Searching for tournaments with β₁ = 1 at n={n}...")

found_b1 = 0
found_b3 = 0
found_both = 0
total = 0

b1_tournaments = []
b3_tournaments = []

t0 = time.time()
for trial in range(500):
    A = random_tournament(n)
    total += 1

    # Quick β₁ check: rank of ∂₂ relative to expected
    # β₁ = dim(ker ∂₁) - rank(∂₂|Ω₂)
    # dim(ker ∂₁) = C(n,2) - n + 1 = 15
    omega2, rank2, _ = compute_omega_and_rank(A, n, 2)

    # For tournaments: Ω₁ = A₁ (all arcs). rank(∂₁) = n-1.
    # β₁ = (C(n,2) - n + 1) - rank(∂₂|Ω₂)
    # Since we computed rank(∂₂|Ω₂) = rank2 (restricted to Omega_2):
    beta1 = max(0, (n*(n-1)//2 - n + 1) - rank2)

    # Quick β₃ check
    omega3, rank3, _ = compute_omega_and_rank(A, n, 3)
    omega4, rank4, _ = compute_omega_and_rank(A, n, 4)
    # β₃ ≈ (omega3 - rank3) - rank4... approximately
    beta3_approx = max(0, (omega3 - rank3) - rank4)

    if beta1 > 0:
        found_b1 += 1
        b1_tournaments.append(A)
    if beta3_approx > 0:
        found_b3 += 1
        b3_tournaments.append(A)
    if beta1 > 0 and beta3_approx > 0:
        found_both += 1

    if (trial+1) % 100 == 0:
        print(f"  [{trial+1}/500] b1>0: {found_b1}, b3>0: {found_b3}, both>0: {found_both}")

print(f"\n  RESULTS: {total} tournaments tested")
print(f"  β₁ > 0: {found_b1} ({found_b1/total*100:.1f}%)")
print(f"  β₃ > 0: {found_b3} ({found_b3/total*100:.1f}%)")
print(f"  BOTH > 0: {found_both} (seesaw {'HOLDS' if found_both == 0 else 'VIOLATED'})")

# ============================================================
# For β₁=1 tournaments: analyze vertex deletion
# ============================================================

if b1_tournaments:
    print(f"\n  === DELETION ANALYSIS of β₁=1 tournaments ===")
    print(f"  For each β₁=1 tournament, delete each vertex and check β₃(T\\v):")

    for idx, A in enumerate(b1_tournaments[:3]):
        print(f"\n  Tournament #{idx}:")
        scores = tuple(sorted(sum(A[i]) for i in range(n)))
        print(f"    Score sequence: {scores}")

        for v in range(n):
            Av = delete_vertex(A, n, v)
            # Compute β₁ and β₃ of T\v
            o2v, r2v, _ = compute_omega_and_rank(Av, n-1, 2)
            b1v = max(0, ((n-1)*(n-2)//2 - (n-1) + 1) - r2v)

            o3v, r3v, _ = compute_omega_and_rank(Av, n-1, 3)
            o4v, r4v, _ = compute_omega_and_rank(Av, n-1, 4)
            b3v = max(0, (o3v - r3v) - r4v)

            out_deg = sum(A[v])
            marker = "  ← 1-hole through v" if b1v == 0 else ""
            print(f"    v={v} (d⁺={out_deg}): β₁(T\\v)={b1v}, β₃(T\\v)={b3v}{marker}")

# ============================================================
# For β₃>0 tournaments: check if any vertex deletion gives β₁>0
# ============================================================

if b3_tournaments:
    print(f"\n  === DELETION ANALYSIS of β₃>0 tournaments ===")

    for idx, A in enumerate(b3_tournaments[:3]):
        print(f"\n  Tournament #{idx} (β₃>0):")
        scores = tuple(sorted(sum(A[i]) for i in range(n)))
        print(f"    Score sequence: {scores}")

        for v in range(n):
            Av = delete_vertex(A, n, v)
            o2v, r2v, _ = compute_omega_and_rank(Av, n-1, 2)
            b1v = max(0, ((n-1)*(n-2)//2 - (n-1) + 1) - r2v)

            o3v, r3v, _ = compute_omega_and_rank(Av, n-1, 3)
            o4v, r4v, _ = compute_omega_and_rank(Av, n-1, 4)
            b3v = max(0, (o3v - r3v) - r4v)

            out_deg = sum(A[v])
            marker = "  ← 3-hole removed" if b3v == 0 else ""
            print(f"    v={v} (d⁺={out_deg}): β₁(T\\v)={b1v}, β₃(T\\v)={b3v}{marker}")

print(f"""
  THE PROOF IDEA:
  If β₁(T) = 1, there exists a vertex v where β₁(T\\v) = 0.
  (The 1-hole "goes through" v — removing v kills it.)

  CONJECTURE: For this specific v, β₃(T\\v) = 0.
  (The vertex that supports the 1-hole is "critical" enough to prevent 3-holes.)

  Then by the LES of (T, T\\v):
    β₃(T) ≤ β₃(T\\v) + dim H₃(T, T\\v) = 0 + dim H₃(T, T\\v)

  If dim H₃(T, T\\v) = 0 (no new 3-holes from inserting v):
    β₃(T) = 0. Combined with β₁(T) = 1: seesaw holds.

  The key question: when a vertex v supports the 1-hole (β₁(T\\v)=0),
  does it also kill any potential 3-holes?
""")
