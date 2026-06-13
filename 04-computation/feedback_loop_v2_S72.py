#!/usr/bin/env python3
"""
feedback_loop_v2_S72.py — opus-2026-03-15-S72

Faster, focused version. Tests refined invariant combinations.

KEY QUESTION: What minimal set of invariants determines Ω_k?
From v1: score sequence is NOT enough at n=5.
Test: does (score_seq, t3) work? What about (score_seq, β₁)?

Also: test H = n · H(T) deletion identity and rank formulas.
"""

import numpy as np
from math import comb, factorial
from collections import Counter, defaultdict
import sys, time
sys.stdout.reconfigure(line_buffering=True)

PRIME = 2**31 - 1

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=np.int8)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def score_seq(A, n):
    return tuple(sorted(int(np.sum(A[i])) for i in range(n)))

def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i,j] and A[j,k] and A[k,i]) or \
                   (A[i,k] and A[k,j] and A[j,i]):
                    c3 += 1
    return c3

def ham_path_count(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v, u] == 1:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def enum_allowed(A, n, length):
    if length == 0:
        return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                adj[i].append(j)
    paths = []
    for start in range(n):
        stack = [([start], 1 << start)]
        while stack:
            path, visited = stack.pop()
            if len(path) == length + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def gauss_rank_np(M, prime):
    M = M.copy() % prime
    nrows, ncols = M.shape
    rank = 0
    for col in range(ncols):
        nonzero = np.where(M[rank:, col] != 0)[0]
        if len(nonzero) == 0:
            continue
        pivot = nonzero[0] + rank
        if pivot != rank:
            M[[rank, pivot]] = M[[pivot, rank]]
        inv = pow(int(M[rank, col]), prime - 2, prime)
        M[rank] = M[rank] * inv % prime
        factors = M[:, col].copy()
        factors[rank] = 0
        nzr = np.where(factors != 0)[0]
        if len(nzr) > 0:
            M[nzr] = (M[nzr] - np.outer(factors[nzr], M[rank])) % prime
        rank += 1
    return rank

def compute_omega_and_ranks(A, n, max_p=None):
    """Compute Ω dims and boundary ranks. Returns (omega_dims, boundary_ranks, betti)."""
    if max_p is None:
        max_p = n - 1
    max_p = min(max_p, n - 1)

    allowed = {}
    for p in range(max_p + 2):
        allowed[p] = enum_allowed(A, n, p) if p <= n - 1 else []

    omega_dims = {}
    omega_bases = {}
    for p in range(max_p + 1):
        paths = allowed.get(p, [])
        if not paths:
            omega_dims[p] = 0; omega_bases[p] = None; continue
        if p == 0:
            omega_dims[p] = len(paths); omega_bases[p] = None; continue

        set_prev = set(allowed.get(p-1, []))
        non_allowed = {}; na_count = 0
        for path in paths:
            for i in range(len(path)):
                face = path[:i] + path[i+1:]
                if len(set(face)) == len(face) and face not in set_prev:
                    if face not in non_allowed:
                        non_allowed[face] = na_count; na_count += 1
        if na_count == 0:
            omega_dims[p] = len(paths); omega_bases[p] = None; continue

        P = np.zeros((na_count, len(paths)), dtype=np.int64)
        for j, path in enumerate(paths):
            for i in range(len(path)):
                face = path[:i] + path[i+1:]
                if face in non_allowed:
                    P[non_allowed[face], j] = (P[non_allowed[face], j] + (-1)**i) % PRIME

        # Just compute rank for Ω dim (no need for full basis for most tests)
        rank_P = gauss_rank_np(P % PRIME, PRIME)
        omega_dims[p] = len(paths) - rank_P
        omega_bases[p] = None  # Skip basis for speed

    # Skip boundary rank computation (expensive) — just return Ω dims
    betti = {}
    return omega_dims, {}, betti

def compute_beta1_fast(A, n):
    """Fast β₁ = C(n,2) - (n-1) - rank(∂₂ restricted to Ω₂)."""
    edges = []
    edge_idx = {}
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                edges.append((i, j))
                edge_idx[(i, j)] = len(edges) - 1

    triples = []
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0: continue
            for k in range(n):
                if k in (i, j) or A[j][k] == 0: continue
                if A[i][k] == 1:
                    triples.append((i, j, k))

    if not triples:
        return comb(n, 2) - (n - 1)

    bd = np.zeros((len(edges), len(triples)), dtype=np.int64)
    for t, (i, j, k) in enumerate(triples):
        if (j, k) in edge_idx:
            bd[edge_idx[(j, k)], t] = (bd[edge_idx[(j, k)], t] + 1) % PRIME
        if (i, k) in edge_idx:
            bd[edge_idx[(i, k)], t] = (bd[edge_idx[(i, k)], t] - 1) % PRIME
        if (i, j) in edge_idx:
            bd[edge_idx[(i, j)], t] = (bd[edge_idx[(i, j)], t] + 1) % PRIME

    rank_d2 = gauss_rank_np(bd % PRIME, PRIME)
    return max(0, comb(n, 2) - (n - 1) - rank_d2)

print("=" * 72)
print("FEEDBACK LOOP v2 — REFINED INVARIANT SEARCH")
print("opus-2026-03-15-S72")
print("=" * 72)

# ============================================================
# TEST 1: Does (score_seq, t3) determine Ω_k?
# ============================================================

print("\n" + "=" * 72)
print("TEST 1: Ω_k from (score_seq, t3)?")
print("=" * 72)

for n in [5, 6]:
    print(f"\n  n={n}:")
    m = n * (n - 1) // 2
    total = 2**m

    t0 = time.time()
    groups = defaultdict(list)

    for bits in range(total):
        A = bits_to_adj(bits, n)
        ss = score_seq(A, n)
        t3 = count_3cycles(A, n)
        omega_dims, _, _ = compute_omega_and_ranks(A, n, max_p=min(n-1, 4))
        omega_tuple = tuple(omega_dims.get(k, 0) for k in range(min(n, 5)))
        groups[(ss, t3)].append(omega_tuple)

    violations = 0
    for key, omegas in sorted(groups.items()):
        unique = set(omegas)
        if len(unique) > 1:
            violations += 1
            ss, t3 = key
            if violations <= 5:
                print(f"    VIOLATION: score_seq={ss}, t3={t3} → {len(unique)} Ω vectors")
                for u in sorted(unique):
                    print(f"      Ω={list(u)}: {omegas.count(u)}")

    elapsed = time.time() - t0
    if violations == 0:
        print(f"    ★ Ω IS determined by (score_seq, t3) at n={n}! ({elapsed:.1f}s)")
    else:
        print(f"    ✗ {violations} violations ({elapsed:.1f}s)")

# ============================================================
# TEST 2: Does (score_seq, t3, β₁) determine Ω_k?
# ============================================================

print("\n" + "=" * 72)
print("TEST 2: Ω_k from (score_seq, t3, β₁)?")
print("=" * 72)

for n in [5]:
    print(f"\n  n={n}:")
    m = n * (n - 1) // 2
    total = 2**m

    groups = defaultdict(list)

    for bits in range(total):
        A = bits_to_adj(bits, n)
        ss = score_seq(A, n)
        t3 = count_3cycles(A, n)
        b1 = compute_beta1_fast(A, n)
        omega_dims, _, _ = compute_omega_and_ranks(A, n, max_p=min(n-1, 4))
        omega_tuple = tuple(omega_dims.get(k, 0) for k in range(min(n, 5)))
        groups[(ss, t3, b1)].append(omega_tuple)

    violations = 0
    for key, omegas in sorted(groups.items()):
        unique = set(omegas)
        if len(unique) > 1:
            violations += 1
            ss, t3, b1 = key
            print(f"    VIOLATION: {ss}, t3={t3}, β₁={b1} → {len(unique)} Ω vectors")

    if violations == 0:
        print(f"    ★ Ω IS determined by (score_seq, t3, β₁) at n={n}!")
        # Print the mapping
        for key, omegas in sorted(groups.items()):
            ss, t3, b1 = key
            print(f"    {ss}, t3={t3}, β₁={b1} → Ω={list(omegas[0])}")
    else:
        print(f"    ✗ {violations} violations")

# ============================================================
# TEST 3: H deletion identity
# ============================================================

print("\n" + "=" * 72)
print("TEST 3: H(T) = Σ H(T\\v) / n ?  Or H = Σᵥ H_v / (n-1)?")
print("If a universal formula exists → O(n!) → O(n · (n-1)!)")
print("=" * 72)

for n in [4, 5, 6]:
    print(f"\n  n={n}:")
    m = n * (n - 1) // 2
    total = 2**m
    if n == 6:
        # Sample for speed
        np.random.seed(42)
        sample = [np.random.randint(0, total) for _ in range(1000)]
    else:
        sample = range(total)

    sum_match = 0
    ratios = set()

    for bits in sample:
        A = bits_to_adj(bits, n)
        H = ham_path_count(A, n)

        sum_del = 0
        for v in range(n):
            verts = [i for i in range(n) if i != v]
            A_del = np.zeros((n-1, n-1), dtype=np.int8)
            for i, vi in enumerate(verts):
                for j, vj in enumerate(verts):
                    A_del[i, j] = A[vi, vj]
            sum_del += ham_path_count(A_del, n-1)

        if H > 0:
            ratios.add(sum_del)
            if sum_del == n * H:
                sum_match += 1

    total_tested = len(list(sample))
    print(f"    Σ H(T\\v) = n · H(T): {sum_match}/{total_tested}")

    # Maybe H(T) = Σ H(T\\v) exactly?
    exact_match = 0
    for bits in (range(total) if n < 6 else [np.random.randint(0, total) for _ in range(200)]):
        A = bits_to_adj(bits, n)
        H = ham_path_count(A, n)
        sum_del = 0
        for v in range(n):
            verts = [i for i in range(n) if i != v]
            A_del = np.zeros((n-1, n-1), dtype=np.int8)
            for i, vi in enumerate(verts):
                for j, vj in enumerate(verts):
                    A_del[i, j] = A[vi, vj]
            sum_del += ham_path_count(A_del, n-1)
        if sum_del == H:
            exact_match += 1

    tested2 = total if n < 6 else 200
    print(f"    Σ H(T\\v) = H(T): {exact_match}/{tested2}")

    # What IS the relationship?
    # Check: each HP goes through exactly n vertices, and deleting the LAST
    # vertex v gives an HP in T\v. So H(T) = Σᵥ (# HP ending at v in T)
    # And H(T\v) counts all HP in T\v, not just those related to T.
    # Actually: Σ H(T\v) counts each (n-1)-path once for each vertex not in it.
    # But in a TOURNAMENT, every (n-1)-subset spans a tournament.
    # There's no simple universal relationship because H(T\v) depends on ALL
    # arcs among V\{v}, not just those in HP through v.

# ============================================================
# TEST 4: Ω₂ formula (most impactful)
# ============================================================

print("\n" + "=" * 72)
print("TEST 4: Ω₂ closed form from (n, t3, score_seq)")
print("This is the BIGGEST win — Ω₂ computation dominates runtime")
print("=" * 72)

for n in [4, 5, 6]:
    print(f"\n  n={n}:")
    m = n * (n - 1) // 2
    total = 2**m

    omega2_data = []
    t0 = time.time()

    for bits in range(total):
        A = bits_to_adj(bits, n)
        t3 = count_3cycles(A, n)

        # Compute |A_2| (allowed 2-paths = transitive triples)
        a2 = 0
        for i in range(n):
            for j in range(n):
                if i == j or A[i][j] == 0: continue
                for k in range(n):
                    if k in (i,j) or A[j][k] == 0: continue
                    if A[i][k] == 1:
                        a2 += 1

        # Compute Ω₂
        omega_dims, _, _ = compute_omega_and_ranks(A, n, max_p=2)
        omega2 = omega_dims.get(2, 0)

        # The constraint dim: # non-allowed faces of 2-paths
        # Each 2-path (i,j,k) has 3 faces: (j,k), (i,k), (i,j)
        # Non-allowed face: a pair (a,b) where A[a][b]=0 and A[b][a]=0 (impossible in tournament!)
        # OR: a pair (a,b) where b->a (reverse direction)
        # In a tournament: EVERY pair is an edge in one direction
        # So the only non-allowed faces are those pointing the wrong way
        # But wait — (i,k) is a face of (i,j,k). If k->i (not i->k), this face is NOT in A_1.
        # But A_1 = all directed edges. (i,k) is in A_1 iff A[i][k]=1.
        # For transitive triple (i,j,k): we have i->j, j->k, i->k. So face (i,k) IS allowed.
        # For non-transitive but allowed: (i,j,k) with i->j, j->k. Face (i,k): i->k or k->i.
        # If i->k: all 3 faces allowed → no constraint
        # If k->i: face (i,k) is not an allowed 1-path direction...
        # Wait, A_1 = {(a,b) : a->b}. Face (i,k) needs A[i][k]=1.
        # In transitive triple: A[i][k]=1 → face allowed.
        # But (i,j,k) is only a valid 2-path if i->j->k. Does NOT require i->k.
        # So Ω₂ constraint: for each 2-path (i,j,k) where k->i, the face (i,k) is not in A_1.
        # This gives a constraint on Ω₂.

        # Count non-transitive 2-paths (allowed path i->j->k where k->i)
        nt2 = 0
        for i in range(n):
            for j in range(n):
                if i == j or A[i][j] == 0: continue
                for k in range(n):
                    if k in (i,j) or A[j][k] == 0: continue
                    if A[k][i] == 1:  # NOT i->k, but k->i
                        nt2 += 1

        omega2_data.append((t3, a2, nt2, omega2))

    elapsed = time.time() - t0

    # Check: is Ω₂ = A₂ - nt2? (non-transitive paths give constraints)
    # Actually Ω₂ = dim(ker of constraint matrix)
    # Constraint matrix has nt2 unique constraints (one per non-allowed face)

    # Count unique non-allowed faces per tournament
    # Each non-transitive 2-path (i,j,k) with k->i has face (i,k) non-allowed.
    # Multiple paths can share the same non-allowed face.
    # The constraint count is the number of UNIQUE non-allowed (p-1)-paths.

    # Test formulas
    formulas_tested = {}

    # Formula 1: Ω₂ = |A₂| - rank(constraints) where rank = some function of nt2
    # We know from THM-119: constraint matrix is always full rank
    # So Ω₂ = |A₂| - (# unique non-allowed faces)
    # # unique non-allowed faces for degree 2 = # edges (a,b) that appear as
    # face of some 2-path but NOT in A₁
    # In tournament: all edges are in A₁ (just pick direction), so... hmm

    # Actually, the non-allowed face is an ORDERED pair where the edge goes the wrong way.
    # (i,k) is a face of (i,j,k) and is non-allowed iff A[i][k]=0 (i.e., k->i).
    # Each such pair (i,k) with k->i can appear as face of MULTIPLE 2-paths.
    # # unique non-allowed faces = # directed pairs (a,b) where b->a AND ∃ middle vertex j
    # with a->j->b. In a tournament, b->a means (a,b) is NOT an edge,
    # but it appears as face of a 2-path a->j->b where a->j and j->b.
    # The # such pairs = # pairs (a,b) with b->a such that some j has a->j->b.

    # Let me just test: Ω₂ = |A₂| - (# non-allowed face pairs)?
    # This requires computing # unique non-allowed faces per tournament.

    # For now, just check if (t3, |A₂|) determines Ω₂
    groups_t3_a2 = defaultdict(set)
    for t3, a2, nt2, omega2 in omega2_data:
        groups_t3_a2[(t3, a2)].add(omega2)

    unique_keys = sum(1 for v in groups_t3_a2.values() if len(v) == 1)
    total_keys = len(groups_t3_a2)
    print(f"    (t3, |A₂|) determines Ω₂: {unique_keys}/{total_keys} unique")

    # Does |A₂| alone determine Ω₂?
    groups_a2 = defaultdict(set)
    for t3, a2, nt2, omega2 in omega2_data:
        groups_a2[a2].add(omega2)
    a2_unique = sum(1 for v in groups_a2.values() if len(v) == 1)
    print(f"    |A₂| alone determines Ω₂: {a2_unique}/{len(groups_a2)} unique")

    # Does t3 alone determine Ω₂?
    groups_t3 = defaultdict(set)
    for t3, a2, nt2, omega2 in omega2_data:
        groups_t3[t3].add(omega2)
    t3_unique = sum(1 for v in groups_t3.values() if len(v) == 1)
    print(f"    t3 alone determines Ω₂: {t3_unique}/{len(groups_t3)} unique")
    for t3, o2_set in sorted(groups_t3.items()):
        print(f"      t3={t3}: Ω₂ ∈ {sorted(o2_set)}")

    # What about Ω₂ = |A₂| - f(non-transitive) ?
    groups_formula = defaultdict(set)
    for t3, a2, nt2, omega2 in omega2_data:
        groups_formula[(a2, nt2)].add(omega2)
    f_unique = sum(1 for v in groups_formula.values() if len(v) == 1)
    print(f"    (|A₂|, nt₂) determines Ω₂: {f_unique}/{len(groups_formula)} unique")

    # Check if Ω₂ = |A₂| - something?
    for t3, a2, nt2, omega2 in omega2_data[:3]:
        print(f"      t3={t3}, |A₂|={a2}, nt₂={nt2}, Ω₂={omega2}, gap={a2-omega2}")

    # Test: Ω₂ = f(t3)?
    if len(groups_t3) > 1 and t3_unique == len(groups_t3):
        t3_list = sorted(groups_t3.keys())
        o2_list = [list(groups_t3[t])[0] for t in t3_list]
        # Try linear: Ω₂ = a + b*t3
        if len(t3_list) >= 2:
            coeffs = np.polyfit(t3_list, o2_list, 1)
            pred = [round(np.polyval(coeffs, t)) for t in t3_list]
            if pred == o2_list:
                print(f"    ★ LINEAR FORMULA: Ω₂ = {coeffs[1]:.1f} + {coeffs[0]:.1f} * t3")

    print(f"    ({elapsed:.1f}s)")

# ============================================================
# TEST 5: β₁ formula from score sequence + t3
# ============================================================

print("\n" + "=" * 72)
print("TEST 5: β₁ from (score_seq, t3)?")
print("=" * 72)

for n in [5, 6]:
    print(f"\n  n={n}:")
    m = n * (n - 1) // 2
    total = 2**m

    groups = defaultdict(set)
    t0 = time.time()

    for bits in range(total):
        A = bits_to_adj(bits, n)
        ss = score_seq(A, n)
        t3 = count_3cycles(A, n)
        b1 = compute_beta1_fast(A, n)
        groups[(ss, t3)].add(b1)

    violations = sum(1 for v in groups.values() if len(v) > 1)
    elapsed = time.time() - t0

    if violations == 0:
        print(f"    ★ β₁ IS determined by (score_seq, t3) at n={n}! ({elapsed:.1f}s)")
        # Print mapping
        for key, b1_set in sorted(groups.items()):
            ss, t3 = key
            print(f"    {ss}, t3={t3} → β₁={list(b1_set)[0]}")
    else:
        print(f"    ✗ {violations} violations ({elapsed:.1f}s)")
        # Show violations
        count = 0
        for key, b1_set in sorted(groups.items()):
            if len(b1_set) > 1:
                ss, t3 = key
                print(f"    VIOLATION: {ss}, t3={t3} → β₁ ∈ {sorted(b1_set)}")
                count += 1
                if count >= 5:
                    break

# ============================================================
# TEST 6: |A_k| formula (allowed path counts)
# ============================================================

print("\n" + "=" * 72)
print("TEST 6: Is |A_k| (allowed k-path count) determined by t3?")
print("If YES → we know the PATH SPACE size from cycle count alone")
print("=" * 72)

for n in [5, 6]:
    print(f"\n  n={n}:")
    m = n * (n - 1) // 2
    total = 2**m

    groups = defaultdict(set)

    for bits in range(total):
        A = bits_to_adj(bits, n)
        t3 = count_3cycles(A, n)

        ak_tuple = tuple(len(enum_allowed(A, n, k)) for k in range(min(n, 5)))
        groups[t3].add(ak_tuple)

    for t3, ak_set in sorted(groups.items()):
        if len(ak_set) == 1:
            print(f"    t3={t3}: |A_k| = {list(ak_set)[0]} (unique)")
        else:
            print(f"    t3={t3}: {len(ak_set)} different |A_k| vectors")

print("\n" + "=" * 72)
print("SYNTHESIS")
print("=" * 72)
print("""
PROOF → SPEEDUP TARGETS (ranked by impact):

1. If Ω₂ = f(t3): skip the LARGEST constraint matrix → O(n³) → O(1)
   (Test 4 results show if this works)

2. If β₁ = g(score_seq, t3): skip rank computation → O(n³) → O(n²)
   (Test 5 results show if this works)

3. If (score_seq, t3, β₁) determines all Ω_k: universal formula → skip ALL constraint matrices
   (Test 2 results show if this works)

4. If |A_k| = h(t3): skip path enumeration → O(n!) → O(n²)
   (Test 6 results show if this works)

The FEEDBACK LOOP:
  PROVE any of these → SPEED UP by 100-1000× → COMPUTE at n=15-20
  → DISCOVER new patterns → PROVE more formulas → ...
""")

print("\nDONE — feedback_loop_v2_S72.py complete")
