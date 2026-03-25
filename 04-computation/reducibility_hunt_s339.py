#!/usr/bin/env python3
"""
reducibility_hunt_s339.py — Hunt for the missing invariant
opus-2026-03-25-S339

Score + c3 + H capture ~54% of the multiway structure.
What captures the remaining ~46%?

IMPLICIT ASSUMPTIONS WE'RE BREAKING:
  A1: Only global invariants → try LOCAL (per-vertex, per-edge)
  A2: Only S_n-symmetric → try PARTIALLY symmetric (degree-sorted)
  A3: Only the tournament itself → try SUB-TOURNAMENTS
  A4: Only numbers → try STRUCTURES (graphs, polynomials, sequences)
  A5: Only undirected invariants → try DIRECTED structure
  A6: Only cycle counts → try PATH structure
  A7: Only algebraic → try TOPOLOGICAL (Betti numbers, homology)

20 CANDIDATE INVARIANTS TO TEST:
  1.  Degree-sorted adjacency matrix (canonical form "fingerprint")
  2.  Sorted neighbor-degree sequence (each vertex: degrees of out-neighbors)
  3.  c5 (5-cycle count)
  4.  bc33 (number of disjoint 3-cycle pairs)
  5.  Eigenvalue multiset of adjacency matrix
  6.  Eigenvalue multiset of SKEW adjacency (A - A^T)
  7.  Number of sources and sinks (vertices with extreme scores)
  8.  Degree sequence of the 3-cycle graph (which vertices are in 3-cycles)
  9.  The feedback arc set size (minimum arcs to remove for acyclicity)
  10. The independence number of the conflict graph Ω(T)
  11. Transitivity ratio (fraction of transitive triples)
  12. The number of Hamiltonian CYCLES (not just paths)
  13. β₁ (first Betti number from GLMY path homology)
  14. The automorphism group ORDER |Aut(T)|
  15. The number of FIXED POINTS under complement (self-comp detection)
  16. The parity vector (which arcs are "forward" in the base path)
  17. The sorted row-sum vector of A² (2-path counts)
  18. The diameter of the tournament (longest shortest directed path)
  19. The number of strongly connected components of T\{v} (vertex removal)
  20. The RANK of the adjacency matrix over GF(2)

For each: compute it for all tilings at n=5,6 and measure the
ADDITIONAL entropy reduction beyond score+c3+H.
"""

import sys
import numpy as np
from math import comb, factorial, log2
from collections import defaultdict, Counter
from itertools import combinations

sys.stdout.reconfigure(line_buffering=True)

def build_tilings(n):
    m = comb(n-1, 2)
    P = [(i,j) for i in range(n) for j in range(i+1, n)]
    base = [(i, i+1) for i in range(n-1)]
    non_base = [(i,j) for i,j in P if (i,j) not in base]
    N = 2**m

    tilings = []
    for tb in range(N):
        A = [[0]*n for _ in range(n)]
        for k in range(n-1): A[k+1][k] = 1
        for k, (i,j) in enumerate(non_base):
            if tb & (1 << k): A[j][i] = 1
            else: A[i][j] = 1

        s = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        c3 = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if A[i][j] and A[j][k] and A[k][i]: c3 += 1
                    if A[i][k] and A[k][j] and A[j][i]: c3 += 1

        dp = {}
        for v in range(n): dp[(1<<v, v)] = 1
        for S in range(1, 1<<n):
            for v in range(n):
                if not (S & (1<<v)): continue
                val = dp.get((S, v), 0)
                if val == 0: continue
                for u in range(n):
                    if S & (1<<u): continue
                    if A[v][u]: dp[(S|(1<<u), u)] = dp.get((S|(1<<u), u), 0) + val
        H = sum(dp.get(((1<<n)-1, v), 0) for v in range(n))

        tilings.append({'tb': tb, 'A': A, 'score': s, 'c3': c3, 'H': H, 'n': n})

    return tilings, m

def partition_entropy(values):
    """Entropy of a partition induced by a list of values."""
    counts = Counter(values)
    total = len(values)
    if total == 0: return 0
    h = 0
    for c in counts.values():
        p = c / total
        if p > 0: h -= p * log2(p)
    return h

def additional_reduction(base_keys, new_keys):
    """How much additional entropy reduction does new_key provide over base_keys?"""
    base_ent = partition_entropy(base_keys)
    combined = partition_entropy([str(b) + '|' + str(n) for b, n in zip(base_keys, new_keys)])
    return base_ent - combined

# ============================================================================
# INVARIANT COMPUTERS
# ============================================================================

def inv_score(t): return t['score']
def inv_c3(t): return t['c3']
def inv_H(t): return t['H']
def inv_base(t): return (t['score'], t['c3'], t['H'])

def inv_c5(t):
    """5-cycle count."""
    A = t['A']; n = t['n']; c = 0
    for path in combinations(range(n), 5):
        for perm in [(0,1,2,3,4),(0,1,2,4,3),(0,1,3,2,4),(0,1,3,4,2),(0,1,4,2,3),(0,1,4,3,2),
                     (0,2,1,3,4),(0,2,1,4,3),(0,2,3,1,4),(0,2,3,4,1),(0,2,4,1,3),(0,2,4,3,1),
                     (0,3,1,2,4),(0,3,1,4,2),(0,3,2,1,4),(0,3,2,4,1),(0,3,4,1,2),(0,3,4,2,1),
                     (0,4,1,2,3),(0,4,1,3,2),(0,4,2,1,3),(0,4,2,3,1),(0,4,3,1,2),(0,4,3,2,1)]:
            verts = [path[p] for p in perm]
            if all(A[verts[i]][verts[(i+1)%5]] for i in range(5)):
                c += 1
    return c // 10  # each 5-cycle counted 10 times (5 rotations × 2 directions)

def inv_bc33(t):
    """Disjoint 3-cycle pairs."""
    A = t['A']; n = t['n']
    cycles3 = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles3.append(frozenset([i,j,k]))
                if A[i][k] and A[k][j] and A[j][i]:
                    cycles3.append(frozenset([i,j,k]))
    # Deduplicate (each 3-set appears at most twice: CW and CCW)
    unique = list(set(cycles3))
    bc = 0
    for a in range(len(unique)):
        for b in range(a+1, len(unique)):
            if not unique[a] & unique[b]:
                bc += 1
    return bc

def inv_eigenvalues(t):
    """Sorted eigenvalue magnitudes of adjacency matrix."""
    A = np.array(t['A'], dtype=float)
    eigs = np.sort(np.abs(np.linalg.eigvals(A)))[::-1]
    return tuple(round(e, 4) for e in eigs)

def inv_skew_eigenvalues(t):
    """Sorted eigenvalues of skew matrix A - A^T."""
    A = np.array(t['A'], dtype=float)
    S = A - A.T
    eigs = np.sort(np.abs(np.linalg.eigvals(S)))[::-1]
    return tuple(round(e, 4) for e in eigs)

def inv_neighbor_degrees(t):
    """For each vertex, sorted degrees of out-neighbors. Then sort the list."""
    A = t['A']; n = t['n']
    out_degs = [sum(A[i]) for i in range(n)]
    profiles = []
    for i in range(n):
        out_nbrs = sorted([out_degs[j] for j in range(n) if A[i][j]])
        profiles.append(tuple(out_nbrs))
    return tuple(sorted(profiles))

def inv_row_sums_A2(t):
    """Sorted row sums of A² (2-path counts per vertex)."""
    A = np.array(t['A'], dtype=int)
    A2 = A @ A
    return tuple(sorted(A2.sum(axis=1)))

def inv_row_sums_A3(t):
    """Sorted row sums of A³ (3-path counts per vertex)."""
    A = np.array(t['A'], dtype=int)
    A3 = A @ A @ A
    return tuple(sorted(A3.sum(axis=1)))

def inv_trace_A2(t):
    """tr(A²) = number of 2-cycles = 0 for tournaments (but tr(A²) counts directed 2-paths)."""
    A = np.array(t['A'], dtype=int)
    return int(np.trace(A @ A))

def inv_trace_A3(t):
    """tr(A³) = 6 × c3 (number of directed triangles)."""
    A = np.array(t['A'], dtype=int)
    return int(np.trace(A @ A @ A))

def inv_gf2_rank(t):
    """Rank of adjacency matrix over GF(2)."""
    A = np.array(t['A'], dtype=int) % 2
    n = t['n']
    M = A.copy()
    rank = 0
    for col in range(n):
        pivot = None
        for row in range(rank, n):
            if M[row, col] == 1:
                pivot = row; break
        if pivot is None: continue
        M[[rank, pivot]] = M[[pivot, rank]]
        for row in range(n):
            if row != rank and M[row, col] == 1:
                M[row] = (M[row] + M[rank]) % 2
        rank += 1
    return rank

def inv_aut_size(t):
    """Automorphism group size (brute force for small n)."""
    A = t['A']; n = t['n']
    from itertools import permutations
    count = 0
    for p in permutations(range(n)):
        is_aut = True
        for i in range(n):
            for j in range(n):
                if A[i][j] != A[p[i]][p[j]]:
                    is_aut = False; break
            if not is_aut: break
        if is_aut: count += 1
    return count

def inv_feedback_arc(t):
    """Size of minimum feedback arc set (arcs to remove for acyclicity).
    For tournament: = C(n,2) - longest path = C(n,2) - (n-1) for transitive.
    Actually: = number of "backward" arcs in a topological ordering.
    Use greedy: topological sort by score, count backward arcs."""
    A = t['A']; n = t['n']
    scores = [sum(A[i]) for i in range(n)]
    order = sorted(range(n), key=lambda v: -scores[v])
    backward = 0
    for idx_i in range(n):
        for idx_j in range(idx_i+1, n):
            i, j = order[idx_i], order[idx_j]
            if A[j][i]:  # j beats i but j comes later in the ordering
                backward += 1
    return backward

def inv_hamcycles(t):
    """Number of Hamiltonian cycles."""
    A = t['A']; n = t['n']
    # Fix vertex 0 as start, count cycles returning to 0
    dp = {}
    dp[(1 << 0, 0)] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)): continue
            val = dp.get((S, v), 0)
            if val == 0: continue
            for u in range(n):
                if S & (1 << u): continue
                if A[v][u]:
                    dp[(S | (1 << u), u)] = dp.get((S | (1 << u), u), 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n) if A[v][0])

# ============================================================================
# MAIN TEST
# ============================================================================

print("=" * 72)
print("  HUNT FOR THE MISSING INVARIANT")
print("  opus-2026-03-25-S339")
print("=" * 72)

for n in [5, 6]:
    print(f"\n{'#'*72}")
    print(f"  n = {n}")
    print(f"{'#'*72}")

    tilings, m = build_tilings(n)
    N = len(tilings)
    full_ent = log2(N)

    # Base: score + c3 + H
    base_keys = [inv_base(t) for t in tilings]
    base_ent = partition_entropy(base_keys)
    base_classes = len(set(base_keys))

    print(f"  {N} tilings, {full_ent:.1f} bits full entropy")
    print(f"  Base (score+c3+H): {base_ent:.2f} bits, {base_classes} classes")
    print(f"  Reducibility: {(1-base_ent/full_ent)*100:.1f}%")
    print(f"  Irreducible: {base_ent/full_ent*100:.1f}% ({base_ent:.2f} bits)")

    # Test each candidate invariant
    invariants = [
        ("c5",             inv_c5),
        ("bc33",           inv_bc33),
        ("eigenvalues",    inv_eigenvalues),
        ("skew_eig",       inv_skew_eigenvalues),
        ("nbr_degrees",    inv_neighbor_degrees),
        ("row_A2",         inv_row_sums_A2),
        ("row_A3",         inv_row_sums_A3),
        ("tr_A2",          inv_trace_A2),
        ("tr_A3",          inv_trace_A3),
        ("gf2_rank",       inv_gf2_rank),
        ("feedback_arc",   inv_feedback_arc),
        ("ham_cycles",     inv_hamcycles),
    ]

    if n <= 5:
        invariants.append(("aut_size", inv_aut_size))

    print(f"\n  ADDITIONAL REDUCTION from each invariant (beyond score+c3+H):")
    print(f"  {'Invariant':>15} {'Classes':>8} {'Entropy':>8} {'Reduction':>10} {'New info':>10}")

    results = []
    for name, fn in invariants:
        try:
            values = [fn(t) for t in tilings]
            # Combined with base
            combined_keys = [str(base_keys[i]) + '|' + str(values[i]) for i in range(N)]
            combined_ent = partition_entropy(combined_keys)
            combined_classes = len(set(combined_keys))
            reduction = base_ent - combined_ent
            pct = reduction / full_ent * 100

            # Standalone entropy
            standalone_ent = partition_entropy(values)
            standalone_classes = len(set(values))

            results.append((name, combined_classes, combined_ent, reduction, pct, standalone_classes))
            print(f"  {name:>15} {combined_classes:>8} {combined_ent:>8.2f} {reduction:>9.2f}b {pct:>9.1f}%"
                  f"  ({standalone_classes} standalone)")
        except Exception as e:
            print(f"  {name:>15} ERROR: {e}")

    # Sort by additional reduction
    results.sort(key=lambda x: -x[3])

    print(f"\n  RANKED BY ADDITIONAL INFORMATION:")
    for name, cls, ent, red, pct, sc in results[:5]:
        print(f"    {name:>15}: +{red:.2f} bits ({pct:.1f}%), total {cls} classes")

    # What's the BEST combination of base + top 3 invariants?
    if results:
        top3 = [r[0] for r in results[:3]]
        top3_fns = {name: fn for name, fn in invariants}
        combo_keys = base_keys[:]
        for name in top3:
            if name in top3_fns:
                vals = [top3_fns[name](t) for t in tilings]
                combo_keys = [str(combo_keys[i]) + '|' + str(vals[i]) for i in range(N)]
        combo_ent = partition_entropy(combo_keys)
        combo_cls = len(set(combo_keys))
        total_red = (1 - combo_ent / full_ent) * 100
        print(f"\n  BEST COMBO (base + {top3}):")
        print(f"    Classes: {combo_cls}, Entropy: {combo_ent:.2f}, Reducibility: {total_red:.1f}%")

print(f"\n{'='*72}")
print(f"  WHAT WE LEARNED")
print(f"{'='*72}")

print("DONE.")
