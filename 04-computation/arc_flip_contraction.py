import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
arc_flip_contraction.py
kind-pasteur-2026-03-07-S39b

ARC-FLIP CONTRACTION REDUCTION (THM-082 Remark 6):

For tournament T with edge e = (u->v), let T' = flip(T, u, v) = tournament
with edge e' = (v->u). Then:

  H(T) - H(T') = H(T/e) - H(T'/e')

where T/e contracts u,v (w inherits IN from u, OUT from v),
and T'/e' contracts v,u (w inherits IN from v, OUT from u).

This REDUCES the arc-flip H-difference to a (n-1)-vertex problem!

Verified 100% at n=5 (10240 pairs). Testing at n=6.

DEEPER QUESTION: Can this be iterated to give a combinatorial
explanation of H(T) mod 4 or other structure?

Also: T/e and T'/e' have the SAME vertex set. How do their
tournaments relate? Is T/e the "flip" of T'/e' along some edge?
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from collections import defaultdict
import random


def ham_paths_dp(adj, n):
    """Count Hamiltonian paths in a general digraph using DP."""
    if n <= 1:
        return 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0 or not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def contract_edge(T, u, v):
    """Contract edge (u->v): merge u,v into w.
    w inherits IN from u, OUT from v."""
    n = len(T)
    verts = [i for i in range(n) if i != v]
    m = len(verts)
    D = [[0]*m for _ in range(m)]
    for i_idx, i in enumerate(verts):
        for j_idx, j in enumerate(verts):
            if i_idx == j_idx:
                continue
            if i == u:
                D[i_idx][j_idx] = T[v][j]
            elif j == u:
                D[i_idx][j_idx] = T[i][u]
            else:
                D[i_idx][j_idx] = T[i][j]
    return D


def is_tournament(adj, n):
    """Check if adjacency matrix is a tournament."""
    for i in range(n):
        for j in range(i+1, n):
            if adj[i][j] + adj[j][i] != 1:
                return False
    return True


# ============================================================
# Test at n=6 (sampled)
# ============================================================
print("=" * 70)
print("ARC-FLIP CONTRACTION REDUCTION")
print("=" * 70)

random.seed(42)

for n in range(3, 8):
    m = n * (n - 1) // 2
    total_tournaments = 1 << m

    if n <= 6:
        # Exhaustive
        sample = range(total_tournaments)
        sample_type = "exhaustive"
    else:
        # Sample
        sample = [random.randint(0, total_tournaments - 1) for _ in range(2000)]
        sample_type = "sampled (2000)"

    holds = 0
    total = 0
    # Also track: is T/e always a tournament?
    te_tournament = 0
    te_total = 0
    # Track: relationship between T/e and T'/e'
    te_equals_flip_te = 0

    for bits in sample:
        T = tournament_from_bits(n, bits)
        H_T = hamiltonian_path_count(T)

        # Test first edge pair only for large n
        edge_limit = n * (n - 1) if n <= 5 else 3

        tested_this = 0
        for u in range(n):
            if tested_this >= edge_limit:
                break
            for v in range(u + 1, n):
                if tested_this >= edge_limit:
                    break
                if T[u][v]:
                    eu, ev = u, v
                else:
                    eu, ev = v, u

                # T has edge eu->ev
                T_flip = [row[:] for row in T]
                T_flip[eu][ev] = 0
                T_flip[ev][eu] = 1
                H_flip = hamiltonian_path_count(T_flip)

                # Contractions
                T_con = contract_edge(T, eu, ev)
                H_con = ham_paths_dp(T_con, n - 1)

                T_flip_con = contract_edge(T_flip, ev, eu)
                H_flip_con = ham_paths_dp(T_flip_con, n - 1)

                total += 1
                if H_T - H_flip == H_con - H_flip_con:
                    holds += 1

                # Check if contractions are tournaments
                te_total += 2
                if is_tournament(T_con, n - 1):
                    te_tournament += 1
                if is_tournament(T_flip_con, n - 1):
                    te_tournament += 1

                # Check if T/e and T'/e' are the same tournament
                if T_con == T_flip_con:
                    te_equals_flip_te += 1

                tested_this += 1

    print(f"\nn={n} ({sample_type}):")
    print(f"  Arc-flip reduction: {holds}/{total} = {100*holds/total:.1f}%")
    print(f"  T/e is tournament: {te_tournament}/{te_total} = {100*te_tournament/te_total:.1f}%")
    print(f"  T/e == T'/e': {te_equals_flip_te}/{total} = {100*te_equals_flip_te/total:.1f}%")

# ============================================================
# Deeper: what IS the relationship between T/e and T'/e'?
# ============================================================
print("\n" + "=" * 70)
print("STRUCTURE OF T/e vs T'/e'")
print("=" * 70)

n = 5
T = [[0]*5 for _ in range(5)]
for i in range(5):
    for d in [1, 2]:
        T[i][(i + d) % 5] = 1

print(f"\nCyclic T_5:")
for u in range(5):
    for v in range(u + 1, 5):
        if T[u][v]:
            eu, ev = u, v
        else:
            eu, ev = v, u

        T_con = contract_edge(T, eu, ev)
        T_flip = [row[:] for row in T]
        T_flip[eu][ev] = 0
        T_flip[ev][eu] = 1
        T_flip_con = contract_edge(T_flip, ev, eu)

        # Compare the two n-1 tournaments
        diff_edges = []
        for a in range(4):
            for b in range(4):
                if a != b and T_con[a][b] != T_flip_con[a][b]:
                    diff_edges.append((a, b))

        H_con = ham_paths_dp(T_con, 4)
        H_fcon = ham_paths_dp(T_flip_con, 4)

        print(f"  edge ({eu}->{ev}): H(T/e)={H_con}, H(T'/e')={H_fcon}, "
              f"diff_H={H_con - H_fcon}, diff_edges={len(diff_edges)//2} pairs")

# ============================================================
# Even deeper: iterated contraction reduction
# ============================================================
print("\n" + "=" * 70)
print("ITERATED CONTRACTION: reduce to n=3")
print("=" * 70)

# Can we iterate? H(T) - H(T') = H(T/e) - H(T'/e')
# Both T/e and T'/e' are (n-1)-vertex tournaments (when they are tournaments).
# If we flip another edge in T/e, do we get a further reduction?

n = 5
T = [[0]*5 for _ in range(5)]
for i in range(5):
    for d in [1, 2]:
        T[i][(i + d) % 5] = 1

# Level 0: T_5 (cyclic), H=15
eu, ev = 0, 1  # edge to flip
T_flip = [row[:] for row in T]
T_flip[eu][ev] = 0
T_flip[ev][eu] = 1
H_T = hamiltonian_path_count(T)
H_flip = hamiltonian_path_count(T_flip)
print(f"Level 0: H(T)={H_T}, H(T')={H_flip}, diff={H_T - H_flip}")

# Contract to level 1
T1a = contract_edge(T, eu, ev)
T1b = contract_edge(T_flip, ev, eu)
H1a = ham_paths_dp(T1a, 4)
H1b = ham_paths_dp(T1b, 4)
print(f"Level 1: H(T/e)={H1a}, H(T'/e')={H1b}, diff={H1a - H1b}")

# Now flip an edge in T1a, and contract again
if is_tournament(T1a, 4) and is_tournament(T1b, 4):
    # Find an edge that differs between T1a and T1b
    for a in range(4):
        for b in range(a + 1, 4):
            if T1a[a][b] != T1b[a][b]:
                print(f"  T1a and T1b differ at undirected pair ({a},{b})")

    # Pick first differing edge for further contraction
    for a in range(4):
        for b in range(a + 1, 4):
            if T1a[a][b] != T1b[a][b]:
                if T1a[a][b]:
                    ea, eb = a, b
                else:
                    ea, eb = b, a

                T2a = contract_edge(T1a, ea, eb)
                # For T1b, the edge goes the other way
                T2b = contract_edge(T1b, eb, ea)
                H2a = ham_paths_dp(T2a, 3)
                H2b = ham_paths_dp(T2b, 3)
                print(f"Level 2 (via ({ea},{eb})): H(T1a/e)={H2a}, H(T1b/e')={H2b}, diff={H2a - H2b}")

                # Check reduction still works
                # H(T1a) - H(T1b) should = H(T1a\f) + ... but that's deletion
                T1a_del = [row[:] for row in T1a]
                T1a_del[ea][eb] = 0
                T1b_del = [row[:] for row in T1b]
                T1b_del[eb][ea] = 0

                H1a_del = ham_paths_dp(T1a_del, 4)
                H1b_del = ham_paths_dp(T1b_del, 4)
                print(f"  H(T1a\\f)={H1a_del}, H(T1b\\f')={H1b_del}")
                print(f"  Check: H1a - H1b = {H1a - H1b}, (H1a\\f - H1b\\f') + (T2a - T2b) = {H1a_del - H1b_del} + {H2a - H2b} = {H1a_del - H1b_del + H2a - H2b}")
                break
        else:
            continue
        break
else:
    print("  T1a or T1b is not a tournament — cannot iterate")

# ============================================================
# CRITICAL TEST: Is T/e always a tournament when e is a T-edge?
# ============================================================
print("\n" + "=" * 70)
print("WHEN IS T/e A TOURNAMENT?")
print("=" * 70)

n = 4
m = n * (n - 1) // 2
tournament_count = 0
not_tournament_count = 0
details = defaultdict(int)

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    for u in range(n):
        for v in range(n):
            if u == v or not T[u][v]:
                continue
            Tc = contract_edge(T, u, v)
            t = is_tournament(Tc, n - 1)
            if t:
                tournament_count += 1
            else:
                not_tournament_count += 1
                # Why not? Find the problematic pair
                for a in range(n - 1):
                    for b in range(a + 1, n - 1):
                        if Tc[a][b] + Tc[b][a] != 1:
                            details[(Tc[a][b], Tc[b][a])] += 1

print(f"n={n}: Tournament={tournament_count}, Not={not_tournament_count}")
print(f"Failure types: {dict(details)}")
print(f"  (0,0) means BOTH edges missing, (1,1) means BOTH edges present")

# Now check: for which edges (u,v) is T/e a tournament?
# The merged vertex w has: (x,w) iff (x,u), and (w,x) iff (v,x)
# For pair (x, y) where x,y != u,v: both original edges, fine.
# For pair (w, x): (w,x) = T[v][x], (x,w) = T[x][u].
#   Tournament iff T[v][x] + T[x][u] = 1.
#   T[v][x] + T[x][u] = T[v][x] + (1 - T[u][x]) = 1 + T[v][x] - T[u][x]
#   This equals 1 iff T[v][x] = T[u][x].
#   So T/e is a tournament iff u and v have identical out-neighborhoods
#   (restricted to other vertices).

print(f"\n  T/e is tournament iff u,v have SAME out-neighborhood on V\\{{u,v}}")
print(f"  i.e., for all x not in {{u,v}}: T[u][x] = T[v][x]")
print(f"  = u and v are 'equivalent' modulo the edge between them")

# Verify this characterization
n = 5
m = n * (n - 1) // 2
correct = 0
incorrect = 0

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    for u in range(n):
        for v in range(n):
            if u == v or not T[u][v]:
                continue
            # Predicted: T/e is tournament iff u,v have same out-nbrs
            same_out = all(T[u][x] == T[v][x] for x in range(n) if x != u and x != v)
            Tc = contract_edge(T, u, v)
            actual = is_tournament(Tc, n - 1)
            if same_out == actual:
                correct += 1
            else:
                incorrect += 1

print(f"  Verification at n={n}: correct={correct}, incorrect={incorrect}")

print("\nDone.")
