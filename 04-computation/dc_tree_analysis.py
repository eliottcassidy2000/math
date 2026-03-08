import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
dc_tree_analysis.py
kind-pasteur-2026-03-07-S39b

DELETION-CONTRACTION TREE ANALYSIS:
How does the DC decomposition H(T) = H(T\e) + H(T/e) relate to OCF?

Key insight: T/e IS a tournament (n-1 vertices), so we can recurse.
After contracting (n-1) edges, we reach a 1-vertex tournament with H=1.
After n-1 contractions, we're at 1 vertex.

The DC tree has 2^m leaves (m = n(n-1)/2), each a highly contracted
or deleted minor. The sum of H values at leaves = H(T).

QUESTION 1: Is the contraction T/e always a tournament?
QUESTION 2: If we ONLY contract (never delete), what do we get?
QUESTION 3: How does the OCF decomposition relate to the DC tree?

For the contraction-only branch: after contracting edges e_1, ..., e_k,
we get a tournament on n-k vertices. H(T) = H(T\e) + H(T/e), and
iterating: H(T) = sum over subsets S of "contraction paths" contributions.

FINDING (from THM-082): The contraction convention is:
  w inherits IN-edges from u (tail), OUT-edges from v (head).
  So T/e IS a tournament on n-1 vertices.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count, find_odd_cycles
from itertools import combinations, permutations
from collections import defaultdict, Counter


def ham_paths_dp(adj, n):
    """Count Hamiltonian paths in a general digraph using DP."""
    if n == 0:
        return 1
    if n == 1:
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
    w inherits IN from u, OUT from v.
    Result IS a tournament on n-1 vertices."""
    n = len(T)
    verts = [i for i in range(n) if i != v]
    m = len(verts)
    D = [[0]*m for _ in range(m)]
    for i_idx, i in enumerate(verts):
        for j_idx, j in enumerate(verts):
            if i_idx == j_idx:
                continue
            if i == u:
                # merged -> j: use v's out-edges
                D[i_idx][j_idx] = T[v][j]
            elif j == u:
                # i -> merged: use u's in-edges
                D[i_idx][j_idx] = T[i][u]
            else:
                D[i_idx][j_idx] = T[i][j]
    return D


def delete_edge(T, u, v):
    """Delete edge (u->v). Result is NOT a tournament."""
    n = len(T)
    D = [row[:] for row in T]
    D[u][v] = 0
    return D


def is_tournament(adj, n):
    """Check if adjacency matrix is a tournament."""
    for i in range(n):
        for j in range(i+1, n):
            if adj[i][j] + adj[j][i] != 1:
                return False
    return True


def get_edges(T):
    """Get all directed edges of tournament T."""
    n = len(T)
    return [(i, j) for i in range(n) for j in range(n) if i != j and T[i][j]]


# ============================================================
# Q1: Is T/e always a tournament?
# ============================================================
print("=" * 70)
print("Q1: Is T/e always a tournament?")
print("=" * 70)

for n in range(3, 7):
    m = n * (n - 1) // 2
    all_tournament = True
    tested = 0
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        edges = get_edges(T)
        for u, v in edges:
            Tc = contract_edge(T, u, v)
            tested += 1
            if not is_tournament(Tc, n - 1):
                all_tournament = False
                print(f"  FAIL at n={n}, bits={bits}, e=({u},{v})")
                break
        if not all_tournament:
            break
    if all_tournament:
        print(f"  n={n}: T/e is ALWAYS a tournament ({tested} tests)")

# ============================================================
# Q2: Contraction-only path analysis
# ============================================================
print("\n" + "=" * 70)
print("Q2: Contraction-only analysis")
print("=" * 70)

# For a tournament T, the DC formula is H(T) = H(T\e) + H(T/e).
# If we keep choosing to contract (right branch), we get a chain:
# T -> T/e1 -> T/e1/e2 -> ... -> T1 (single vertex, H=1)
# But which edges do we contract? And how does H change?

# For n=5, examine the contraction tree for one tournament
n = 5
T5c = [[0]*5 for _ in range(5)]
for i in range(5):
    for d in [1, 2]:
        T5c[i][(i + d) % 5] = 1

H5 = hamiltonian_path_count(T5c)
print(f"\nCyclic T_5: H = {H5}")
print(f"Testing contraction along first available edge at each step:")

T = [row[:] for row in T5c]
curr_n = 5
step = 0
while curr_n > 1:
    edges = get_edges(T)
    u, v = edges[0]  # always pick first edge
    T_del = delete_edge(T, u, v)
    T_con = contract_edge(T, u, v)
    H_del = ham_paths_dp(T_del, curr_n)
    H_con = ham_paths_dp(T_con, curr_n - 1)
    H_cur = ham_paths_dp(T, curr_n)
    print(f"  Step {step}: n={curr_n}, e=({u},{v}), H={H_cur} = H(T\\e)={H_del} + H(T/e)={H_con}")
    T = T_con
    curr_n -= 1
    step += 1

# ============================================================
# Q3: Full DC tree — which leaves contribute?
# ============================================================
print("\n" + "=" * 70)
print("Q3: Full DC tree analysis at n=4")
print("=" * 70)

# At n=4, we have m=6 edges. Full DC tree has 2^6 = 64 leaves.
# H(T) = sum over leaves (weighted by +-1) of H(minor).
# But since H = H(del) + H(con), the sum is positive throughout.

n = 4
for bits in range(1 << 6):
    T = tournament_from_bits(n, bits)
    H_T = hamiltonian_path_count(T)

    # Enumerate all paths through the DC tree
    # Edge ordering: list directed edges
    edges = get_edges(T)
    m_edges = len(edges)  # This is n(n-1), not n(n-1)/2

    # For DC tree, we process UNDIRECTED edge positions
    # At each step, pick the first available directed edge
    # Actually, let's use undirected edges for the tree structure
    undirected = []
    seen_pairs = set()
    for u, v in edges:
        pair = (min(u, v), max(u, v))
        if pair not in seen_pairs:
            seen_pairs.add(pair)
            undirected.append((u, v) if T[u][v] else (v, u))  # directed version

    # Build the DC tree: process edges in order
    # At each node, we either delete or contract the next edge
    # A "path" through the tree is a binary string of length m (=6)
    # 0 = delete, 1 = contract

    leaf_h_sum = 0
    nonzero_leaves = 0
    leaf_contributions = defaultdict(int)

    for path_bits in range(1 << len(undirected)):
        # Apply operations in order
        current_adj = [row[:] for row in T]
        current_n = n
        valid = True
        contracted = 0

        # We need to track vertex relabeling through contractions
        # This is complex. Let's use a simpler approach:
        # Build a sequence of operations and apply them.
        # Actually, the edges change after each contraction, so this is tricky.
        # Instead, let's use the identity recursively.
        pass

    if bits < 2:
        print(f"\n  bits={bits}, H(T)={H_T}")
        # Show just the first level DC
        for u, v in undirected[:3]:
            T_del = delete_edge(T, u, v)
            T_con = contract_edge(T, u, v)
            H_del = ham_paths_dp(T_del, n)
            H_con = ham_paths_dp(T_con, n - 1)
            print(f"    e=({u},{v}): H(T\\e)={H_del} + H(T/e)={H_con} = {H_del + H_con}")

# ============================================================
# Q4: OCF on T/e — does contraction preserve OCF?
# ============================================================
print("\n" + "=" * 70)
print("Q4: Does OCF hold for T/e (contracted tournament)?")
print("=" * 70)

for n in range(3, 7):
    m = n * (n - 1) // 2
    ocf_holds = 0
    ocf_total = 0

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        edges_dir = [(u, v) for u in range(n) for v in range(n) if u != v and T[u][v]]

        # Pick first edge
        u, v = edges_dir[0]
        Tc = contract_edge(T, u, v)

        # Check OCF for Tc
        H_con = hamiltonian_path_count(Tc)  # This uses the tournament_lib DP

        # Compute I(Omega(Tc), 2)
        cycles = find_odd_cycles(Tc)
        cycle_vsets = [frozenset(c) for c in cycles]
        m_cyc = len(cycle_vsets)

        # Build conflict graph adjacency
        adj_bits = [0] * m_cyc
        for a in range(m_cyc):
            for b in range(a + 1, m_cyc):
                if cycle_vsets[a] & cycle_vsets[b]:
                    adj_bits[a] |= 1 << b
                    adj_bits[b] |= 1 << a

        # Evaluate I(Omega, 2)
        I_val = 1  # empty set
        for mask in range(1, 1 << m_cyc):
            verts_in = [i for i in range(m_cyc) if mask & (1 << i)]
            is_indep = True
            for a in range(len(verts_in)):
                for b in range(a + 1, len(verts_in)):
                    if adj_bits[verts_in[a]] & (1 << verts_in[b]):
                        is_indep = False
                        break
                if not is_indep:
                    break
            if is_indep:
                I_val += 2 ** len(verts_in)

        ocf_total += 1
        if H_con == I_val:
            ocf_holds += 1

    print(f"  n={n}: OCF holds for T/e in {ocf_holds}/{ocf_total} = {100*ocf_holds/ocf_total:.1f}%")

# ============================================================
# Q5: Arc-flip reduction via DC
# ============================================================
print("\n" + "=" * 70)
print("Q5: Arc-flip reduction: H(T) - H(T') = H(T/e) - H(T'/e')")
print("=" * 70)
print("(From THM-082 Remark 6)")

n = 5
m = n * (n - 1) // 2
flip_reduces = 0
flip_total = 0

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    H_T = hamiltonian_path_count(T)

    # For each edge, flip it and check the contraction reduction
    for u in range(n):
        for v in range(u + 1, n):
            if T[u][v]:
                eu, ev = u, v
            else:
                eu, ev = v, u

            # T has edge eu->ev
            # Flip: T' has edge ev->eu
            T_flip = [row[:] for row in T]
            T_flip[eu][ev] = 0
            T_flip[ev][eu] = 1
            H_flip = hamiltonian_path_count(T_flip)

            # T/e: contract eu->ev
            T_con = contract_edge(T, eu, ev)
            H_con = ham_paths_dp(T_con, n - 1)

            # T'/e': contract ev->eu
            T_flip_con = contract_edge(T_flip, ev, eu)
            H_flip_con = ham_paths_dp(T_flip_con, n - 1)

            # Check: H(T) - H(T') = H(T/e) - H(T'/e')
            lhs = H_T - H_flip
            rhs = H_con - H_flip_con

            flip_total += 1
            if lhs == rhs:
                flip_reduces += 1

print(f"  n={n}: Arc-flip reduction holds in {flip_reduces}/{flip_total} = {100*flip_reduces/flip_total:.1f}%")

# ============================================================
# Q6: What is T/e in terms of the original T?
# ============================================================
print("\n" + "=" * 70)
print("Q6: Structure of T/e — how does score sequence change?")
print("=" * 70)

n = 5
T5c = [[0]*5 for _ in range(5)]
for i in range(5):
    for d in [1, 2]:
        T5c[i][(i + d) % 5] = 1

scores = sorted(sum(T5c[i]) for i in range(5))
print(f"T_5 scores: {scores}")

edges_dir = [(u, v) for u in range(5) for v in range(5) if u != v and T5c[u][v]]
print(f"\nContracting each edge of T_5:")
for u, v in edges_dir[:5]:  # first 5 directed edges
    Tc = contract_edge(T5c, u, v)
    scores_c = sorted(sum(Tc[i]) for i in range(4))
    H_c = hamiltonian_path_count(Tc)
    c3_c = sum(1 for a in range(4) for b in range(a+1,4) for c in range(b+1,4)
              if (Tc[a][b] and Tc[b][c] and Tc[c][a]) or
                 (Tc[a][c] and Tc[c][b] and Tc[b][a]))
    print(f"  Contract ({u}->{v}): scores={scores_c}, H={H_c}, c3={c3_c}")

# ============================================================
# Q7: OCF difference: how does I(Omega, 2) change under contraction?
# ============================================================
print("\n" + "=" * 70)
print("Q7: Cycle analysis under contraction")
print("=" * 70)

n = 5
T = T5c
cycles_T = find_odd_cycles(T)
print(f"T_5 cycles: {len(cycles_T)} odd cycles")
for c in sorted(set(tuple(sorted(c)) for c in cycles_T)):
    count = sum(1 for cc in cycles_T if tuple(sorted(cc)) == c)
    print(f"  vertex set {c}: {count} directed cycle(s)")

u, v = 0, 1  # Contract edge 0->1
Tc = contract_edge(T, u, v)
cycles_Tc = find_odd_cycles(Tc)
print(f"\nT_5 / (0->1) cycles: {len(cycles_Tc)} odd cycles")
for c in sorted(set(tuple(sorted(c)) for c in cycles_Tc)):
    count = sum(1 for cc in cycles_Tc if tuple(sorted(cc)) == c)
    print(f"  vertex set {c}: {count} directed cycle(s)")

# Show the relationship
H_T = hamiltonian_path_count(T)
H_del = ham_paths_dp(delete_edge(T, u, v), n)
H_con = hamiltonian_path_count(Tc)
print(f"\nH(T) = {H_T} = H(T\\e) + H(T/e) = {H_del} + {H_con}")

print("\nDone.")
