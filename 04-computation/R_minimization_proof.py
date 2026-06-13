#!/usr/bin/env python3
"""
R-Minimization Proof Attempt.

We want to prove: T maximizes H(T) => T minimizes R(T) = sum_v H(T-v)/H(T).

Using OCF: H(T) = I(Omega(T), 2) = sum_k alpha_k * 2^k.

Key identity: sum_v H(T-v) = sum_v I(Omega(T-v), 2).

For each vertex deletion v, Omega(T-v) is the conflict graph on odd cycles
of T-v. These include:
  (a) Cycles of T not through v (survive in T-v)
  (b) Possibly NEW cycles created by deleting v (no, can't happen: T-v is
      a subtournament, its cycles are exactly T's cycles not through v)

Wait: deleting v from T cannot CREATE new cycles. A directed cycle in T-v
must use arcs that all exist in T. So cycles of T-v = cycles of T not
passing through v.

Therefore: Omega(T-v) = Omega(T) restricted to cycles not through v.
This is the INDUCED SUBGRAPH of Omega(T) on vertices (cycles) not through v.

So: H(T-v) = I(Omega(T)[cycles not through v], 2)

And: sum_v H(T-v) = sum_v I(Omega(T) restricted to cycles avoiding v, 2)

Each cycle C avoids n - |V(C)| vertices. So cycle C appears in
exactly n - |V(C)| of the n restricted graphs.

For independent sets: an independent set S in Omega(T) avoids vertex v
iff EVERY cycle in S avoids v. S avoids v iff v not in union(V(C) for C in S).

Let U(S) = union of vertex sets of cycles in S. |U(S)| = total vertices
used by S. S appears in exactly n - |U(S)| restricted graphs
(the deletions of vertices outside U(S)... wait, S appears in the
restricted graph for v iff ALL cycles in S avoid v, i.e. v not in U(S)).

Actually, S contributes to H(T-v) iff S is an independent set of
Omega(T-v). This means all cycles in S exist in T-v (they avoid v)
AND they're pairwise vertex-disjoint in T-v (same as in T since no
new adjacencies are created).

So S contributes to sum_v H(T-v) for exactly (n - |U(S)|) values of v.

Therefore:
sum_v H(T-v) = sum_{indep sets S} 2^{|S|} * (n - |U(S)|)

And H(T) = sum_{indep sets S} 2^{|S|}

So R(T) = sum_v H(T-v) / H(T)
       = [sum_S 2^{|S|} * (n - |U(S)|)] / [sum_S 2^{|S|}]
       = n - [sum_S 2^{|S|} * |U(S)|] / [sum_S 2^{|S|}]
       = n - E[|U(S)|]

where E[|U(S)|] is the WEIGHTED AVERAGE of |U(S)| over independent
sets S, weighted by 2^{|S|}.

So: R(T) = n - E_weighted[|U(S)|]

Minimizing R <=> MAXIMIZING E_weighted[|U(S)|].

The maximizer of H has the largest sum_S 2^{|S|} (denominator) AND
the largest E[|U|] (which means its independent sets use many vertices
on average).

This is a CLEAN FORMULA! Let me verify it computationally.

kind-pasteur-2026-03-06-S18g
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count, find_odd_cycles, conflict_graph

MAX_H = {1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189}

def score_seq(T):
    return tuple(sorted(sum(T[i]) for i in range(len(T))))

def delete_vertex(T, v):
    n = len(T)
    verts = [i for i in range(n) if i != v]
    return [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]

def all_independent_sets(adj):
    """Return list of all independent sets as frozensets of indices."""
    m = len(adj)
    if m == 0:
        return [frozenset()]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
    result = []
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            s = frozenset()
            temp = mask
            while temp:
                v = (temp & -temp).bit_length() - 1
                s = s | {v}
                temp &= temp - 1
            result.append(s)
    return result

# ============================================================
# Verify R = n - E_weighted[|U(S)|]
# ============================================================
print("=" * 70)
print("R-MINIMIZATION: R(T) = n - E_weighted[|U(S)|]")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    max_h = MAX_H[n]

    print(f"\nn={n}:")

    # Check a few tournaments
    checked = 0
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h == 0:
            continue

        cycles = find_odd_cycles(T)
        if not cycles:
            # H should be 1 (no odd cycles)
            ds = sum(hamiltonian_path_count(delete_vertex(T, v)) for v in range(n))
            r_actual = ds / h
            r_formula = n - 0  # E[|U|] = 0 since only empty set
            if abs(r_actual - r_formula) > 0.001 and checked < 3:
                print(f"  bits={bits}: H={h}, R_actual={r_actual:.4f}, R_formula={r_formula:.4f} {'OK' if abs(r_actual - r_formula) < 0.001 else 'MISMATCH'}")
            checked += 1
            continue

        cg = conflict_graph(cycles)
        indep_sets = all_independent_sets(cg)

        # Compute E_weighted[|U(S)|]
        sum_2k_U = 0
        sum_2k = 0
        for S in indep_sets:
            k = len(S)
            # U(S) = union of vertex sets of cycles in S
            U = set()
            for idx in S:
                U.update(cycles[idx])
            sum_2k += 2**k
            sum_2k_U += 2**k * len(U)

        E_U = sum_2k_U / sum_2k
        r_formula = n - E_U

        # Compute R directly
        ds = sum(hamiltonian_path_count(delete_vertex(T, v)) for v in range(n))
        r_actual = ds / h

        match = abs(r_actual - r_formula) < 0.001
        if not match or (h == max_h and checked < 10):
            marker = " *** MAX" if h == max_h else ""
            status = "OK" if match else "MISMATCH"
            print(f"  bits={bits}: H={h}, R_actual={r_actual:.4f}, R_formula={r_formula:.4f} {status}{marker}")
            if not match:
                print(f"    sum_2k={sum_2k}, sum_2k_U={sum_2k_U}, h={h}")

        checked += 1
        if checked >= 200 and h != max_h:
            break

    # Now for ALL maximizers, compute E[|U|]
    print(f"\n  All maximizers E[|U|]:")
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h != max_h:
            continue

        cycles = find_odd_cycles(T)
        cg = conflict_graph(cycles)
        indep_sets = all_independent_sets(cg)

        sum_2k_U = 0
        sum_2k = 0
        for S in indep_sets:
            k = len(S)
            U = set()
            for idx in S:
                U.update(cycles[idx])
            sum_2k += 2**k
            sum_2k_U += 2**k * len(U)

        E_U = sum_2k_U / sum_2k
        r = n - E_U
        s = score_seq(T)
        print(f"    bits={bits}: score={s}, E[|U|]={E_U:.4f}, R={r:.4f}")
        break  # Just one example per score

# ============================================================
# The formula suggests: maximizing H <=> maximizing E_weighted[|U|]
# But that's NOT obvious! H = sum_S 2^{|S|} is the total weight,
# while E[|U|] is the average vertex coverage.
#
# Can we relate these? For the maximizer, both sum_S 2^{|S|} and
# sum_S 2^{|S|} * |U(S)| are maximized.
# ============================================================
print(f"\n{'='*70}")
print("DECOMPOSITION: H = sum_S 2^|S|, R-numerator = sum_S 2^|S| * (n-|U|)")
print("=" * 70)

n = 5
m = n * (n - 1) // 2
max_h = MAX_H[n]

# For each H value, compute the sum_S 2^|S| * |U|
data = []
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    if h == 0:
        continue

    cycles = find_odd_cycles(T)
    if not cycles:
        data.append((h, 0, n * h))
        continue

    cg = conflict_graph(cycles)
    indep_sets = all_independent_sets(cg)

    sum_2k_U = 0
    for S in indep_sets:
        k = len(S)
        U = set()
        for idx in S:
            U.update(cycles[idx])
        sum_2k_U += 2**k * len(U)

    data.append((h, sum_2k_U, n * h - sum_2k_U))

# Group by H, show sum_2k_U range
from collections import defaultdict
by_h = defaultdict(list)
for h, u_sum, r_num in data:
    by_h[h].append((u_sum, r_num))

print(f"\nn={n}: H vs sum_S 2^|S|*|U(S)| and R-numerator=sum_S 2^|S|*(n-|U|)")
for h in sorted(by_h.keys(), reverse=True)[:8]:
    items = by_h[h]
    u_sums = [x[0] for x in items]
    r_nums = [x[1] for x in items]
    marker = " *** MAX" if h == max_h else ""
    print(f"  H={h}: U_sum in [{min(u_sums)},{max(u_sums)}], "
          f"R_num in [{min(r_nums)},{max(r_nums)}]{marker}")

# ============================================================
# KEY: The R-numerator = sum_v H(T-v) = sum_S 2^|S| * (n - |U(S)|)
# For the maximizer, this is the DELETION SUM.
# ============================================================

# Check if U_sum is also maximized by H-maximizer
print(f"\n{'='*70}")
print("IS U_sum = sum_S 2^|S|*|U(S)| ALSO MAXIMIZED BY H-MAXIMIZER?")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    max_h = MAX_H[n]

    max_u_sum = 0
    max_u_h = 0
    max_h_u = 0

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h == 0:
            continue

        cycles = find_odd_cycles(T)
        if not cycles:
            continue

        cg = conflict_graph(cycles)
        indep_sets = all_independent_sets(cg)

        sum_2k_U = 0
        for S in indep_sets:
            k = len(S)
            U = set()
            for idx in S:
                U.update(cycles[idx])
            sum_2k_U += 2**k * len(U)

        if sum_2k_U > max_u_sum:
            max_u_sum = sum_2k_U
            max_u_h = h
        if h == max_h:
            max_h_u = max(max_h_u, sum_2k_U)

    print(f"\nn={n}: max U_sum = {max_u_sum} (at H={max_u_h}), "
          f"H-maximizer U_sum = {max_h_u}")
    print(f"  H-maximizer maximizes U_sum? {max_h_u >= max_u_sum}")

# ============================================================
# So the question reduces to: does the maximizer of I(G,2) also
# maximize sum_S 2^|S| * |U(S)| / I(G,2) ?
#
# Write U_sum = sum_S 2^|S| * |U(S)| = sum_v sum_{S: v in U(S)} 2^|S|
#            = sum_v [H(T) - H(T-v)]
#            = n*H(T) - sum_v H(T-v)
#
# Wait! U_sum = n*H - sum_v H(T-v) = n*H - R_num
# And R = R_num / H = n - U_sum/H
# So U_sum/H = n - R
# Maximizing E[|U|] = U_sum/H <=> minimizing R. Already knew this.
#
# But U_sum = n*H - sum_v H(T-v). And:
# R = sum_v H(T-v) / H = n - U_sum/H
#
# The REAL question: is U_sum = n*H - sum_v H(T-v) maximized by H-max?
# U_sum = sum_v [H(T) - H(T-v)] = sum_v delta_v
# where delta_v = H(T) - H(T-v) = 2 * sum_{C through v} mu(C) by Claim A.
#
# So U_sum = 2 * sum_v sum_{C through v} mu(C)
#          = 2 * sum_C |V(C)| * mu(C)
#          (each cycle C is counted |V(C)| times, once per vertex)
#
# And H(T) = 1 + sum_C 2*mu(C) (by Claim A recursion... actually this
# isn't quite right, the recursion is more complex)
#
# Actually, from OCF: H(T) = sum_S 2^|S| where S ranges over
# independent sets in Omega(T).
# And U_sum = sum_S 2^|S| * |U(S)|.
#
# For |S|=0: contribution to H is 1, to U_sum is 0.
# For |S|=1: contribution to H is 2*m (m cycles), to U_sum is 2*sum_C |V(C)|
# For |S|=2: contribution to H is 4*alpha_2, to U_sum is 4*sum_{disjoint pairs} |U(C1,C2)|
#
# The ratio U_sum/H = E[|U|] weights larger independent sets (with 2^|S|)
# more heavily. The maximizer has more large independent sets, which
# use more vertices, so E[|U|] is large.
# ============================================================

print(f"\n{'='*70}")
print("DELTA_v = H(T) - H(T-v) = 2*sum_C mu(C) ANALYSIS")
print("=" * 70)

for n in [5, 6, 7]:
    m = n * (n - 1) // 2
    max_h = MAX_H[n]

    found = False
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h != max_h:
            continue

        del_hs = [hamiltonian_path_count(delete_vertex(T, v)) for v in range(n)]
        deltas = [h - dh for dh in del_hs]
        u_sum = sum(deltas)
        r = sum(del_hs) / h

        print(f"\nn={n}: H={h}, deltas={deltas}")
        print(f"  sum deltas = {u_sum} = n*H - sum H(T-v) = {n*h} - {sum(del_hs)}")
        print(f"  E[|U|] = {u_sum/h:.4f}, R = {r:.4f}")
        print(f"  Check: E[|U|] + R = {u_sum/h + r:.4f} (should be {n})")

        found = True
        break

print("\nDone.")
