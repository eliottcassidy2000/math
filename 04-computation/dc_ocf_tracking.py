import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
dc_ocf_tracking.py
kind-pasteur-2026-03-07-S39b

DELETION-CONTRACTION vs OCF TRACKING

For a tournament T with edge e = (u->v):
  H(T) = H(T\e) + H(T/e)    [THM-082, PROVED]

OCF says H(T) = I(Omega(T), 2). Question: how do cycles in Omega change?

Key decomposition:
  Omega(T) has odd cycles. Under deletion of edge e:
  - Cycles using e lose one edge; some may no longer be cycles
  - Cycles not using e are unchanged

Under contraction of edge e:
  - The merged digraph T/e has different cycle structure
  - But T/e IS a tournament (when u,v have same neighborhoods)

APPROACH: Track I(Omega(T),2), I(Omega(T\e),2), I(Omega(T/e),2)
and see if there's a clean relationship.

Note: T\e is NOT a tournament, so OCF may not hold for it.
But we can still compute I(Omega(T\e), 2) and H(T\e) separately.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count, find_odd_cycles
from itertools import combinations, permutations
from collections import defaultdict


def ham_paths_dp(adj, n):
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


def find_odd_cycles_general(adj, n):
    """Find all directed odd cycles in a general digraph.
    Returns list of tuples (one per directed cycle, vertex-set canonical)."""
    results = []
    for k in range(3, n + 1, 2):
        for verts in combinations(range(n), k):
            first = verts[0]
            for perm in permutations(verts[1:]):
                path = (first,) + perm
                if all(adj[path[i]][path[(i+1) % k]] for i in range(k)):
                    results.append(frozenset(verts))
    return results


def independence_poly_at_2(cycle_vsets):
    """Compute I(Omega, 2) from list of cycle vertex sets."""
    if not cycle_vsets:
        return 1
    m = len(cycle_vsets)
    # Build adjacency
    adj = [0] * m
    for a in range(m):
        for b in range(a+1, m):
            if cycle_vsets[a] & cycle_vsets[b]:
                adj[a] |= 1 << b
                adj[b] |= 1 << a
    # Enumerate independent sets
    total = 1
    for mask in range(1, 1 << m):
        bits = []
        temp = mask
        while temp:
            b = temp & (-temp)
            bits.append(b.bit_length() - 1)
            temp ^= b
        is_indep = True
        for i in range(len(bits)):
            for j in range(i+1, len(bits)):
                if adj[bits[i]] & (1 << bits[j]):
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            total += 2 ** len(bits)
    return total


def contract_edge(T, u, v):
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


def delete_edge(T, u, v):
    n = len(T)
    D = [row[:] for row in T]
    D[u][v] = 0
    return D


# ============================================================
# Main analysis
# ============================================================
print("=" * 70)
print("DELETION-CONTRACTION vs OCF TRACKING")
print("=" * 70)

for n in range(3, 7):
    m = n * (n - 1) // 2
    print(f"\nn={n} (exhaustive, {1 << m} tournaments)")

    # Track the relationship between I values
    relationships = defaultdict(int)
    ocf_holds_del = 0
    ocf_fails_del = 0
    ocf_holds_con = 0
    ocf_fails_con = 0
    total_edges = 0

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        H_T = hamiltonian_path_count(T)

        # OCF for T
        cycles_T = [frozenset(c) for c in find_odd_cycles(T)]
        I_T = independence_poly_at_2(cycles_T)
        assert I_T == H_T, f"OCF fails for tournament bits={bits}"

        # Pick first directed edge
        for u in range(n):
            for v in range(n):
                if u != v and T[u][v]:
                    eu, ev = u, v
                    break
            else:
                continue
            break

        # Deletion
        T_del = delete_edge(T, eu, ev)
        H_del = ham_paths_dp(T_del, n)
        cycles_del = find_odd_cycles_general(T_del, n)
        I_del = independence_poly_at_2(cycles_del)

        if H_del == I_del:
            ocf_holds_del += 1
        else:
            ocf_fails_del += 1

        # Contraction
        T_con = contract_edge(T, eu, ev)
        H_con = ham_paths_dp(T_con, n - 1)
        cycles_con = find_odd_cycles_general(T_con, n - 1)
        I_con = independence_poly_at_2(cycles_con)

        if H_con == I_con:
            ocf_holds_con += 1
        else:
            ocf_fails_con += 1

        # Track relationship
        # H(T) = H(T\e) + H(T/e) always
        # Does I(T) = I(T\e) + I(T/e)?
        I_sum = I_del + I_con
        diff = I_T - I_sum
        relationships[diff] += 1

        total_edges += 1

    print(f"  OCF for T\\e: {ocf_holds_del}/{total_edges} ({100*ocf_holds_del/total_edges:.1f}%)")
    print(f"  OCF for T/e: {ocf_holds_con}/{total_edges} ({100*ocf_holds_con/total_edges:.1f}%)")
    print(f"  I(T) = I(T\\e) + I(T/e): {relationships.get(0,0)}/{total_edges} ({100*relationships.get(0,0)/total_edges:.1f}%)")
    print(f"  I(T) - [I(T\\e) + I(T/e)] distribution: {dict(sorted(relationships.items()))}")

# ============================================================
# Deeper: when OCF fails for T\e, what's the error?
# ============================================================
print("\n" + "=" * 70)
print("WHEN OCF FAILS FOR T\\e: error analysis")
print("=" * 70)

n = 4
m = n * (n - 1) // 2
errors_del = []
errors_con = []

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)

    for u in range(n):
        for v in range(n):
            if u != v and T[u][v]:
                eu, ev = u, v
                break
        else:
            continue
        break

    T_del = delete_edge(T, eu, ev)
    H_del = ham_paths_dp(T_del, n)
    cycles_del = find_odd_cycles_general(T_del, n)
    I_del = independence_poly_at_2(cycles_del)

    if H_del != I_del:
        errors_del.append({
            'bits': bits, 'e': (eu, ev),
            'H': H_del, 'I': I_del,
            'err': H_del - I_del,
            'num_cycles': len(cycles_del)
        })

    T_con = contract_edge(T, eu, ev)
    H_con = ham_paths_dp(T_con, n - 1)
    cycles_con = find_odd_cycles_general(T_con, n - 1)
    I_con = independence_poly_at_2(cycles_con)

    if H_con != I_con:
        errors_con.append({
            'bits': bits, 'e': (eu, ev),
            'H': H_con, 'I': I_con,
            'err': H_con - I_con,
            'num_cycles': len(cycles_con)
        })

print(f"\nn={n}: OCF errors for T\\e:")
for err in errors_del[:5]:
    print(f"  bits={err['bits']}, e={err['e']}: H={err['H']}, I={err['I']}, "
          f"diff={err['err']}, #cycles={err['num_cycles']}")

err_vals = [e['err'] for e in errors_del]
print(f"  Error values: {sorted(set(err_vals))}")
print(f"  Error always even? {all(e % 2 == 0 for e in err_vals)}")

print(f"\nn={n}: OCF errors for T/e:")
for err in errors_con[:5]:
    print(f"  bits={err['bits']}, e={err['e']}: H={err['H']}, I={err['I']}, "
          f"diff={err['err']}, #cycles={err['num_cycles']}")

err_vals_con = [e['err'] for e in errors_con]
if err_vals_con:
    print(f"  Error values: {sorted(set(err_vals_con))}")
else:
    print(f"  (no errors)")

# ============================================================
# KEY QUESTION: Does OCF hold for ALL digraphs, or just tournaments?
# ============================================================
print("\n" + "=" * 70)
print("OCF FOR GENERAL DIGRAPHS — does H(D) = I(Omega(D), 2)?")
print("=" * 70)

# Grinberg-Stanley proves it for all digraphs D with D-bar = D^op.
# For tournaments, D-bar = D^op. But T\e does NOT have this property!
# T\e has one missing edge, so (T\e)-bar != (T\e)^op.

# Specifically: for T\e, the pair (u,v) has NO edge in either direction.
# In (T\e)-bar (complement): (u,v) has BOTH edges.
# In (T\e)^op (converse): (u,v) has NO edge.
# So (T\e)-bar != (T\e)^op. GS does not apply.

# For T/e: it IS a tournament (when applicable), so GS DOES apply.
# Let's verify this directly.

n = 4
m = n * (n - 1) // 2
te_tournament = 0
te_not_tournament = 0
ocf_when_tournament = 0
ocf_when_not = 0

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    for u in range(n):
        for v in range(n):
            if u != v and T[u][v]:
                T_con = contract_edge(T, u, v)
                H_con = ham_paths_dp(T_con, n - 1)
                cycles_con = find_odd_cycles_general(T_con, n - 1)
                I_con = independence_poly_at_2(cycles_con)

                is_tourn = all(
                    T_con[a][b] + T_con[b][a] == 1
                    for a in range(n-1) for b in range(a+1, n-1)
                )
                if is_tourn:
                    te_tournament += 1
                    if H_con == I_con:
                        ocf_when_tournament += 1
                else:
                    te_not_tournament += 1
                    if H_con == I_con:
                        ocf_when_not += 1
                break
        break

print(f"\nn={n}: T/e analysis:")
print(f"  T/e is tournament: {te_tournament}, OCF holds: {ocf_when_tournament}/{te_tournament}")
print(f"  T/e NOT tournament: {te_not_tournament}, OCF holds: {ocf_when_not}/{te_not_tournament}")

# When T/e IS a tournament, OCF always holds (by GS theorem).
# When T/e is NOT a tournament, OCF may or may not hold.

# So the DC approach for OCF has this structure:
# H(T) = H(T\e) + H(T/e)
# I(Omega(T), 2) = I(Omega(T), 2)  (tautology)
# We want: H(T) = I(Omega(T), 2)
# By induction (if T/e is tournament): H(T/e) = I(Omega(T/e), 2)
# So need: H(T\e) = I(Omega(T), 2) - I(Omega(T/e), 2)

# The question becomes: what is I(Omega(T), 2) - I(Omega(T/e), 2)?

print("\n" + "=" * 70)
print("KEY IDENTITY: I(Omega(T),2) - I(Omega(T/e),2) vs H(T\\e)")
print("=" * 70)

for n in range(3, 7):
    m_bits = n * (n - 1) // 2
    matches = 0
    total = 0
    diffs = defaultdict(int)

    for bits in range(1 << m_bits):
        T = tournament_from_bits(n, bits)

        for u in range(n):
            for v in range(n):
                if u != v and T[u][v]:
                    eu, ev = u, v
                    break
            else:
                continue
            break

        H_T = hamiltonian_path_count(T)
        cycles_T = [frozenset(c) for c in find_odd_cycles(T)]
        I_T = independence_poly_at_2(cycles_T)

        T_del = delete_edge(T, eu, ev)
        H_del = ham_paths_dp(T_del, n)

        T_con = contract_edge(T, eu, ev)
        cycles_con = find_odd_cycles_general(T_con, n - 1)
        I_con = independence_poly_at_2(cycles_con)

        # Key identity test: H(T\e) = I(Omega(T),2) - I(Omega(T/e),2)?
        target = I_T - I_con
        diff = H_del - target
        diffs[diff] += 1
        total += 1
        if diff == 0:
            matches += 1

    print(f"\n  n={n}: H(T\\e) = I(T) - I(T/e): {matches}/{total} ({100*matches/total:.1f}%)")
    print(f"    Diff distribution: {dict(sorted(diffs.items()))}")

print("\nDone.")
