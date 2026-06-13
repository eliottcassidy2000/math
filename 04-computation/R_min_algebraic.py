#!/usr/bin/env python3
"""
Algebraic analysis of R-minimization.

We proved: R(T) = n - (sum_S 2^|S| * |U(S)|) / (sum_S 2^|S|)
         = n - U_sum / H

where U_sum = sum_S 2^|S| * |U(S)| = sum_v delta_v = sum_v [H(T) - H(T-v)].

By Claim A: delta_v = H(T) - H(T-v) = 2 * sum_{C through v} mu(C).

So U_sum = 2 * sum_v sum_{C through v} mu(C)
         = 2 * sum_C |V(C)| * mu(C)    [each C counted |V(C)| times]

And R(T) = n - (2/H) * sum_C |V(C)| * mu(C)

Minimizing R <=> maximizing (1/H) * sum_C |V(C)| * mu(C)
<=> maximizing sum_C |V(C)| * mu(C) / H.

For the empty-cycle tournament (H=1, no cycles): R = n.
For tournaments with many cycles: R → small.

Can we express sum_C |V(C)| * mu(C) in terms of the independence polynomial?

H = I(Omega, 2) = sum_k alpha_k * 2^k
U_sum = sum_S 2^|S| * |U(S)|

For |S|=0: contributes 1 to H, 0 to U_sum.
For |S|=1 (single cycle C_i): contributes 2 to H, 2*|V(C_i)| to U_sum.
For |S|=2 (pair {C_i, C_j}, disjoint): contributes 4 to H, 4*|V(C_i) ∪ V(C_j)| to U_sum.

Since cycles in an independent set are vertex-disjoint:
|U(S)| = sum_{C in S} |V(C)|.

So U_sum = sum_S 2^|S| * sum_{C in S} |V(C)|
         = sum_C sum_{S containing C, independent} 2^|S| * |V(C)|
         = sum_C |V(C)| * sum_{S containing C, independent} 2^|S|

Let f(C) = sum_{S containing C, S independent} 2^|S|.
This is 2 * I(Omega[N̄(C)], 2) where N̄(C) = vertices not adjacent to C in Omega
(i.e., cycles vertex-disjoint from C). The factor 2 comes from the C-contribution
to |S|.

Actually: f(C) = 2 * I(Omega - N[C], 2) where N[C] = {C} ∪ neighbors of C in Omega.
Omega - N[C] is the subgraph of cycles disjoint from C (both vertex-disjoint AND
non-conflicting in Omega, but for Omega the adjacency IS vertex-sharing, so
non-adjacent means vertex-disjoint).

So f(C) = 2 * H(T[V \ V(C)]) = 2 * mu(C).

Therefore: U_sum = sum_C |V(C)| * 2 * mu(C) = 2 * sum_C |V(C)| * mu(C).

This confirms U_sum = 2 * sum_C |V(C)| * mu(C) from both the Claim A route
and the independence polynomial route! Consistency check.

Now: R = n - 2*sum_C |V(C)| * mu(C) / H(T).

By Claim A: delta_v = 2 * sum_{C through v} mu(C).
By OCF: H(T) = 1 + sum_{v's contribution from cycles through v via recursion}

The question: can we prove sum_C |V(C)| * mu(C) / H(T) is maximized when H is maximal?

Let's look at this ratio for various tournaments.

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

# ============================================================
# For n=5,6: compute sum_C |V(C)|*mu(C) and compare with H
# ============================================================
print("=" * 70)
print("CYCLE-WEIGHTED MU SUM vs H")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    max_h = MAX_H[n]

    from collections import defaultdict
    by_h = defaultdict(list)

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h == 0:
            continue

        cycles = find_odd_cycles(T)
        if not cycles:
            by_h[h].append((0, 0))
            continue

        # Compute mu(C) for each cycle
        weighted_sum = 0
        total_mu = 0
        for c in cycles:
            # mu(C) = H(T[V \ V(C)])
            remaining = [v for v in range(n) if v not in c]
            if not remaining:
                mu_c = 1
            else:
                sub = [[T[remaining[i]][remaining[j]]
                        for j in range(len(remaining))]
                       for i in range(len(remaining))]
                mu_c = hamiltonian_path_count(sub)
            weighted_sum += len(c) * mu_c
            total_mu += mu_c

        by_h[h].append((weighted_sum, total_mu))

    print(f"\nn={n}:")
    print(f"  H | count | avg sum_C|V|*mu | avg sum_mu | avg sum_C|V|*mu/H")
    for h in sorted(by_h.keys(), reverse=True)[:8]:
        items = by_h[h]
        avg_w = sum(w for w, m in items) / len(items)
        avg_m = sum(m for w, m in items) / len(items)
        avg_ratio = avg_w / h
        marker = " *** MAX" if h == max_h else ""
        print(f"  {h:>3} | {len(items):>5} | {avg_w:>14.1f} | {avg_m:>10.1f} | {avg_ratio:>8.4f}{marker}")

# ============================================================
# INDEPENDENCE POLYNOMIAL STRUCTURE of maximizers
# ============================================================
print(f"\n{'='*70}")
print("INDEPENDENCE POLYNOMIAL COEFFICIENTS OF MAXIMIZERS")
print("=" * 70)

def indep_poly_coeffs(adj):
    m = len(adj)
    if m == 0:
        return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
    coeffs = [0] * (m + 1)
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
            coeffs[bin(mask).count('1')] += 1
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

for n in [5, 6, 7]:
    m = n * (n - 1) // 2
    max_h = MAX_H[n]

    seen_ips = set()
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h != max_h:
            continue

        cycles = find_odd_cycles(T)
        cg = conflict_graph(cycles)
        ip = tuple(indep_poly_coeffs(cg))

        if ip not in seen_ips:
            seen_ips.add(ip)
            s = score_seq(T)
            # Verify H = I(Omega, 2)
            h_check = sum(c * (2**k) for k, c in enumerate(ip))
            print(f"  n={n}: score={s}, IP={list(ip)}, I(Omega,2)={h_check}")

            # Compute weighted sum from IP
            # U_sum/H = n - R
            # But can we compute U_sum from IP alone?
            # U_sum = sum_S 2^|S| * |U(S)|
            # This depends on the SIZES of cycles in each independent set,
            # not just the number of independent sets. So IP alone isn't enough.

# ============================================================
# THE WEIGHTED IP: sum_S 2^|S| * |U(S)| by independent set size
# ============================================================
print(f"\n{'='*70}")
print("WEIGHTED IP BREAKDOWN (by |S| size)")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    max_h = MAX_H[n]

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h != max_h:
            continue

        cycles = find_odd_cycles(T)
        cg = conflict_graph(cycles)

        # Enumerate all independent sets and compute |U(S)|
        adj = cg
        m_cycles = len(adj)
        nbr = [0] * m_cycles
        for i in range(m_cycles):
            for j in range(m_cycles):
                if adj[i][j]:
                    nbr[i] |= 1 << j

        # Group by |S|: total 2^|S| and total 2^|S|*|U(S)|
        from collections import defaultdict
        h_by_k = defaultdict(int)
        u_by_k = defaultdict(int)

        for mask in range(1 << m_cycles):
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
            if not ok:
                continue

            k = bin(mask).count('1')
            U = set()
            temp = mask
            while temp:
                v = (temp & -temp).bit_length() - 1
                U.update(cycles[v])
                temp &= temp - 1

            h_by_k[k] += 2**k
            u_by_k[k] += 2**k * len(U)

        s = score_seq(T)
        print(f"\nn={n}, score={s}, H={h}:")
        print(f"  |S| | count*2^|S| | count*2^|S|*|U| | avg |U|")
        for k in sorted(h_by_k.keys()):
            hk = h_by_k[k]
            uk = u_by_k[k]
            alpha_k = hk // (2**k) if k > 0 else 1
            avg_u = uk / hk if hk > 0 else 0
            print(f"  {k:>3} | {hk:>11} | {uk:>15} | {avg_u:>7.3f}")

        print(f"  Total H = {sum(h_by_k.values())}")
        print(f"  Total U_sum = {sum(u_by_k.values())}")
        print(f"  R = {n - sum(u_by_k.values())/sum(h_by_k.values()):.4f}")

        break  # Just one per n

print("\nDone.")
