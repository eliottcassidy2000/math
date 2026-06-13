#!/usr/bin/env python3
"""sc_maximizer_orbit_fixed.py -- Fixed complementary pair analysis.

Session: kind-pasteur-2026-03-20-S1

BUG FIX: frozenset comparison s < comp uses SUBSET ordering, not lexicographic.
Two disjoint sets of equal size are never subsets, so s < comp always False.
Fix: use tuple(sorted(s)) < tuple(sorted(comp)) for lexicographic comparison.
"""

from itertools import combinations, permutations
from collections import Counter, defaultdict

def adj_from_bits(bits, n):
    adj = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    return adj

def held_karp_H(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))

def find_involutory_anti_auts(adj, n):
    result = []
    for perm in permutations(range(n)):
        if not all(perm[perm[i]] == i for i in range(n)):
            continue
        is_anti = True
        for i in range(n):
            for j in range(i+1, n):
                if adj[perm[i]][perm[j]] != (1 - adj[i][j]):
                    is_anti = False
                    break
            if not is_anti:
                break
        if is_anti:
            result.append(perm)
    return result

def count_all_directed_odd_cycles(adj, n):
    """Return list of (vertex_frozenset, cycle_tuple) for ALL directed odd cycles."""
    cycles = []
    for size in range(3, n+1, 2):
        for verts in combinations(range(n), size):
            sub_adj = [[adj[verts[i]][verts[j]] for j in range(size)] for i in range(size)]
            for perm in permutations(range(1, size)):
                path = [0] + list(perm)
                is_cycle = True
                for idx in range(size):
                    if not sub_adj[path[idx]][path[(idx+1) % size]]:
                        is_cycle = False
                        break
                if is_cycle:
                    actual = tuple(verts[path[i]] for i in range(size))
                    cycles.append((frozenset(verts), actual))
    return cycles

def independence_polynomial(adj, n):
    """Compute independence polynomial of Omega(T) correctly."""
    all_cycles = count_all_directed_odd_cycles(adj, n)
    num_cycles = len(all_cycles)
    if num_cycles == 0:
        return [1]

    vsets = [c[0] for c in all_cycles]
    conflict = [[False]*num_cycles for _ in range(num_cycles)]
    for i in range(num_cycles):
        for j in range(i+1, num_cycles):
            if vsets[i] & vsets[j]:
                conflict[i][j] = True
                conflict[j][i] = True

    alpha = [1]
    for k in range(1, num_cycles + 1):
        count = 0
        for subset in combinations(range(num_cycles), k):
            independent = True
            for a in range(len(subset)):
                for b in range(a+1, len(subset)):
                    if conflict[subset[a]][subset[b]]:
                        independent = False
                        break
                if not independent:
                    break
            if independent:
                count += 1
        if count == 0:
            break
        alpha.append(count)
    return alpha

def eval_ip(alpha, x):
    return sum(a * x**k for k, a in enumerate(alpha))

def main():
    print("=" * 70)
    print("SC MAXIMIZER — FIXED COMPLEMENTARY PAIR ANALYSIS (n=6)")
    print("=" * 70)

    n = 6
    m = n * (n - 1) // 2
    target = (2, 2, 2, 3, 3, 3)

    all_data = []
    count = 0
    for bits in range(1 << m):
        adj = adj_from_bits(bits, n)
        sc = score_seq(adj, n)
        if sc != target:
            continue

        H = held_karp_H(adj, n)

        # Find cyclic 3-subsets
        cyclic = set()
        for triple in combinations(range(n), 3):
            i, j, k = triple
            if (adj[i][j] and adj[j][k] and adj[k][i]) or \
               (adj[i][k] and adj[k][j] and adj[j][i]):
                cyclic.add(frozenset(triple))

        # Count complementary pairs (FIXED: use tuple comparison)
        comp_both = 0
        for triple in combinations(range(n), 3):
            s = frozenset(triple)
            comp = frozenset(range(n)) - s
            # Use lexicographic comparison on sorted tuples to avoid double counting
            if tuple(sorted(s)) < tuple(sorted(comp)):
                if s in cyclic and comp in cyclic:
                    comp_both += 1

        # Find involutory anti-automorphisms
        inv_auts = find_involutory_anti_auts(adj, n)

        # Full independence polynomial
        ip = independence_polynomial(adj, n)

        all_data.append({
            'bits': bits, 'H': H, 'c3': len(cyclic),
            'comp_both': comp_both,
            'has_anti_aut': len(inv_auts) > 0,
            'num_inv_aut': len(inv_auts),
            'IP': tuple(ip),
            'alpha_2': ip[2] if len(ip) > 2 else 0,
        })
        count += 1

    print(f"\n  Total regular n=6 tournaments: {count}")

    # Group by H
    for H_val in sorted(set(d['H'] for d in all_data)):
        entries = [d for d in all_data if d['H'] == H_val]
        has_aa = sum(1 for d in entries if d['has_anti_aut'])
        comp_vals = Counter(d['comp_both'] for d in entries)
        alpha2_vals = Counter(d['alpha_2'] for d in entries)
        ip_vals = Counter(d['IP'] for d in entries)

        print(f"\n  H={H_val}: {len(entries)} tournaments")
        print(f"    Has invol anti-aut: {has_aa}/{len(entries)}")
        print(f"    comp_both (complement pairs both cyclic): {dict(sorted(comp_vals.items()))}")
        print(f"    alpha_2 (from full I.P.): {dict(sorted(alpha2_vals.items()))}")
        print(f"    Full IP: {dict(sorted(ip_vals.items()))}")

    # THEOREM CHECK
    print(f"\n  ===== THEOREM STRUCTURE =====")
    print(f"  comp_both -> alpha_2 -> H mapping:")
    for d in all_data[:1]:  # Just show structure
        pass

    by_comp = defaultdict(list)
    for d in all_data:
        by_comp[d['comp_both']].append(d)

    for cb in sorted(by_comp.keys()):
        entries = by_comp[cb]
        H_vals = Counter(d['H'] for d in entries)
        a2_vals = Counter(d['alpha_2'] for d in entries)
        aa = sum(1 for d in entries if d['has_anti_aut'])
        print(f"  comp_both={cb}: {len(entries)} tours, H={dict(H_vals)}, "
              f"alpha_2={dict(a2_vals)}, anti_aut={aa}/{len(entries)}")

    # KEY QUESTION: Does comp_both determine alpha_2?
    print(f"\n  Does comp_both determine alpha_2?")
    for d in all_data:
        if d['comp_both'] > 0:
            print(f"  bits={d['bits']}: comp_both={d['comp_both']}, alpha_2={d['alpha_2']}, "
                  f"H={d['H']}, IP={d['IP']}")

    # Check: alpha_2 vs comp_both for a sample
    print(f"\n  Sample of (comp_both, alpha_2, H) tuples:")
    for H_val in sorted(set(d['H'] for d in all_data)):
        entries = [d for d in all_data if d['H'] == H_val]
        sample = entries[:3]
        for d in sample:
            print(f"  H={d['H']}: comp_both={d['comp_both']}, alpha_2={d['alpha_2']}")

    # Final analysis: what IS alpha_2 counting?
    # For H=45 tournaments, show explicit independent set structure
    print(f"\n  ===== EXPLICIT ANALYSIS OF H=45 TOURNAMENTS =====")
    for d in all_data:
        if d['H'] == 45:
            adj = adj_from_bits(d['bits'], n)
            all_cycles = count_all_directed_odd_cycles(adj, n)
            by_len = Counter(len(c[0]) for c in all_cycles)
            print(f"\n  bits={d['bits']}, H=45:")
            print(f"    Directed cycles by length: {dict(sorted(by_len.items()))}")
            print(f"    Total directed cycles: {len(all_cycles)}")

            # Find vertex-disjoint pairs
            vsets = [c[0] for c in all_cycles]
            disjoint_pairs = []
            for i in range(len(vsets)):
                for j in range(i+1, len(vsets)):
                    if not (vsets[i] & vsets[j]):
                        disjoint_pairs.append((i, j, vsets[i], vsets[j],
                                               len(vsets[i]), len(vsets[j])))
            print(f"    Vertex-disjoint cycle pairs: {len(disjoint_pairs)}")
            for p in disjoint_pairs[:10]:
                print(f"      {set(p[2])} ({p[4]}-cycle) + {set(p[3])} ({p[5]}-cycle)")
            break  # Just first example

    # Same for H=43 and H=41
    for H_target in [43, 41]:
        for d in all_data:
            if d['H'] == H_target:
                adj = adj_from_bits(d['bits'], n)
                all_cycles = count_all_directed_odd_cycles(adj, n)
                by_len = Counter(len(c[0]) for c in all_cycles)
                vsets = [c[0] for c in all_cycles]
                dp_count = sum(1 for i in range(len(vsets))
                              for j in range(i+1, len(vsets))
                              if not (vsets[i] & vsets[j]))
                print(f"\n  bits={d['bits']}, H={H_target}:")
                print(f"    Directed cycles by length: {dict(sorted(by_len.items()))}")
                print(f"    Total: {len(all_cycles)}, disjoint pairs: {dp_count}")
                break

if __name__ == "__main__":
    main()
