#!/usr/bin/env python3
"""sc_maximizer_orbit_theorem.py -- Prove the SC Maximizer for n=6 regular.

Session: kind-pasteur-2026-03-20-S1

THEOREM ATTEMPT: For the regular score class (2,2,2,3,3,3) at n=6,
every SC tournament T with alpha_2(T) = 4 (all complementary 3-set
pairs are both cyclic) achieves H(T) = 45 = max H within this score class.

PROOF STRATEGY:
1. Show that SC tournaments with involutory anti-aut sigma have a
   natural pairing of 3-cycles into vertex-disjoint pairs.
2. For regular score at n=6, c3=8 always. The 8 cyclic 3-subsets
   can form at most 4 complementary pairs.
3. SC achieves this maximum (4 comp pairs) iff all 4 orbit-complementary
   selections are both cyclic.
4. This gives alpha_2 >= 4, hence H >= 1 + 2*alpha_1 + 4*4 = 1 + 2*alpha_1 + 16.
5. Need to show alpha_1 >= 14 when all comp pairs are cyclic.

ALTERNATIVE: Direct computation showing H = I(Omega,2) determines
the structure completely.
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
    """Find all involutory anti-automorphisms (sigma^2 = id)."""
    result = []
    for perm in permutations(range(n)):
        # Check if sigma^2 = id
        if not all(perm[perm[i]] == i for i in range(n)):
            continue
        # Check anti-automorphism: adj[sigma(i)][sigma(j)] = 1 - adj[i][j]
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

def analyze_orbit_structure(sigma, n):
    """Analyze the orbit structure of an involution."""
    orbits = []
    used = set()
    for i in range(n):
        if i not in used:
            if sigma[i] == i:
                orbits.append((i,))
            else:
                orbits.append((i, sigma[i]))
            used.add(i)
            used.add(sigma[i])
    return orbits

def main():
    print("=" * 70)
    print("SC MAXIMIZER ORBIT THEOREM — n=6 REGULAR ANALYSIS")
    print("=" * 70)

    n = 6
    m = n * (n - 1) // 2
    target = (2, 2, 2, 3, 3, 3)

    # Exhaustively analyze ALL regular n=6 tournaments
    all_data = []
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

        # Find involutory anti-automorphisms
        inv_auts = find_involutory_anti_auts(adj, n)

        # Count complementary pairs (both halves cyclic)
        comp_both = 0
        for triple in combinations(range(n), 3):
            s = frozenset(triple)
            comp = frozenset(range(n)) - s
            if s < comp and s in cyclic and comp in cyclic:
                comp_both += 1

        all_data.append({
            'bits': bits, 'H': H, 'c3': len(cyclic),
            'comp_both': comp_both, 'num_inv_aut': len(inv_auts),
            'has_anti_aut': len(inv_auts) > 0,
        })

    # Group by (H, comp_both, has_anti_aut)
    print(f"\n  Total regular n=6 tournaments: {len(all_data)}")

    for H_val in sorted(set(d['H'] for d in all_data)):
        entries = [d for d in all_data if d['H'] == H_val]
        has_aa = sum(1 for d in entries if d['has_anti_aut'])
        no_aa = sum(1 for d in entries if not d['has_anti_aut'])
        comp_vals = Counter(d['comp_both'] for d in entries)
        print(f"\n  H={H_val}: {len(entries)} tournaments")
        print(f"    Has invol anti-aut: {has_aa}, No anti-aut: {no_aa}")
        print(f"    comp_both values: {dict(sorted(comp_vals.items()))}")

    # THE KEY THEOREM ATTEMPT:
    # Show that comp_both = 4 <=> H = 45 (for regular n=6)
    print(f"\n  --- THEOREM CHECK ---")
    print(f"  comp_both=4 implies H=45?")
    comp4 = [d for d in all_data if d['comp_both'] == 4]
    print(f"    {len(comp4)} tournaments with comp_both=4")
    print(f"    H values: {Counter(d['H'] for d in comp4)}")

    print(f"  H=45 implies comp_both=4?")
    h45 = [d for d in all_data if d['H'] == 45]
    print(f"    {len(h45)} tournaments with H=45")
    print(f"    comp_both values: {Counter(d['comp_both'] for d in h45)}")

    # Check: is comp_both=4 equivalent to having invol anti-aut with specific orbit structure?
    print(f"\n  --- ORBIT STRUCTURE OF ANTI-AUTOMORPHISMS ---")
    for d in all_data[:5]:
        if d['has_anti_aut'] and d['H'] == 45:
            adj = adj_from_bits(d['bits'], n)
            inv_auts = find_involutory_anti_auts(adj, n)
            print(f"\n  bits={d['bits']}, H={d['H']}:")
            for sigma in inv_auts[:2]:
                orbits = analyze_orbit_structure(sigma, n)
                print(f"    sigma={sigma}, orbits={orbits}")

                # Check: which 3-cycle vertex sets come from orbit selections?
                pairs = [o for o in orbits if len(o) == 2]
                if len(pairs) == 3:
                    # Type (a): one from each pair = 2^3 = 8 selections
                    orbit_cyclic = 0
                    for sel in range(8):
                        verts = [pairs[i][1 if (sel >> i) & 1 else 0] for i in range(3)]
                        s = frozenset(verts)
                        i, j, k = verts
                        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                           (adj[i][k] and adj[k][j] and adj[j][i]):
                            orbit_cyclic += 1
                    print(f"    Orbit selections cyclic: {orbit_cyclic}/8")

    # DEEP ANALYSIS: For comp_both=4, are ALL 8 three-cycles orbit-generated?
    print(f"\n  --- DEEP: ALL 8 CYCLES FROM ORBITS? ---")
    for d in all_data:
        if d['comp_both'] == 4 and d['has_anti_aut']:
            adj = adj_from_bits(d['bits'], n)
            inv_auts = find_involutory_anti_auts(adj, n)
            sigma = inv_auts[0]
            orbits = analyze_orbit_structure(sigma, n)
            pairs = [o for o in orbits if len(o) == 2]

            if len(pairs) != 3:
                continue

            # Count orbit-generated 3-cycles
            orbit_cycles = set()
            for sel in range(8):
                verts = tuple(sorted(pairs[i][1 if (sel >> i) & 1 else 0] for i in range(3)))
                i, j, k = verts
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    orbit_cycles.add(frozenset(verts))

            # Count non-orbit 3-cycles (type (b): 2 from one pair, 1 from another)
            non_orbit_cycles = set()
            cyclic = set()
            for triple in combinations(range(n), 3):
                i, j, k = triple
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    cyclic.add(frozenset(triple))
            non_orbit = cyclic - orbit_cycles

            print(f"  bits={d['bits']}: orbit_cycles={len(orbit_cycles)}, "
                  f"non_orbit={len(non_orbit)}, total_c3={len(cyclic)}")
            break  # Just show first example

    # FINAL THEOREM SUMMARY
    print(f"\n  ===== THEOREM SUMMARY =====")
    print(f"  For regular n=6 tournaments (score 2,2,2,3,3,3):")
    print(f"  - ALL tournaments have c3=8 (verified: {set(d['c3'] for d in all_data)})")
    print(f"  - H takes values: {sorted(set(d['H'] for d in all_data))}")
    print(f"  - comp_both (# complementary pairs both cyclic):")

    for cb in sorted(set(d['comp_both'] for d in all_data)):
        entries = [d for d in all_data if d['comp_both'] == cb]
        H_vals = sorted(set(d['H'] for d in entries))
        aa = sum(1 for d in entries if d['has_anti_aut'])
        print(f"    comp_both={cb}: count={len(entries)}, H_values={H_vals}, "
              f"has_anti_aut={aa}/{len(entries)}")


if __name__ == "__main__":
    main()
