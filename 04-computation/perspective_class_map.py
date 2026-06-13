#!/usr/bin/env python3
"""perspective_class_map.py -- Map perspectives at n to classes at n+1.

Session: kind-pasteur-2026-03-20-S4

P(n) = T(n+1) for n=1,2,3,4 but P(5) = 48 < T(6) = 56.
Deficit = 8. Study which perspectives share classes and why.

APPROACH:
1. Enumerate all 48 rooted tournaments at n=5
2. For each, extend to n=6 by adding vertex with all 2^5 orientations
3. Determine which of the 56 n=6 classes each extension belongs to
4. Build the bipartite map: perspectives -> classes
5. Find which classes have multiple parent perspectives
6. Find which perspectives reach the most classes
"""

from itertools import permutations
from collections import defaultdict, Counter
from math import factorial

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

def count_3cycles(adj, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    c3 += 1
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    c3 += 1
    return c3

def are_isomorphic(adj1, adj2, n):
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(i+1, n):
                if adj1[perm[i]][perm[j]] != adj2[i][j]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return True
    return False

def find_aut_group(adj, n):
    auts = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(i+1, n):
                if adj[perm[i]][perm[j]] != adj[i][j]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            auts.append(perm)
    return auts

def vertex_orbits(adj, n):
    auts = find_aut_group(adj, n)
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry
    for aut in auts:
        for v in range(n):
            union(v, aut[v])
    orbits = defaultdict(set)
    for v in range(n):
        orbits[find(v)].add(v)
    return list(orbits.values())


def main():
    # ============================================================
    # Step 1: Find all rooted tournament types at n=5
    # ============================================================
    n = 5
    m = n * (n-1) // 2
    N = 1 << m

    print("=" * 70)
    print(f"PERSPECTIVE -> CLASS MAP: n={n} -> n+1={n+1}")
    print("=" * 70)

    # Find iso classes at n=5
    classes_5 = []
    seen_5 = {}

    for bits in range(N):
        adj = adj_from_bits(bits, n)
        sc = score_seq(adj, n)
        H = held_karp_H(adj, n)
        key = (sc, H)

        is_new = True
        for cls in classes_5:
            if cls['key'] == key:
                rep_adj = adj_from_bits(cls['bits'], n)
                if are_isomorphic(adj, rep_adj, n):
                    is_new = False
                    break

        if is_new:
            orbits = vertex_orbits(adj, n)
            classes_5.append({
                'bits': bits, 'key': key, 'score': sc, 'H': H,
                'orbits': orbits,
                'aut_size': len(find_aut_group(adj, n)),
            })

    # Build rooted tournament list
    rooted_5 = []
    for cls_idx, cls in enumerate(classes_5):
        adj = adj_from_bits(cls['bits'], n)
        for orb in cls['orbits']:
            v_rep = min(orb)
            sv = sum(adj[v_rep][j] for j in range(n) if j != v_rep)
            rooted_5.append({
                'class_idx': cls_idx,
                'bits': cls['bits'],
                'vertex': v_rep,
                'score_v': sv,
                'orbit_size': len(orb),
                'class_score': cls['score'],
                'class_H': cls['H'],
            })

    print(f"\n  n={n}: {len(classes_5)} classes, {len(rooted_5)} perspectives (rooted types)")

    # ============================================================
    # Step 2: Find all iso classes at n+1=6
    # ============================================================
    n1 = n + 1
    m1 = n1 * (n1-1) // 2
    N1 = 1 << m1

    print(f"  Building n={n1} class table...")
    classes_6 = []
    for bits in range(N1):
        adj = adj_from_bits(bits, n1)
        sc = score_seq(adj, n1)
        H = held_karp_H(adj, n1)
        c3 = count_3cycles(adj, n1)
        key = (sc, H, c3)

        is_new = True
        for cls in classes_6:
            if cls['key'] == key:
                rep_adj = adj_from_bits(cls['bits'], n1)
                if are_isomorphic(adj, rep_adj, n1):
                    is_new = False
                    break

        if is_new:
            classes_6.append({
                'bits': bits, 'key': key, 'score': sc, 'H': H, 'c3': c3,
            })

    print(f"  n={n1}: {len(classes_6)} classes")

    # ============================================================
    # Step 3: For each rooted tournament, find reachable classes at n+1
    # ============================================================
    print(f"\n  Computing extension map...")

    # Function to find which n=6 class a tournament belongs to
    def find_class_6(adj_ext):
        sc = score_seq(adj_ext, n1)
        H = held_karp_H(adj_ext, n1)
        c3 = count_3cycles(adj_ext, n1)
        key = (sc, H, c3)
        for idx, cls in enumerate(classes_6):
            if cls['key'] == key:
                rep_adj = adj_from_bits(cls['bits'], n1)
                if are_isomorphic(adj_ext, rep_adj, n1):
                    return idx
        return -1  # should not happen

    ext_map = defaultdict(set)  # rooted_idx -> set of class_6 indices
    class_parents = defaultdict(set)  # class_6_idx -> set of rooted_idx

    for rt_idx, rt in enumerate(rooted_5):
        adj = adj_from_bits(rt['bits'], n)
        v = rt['vertex']

        for new_arcs in range(1 << n):
            # Build extension: add vertex w=n to the tournament
            adj_ext = [row[:] + [0] for row in adj]
            adj_ext.append([0] * n1)
            w = n
            for i in range(n):
                if (new_arcs >> i) & 1:
                    adj_ext[w][i] = 1
                else:
                    adj_ext[i][w] = 1

            cls_idx = find_class_6(adj_ext)
            ext_map[rt_idx].add(cls_idx)
            class_parents[cls_idx].add(rt_idx)

    # ============================================================
    # Step 4: Analyze the map
    # ============================================================
    print(f"\n  {'='*60}")
    print(f"  EXTENSION MAP ANALYSIS")
    print(f"  {'='*60}")

    # How many classes does each perspective reach?
    print(f"\n  Perspectives and their reachable classes:")
    for rt_idx, rt in enumerate(rooted_5):
        classes_reached = ext_map[rt_idx]
        class_descs = []
        for c_idx in sorted(classes_reached):
            cls = classes_6[c_idx]
            class_descs.append(f"H={cls['H']}")
        print(f"    RT[{rt_idx:2d}] {str(rt['class_score']):>20} v_score={rt['score_v']} "
              f"| reaches {len(classes_reached)} classes: {', '.join(class_descs)}")

    # How many parent perspectives does each class have?
    print(f"\n  Classes and their parent perspectives:")
    for c_idx, cls in enumerate(classes_6):
        parents = class_parents[c_idx]
        parent_descs = []
        for rt_idx in sorted(parents):
            rt = rooted_5[rt_idx]
            parent_descs.append(f"({rt['class_score']},v={rt['score_v']})")
        num_parents = len(parents)
        marker = " <<<" if num_parents > 1 else ""
        print(f"    CLS[{c_idx:2d}] {str(cls['score']):>24} H={cls['H']:>3} "
              f"| {num_parents} parents{marker}")

    # Summary stats
    parent_counts = [len(class_parents[i]) for i in range(len(classes_6))]
    reach_counts = [len(ext_map[i]) for i in range(len(rooted_5))]

    print(f"\n  Summary:")
    print(f"    Perspectives: {len(rooted_5)}")
    print(f"    Classes: {len(classes_6)}")
    print(f"    Deficit: {len(classes_6) - len(rooted_5)}")
    print(f"    Classes with 1 parent: {sum(1 for c in parent_counts if c == 1)}")
    print(f"    Classes with 2+ parents: {sum(1 for c in parent_counts if c >= 2)}")
    print(f"    Classes with 0 parents: {sum(1 for c in parent_counts if c == 0)}")
    print(f"    Parent count distribution: {Counter(parent_counts)}")
    print(f"    Reach count distribution: {Counter(reach_counts)}")


if __name__ == "__main__":
    main()
