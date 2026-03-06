#!/usr/bin/env python3
"""
Exhaustive anti-automorphism analysis for SC tournaments.

Questions investigated:
1. Is every anti-automorphism an involution? (sigma^2 = id)
2. What is the orbit structure of involutory anti-automorphisms?
3. How does the involution structure relate to H(T)?
4. Do ALL SC tournaments have an involutory anti-aut? (always yes for finite groups)

Runs exhaustively at n=4,5,6,7. n=7 takes ~30-60min.

Instance: opus-2026-03-06-S4
"""
import sys
import time
from itertools import permutations
from collections import defaultdict, Counter

def build_arcs(n):
    return [(i, j) for i in range(n) for j in range(i + 2, n)]

def build_adj_flat(n, bits, arcs):
    adj = bytearray(n * n)
    for i in range(n - 1):
        adj[i * n + (i + 1)] = 1
    for k, (i, j) in enumerate(arcs):
        if bits & (1 << k):
            adj[j * n + i] = 1
        else:
            adj[i * n + j] = 1
    return adj

def adj_to_matrix(adj, n):
    return [[adj[i * n + j] for j in range(n)] for i in range(n)]

def count_ham_dp(adj, n):
    full = (1 << n) - 1
    dp = {}
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get(mask * n + v, 0)
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v * n + u]:
                    nk = (mask | (1 << u)) * n + u
                    dp[nk] = dp.get(nk, 0) + c
    return sum(dp.get(full * n + v, 0) for v in range(n))

def canon_full(adj, n):
    best = None
    for perm in permutations(range(n)):
        s = tuple(adj[perm[i] * n + perm[j]] for i in range(n) for j in range(n))
        if best is None or s < best:
            best = s
    return best

def find_anti_auts(adj, n):
    """Find all anti-automorphisms: sigma with adj[sigma(i)][sigma(j)] = 1 - adj[i][j] for i!=j."""
    auts = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            if not ok:
                break
            for j in range(n):
                if i == j:
                    continue
                if adj[perm[i] * n + perm[j]] != (1 - adj[i * n + j]):
                    ok = False
                    break
        if ok:
            auts.append(perm)
    return auts

def find_auts(adj, n):
    """Find all automorphisms."""
    auts = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            if not ok:
                break
            for j in range(n):
                if i == j:
                    continue
                if adj[perm[i] * n + perm[j]] != adj[i * n + j]:
                    ok = False
                    break
        if ok:
            auts.append(perm)
    return auts

def compose(sigma, tau):
    return tuple(sigma[tau[i]] for i in range(len(sigma)))

def perm_order(sigma):
    n = len(sigma)
    x = list(range(n))
    for k in range(1, 2 * n + 1):
        x = [sigma[xi] for xi in x]
        if all(x[i] == i for i in range(n)):
            return k
    return None

def cycle_type(sigma):
    n = len(sigma)
    visited = [False] * n
    cycles = []
    for i in range(n):
        if visited[i]:
            continue
        cycle = []
        j = i
        while not visited[j]:
            visited[j] = True
            cycle.append(j)
            j = sigma[j]
        cycles.append(len(cycle))
    return tuple(sorted(cycles, reverse=True))

def find_odd_cycles(adj, n):
    """Find all directed odd cycles in the tournament."""
    cycles = []
    for length in range(3, n + 1, 2):  # odd lengths only
        # Enumerate all vertex subsets of this size
        for mask in range(1, 1 << n):
            if bin(mask).count('1') != length:
                continue
            verts = [v for v in range(n) if mask & (1 << v)]
            # Try all orderings starting with smallest vertex (avoid counting each cycle length times)
            # Actually, find all directed cycles on this vertex set
            # Use DFS from verts[0]
            start = verts[0]
            def dfs(path, visited_set):
                last = path[-1]
                if len(path) == length:
                    # Check if there's an edge back to start
                    if adj[last * n + start]:
                        cycles.append(tuple(path))
                    return
                for v in verts:
                    if v in visited_set:
                        continue
                    if adj[last * n + v]:
                        visited_set.add(v)
                        path.append(v)
                        dfs(path, visited_set)
                        path.pop()
                        visited_set.remove(v)
            dfs([start], {start})
    # Each cycle of length L is found L times (once for each starting point)
    # But we fixed start = verts[0], so we find each cycle exactly once
    # Actually no - we find each directed cycle that starts at verts[0]
    # A directed cycle on L vertices has L rotations, all starting at different vertices
    # By fixing start = min vertex, we get exactly 1 representative per cycle
    return cycles

def count_disjoint_pairs(cycles, n):
    """Count pairs of vertex-disjoint cycles."""
    count = 0
    for i in range(len(cycles)):
        si = set(cycles[i])
        for j in range(i + 1, len(cycles)):
            sj = set(cycles[j])
            if si.isdisjoint(sj):
                count += 1
    return count

def independence_poly_at_2(cycles, n):
    """Compute I(Omega(T), 2) where Omega is the odd-cycle conflict graph."""
    # Build conflict graph
    nc = len(cycles)
    if nc == 0:
        return 1
    cycle_sets = [set(c) for c in cycles]
    # adj_omega[i][j] = 1 if cycles i,j share a vertex
    adj_omega = [[0] * nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i + 1, nc):
            if not cycle_sets[i].isdisjoint(cycle_sets[j]):
                adj_omega[i][j] = 1
                adj_omega[j][i] = 1

    # Enumerate independent sets
    total = 0
    for mask in range(1 << nc):
        bits = []
        m = mask
        while m:
            bits.append(m & -m)
            m &= m - 1
        indices = []
        m = mask
        idx = 0
        while m:
            if m & 1:
                indices.append(idx)
            m >>= 1
            idx += 1

        # Check independence
        ok = True
        for a in range(len(indices)):
            for b in range(a + 1, len(indices)):
                if adj_omega[indices[a]][indices[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            total += (1 << len(indices))  # 2^|S|
    return total


def main():
    target = int(sys.argv[1]) if len(sys.argv) > 1 else 7

    for n in range(4, target + 1):
        arcs = build_arcs(n)
        m = len(arcs)
        total = 1 << m

        print(f"\n{'=' * 70}")
        print(f"n={n}: {m} arcs, {total} tilings, n!={len(list(permutations(range(n))))} perms")
        print(f"{'=' * 70}")

        t0 = time.time()
        seen = {}
        sc_classes = []
        nsc_classes = []

        for bits in range(total):
            if bits % 5000 == 0 and bits > 0:
                elapsed = time.time() - t0
                rate = bits / elapsed
                eta = (total - bits) / rate
                print(f"  Scanning {bits}/{total} ({100*bits/total:.1f}%), "
                      f"{rate:.0f}/s, ETA {eta:.0f}s", flush=True)

            adj = build_adj_flat(n, bits, arcs)
            cf = canon_full(adj, n)
            if cf in seen:
                continue
            seen[cf] = True

            # Check SC: is converse isomorphic to T?
            # Converse: reverse all arcs
            conv_adj = bytearray(n * n)
            for i in range(n):
                for j in range(n):
                    if i != j:
                        conv_adj[i * n + j] = 1 - adj[i * n + j]
            conv_cf = canon_full(conv_adj, n)

            h = count_ham_dp(adj, n)

            if cf == conv_cf:
                sc_classes.append((bits, adj, h, cf))
            else:
                nsc_classes.append((bits, adj, h, cf))

        elapsed = time.time() - t0
        print(f"\n  Found {len(sc_classes)} SC classes, {len(nsc_classes)} NSC classes in {elapsed:.1f}s")

        # Analyze anti-automorphisms of SC tournaments
        print(f"\n  ANTI-AUTOMORPHISM ANALYSIS:")

        all_involutory = 0
        has_non_involutory = 0
        cycle_type_dist = Counter()
        involution_count_dist = Counter()

        h_by_aut_size = defaultdict(list)
        h_by_num_involutions = defaultdict(list)

        for bits, adj, h, cf in sc_classes:
            anti_auts = find_anti_auts(adj, n)
            auts = find_auts(adj, n)

            # Classify anti-auts
            involutions = [s for s in anti_auts if compose(s, s) == tuple(range(n))]
            non_inv = [s for s in anti_auts if compose(s, s) != tuple(range(n))]

            h_by_aut_size[len(auts)].append(h)
            h_by_num_involutions[len(involutions)].append(h)

            if non_inv:
                has_non_involutory += 1
                print(f"\n    NON-INVOLUTORY in SC class (H={h}, |Aut|={len(auts)}, |AntiAut|={len(anti_auts)}):")
                for s in non_inv[:3]:
                    s2 = compose(s, s)
                    ct = cycle_type(s)
                    print(f"      sigma={s}, sigma^2={s2}, order={perm_order(s)}, cycle_type={ct}")
            else:
                all_involutory += 1

            # Record cycle types of involutory anti-auts
            for s in involutions:
                ct = cycle_type(s)
                cycle_type_dist[ct] += 1

            involution_count_dist[len(involutions)] += 1

        print(f"\n    All-involutory SC classes: {all_involutory}/{len(sc_classes)}")
        print(f"    Has non-involutory anti-aut: {has_non_involutory}/{len(sc_classes)}")
        print(f"\n    Involutory anti-aut cycle types:")
        for ct in sorted(cycle_type_dist.keys()):
            print(f"      {ct}: {cycle_type_dist[ct]} occurrences")

        print(f"\n    Number of involutory anti-auts per SC class:")
        for k in sorted(involution_count_dist.keys()):
            print(f"      {k} involutions: {involution_count_dist[k]} classes")

        # H vs automorphism group size
        print(f"\n    H by |Aut(T)|:")
        for aut_size in sorted(h_by_aut_size.keys()):
            hs = sorted(h_by_aut_size[aut_size])
            print(f"      |Aut|={aut_size}: H={hs}, max={max(hs)}")

        # H vs number of involutory anti-auts
        print(f"\n    H by #involutory anti-auts:")
        for k in sorted(h_by_num_involutions.keys()):
            hs = sorted(h_by_num_involutions[k])
            print(f"      #inv_antiaut={k}: H={hs}, max={max(hs)}")

        # Clique analysis of Omega for SC vs NSC
        if n <= 6:  # Omega analysis only feasible for small n
            print(f"\n  OMEGA CLIQUE ANALYSIS (SC vs NSC):")

            for label, classes in [("SC", sc_classes), ("NSC", nsc_classes[:20])]:
                if not classes:
                    continue
                print(f"\n    {label} tournaments:")
                for bits, adj, h, cf in classes[:10]:
                    cycles = find_odd_cycles(adj, n)
                    nc = len(cycles)
                    dp = count_disjoint_pairs(cycles, n)

                    # Compute independence polynomial coefficients
                    cycle_sets = [set(c) for c in cycles]
                    # alpha_k = number of independent sets of size k
                    alphas = [0] * (nc + 1)
                    for mask in range(1 << nc):
                        indices = [i for i in range(nc) if mask & (1 << i)]
                        ok = True
                        for a in range(len(indices)):
                            for b in range(a + 1, len(indices)):
                                if not cycle_sets[indices[a]].isdisjoint(cycle_sets[indices[b]]):
                                    ok = False
                                    break
                            if not ok:
                                break
                        if ok:
                            alphas[len(indices)] += 1

                    # Clique number
                    max_clique = 0
                    for mask in range(1 << nc):
                        indices = [i for i in range(nc) if mask & (1 << i)]
                        all_adj = True
                        for a in range(len(indices)):
                            for b in range(a + 1, len(indices)):
                                if cycle_sets[indices[a]].isdisjoint(cycle_sets[indices[b]]):
                                    all_adj = False
                                    break
                            if not all_adj:
                                break
                        if all_adj and len(indices) > max_clique:
                            max_clique = len(indices)

                    alpha_str = [alphas[k] for k in range(min(4, nc + 1))]
                    scores = tuple(sorted(sum(adj[i * n + j] for j in range(n)) for i in range(n)))
                    print(f"      H={h}, #cycles={nc}, disjoint_pairs={dp}, "
                          f"clique#={max_clique}, alpha={alpha_str}, scores={scores}")

        print(f"\n  Total time for n={n}: {time.time() - t0:.1f}s")


if __name__ == '__main__':
    main()
