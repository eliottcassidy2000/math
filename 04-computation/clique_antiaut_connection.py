#!/usr/bin/env python3
"""
Investigate connections between:
1. Cliques in Omega(T) (odd-cycle conflict graph)
2. Anti-automorphisms of T
3. The SC maximizer phenomenon

Key questions:
- How does the involutory anti-aut sigma create structure in Omega(T)?
- Does sigma induce an automorphism of Omega(T)?
- Do sigma-orbits of cycles explain the clique structure?
- How do clique covers of Omega relate to H(T)?

Instance: opus-2026-03-06-S4
"""
import sys
import time
from itertools import permutations, combinations
from collections import defaultdict, Counter

def build_adj(n, bits):
    """Build adjacency matrix from tiling bits (path arcs i->i+1 implicit)."""
    arcs = [(i, j) for i in range(n) for j in range(i + 2, n)]
    adj = [[0] * n for _ in range(n)]
    for i in range(n - 1):
        adj[i][i + 1] = 1
    for k, (i, j) in enumerate(arcs):
        if bits & (1 << k):
            adj[j][i] = 1
        else:
            adj[i][j] = 1
    return adj

def canon(adj, n):
    best = None
    for perm in permutations(range(n)):
        s = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or s < best:
            best = s
    return best

def count_ham(adj, n):
    full = (1 << n) - 1
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    nk = (mask | (1 << u), u)
                    dp[nk] = dp.get(nk, 0) + c
    return sum(dp.get((full, v), 0) for v in range(n))

def find_directed_odd_cycles(adj, n):
    """Find all directed odd cycles, each represented by sorted vertex set + direction."""
    cycles = []
    for length in range(3, n + 1, 2):
        for verts in combinations(range(n), length):
            verts_list = list(verts)
            start = verts_list[0]
            # Find all directed Hamiltonian cycles on this vertex set starting at start
            def dfs(path, used):
                last = path[-1]
                if len(path) == length:
                    if adj[last][start]:
                        cycles.append(tuple(path))
                    return
                for v in verts_list:
                    if v in used:
                        continue
                    if adj[last][v]:
                        used.add(v)
                        path.append(v)
                        dfs(path, used)
                        path.pop()
                        used.remove(v)
            dfs([start], {start})
    return cycles

def build_omega(cycles, n):
    """Build conflict graph: cycles adjacent iff sharing a vertex."""
    nc = len(cycles)
    cycle_sets = [frozenset(c) for c in cycles]
    adj = [[0] * nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i + 1, nc):
            if cycle_sets[i] & cycle_sets[j]:
                adj[i][j] = 1
                adj[j][i] = 1
    return adj, cycle_sets

def find_anti_auts(adj, n):
    auts = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            if not ok:
                break
            for j in range(n):
                if i == j:
                    continue
                if adj[perm[i]][perm[j]] != (1 - adj[i][j]):
                    ok = False
                    break
        if ok:
            auts.append(perm)
    return auts

def apply_perm_to_cycle(sigma, cycle):
    """Apply permutation sigma to a directed cycle."""
    return tuple(sigma[v] for v in cycle)

def normalize_cycle(cycle):
    """Normalize a directed cycle to start at min vertex."""
    n = len(cycle)
    min_idx = cycle.index(min(cycle))
    return tuple(cycle[(min_idx + i) % n] for i in range(n))

def reverse_cycle(cycle):
    """Reverse direction of a cycle."""
    return (cycle[0],) + tuple(reversed(cycle[1:]))

def compose(sigma, tau):
    return tuple(sigma[tau[i]] for i in range(len(sigma)))

def cycle_type(sigma):
    n = len(sigma)
    visited = [False] * n
    cycles = []
    for i in range(n):
        if visited[i]:
            continue
        cycle_len = 0
        j = i
        while not visited[j]:
            visited[j] = True
            j = sigma[j]
            cycle_len += 1
        cycles.append(cycle_len)
    return tuple(sorted(cycles, reverse=True))


def analyze_sigma_on_omega(adj, n, sigma, cycles, omega_adj, cycle_sets):
    """Analyze how anti-automorphism sigma acts on Omega(T)."""
    nc = len(cycles)

    # sigma reverses all arcs, so sigma maps cycle C to reversed cycle on sigma(V(C))
    # The reversed cycle is also a directed cycle (just opposite direction)
    # For odd cycles, reverse = different cycle (not same up to rotation)

    # Build mapping: for each cycle, find its sigma-image
    cycle_to_idx = {}
    for i, c in enumerate(cycles):
        nc_norm = normalize_cycle(c)
        cycle_to_idx[nc_norm] = i
        # Also store reversed
        rc = reverse_cycle(c)
        rc_norm = normalize_cycle(rc)
        # Note: reversed cycle might or might not be in our list

    sigma_image = [None] * nc
    sigma_pairs = []  # (i, j) where sigma maps cycle i to cycle j
    sigma_fixed = []  # cycles fixed by sigma (up to rotation/direction)

    for i, c in enumerate(cycles):
        # Apply sigma to each vertex
        img = apply_perm_to_cycle(sigma, c)
        # sigma reverses arcs, so the actual directed cycle is the REVERSE
        img_rev = reverse_cycle(img)
        img_norm = normalize_cycle(img_rev)

        j = cycle_to_idx.get(img_norm)
        if j is not None:
            sigma_image[i] = j
            if i == j:
                sigma_fixed.append(i)
            elif i < j:
                sigma_pairs.append((i, j))
        else:
            # Try without reversing (sigma might preserve direction for some cycles)
            img_norm2 = normalize_cycle(img)
            j2 = cycle_to_idx.get(img_norm2)
            if j2 is not None:
                sigma_image[i] = j2
                if i == j2:
                    sigma_fixed.append(i)
                elif i < j2:
                    sigma_pairs.append((i, j2))

    # Check: does sigma induce an automorphism of Omega?
    omega_aut = True
    for i in range(nc):
        if sigma_image[i] is None:
            omega_aut = False
            break
        for j in range(i + 1, nc):
            if sigma_image[j] is None:
                omega_aut = False
                break
            si, sj = sigma_image[i], sigma_image[j]
            if omega_adj[i][j] != omega_adj[si][sj]:
                omega_aut = False
                break
        if not omega_aut:
            break

    # Vertex-disjoint pairs from sigma
    disjoint_sigma_pairs = [(i, j) for i, j in sigma_pairs
                            if not (cycle_sets[i] & cycle_sets[j])]

    return {
        'omega_aut': omega_aut,
        'fixed': len(sigma_fixed),
        'pairs': len(sigma_pairs),
        'disjoint_pairs': len(disjoint_sigma_pairs),
        'unmapped': sum(1 for x in sigma_image if x is None),
    }


def main():
    target = int(sys.argv[1]) if len(sys.argv) > 1 else 6

    for n in range(4, target + 1):
        arcs = [(i, j) for i in range(n) for j in range(i + 2, n)]
        m = len(arcs)
        total = 1 << m

        print(f"\n{'=' * 70}")
        print(f"n={n}: Clique-AntiAut Connection Analysis")
        print(f"{'=' * 70}")

        t0 = time.time()
        seen = {}

        sc_data = []
        nsc_data = []

        for bits in range(total):
            if bits % 5000 == 0 and bits > 0:
                elapsed = time.time() - t0
                rate = bits / elapsed
                print(f"  {bits}/{total} ({100*bits/total:.1f}%), ETA {(total-bits)/rate:.0f}s",
                      flush=True)

            adj = build_adj(n, bits)
            cf = canon(adj, n)
            if cf in seen:
                continue
            seen[cf] = True

            h = count_ham(adj, n)
            scores = tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

            # Check SC
            conv = [[1 - adj[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
            conv_cf = canon(conv, n)
            is_sc = (cf == conv_cf)

            # Find odd cycles and build Omega
            cycles = find_directed_odd_cycles(adj, n)
            nc = len(cycles)
            omega_adj, cycle_sets = build_omega(cycles, n)

            # Independence polynomial coefficients
            alphas = []
            for k in range(nc + 1):
                count = 0
                for combo in combinations(range(nc), k):
                    ok = True
                    for a in range(len(combo)):
                        for b in range(a + 1, len(combo)):
                            if omega_adj[combo[a]][combo[b]]:
                                ok = False
                                break
                        if not ok:
                            break
                    if ok:
                        count += 1
                alphas.append(count)
                if count == 0:
                    break

            # Clique number of Omega
            max_clique = 0
            for k in range(1, nc + 1):
                found = False
                for combo in combinations(range(nc), k):
                    all_adj = True
                    for a in range(len(combo)):
                        for b in range(a + 1, len(combo)):
                            if not omega_adj[combo[a]][combo[b]]:
                                all_adj = False
                                break
                        if not all_adj:
                            break
                    if all_adj:
                        found = True
                        max_clique = k
                        break
                if not found:
                    break

            # Chromatic number = clique cover number (Omega is perfect for n<=7)
            # For now just record clique number

            entry = {
                'bits': bits, 'h': h, 'scores': scores, 'is_sc': is_sc,
                'nc': nc, 'alphas': alphas, 'clique_num': max_clique,
            }

            if is_sc:
                # Anti-automorphism analysis
                anti_auts = find_anti_auts(adj, n)
                involutions = [s for s in anti_auts if compose(s, s) == tuple(range(n))]
                entry['num_anti_auts'] = len(anti_auts)
                entry['num_involutions'] = len(involutions)

                # Analyze best involution's action on Omega
                if involutions and cycles:
                    best_result = None
                    for sigma in involutions:
                        result = analyze_sigma_on_omega(adj, n, sigma, cycles, omega_adj, cycle_sets)
                        if best_result is None or result['disjoint_pairs'] > best_result['disjoint_pairs']:
                            best_result = result
                            best_sigma = sigma
                    entry['sigma_analysis'] = best_result
                    entry['sigma_cycle_type'] = cycle_type(best_sigma)

                sc_data.append(entry)
            else:
                nsc_data.append(entry)

        elapsed = time.time() - t0
        print(f"\n  {len(sc_data)} SC classes, {len(nsc_data)} NSC classes ({elapsed:.1f}s)")

        # Report: SC tournaments
        print(f"\n  SC TOURNAMENT OMEGA STRUCTURE:")
        print(f"  {'H':>5} {'#cyc':>5} {'cliq#':>5} {'alpha':>15} {'#aa':>4} {'#inv':>4} "
              f"{'omega_aut':>9} {'sigma_dpairs':>12} {'sigma_ct':>15} {'scores'}")
        for d in sorted(sc_data, key=lambda x: x['h']):
            sa = d.get('sigma_analysis', {})
            aa_str = f"{d.get('num_anti_auts', '?')}"
            inv_str = f"{d.get('num_involutions', '?')}"
            oa_str = str(sa.get('omega_aut', '?'))
            dp_str = str(sa.get('disjoint_pairs', '?'))
            ct_str = str(d.get('sigma_cycle_type', '?'))
            print(f"  {d['h']:5d} {d['nc']:5d} {d['clique_num']:5d} {str(d['alphas']):>15} "
                  f"{aa_str:>4} {inv_str:>4} {oa_str:>9} {dp_str:>12} {ct_str:>15} {d['scores']}")

        # Report: NSC tournaments (sample)
        if nsc_data:
            print(f"\n  NSC TOURNAMENT OMEGA STRUCTURE (sample):")
            print(f"  {'H':>5} {'#cyc':>5} {'cliq#':>5} {'alpha':>15} {'scores'}")
            for d in sorted(nsc_data, key=lambda x: x['h'])[:20]:
                print(f"  {d['h']:5d} {d['nc']:5d} {d['clique_num']:5d} {str(d['alphas']):>15} {d['scores']}")

        # Compare SC vs NSC within same score sequence
        print(f"\n  SC vs NSC WITHIN SCORE CLASSES:")
        score_groups = defaultdict(lambda: {'sc': [], 'nsc': []})
        for d in sc_data:
            score_groups[d['scores']]['sc'].append(d)
        for d in nsc_data:
            score_groups[d['scores']]['nsc'].append(d)

        for scores in sorted(score_groups.keys()):
            g = score_groups[scores]
            if g['sc'] and g['nsc']:
                sc_h = [d['h'] for d in g['sc']]
                nsc_h = [d['h'] for d in g['nsc']]
                sc_alpha2 = [d['alphas'][2] if len(d['alphas']) > 2 else 0 for d in g['sc']]
                nsc_alpha2 = [d['alphas'][2] if len(d['alphas']) > 2 else 0 for d in g['nsc']]
                sc_dp = [d.get('sigma_analysis', {}).get('disjoint_pairs', '?') for d in g['sc']]

                print(f"    scores={scores}:")
                print(f"      SC:  H={sc_h}, alpha2={sc_alpha2}, sigma_dpairs={sc_dp}")
                print(f"      NSC: H={nsc_h}, alpha2={nsc_alpha2}")
                print(f"      SC max H >= NSC max H: {max(sc_h) >= max(nsc_h)}")

        # Key statistics
        print(f"\n  KEY FINDINGS:")
        if sc_data:
            all_omega_aut = all(d.get('sigma_analysis', {}).get('omega_aut', False) for d in sc_data if d.get('sigma_analysis'))
            print(f"    sigma induces Omega automorphism: {'ALL' if all_omega_aut else 'NOT ALL'}")

            max_h_sc = max(d['h'] for d in sc_data)
            max_h_all = max(max_h_sc, max((d['h'] for d in nsc_data), default=0))
            print(f"    Global max H = {max_h_all}, achieved by SC: {max_h_sc == max_h_all}")

        print(f"\n  Time: {time.time() - t0:.1f}s")


if __name__ == '__main__':
    main()
