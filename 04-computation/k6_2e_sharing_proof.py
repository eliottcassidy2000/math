#!/usr/bin/env python3
"""
K_6-2e cannot be realized as Omega(T): computational proof.

THEOREM: No tournament T has Omega(T) containing a connected component
isomorphic to K_6 minus 2 edges (K_6-2e), where Omega(T) is the conflict
graph on ALL directed odd cycles (3-cycles, 5-cycles, 7-cycles, ...).

CONTEXT (THM-079):
- H(T) = I(Omega(T), 2) by OCF (Grinberg-Stanley)
- I(K_6-2e, 2) = 1 + 12 + 8 = 21
- All other routes to H=21 are ruled out
- If K_6-2e is impossible as an Omega component, then H=21 is a permanent gap

APPROACH:
1. Verify Sharing Lemma: two 3-cycles sharing a vertex on 5 vertices
   always force additional cycles
2. For small n (5,6,7): exhaustively enumerate all tournaments and check
   if any Omega(T) has a K_6-2e connected component
3. For n=8: sampling check
4. Combinatorial analysis of triple/cycle arrangements

IMPORTANT: Omega(T) uses ALL directed odd cycles, not just 3-cycles.
Each directed cycle (up to rotation) is a distinct vertex of Omega.
Two cycles are adjacent iff they share at least one vertex.

Instance: opus-2026-03-07
"""

import itertools
import sys
import time
from collections import Counter, defaultdict

# =====================================================================
# Core utilities
# =====================================================================

def find_all_directed_odd_cycles(adj, n):
    """Find all directed odd cycles in a tournament.
    Returns list of tuples, each a canonical representation of a directed cycle.
    Canonical form: starts with min vertex, second element < last element."""
    cycles = []
    for size in range(3, n+1, 2):
        for verts in itertools.combinations(range(n), size):
            v0 = verts[0]
            for perm in itertools.permutations(verts[1:]):
                full = (v0,) + perm
                is_cycle = True
                for idx in range(size):
                    if not adj[full[idx]][full[(idx+1)%size]]:
                        is_cycle = False
                        break
                if is_cycle:
                    # Normalize
                    min_idx = full.index(min(full))
                    rotated = full[min_idx:] + full[:min_idx]
                    if rotated[1] > rotated[-1]:
                        rotated = (rotated[0],) + tuple(reversed(rotated[1:]))
                    cycles.append(rotated)
    return list(set(cycles))

def cycle_vertex_set(c):
    """Get the vertex set of a directed cycle."""
    return frozenset(c)

def build_omega(cycles):
    """Build conflict graph adjacency matrix.
    Cycles adjacent iff they share a vertex."""
    m = len(cycles)
    adj = [[False]*m for _ in range(m)]
    vsets = [cycle_vertex_set(c) for c in cycles]
    for a in range(m):
        for b in range(a+1, m):
            if vsets[a] & vsets[b]:
                adj[a][b] = True
                adj[b][a] = True
    return adj, vsets

def connected_components(adj, n):
    """Find connected components."""
    visited = [False]*n
    components = []
    for start in range(n):
        if visited[start]:
            continue
        comp = []
        stack = [start]
        while stack:
            v = stack.pop()
            if visited[v]:
                continue
            visited[v] = True
            comp.append(v)
            for u in range(n):
                if adj[v][u] and not visited[u]:
                    stack.append(u)
        components.append(comp)
    return components

def ham_count(adj, n):
    """Count Hamiltonian paths in tournament."""
    count = 0
    for perm in itertools.permutations(range(n)):
        ok = True
        for idx in range(n-1):
            if not adj[perm[idx]][perm[idx+1]]:
                ok = False
                break
        if ok:
            count += 1
    return count

def indep_poly_at_2(adj, m):
    """Compute I(G, 2) for graph G given as adjacency matrix."""
    total = 0
    for mask in range(2**m):
        verts = [i for i in range(m) if (mask >> i) & 1]
        indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            total += 2**len(verts)
    return total

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(2**len(edges)):
        adj = [[False]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
        yield adj

# =====================================================================
# Part 1: Verify OCF at n=5 (sanity check)
# =====================================================================

def verify_ocf(n):
    """Verify H(T) = I(Omega(T), 2) for all tournaments on n vertices."""
    print(f"Verifying OCF at n={n}...")
    mismatches = 0
    total = 0
    for adj in all_tournaments(n):
        total += 1
        cycles = find_all_directed_odd_cycles(adj, n)
        h = ham_count(adj, n)

        m = len(cycles)
        if m == 0:
            ip = 1
        else:
            omega, vsets = build_omega(cycles)
            ip = indep_poly_at_2(omega, m)

        if h != ip:
            mismatches += 1
            if mismatches <= 2:
                print(f"  MISMATCH: H={h}, I(Omega,2)={ip}, #cycles={m}")

    if mismatches == 0:
        print(f"  VERIFIED: All {total} tournaments match.")
    else:
        print(f"  FAILED: {mismatches}/{total} mismatches")
    return mismatches == 0

# =====================================================================
# Part 2: Verify Sharing Lemma
# =====================================================================

def verify_sharing_lemma():
    """
    Verify: Two 3-cycles sharing a vertex with |union| = 5
    always force at least one additional odd cycle on a subset of those 5 vertices.
    """
    print("\n" + "=" * 70)
    print("SHARING LEMMA VERIFICATION")
    print("=" * 70)

    c1 = frozenset([0, 1, 2])
    c2 = frozenset([2, 3, 4])

    cases = 0
    min_extra = float('inf')

    for adj in all_tournaments(5):
        # Check if {0,1,2} and {2,3,4} are both 3-cycles
        has_c1 = ((adj[0][1] and adj[1][2] and adj[2][0]) or
                  (adj[0][2] and adj[2][1] and adj[1][0]))
        has_c2 = ((adj[2][3] and adj[3][4] and adj[4][2]) or
                  (adj[2][4] and adj[4][3] and adj[3][2]))

        if not (has_c1 and has_c2):
            continue

        cases += 1
        all_cycles = find_all_directed_odd_cycles(adj, 5)
        # Count odd cycles whose vertex set is a subset of {0,1,2,3,4}
        # (which is all of them, since n=5)
        cycle_vsets = [cycle_vertex_set(c) for c in all_cycles]

        # Exclude the two original 3-cycles (as vertex sets)
        extra = [c for c in cycle_vsets if c != c1 and c != c2]
        # But we need to count distinct vertex sets, not directed cycles
        extra_vsets = set(extra)
        min_extra = min(min_extra, len(extra_vsets))

    print(f"  Cases with C1={{0,1,2}} and C2={{2,3,4}}: {cases}")
    print(f"  Minimum extra odd cycle vertex sets: {min_extra}")
    print(f"  Sharing Lemma VERIFIED: {min_extra >= 1}")

    # Also check with 5-cycle cycles: the extra cycle might be a 5-cycle
    print("\n  Breakdown: what types of extra cycles appear?")
    min_extra_3 = float('inf')
    min_extra_5 = float('inf')

    for adj in all_tournaments(5):
        has_c1 = ((adj[0][1] and adj[1][2] and adj[2][0]) or
                  (adj[0][2] and adj[2][1] and adj[1][0]))
        has_c2 = ((adj[2][3] and adj[3][4] and adj[4][2]) or
                  (adj[2][4] and adj[4][3] and adj[3][2]))

        if not (has_c1 and has_c2):
            continue

        all_cycles = find_all_directed_odd_cycles(adj, 5)
        extra_3 = set()
        extra_5 = set()
        for c in all_cycles:
            vs = frozenset(c)
            if vs == c1 or vs == c2:
                continue
            if len(c) == 3:
                extra_3.add(vs)
            elif len(c) == 5:
                extra_5.add(vs)

        min_extra_3 = min(min_extra_3, len(extra_3))
        min_extra_5 = min(min_extra_5, len(extra_5))

    print(f"  Min extra 3-cycle vertex sets: {min_extra_3}")
    print(f"  Min extra 5-cycle vertex sets: {min_extra_5}")

    return min_extra >= 1

# =====================================================================
# Part 3: Exhaustive tournament check for K_6-2e component in Omega
# =====================================================================

def check_k6_2e_component(adj, n, cycles=None):
    """Check if Omega(T) has a connected component isomorphic to K_6-2e.
    Returns True if found, with details."""
    if cycles is None:
        cycles = find_all_directed_odd_cycles(adj, n)

    m = len(cycles)
    if m < 6:
        return False, None

    omega, vsets = build_omega(cycles)
    comps = connected_components(omega, m)

    for comp in comps:
        if len(comp) != 6:
            continue
        # Check if induced subgraph is K_6-2e (13 edges, 2 non-edges)
        edge_count = 0
        for i in range(6):
            for j in range(i+1, 6):
                if omega[comp[i]][comp[j]]:
                    edge_count += 1
        if edge_count == 13:
            comp_cycles = [cycles[c] for c in comp]
            comp_vsets = [vsets[c] for c in comp]
            return True, {
                'cycles': comp_cycles,
                'vsets': comp_vsets,
                'total_cycles': m
            }

    return False, None

def exhaustive_tournament_check(n, verbose=True):
    """Check all tournaments on n vertices for K_6-2e in Omega."""
    print(f"\n{'='*70}")
    print(f"EXHAUSTIVE CHECK: n={n}")
    print(f"{'='*70}")

    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    num_edges = len(edges)
    total = 2**num_edges
    print(f"  Total tournaments: {total}")

    found = 0
    count = 0
    t0 = time.time()

    for bits in range(total):
        count += 1
        if count % 500000 == 0 and verbose:
            elapsed = time.time() - t0
            rate = count / elapsed if elapsed > 0 else 0
            print(f"  Progress: {count}/{total} ({100*count/total:.1f}%, "
                  f"{rate:.0f}/s, ETA {(total-count)/rate:.0f}s)", flush=True)

        adj = [[False]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True

        is_k6, info = check_k6_2e_component(adj, n)
        if is_k6:
            found += 1
            if found <= 5:
                print(f"\n  *** FOUND K_6-2e component! ***")
                print(f"  Cycles: {info['cycles']}")
                print(f"  Cycle sizes: {[len(c) for c in info['cycles']]}")
                print(f"  Total Omega vertices: {info['total_cycles']}")
                h = ham_count(adj, n)
                print(f"  H(T) = {h}")

    elapsed = time.time() - t0
    print(f"\n  Done in {elapsed:.1f}s")
    print(f"  K_6-2e components found: {found}")
    if found == 0:
        print(f"  ==> CONFIRMED: No K_6-2e component in Omega(T) for any n={n} tournament")
    return found

# =====================================================================
# Part 4: n=8 sampling
# =====================================================================

def sampled_check(n, num_samples=200000):
    """Sample random tournaments for K_6-2e in Omega."""
    import random

    print(f"\n{'='*70}")
    print(f"SAMPLING CHECK: n={n} ({num_samples} samples)")
    print(f"{'='*70}")

    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    num_edges = len(edges)

    found = 0
    t0 = time.time()

    for s in range(num_samples):
        if s % 100000 == 0 and s > 0:
            print(f"  Progress: {s}/{num_samples}")

        bits = random.randint(0, 2**num_edges - 1)
        adj = [[False]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True

        is_k6, info = check_k6_2e_component(adj, n)
        if is_k6:
            found += 1
            print(f"  *** FOUND K_6-2e component! ***")
            print(f"  Cycles: {info['cycles']}")
            print(f"  Cycle sizes: {[len(c) for c in info['cycles']]}")

    elapsed = time.time() - t0
    print(f"\n  Done in {elapsed:.1f}s. Found: {found}")
    return found

# =====================================================================
# Part 5: Pure 3-cycle analysis with Sharing Lemma
# =====================================================================

def analyze_3cycle_k6_2e():
    """
    Analyze whether 6 three-cycles can form K_6-2e in Omega.
    For each arrangement of 6 triples with K_6-2e intersection pattern,
    check if the Sharing Lemma forces extra cycles.
    """
    print("\n" + "=" * 70)
    print("PURE 3-CYCLE ANALYSIS: K_6-2e via Sharing Lemma")
    print("=" * 70)

    for n in range(5, 9):
        triples = list(itertools.combinations(range(n), 3))
        num_triples = len(triples)
        num_combos = 1
        for i in range(6):
            num_combos = num_combos * (num_triples - i) // (i + 1)

        print(f"\n--- n = {n}: C({num_triples}, 6) = {num_combos} ---")

        if num_combos > 20_000_000:
            print(f"  Skipping (too large for direct enumeration)")
            continue

        total = 0
        killed_by_sharing = 0
        survivors = []

        for six_idx in itertools.combinations(range(num_triples), 6):
            six_sets = [frozenset(triples[i]) for i in six_idx]

            # Check K_6-2e intersection pattern
            disjoint_count = 0
            sharing_pairs = []
            for a in range(6):
                for b in range(a+1, 6):
                    if six_sets[a] & six_sets[b]:
                        sharing_pairs.append((a, b))
                    else:
                        disjoint_count += 1
            if disjoint_count != 2:
                continue

            total += 1

            # For each sharing pair with |union| = 5: check if covered
            has_uncovered_5 = False
            for a, b in sharing_pairs:
                union_ab = six_sets[a] | six_sets[b]
                if len(union_ab) != 5:
                    continue
                # Check: is there another triple fully inside union_ab?
                covered = any(six_sets[k] <= union_ab
                             for k in range(6) if k != a and k != b)
                if not covered:
                    has_uncovered_5 = True
                    break

            if has_uncovered_5:
                killed_by_sharing += 1
            else:
                if len(survivors) < 5:
                    survivors.append([set(s) for s in six_sets])

        print(f"  K_6-2e triple arrangements: {total}")
        print(f"  Killed by Sharing Lemma (uncovered 5-vertex pair): {killed_by_sharing}")
        print(f"  Survivors (all 5-vertex pairs covered): {total - killed_by_sharing}")

        if survivors:
            print(f"\n  Survivor examples:")
            for s in survivors[:3]:
                print(f"    {s}")
                span = set().union(*s)
                print(f"    Vertex span: {len(span)}")

            # For survivors at n=6, check realizability
            if n <= 7:
                print(f"\n  Checking survivor realizability...")
                check_3cycle_survivors(survivors[:5], n)

def check_3cycle_survivors(survivor_list, n):
    """Check if survivor 3-cycle arrangements are realizable AND isolated in Omega."""
    for idx, six_list in enumerate(survivor_list):
        six_sets = [frozenset(s) for s in six_list]
        all_verts = sorted(set().union(*six_sets))
        nv = len(all_verts)
        vert_map = {v: i for i, v in enumerate(all_verts)}
        mapped = [frozenset(vert_map[v] for v in s) for s in six_sets]

        e_list = [(i, j) for i in range(nv) for j in range(i+1, nv)]
        ne = len(e_list)

        realizable_exact = 0  # exactly these 6 as ALL 3-cycles
        realizable_isolated = 0  # these 6 as isolated component in full Omega

        for bits in range(2**ne):
            tadj = [[False]*nv for _ in range(nv)]
            for k, (i, j) in enumerate(e_list):
                if (bits >> k) & 1:
                    tadj[i][j] = True
                else:
                    tadj[j][i] = True

            # Check all 6 triples are 3-cycles
            all_ok = True
            for trip in mapped:
                verts = sorted(trip)
                a, b, c = verts
                is_cyc = ((tadj[a][b] and tadj[b][c] and tadj[c][a]) or
                          (tadj[a][c] and tadj[c][b] and tadj[b][a]))
                if not is_cyc:
                    all_ok = False
                    break
            if not all_ok:
                continue

            # Find ALL odd cycles (not just 3-cycles)
            all_cycles = find_all_directed_odd_cycles(tadj, nv)
            all_vsets = [cycle_vertex_set(c) for c in all_cycles]

            # Check: are exactly these 6 triple vertex-sets the only
            # odd cycle vertex sets that intersect any of them?
            target_vsets = set(mapped)
            target_union = set().union(*mapped)

            # Cycles NOT in our 6 that share vertices with our 6
            intruders = []
            for c, vs in zip(all_cycles, all_vsets):
                if vs not in target_vsets:
                    if vs & target_union:
                        intruders.append((c, vs))

            if not intruders:
                realizable_isolated += 1
                if realizable_isolated <= 2:
                    print(f"\n    Survivor {idx}: REALIZABLE AND ISOLATED!")
                    print(f"    Tournament on {nv} vertices:")
                    for i in range(nv):
                        out = [j for j in range(nv) if tadj[i][j]]
                        print(f"      {i} -> {out}")
                    print(f"    All odd cycles: {len(all_cycles)}")
                    print(f"    Target cycles: {[set(s) for s in mapped]}")
                    h = ham_count(tadj, nv)
                    print(f"    H(T) = {h}")

        if realizable_isolated == 0:
            print(f"    Survivor {idx}: NOT realizable as isolated component")
            print(f"      (every realization has additional cycles sharing vertices)")
        else:
            print(f"    Survivor {idx}: {realizable_isolated} realizations as isolated component")

# =====================================================================
# Part 6: Direct H=21 search
# =====================================================================

def search_h21(n):
    """Search for tournaments with H(T) = 21."""
    print(f"\n{'='*70}")
    print(f"DIRECT SEARCH: H(T) = 21 at n={n}")
    print(f"{'='*70}")

    found = 0
    total = 0
    for adj in all_tournaments(n):
        total += 1
        h = ham_count(adj, n)
        if h == 21:
            found += 1
            if found <= 3:
                cycles = find_all_directed_odd_cycles(adj, n)
                omega, vsets = build_omega(cycles)
                m = len(cycles)
                ip = indep_poly_at_2(omega, m)
                print(f"  H=21 found! #cycles={m}, I(Omega,2)={ip}")
                print(f"  Cycle sizes: {sorted([len(c) for c in cycles])}")

    print(f"  Total tournaments: {total}")
    print(f"  Tournaments with H=21: {found}")
    return found

# =====================================================================
# Main
# =====================================================================

if __name__ == "__main__":
    t0 = time.time()

    # Sanity check: verify OCF at n=4
    print("=" * 70)
    print("SANITY CHECK: OCF verification")
    print("=" * 70)
    verify_ocf(4)

    # Verify Sharing Lemma
    sharing_ok = verify_sharing_lemma()

    # Pure 3-cycle analysis
    analyze_3cycle_k6_2e()

    # Exhaustive checks
    exhaustive_tournament_check(5)
    exhaustive_tournament_check(6)

    # Direct H=21 search at n=5,6
    search_h21(5)
    search_h21(6)

    # n=7: exhaustive (2^21 = ~2M tournaments, but cycle-finding is slow)
    # This is the expensive part
    print("\n" + "=" * 70)
    print("n=7 EXHAUSTIVE CHECK (may take several minutes)")
    print("=" * 70)
    exhaustive_tournament_check(7)

    # n=8: sampling
    sampled_check(8, 100000)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    elapsed = time.time() - t0
    print(f"Total runtime: {elapsed:.1f}s")
