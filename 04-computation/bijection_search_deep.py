#!/usr/bin/env python3
"""
Deep search for the OCF bijection: Ham(T) <-> {(S, f) : S indep set in Omega(T), f: S -> {0,1}}

OCF says H(T) = sum_{k>=0} alpha_k * 2^k = sum_{S independent} 2^{|S|}

So we need to partition Ham(T) into groups indexed by independent sets S,
where the group for S has exactly 2^|S| paths.

At n=3 with one 3-cycle C=(0,1,2): H=3, alpha_0=1, alpha_1=1.
Partition: 1 path for S=emptyset, 2 paths for S={C}.

At n=4: various cases.

Strategy: For each tournament, try to find a CANONICAL assignment
path -> independent set that respects the 2^|S| count.

Key insight: The "transitive path" (path following the transitive ordering)
always exists. This should correspond to S = emptyset (the trivial independent set).
The other paths should be "explained" by the cycles they "use".

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count, find_odd_cycles, conflict_graph
from itertools import permutations
from collections import defaultdict
from math import comb

def get_ham_paths(T):
    """Return all Hamiltonian paths as tuples."""
    n = len(T)
    paths = []
    for perm in permutations(range(n)):
        if all(T[perm[i]][perm[i+1]] for i in range(n-1)):
            paths.append(perm)
    return paths

def get_independent_sets(adj):
    """Return all independent sets of the conflict graph."""
    m = len(adj)
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
    indep_sets = []
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
            indep_sets.append(mask)
    return indep_sets

def cycles_used_by_path(path, cycles):
    """Which cycles have ALL their arcs used by the path?"""
    n = len(path)
    # Build set of arcs used by path
    path_arcs = set()
    for i in range(n - 1):
        path_arcs.add((path[i], path[i+1]))

    used = []
    for ci, cycle in enumerate(cycles):
        L = len(cycle)
        all_in = True
        for j in range(L):
            arc = (cycle[j], cycle[(j+1) % L])
            if arc not in path_arcs:
                all_in = False
                break
        if all_in:
            used.append(ci)
    return used

def cycle_vertices_on_path(path, cycle):
    """Check if the cycle's vertices appear consecutively on the path
    and the cycle arcs agree with the path direction."""
    n = len(path)
    L = len(cycle)
    pos = {v: i for i, v in enumerate(path)}

    # Check all rotations of the cycle
    cycle_verts = set(cycle)
    positions = sorted([pos[v] for v in cycle_verts])

    # Are they consecutive?
    if positions[-1] - positions[0] != L - 1:
        return False, None

    # Get the path's order of these vertices
    start = positions[0]
    path_order = list(path[start:start+L])

    return True, path_order

def reversal_signature(path, cycles):
    """For each cycle, does the path traverse it 'forward' or 'backward'?"""
    sigs = []
    for cycle in cycles:
        verts = set(cycle)
        pos = {v: i for i, v in enumerate(path)}
        positions = sorted([pos[v] for v in verts])
        if positions[-1] - positions[0] != len(cycle) - 1:
            sigs.append(None)  # Not contiguous
            continue
        start = positions[0]
        path_sub = tuple(path[start:start+len(cycle)])
        # Check if path_sub matches cycle orientation
        if path_sub == tuple(cycle):
            sigs.append(0)
        elif path_sub == tuple(reversed(cycle)):
            sigs.append(1)
        else:
            sigs.append(None)
    return tuple(sigs)

print("=" * 70)
print("DEEP BIJECTION SEARCH")
print("=" * 70)

# Analyze small cases in detail
for n in [3, 4, 5]:
    print(f"\n{'='*60}")
    print(f"  n = {n}")
    print(f"{'='*60}")

    m = n * (n - 1) // 2
    interesting_cases = 0

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        paths = get_ham_paths(T)
        H = len(paths)

        cycles = find_odd_cycles(T)
        if not cycles:
            continue  # Transitive tournament, H=1, trivial

        cg = conflict_graph(cycles)
        indep_sets = get_independent_sets(cg)

        # Verify OCF
        ocf_sum = sum(2**bin(mask).count('1') for mask in indep_sets)
        assert ocf_sum == H, f"OCF fails: {ocf_sum} != {H}"

        # For each path, find which cycles it "uses" (all arcs present)
        path_cycle_usage = []
        for path in paths:
            used = cycles_used_by_path(path, cycles)
            path_cycle_usage.append((path, used))

        # Analyze: which independent sets are "used" by paths?
        used_masks = defaultdict(list)
        for path, used in path_cycle_usage:
            mask = 0
            for ci in used:
                mask |= 1 << ci
            # Check if this mask is an independent set
            if mask in indep_sets:
                used_masks[mask].append(path)

        # Group paths by their maximal independent set of used cycles
        if n <= 4 or (interesting_cases < 5 and len(cycles) >= 2):
            interesting_cases += 1
            print(f"\n  bits={bits}, H={H}, cycles={len(cycles)}")
            print(f"  Cycles: {cycles}")
            print(f"  Independent sets: {[bin(m) for m in indep_sets]}")

            # Show cycle-usage pattern
            for path, used in path_cycle_usage:
                mask = 0
                for ci in used:
                    mask |= 1 << ci
                is_indep = mask in indep_sets
                used_names = [f"C{ci}" for ci in used]
                print(f"    Path {path}: uses {used_names}, mask={bin(mask)}, indep={is_indep}")

            # Try a different approach: "signature" based on which arcs are reversals
            # A 3-cycle (a,b,c) has arcs a->b, b->c, c->a.
            # A path either uses all 3 arcs (traverses the cycle) or uses 0-2.
            # For each path, compute a 3-valued signature per cycle:
            #   +1 if traversed forward, -1 if traversed backward, 0 otherwise
            print(f"  Contiguity analysis:")
            for path, used in path_cycle_usage:
                contigs = []
                for ci, cycle in enumerate(cycles):
                    contig, order = cycle_vertices_on_path(path, cycle)
                    if contig:
                        contigs.append(f"C{ci}:{order}")
                    else:
                        contigs.append(f"C{ci}:scattered")
                print(f"    {path}: {', '.join(contigs)}")

# Detailed analysis at n=5: can we find a pattern?
print(f"\n{'='*60}")
print(f"  n=5 STATISTICAL ANALYSIS")
print(f"{'='*60}")

n = 5
m = n * (n - 1) // 2
stats = defaultdict(int)

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    paths = get_ham_paths(T)
    H = len(paths)
    cycles = find_odd_cycles(T)
    if not cycles:
        continue

    # For each path, count how many cycles it fully traverses
    for path in paths:
        used_count = len(cycles_used_by_path(path, cycles))
        stats[used_count] += 1

print(f"  Paths by #cycles fully traversed:")
for k in sorted(stats):
    print(f"    {k} cycles: {stats[k]} paths ({100*stats[k]/sum(stats.values()):.1f}%)")

# Key question: is there a natural way to assign each path to an independent set
# such that paths assigned to S_k form groups of exactly 2^k?
# The "fully traversed cycles" approach may not give independent sets.
# Try: assign path to its MAXIMAL independent set among fully-used cycles.

print(f"\n{'='*60}")
print(f"  MAXIMAL INDEPENDENT SET ASSIGNMENT")
print(f"{'='*60}")

n = 5
success = 0
fail = 0
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    paths = get_ham_paths(T)
    H = len(paths)
    cycles = find_odd_cycles(T)
    if len(cycles) < 2:
        continue

    cg = conflict_graph(cycles)
    indep_sets_list = get_independent_sets(cg)
    indep_set = set(indep_sets_list)

    # For each path, find maximal independent set among used cycles
    assignment = defaultdict(list)
    for path in paths:
        used = cycles_used_by_path(path, cycles)
        # Find maximal independent subset of used cycles
        best_mask = 0
        best_size = 0
        for mask in indep_sets_list:
            # Check if mask is subset of used
            if all((mask >> ci) & 1 == 0 or ci in used for ci in range(len(cycles))):
                size = bin(mask).count('1')
                if size > best_size:
                    best_size = size
                    best_mask = mask
        assignment[best_mask].append(path)

    # Check if this gives the right partition
    ok = all(len(paths_for_s) == 2**bin(s).count('1') for s, paths_for_s in assignment.items())
    # Also check all independent sets are covered
    covered = set(assignment.keys())
    ok = ok and covered == indep_set

    if ok:
        success += 1
    else:
        fail += 1

print(f"  n=5: maximal-independent-set assignment: {success} success, {fail} fail")

# Try another approach: assign path to empty set if it follows transitive order,
# otherwise assign based on the "first reversal"
print(f"\n{'='*60}")
print(f"  TRANSITIVE-DEVIATION APPROACH")
print(f"{'='*60}")

# For each tournament, find the transitive tournament closest to it
# (fewest arc reversals). The paths can be grouped by which arcs
# they "use" that differ from the transitive tournament.

n = 4
m = n * (n - 1) // 2
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    paths = get_ham_paths(T)
    H = len(paths)
    cycles = find_odd_cycles(T)
    if len(cycles) != 1:
        continue

    # Single 3-cycle case
    cycle = cycles[0]
    cg = conflict_graph(cycles)

    # H = 1 + 2*1 = 3 paths (if only one 3-cycle)
    # or more if longer cycles exist
    print(f"\n  bits={bits}, H={H}, cycle={cycle}")

    # Identify which path is the "base" (transitive-like)
    for path in paths:
        used = cycles_used_by_path(path, cycles)
        # Check if path traverses the cycle
        print(f"    Path {path}: uses {len(used)} cycles")

    # The path NOT using the 3-cycle = base path (S=emptyset)
    # The 2 paths using the 3-cycle = S={C} paths (forward and backward traversal?)
    break  # Just show first example

print(f"\n{'='*70}")
print("DONE")
print("=" * 70)
