#!/usr/bin/env python3
"""
burnside_sl_s214.py -- opus-2026-03-23-S214

BURNSIDE FORMULA FOR SL (self-loop arc-orbits)

SL = (1/n!) * sum_{sigma in S_n} |{(T, e) : sigma(T)=T, sigma(e)=e, T_e ≅ T}|

Let's compute this sum term-by-term for each cycle type.
For cycle type lambda with cycle structure (c_1, c_2, ..., c_n):
  Fix_lambda = #{(T, e) : sigma(T)=T, sigma(e)=e, T_e ≅ T}
  for some sigma of type lambda.

The conditions:
1. sigma(T) = T: T is a tournament fixed by sigma
2. sigma(e) = e: arc e = (a,b) is fixed by sigma (both a,b are fixed points of sigma)
3. T_e ≅ T: reversing e gives a tournament isomorphic to T

Condition 2 means a and b are BOTH in 1-cycles of sigma.
So for sigma with k fixed points: the arc e must be among the C(k,2) arcs
between fixed points.

Condition 3 is the hard one: it requires T_e to be in the same iso class as T.

For sigma = identity (k=n): all n(n-1)/2 arcs are candidates.
Fix(id) = #{(T, e) : T_e ≅ T} = total self-reversible (T, e) pairs.

This script computes Fix(sigma) for each cycle type at n=3..7.

Author: opus-2026-03-23-S214
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canonical(A, n):
    if n <= 5:
        best = None
        for perm in permutations(range(n)):
            form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
            if best is None or form < best: best = form
        return best
    else:
        scores = [sum(A[i]) for i in range(n)]
        score_groups = defaultdict(list)
        for v in range(n): score_groups[scores[v]].append(v)
        groups = [score_groups[s] for s in sorted(score_groups.keys())]
        best = None
        def gen_perms(grps):
            if not grps:
                yield []
                return
            for p in permutations(grps[0]):
                for rest in gen_perms(grps[1:]):
                    yield list(p) + rest
        for perm in gen_perms(groups):
            form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
            if best is None or form < best: best = form
        return best

def cycle_type(perm):
    """Return sorted cycle lengths of permutation."""
    n = len(perm)
    visited = [False]*n
    cycles = []
    for i in range(n):
        if visited[i]: continue
        length = 0; j = i
        while not visited[j]:
            visited[j] = True; j = perm[j]; length += 1
        cycles.append(length)
    return tuple(sorted(cycles))

def count_by_cycle_type(n):
    """For each cycle type, compute Fix = #{(T, e) : sigma(T)=T, sigma(e)=e, T_e ~ T}."""
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Build all iso classes
    cmap = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canonical(A, n)
        if c not in cmap:
            cmap[c] = len(cmap)

    # For each permutation sigma:
    # Find Fix(sigma) = #{(T, e) : sigma(T)=T, sigma(e)=e, T_e ~ T}
    fix_by_type = defaultdict(int)
    type_count = defaultdict(int)

    for perm in permutations(range(n)):
        ct = cycle_type(perm)
        type_count[ct] += 1

        # Fixed points of sigma
        fixed_pts = [i for i in range(n) if perm[i] == i]
        k = len(fixed_pts)
        if k < 2:
            # No arc can be fixed by sigma (need 2 fixed points)
            continue

        # Arcs between fixed points
        candidate_arcs = [(a,b) for a in fixed_pts for b in fixed_pts if a < b]

        # Enumerate tournaments fixed by sigma
        # A tournament T is fixed by sigma iff A[sigma(i)][sigma(j)] = A[i][j] for all i,j
        # This constrains arcs: for each orbit of pairs under sigma, all pairs in
        # the orbit must have the same orientation.

        # Orbits of ordered pairs (i,j) with i!=j under sigma:
        # sigma maps (i,j) to (sigma(i), sigma(j))
        pair_orbits = []
        seen_pairs = set()
        for i in range(n):
            for j in range(n):
                if i == j: continue
                if (i,j) in seen_pairs: continue
                orbit = []
                a, b = i, j
                while (a,b) not in seen_pairs:
                    seen_pairs.add((a,b))
                    orbit.append((a,b))
                    a, b = perm[a], perm[b]
                pair_orbits.append(orbit)

        # For T to be fixed by sigma: for each orbit of ordered pairs,
        # all pairs must have the same A value (0 or 1).
        # But: for undirected edge {i,j}, the orbit of (i,j) and (j,i) must
        # have A[i][j] + A[j][i] = 1 (tournament condition).

        # Better: work with unordered pairs. The orbit of {i,j} under sigma
        # is the set of {sigma^k(i), sigma^k(j)} for k.
        # All edges in an orbit must have consistent orientation.

        # Group unordered pairs into orbits under sigma
        edge_orbits = []
        seen_edges = set()
        for i in range(n):
            for j in range(i+1, n):
                if (i,j) in seen_edges: continue
                orbit = []
                a, b = i, j
                while True:
                    edge = (min(a,b), max(a,b))
                    if edge in seen_edges: break
                    seen_edges.add(edge)
                    orbit.append(edge)
                    a, b = perm[a], perm[b]
                edge_orbits.append(orbit)

        # For each edge orbit: the orientation is determined by the first edge.
        # But we need to check consistency: sigma maps (a->b) to (sigma(a)->sigma(b)).
        # If the orbit has length L, then sigma^L maps (a->b) back to itself.
        # The orientation must be consistent: if (a->b) then (sigma(a)->sigma(b)) etc.

        # Actually, for a tournament fixed by sigma:
        # For each edge orbit, there are exactly 2 consistent orientations
        # (unless the orbit contains a "reversal", i.e., sigma maps (i,j) to (j,i)).
        # If the orbit maps (i,j) to (j,i) at some point, only 1 orientation works.

        # Let's enumerate: for each orbit, check if sigma eventually reverses.
        n_free_orbits = 0
        n_forced_orbits = 0
        free_orbits = []

        for orbit in edge_orbits:
            # Check if any edge in the orbit gets reversed by sigma
            # Take first edge (a,b) with a<b. Follow sigma:
            a, b = orbit[0]
            reversed_at_some_point = False
            aa, bb = a, b
            for _ in range(len(orbit)):
                aa, bb = perm[aa], perm[bb]
                if aa > bb:
                    reversed_at_some_point = True
                    break
            if reversed_at_some_point:
                n_forced_orbits += 1
            else:
                n_free_orbits += 1
                free_orbits.append(orbit)

        # Total fixed tournaments = 2^(n_free_orbits) * (1 for each forced orbit)
        # Wait: forced orbits have exactly 1 valid orientation (since the chain
        # eventually reverses, the orientation is constrained to alternate).
        # Free orbits have 2 choices.

        # Actually more carefully: for a free orbit of length L:
        # sigma maps (a,b) consistently (never reverses), so both orientations work.
        # For a forced orbit: the chain a,sigma(a),...,sigma^{L-1}(a) and similarly
        # for b. At some point aa > bb, so (sigma^k(a), sigma^k(b)) has aa > bb,
        # meaning A[aa][bb] = 1-A[bb][aa]. But sigma maps A[a][b] to A[sigma(a)][sigma(b)],
        # which must equal A[a][b]. If at step k we get (bb, aa) with bb < aa, then
        # A[bb][aa] = A[a][b] (from sigma^k), but also A[aa][bb] = 1 - A[bb][aa].
        # For the orbit of length L where reversal happens at step k:
        # A[a][b] is determined. Exactly 1 choice.

        # Hmm, this is getting complicated. For the PURPOSE of this script,
        # let me just ENUMERATE all fixed tournaments directly for small n.

        # Enumerate tournaments fixed by sigma
        # For each edge orbit, choose orientation (2 if free, check both if forced)
        # Total: iterate over 2^n_free choices, for each check forced orbits

        fixed_T = 0
        fixed_T_e = 0  # (T, e) with e fixed by sigma AND T_e ~ T

        # For small number of free orbits, enumerate
        n_total_orbits = len(edge_orbits)
        if n_total_orbits > 20:
            continue  # skip huge cases

        # Represent each orbit as (list_of_edges, is_free)
        orbit_info = []
        for orbit in edge_orbits:
            a, b = orbit[0]
            reversed_flag = False
            aa, bb = a, b
            for _ in range(len(orbit)):
                aa, bb = perm[aa], perm[bb]
                if aa > bb:
                    reversed_flag = True; break
            orbit_info.append((orbit, not reversed_flag))

        # Enumerate over free orbit choices
        free_indices = [i for i, (_, is_free) in enumerate(orbit_info) if is_free]
        forced_indices = [i for i, (_, is_free) in enumerate(orbit_info) if not is_free]

        for bits in range(2**len(free_indices)):
            # Build tournament
            A = [[0]*n for _ in range(n)]

            # Set free orbits
            valid = True
            for bi, fi in enumerate(free_indices):
                orient = (bits >> bi) & 1
                for edge in orbit_info[fi][0]:
                    a, b = edge
                    if orient: A[a][b] = 1
                    else: A[b][a] = 1

            # Set forced orbits (try orient 0, check consistency)
            for fi in forced_indices:
                orbit = orbit_info[fi][0]
                # Try orient = 0 (A[a][b] = 0, A[b][a] = 1 for first edge)
                # Check if sigma-consistent
                a0, b0 = orbit[0]
                # First, set A[b0][a0] = 1
                A[b0][a0] = 1; A[a0][b0] = 0
                # Propagate
                aa, bb = a0, b0
                for _ in range(len(orbit)-1):
                    aa, bb = perm[aa], perm[bb]
                    ea, eb = min(aa,bb), max(aa,bb)
                    if aa < bb:
                        A[bb][aa] = 1; A[aa][bb] = 0
                    else:
                        A[aa][bb] = 1; A[bb][aa] = 0

            # Verify tournament is valid (no contradictions)
            valid = True
            for i in range(n):
                for j in range(i+1, n):
                    if A[i][j] + A[j][i] != 1:
                        valid = False; break
                if not valid: break

            if not valid: continue

            # Verify sigma-fixed
            ok = all(A[i][j] == A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
            if not ok: continue

            fixed_T += 1

            # Check candidate arcs (between fixed points of sigma)
            c_T = canonical(A, n)
            for a, b in candidate_arcs:
                # Arc (a,b) if A[a][b]=1, else (b,a)
                if A[a][b]:
                    arc = (a, b)
                else:
                    arc = (b, a)

                # Reverse this arc
                A2 = [row[:] for row in A]
                A2[arc[0]][arc[1]] = 0
                A2[arc[1]][arc[0]] = 1
                c_T2 = canonical(A2, n)

                if c_T2 == c_T:
                    # Self-reversible! This (T, e) contributes to Fix(sigma)
                    fixed_T_e += 1

        fix_by_type[ct] += fixed_T_e

    return fix_by_type, type_count

# ============================================================
banner("BURNSIDE DECOMPOSITION OF SL")
# ============================================================

for n in [3, 4, 5]:
    t0 = time.time()
    fix_by_type, type_count = count_by_cycle_type(n)
    dt = time.time() - t0

    print("  n = %d (%.1fs):" % (n, dt))

    total_fix = 0
    for ct in sorted(fix_by_type.keys()):
        f = fix_by_type[ct]
        cnt = type_count[ct]
        # Number of sigma with this cycle type
        # (the count is already the sum over all sigma of that type)
        print("    cycle type %s: Fix = %d (from %d permutations)" % (ct, f, cnt))
        total_fix += f

    SL_burnside = total_fix // factorial(n)
    print("    Total Fix = %d, SL = %d/%d = %d" %
          (total_fix, total_fix, factorial(n), SL_burnside))
    SL_expected = {3:2, 4:6, 5:16}[n]
    print("    Expected SL = %d, match = %s" % (SL_expected, SL_burnside == SL_expected))
    print()

# ============================================================
banner("IDENTITY PERMUTATION CONTRIBUTION")
# ============================================================

# For sigma = identity: Fix(id) = #{(T, e) : T_e ~ T}
# = total self-reversible (T, e) pairs across ALL labeled tournaments

for n in [3, 4, 5, 6]:
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Build classes
    cmap = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canonical(A, n)
        if c not in cmap:
            cmap[c] = len(cmap)

    # Count self-reversible pairs
    sr_count = 0
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        cT = canonical(A, n)

        for k in range(m):
            bits2 = bits ^ (1<<k)
            A2 = [[0]*n for _ in range(n)]
            for kk,(ii,jj) in enumerate(pairs):
                if (bits2>>kk)&1: A2[ii][jj]=1
                else: A2[jj][ii]=1
            if canonical(A2, n) == cT:
                sr_count += 1

    print("  n=%d: Fix(id) = %d, Fix(id)/n! = %d" % (n, sr_count, sr_count // factorial(n)))
    print("       (total self-reversible (T, e) pairs)")
    print("       This is the main term in the Burnside sum for SL.")
    print()

print("\nDone. opus-2026-03-23-S214")
