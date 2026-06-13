#!/usr/bin/env python3
r"""beta2_one_free_component.py - Prove at most 1 free component

GOAL: Show that in any tournament, the 3-cycle graph has at most 1
connected component consisting entirely of free (non-dominated) cycles.

KEY THEOREM (THM-106, proved): For a free cycle C = {a,b,c},
every external vertex v creates a bridge cycle sharing an edge with C.

CONSEQUENCE: Every vertex of the tournament is in some cycle in C's
component. So the vertex set of C's component = V(T).

But this doesn't immediately show all CYCLES are in C's component.
We need: for any other cycle C2, C2 shares an edge with some cycle
in C's component.

APPROACH: Show that if C1 is free, any cycle C2 is reachable from C1
via edge-sharing chains. Use the bridge theorem iteratively.

KEY ARGUMENT:
For free C1 = a->b->c->a and cycle C2 = {d,e,f} in different component:
1. Each v in {d,e,f} creates a bridge to C1 (by THM-106)
2. Bridge B_v = {v, u, w} with u,w in C1, shares edge of C1
3. C2 has edges d->e, e->f, f->d. For edge d->e:
   Consider triples {d,e,u} for u in {a,b,c}.
   If any is a 3-cycle, it shares edge (d,e) with C2 AND shares
   vertex pair (d,u) or (e,u) with B_d or B_e (connecting to C1's comp).
4. If ALL triples {d,e,u} are transitive for all u in C1:
   Then d and e have a specific dominance pattern over C1.
   Combined with the bridge patterns of d and e over C1,
   this is highly constrained.

This script explores whether we can prove step 4 leads to contradiction.

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

random.seed(42)


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def is_free_cycle(A, n, a, b, c):
    """Check if 3-cycle {a,b,c} is free (non-dominated)."""
    for d in range(n):
        if d in (a, b, c):
            continue
        if A[d][a] and A[d][b] and A[d][c]:
            return False
        if A[a][d] and A[b][d] and A[c][d]:
            return False
    return True


def get_bridge_pattern(A, v, a, b, c):
    """Get pattern (v->a, v->b, v->c) for vertex v w.r.t. cycle (a,b,c)."""
    return (A[v][a], A[v][b], A[v][c])


def bridge_edge(pattern):
    """Which edge of a->b->c->a does this pattern bridge?
    Returns 'ab', 'bc', or 'ca' (or None if pattern is 000/111)."""
    x, y, z = pattern
    # Bridge from a->b: x=1, y=0
    # Bridge from b->c: y=1, z=0
    # Bridge from c->a: z=1, x=0
    bridges = []
    if x == 1 and y == 0:
        bridges.append('ab')
    if y == 1 and z == 0:
        bridges.append('bc')
    if z == 1 and x == 0:
        bridges.append('ca')
    return bridges


# ============================================================
# Part 1: Detailed bridge analysis for disconnected components
# ============================================================
print("=" * 70)
print("BRIDGE PATTERN ANALYSIS")
print("=" * 70)

# At n=6, find cases with 2+ components and a free component
n = 6
total = 2 ** (n*(n-1)//2)
found = 0

for bits in range(total):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    # Find 3-cycles
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles.append((i, j, k))
                elif A[i][k] and A[k][j] and A[j][i]:
                    cycles.append((i, k, j))

    c3 = len(cycles)
    if c3 < 2:
        continue

    # Cycle graph components
    cycle_edge_sets = []
    for (i, j, k) in cycles:
        es = set()
        for x, y in [(i,j),(j,k),(k,i),(j,i),(k,j),(i,k)]:
            if A[x][y]:
                es.add((x,y))
        cycle_edge_sets.append(es)

    adj = [[] for _ in range(c3)]
    for p in range(c3):
        for q in range(p+1, c3):
            if cycle_edge_sets[p] & cycle_edge_sets[q]:
                adj[p].append(q)
                adj[q].append(p)

    visited = set()
    components = []
    for start in range(c3):
        if start in visited:
            continue
        comp = []
        stack = [start]
        visited.add(start)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v not in visited:
                    visited.add(v)
                    stack.append(v)
        components.append(comp)

    if len(components) < 2:
        continue

    # Check for free components
    free_comps = []
    for ci, comp in enumerate(components):
        all_free = True
        for idx2 in comp:
            a, b, c = cycles[idx2]
            if not is_free_cycle(A, n, a, b, c):
                all_free = False
                break
        if all_free:
            free_comps.append(ci)

    if len(free_comps) == 0:
        continue

    found += 1
    if found <= 3:
        print(f"\nbits={bits}: {len(components)} components, {len(free_comps)} free")
        for ci, comp in enumerate(components):
            is_free = ci in free_comps
            print(f"  Component {ci} ({'FREE' if is_free else 'dominated'}): "
                  f"{[cycles[i] for i in comp]}")

        # Show bridge patterns for the free component
        fc = free_comps[0]
        fc_cycles = [cycles[i] for i in components[fc]]
        # Pick first cycle
        a, b, c = fc_cycles[0]
        print(f"  Free cycle ({a},{b},{c}) = {a}->{b}->{c}->{a}")
        for v in range(n):
            if v in (a, b, c):
                continue
            pat = get_bridge_pattern(A, v, a, b, c)
            br = bridge_edge(pat)
            print(f"    v={v}: pattern={pat}, bridges={br}")

        # Show dominated cycles
        for ci, comp in enumerate(components):
            if ci in free_comps:
                continue
            for idx2 in comp:
                cyc = cycles[idx2]
                x, y, z = cyc
                dominators = []
                for d in range(n):
                    if d in (x, y, z):
                        continue
                    if A[d][x] and A[d][y] and A[d][z]:
                        dominators.append(f"{d} dom")
                    if A[x][d] and A[y][d] and A[z][d]:
                        dominators.append(f"{d} sub")
                print(f"  Dominated cycle {cyc}: {dominators}")

print(f"\nn=6: {found} tournaments with multi-component + free component")


# ============================================================
# Part 2: Why at most 1 free component? Check larger n
# ============================================================
print(f"\n{'=' * 70}")
print("MAX FREE COMPONENTS (SAMPLED)")
print("=" * 70)

for n in [7, 8, 9, 10]:
    random.seed(42)
    trials = {7: 2000, 8: 1000, 9: 500, 10: 200}[n]
    max_fc = 0
    fc_dist = Counter()

    for _ in range(trials):
        A = random_tournament(n)
        cycles = []
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if A[i][j] and A[j][k] and A[k][i]:
                        cycles.append((i, j, k))
                    elif A[i][k] and A[k][j] and A[j][i]:
                        cycles.append((i, k, j))

        c3 = len(cycles)
        if c3 == 0:
            fc_dist[0] += 1
            continue

        # Components
        cycle_edge_sets = []
        for (i, j, k) in cycles:
            es = set()
            for x, y in [(i,j),(j,k),(k,i),(j,i),(k,j),(i,k)]:
                if A[x][y]:
                    es.add((x,y))
            cycle_edge_sets.append(es)

        adj_list = [[] for _ in range(c3)]
        for p in range(c3):
            for q in range(p+1, c3):
                if cycle_edge_sets[p] & cycle_edge_sets[q]:
                    adj_list[p].append(q)
                    adj_list[q].append(p)

        visited = set()
        components = []
        for start in range(c3):
            if start in visited:
                continue
            comp = []
            stack = [start]
            visited.add(start)
            while stack:
                u = stack.pop()
                for v in adj_list[u]:
                    if v not in visited:
                        visited.add(v)
                        stack.append(v)
                comp.append(u)
            components.append(comp)

        fc = 0
        for comp in components:
            all_free = True
            for idx2 in comp:
                a, b, c = cycles[idx2]
                if not is_free_cycle(A, n, a, b, c):
                    all_free = False
                    break
            if all_free:
                fc += 1

        fc_dist[fc] += 1
        max_fc = max(max_fc, fc)

    print(f"n={n}: max_fc={max_fc}, dist={dict(sorted(fc_dist.items()))}")


# ============================================================
# Part 3: Can we prove at most 1 free component via the bridge theorem?
# ============================================================
print(f"\n{'=' * 70}")
print("PROVING AT MOST 1 FREE COMPONENT")
print("=" * 70)

# ARGUMENT:
# Suppose C1 and C2 are free cycles in different components.
# C1 = a->b->c->a, C2 = d->e->f->d.
#
# By THM-106: every external vertex creates a bridge to C1.
# So d, e, f all create bridges to C1. These bridges are in C1's component.
#
# Now: d creates bridge B_d sharing some edge of C1.
# Say B_d shares edge a->b: B_d = d->a->b->d (pattern (1,0,*)).
# So d->a and b->d.
#
# e creates bridge B_e sharing some edge of C1.
# f creates bridge B_f sharing some edge of C1.
#
# We need to show: C2 is connected to C1's component.
# I.e., some cycle sharing an edge with C2 is in C1's component.
#
# C2 has edges d->e, e->f, f->d.
# For edge d->e: if {d,e,u} is a 3-cycle for some u in V(C1),
#   it shares edge (d,e) or (e,d) with C2,
#   and shares vertex pair (d,u) or (e,u) with a bridge cycle.
#   Since bridge cycles share edges with C1, and the new cycle
#   shares a vertex pair with a bridge cycle, the new cycle
#   shares a directed edge with the bridge cycle (hence in C1's component).
#   So C2 is connected to C1's component.
#
# For {d,e,u} NOT a 3-cycle for any u in V(C1):
#   All {d,e,a}, {d,e,b}, {d,e,c} are transitive.
#   With d->e (given):
#   - {d,e,u} transitive means d dom (d->e, d->u) or e dom or u dom
#   - d dom: d->u. e dom: e->d (impossible, d->e). u dom: u->d, u->e.
#   So either d->u or (u->d AND u->e).

# Let me check: for edge d->e of C2 and all u in {a,b,c},
# {d,e,u} transitive implies: d->u OR (u->d AND u->e).

# Similarly for edges e->f and f->d of C2.

# Let me enumerate: at n=6 with free C1 and disconnected C2,
# what are the bridge patterns and transitivity?

print("\nDetailed analysis of disconnected free+dominated cases:")
n = 6
checked = 0
for bits in range(2**(n*(n-1)//2)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    # Find cycles and components
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles.append((i, j, k))
                elif A[i][k] and A[k][j] and A[j][i]:
                    cycles.append((i, k, j))

    c3 = len(cycles)
    if c3 < 2:
        continue

    cycle_edge_sets = []
    for (i, j, k) in cycles:
        es = set()
        for x, y in [(i,j),(j,k),(k,i),(j,i),(k,j),(i,k)]:
            if A[x][y]:
                es.add((x,y))
        cycle_edge_sets.append(es)

    adj_list = [[] for _ in range(c3)]
    for p in range(c3):
        for q in range(p+1, c3):
            if cycle_edge_sets[p] & cycle_edge_sets[q]:
                adj_list[p].append(q)
                adj_list[q].append(p)

    visited = set()
    components = []
    for start in range(c3):
        if start in visited:
            continue
        comp = []
        stack = [start]
        visited.add(start)
        while stack:
            u = stack.pop()
            for v in adj_list[u]:
                if v not in visited:
                    visited.add(v)
                    stack.append(v)
            comp.append(u)
        components.append(comp)

    if len(components) < 2:
        continue

    # Identify free and dominated components
    free_comp_indices = []
    dom_comp_indices = []
    for ci, comp in enumerate(components):
        all_free = True
        for idx2 in comp:
            a, b, c = cycles[idx2]
            if not is_free_cycle(A, n, a, b, c):
                all_free = False
                break
        if all_free:
            free_comp_indices.append(ci)
        else:
            dom_comp_indices.append(ci)

    if not free_comp_indices or not dom_comp_indices:
        continue

    checked += 1

    # For each dominated component, find the dominant vertex
    for dci in dom_comp_indices:
        for idx2 in components[dci]:
            cyc = cycles[idx2]
            x, y, z = cyc
            for d in range(n):
                if d in (x, y, z):
                    continue
                if A[d][x] and A[d][y] and A[d][z]:
                    # d dominates all of cycle. Is d in C1 or external?
                    fci = free_comp_indices[0]
                    fc_verts = set()
                    for fi in components[fci]:
                        fc_verts.update(cycles[fi])
                    if d in fc_verts:
                        pass  # d is in the free cycle
                    break
                if A[x][d] and A[y][d] and A[z][d]:
                    break

if checked > 0:
    print(f"  Checked {checked} tournaments with free+dominated disconnected components")


# ============================================================
# Part 4: Direct structural analysis
# ============================================================
print(f"\n{'=' * 70}")
print("STRUCTURAL: WHO DOMINATES THE ISOLATED CYCLES?")
print("=" * 70)

n = 6
dominator_location = Counter()

for bits in range(2**(n*(n-1)//2)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles.append((i, j, k))
                elif A[i][k] and A[k][j] and A[j][i]:
                    cycles.append((i, k, j))

    c3 = len(cycles)
    if c3 < 2:
        continue

    # Quick check: which cycles are dominated?
    dominated_cycles = []
    free_cycles = []
    for cyc in cycles:
        x, y, z = cyc
        is_dom = False
        for d in range(n):
            if d in (x, y, z):
                continue
            if (A[d][x] and A[d][y] and A[d][z]) or (A[x][d] and A[y][d] and A[z][d]):
                is_dom = True
                break
        if is_dom:
            dominated_cycles.append(cyc)
        else:
            free_cycles.append(cyc)

    if not free_cycles:
        continue

    # For dominated cycles: find dominator, check if dominator is in a free cycle
    for dc in dominated_cycles:
        x, y, z = dc
        for d in range(n):
            if d in (x, y, z):
                continue
            if A[d][x] and A[d][y] and A[d][z]:
                # d dominates. Is d in a free cycle?
                d_in_free = any(d in fc for fc in free_cycles)
                dominator_location[('dom', d_in_free)] += 1
                break
            if A[x][d] and A[y][d] and A[z][d]:
                d_in_free = any(d in fc for fc in free_cycles)
                dominator_location[('sub', d_in_free)] += 1
                break

print(f"n=6: Dominator location for dominated cycles:")
for key, cnt in sorted(dominator_location.items()):
    print(f"  {key}: {cnt}")

# This tells us: is the dominant vertex always in the free cycle(s)?
# If so, the free cycle "absorbs" the dominator, and the dominated
# cycles are isolated precisely because the dominator is "inside" the free component.


# ============================================================
# Part 5: At most 1 free component — proof sketch
# ============================================================
print(f"\n{'=' * 70}")
print("PROOF SKETCH")
print("=" * 70)

print("""
THEOREM: In any tournament T on n >= 4 vertices, b1(T) <= 1.

PROOF SKETCH:

Step 1: b1 = #{free components of 3-cycle graph}
  (Verified computationally n<=8, follows from THM-122 + THM-104 + THM-105)

Step 2: At most 1 free component.
  Suppose for contradiction that FC1 and FC2 are two free components.

  By THM-106 (bridge theorem): for any free cycle C in FC1 and any
  vertex v not in C, there exists a 3-cycle containing v that shares
  an edge with C. This cycle is in FC1's component.

  Therefore: every vertex of the tournament appears in some cycle
  of FC1's component. Similarly for FC2.

  But FC1 and FC2 are DISJOINT components. So FC2's cycles share
  no edges with FC1's cycles. Yet all vertices of FC2's cycles are
  also in FC1's cycles (from step above). This means FC1 and FC2
  share vertices but not edges.

  CLAIM: If two 3-cycles share 2+ vertices, they share a directed edge.
  (Two vertices determine a unique directed edge in a tournament.)

  So FC2's cycles can share at most 1 vertex with any single FC1 cycle.
  But the vertex overlap ensures: each FC2 cycle vertex is in some FC1
  cycle, creating edge-sharing opportunities via intermediate cycles.

  [Remaining: complete the proof that edge-sharing chains always exist]

Step 3: Therefore b1 <= 1. Combined with Sum_v b1(T\v) <= 3 when b1=0,
  good vertex always exists, and beta_2 = 0 by induction.
""")


print("\n\nDone.")
