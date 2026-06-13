#!/usr/bin/env python3
"""
opus-2026-04-05-S24: Study structural properties of Ω(T) (conflict graph).

Questions:
1. Is Ω(T) always perfect? (perfect → real-rooted for some classes)
2. What is the chromatic number?
3. What is the clique cover number?
4. Is it chordal? (chordal → perfect → real-rooted independence poly)
5. What is the maximum degree?
6. Graph invariants that persist across all tournaments.

Focus on n=5,6 where we can compute exhaustively.
"""

from itertools import combinations, permutations
from collections import Counter, defaultdict

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(2**m):
        T = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
        yield T

def tournament_to_canon(T):
    n = len(T)
    best = None
    for perm in permutations(range(n)):
        form = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def find_all_directed_odd_cycles(T):
    n = len(T)
    cycles = set()
    for L in range(3, n+1, 2):
        for verts in combinations(range(n), L):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                path = (v0,) + perm
                ok = True
                for i in range(L):
                    if T[path[i]][path[(i+1) % L]] != 1:
                        ok = False
                        break
                if ok:
                    mi = path.index(min(path))
                    canon = tuple(path[(mi + j) % L] for j in range(L))
                    cycles.add(canon)
                    break
    return list(cycles)

def build_conflict_graph(cycles):
    k = len(cycles)
    adj = [[False]*k for _ in range(k)]
    vsets = [set(c) for c in cycles]
    for i in range(k):
        for j in range(i+1, k):
            if vsets[i] & vsets[j]:
                adj[i][j] = adj[j][i] = True
    return adj, vsets

def graph_properties(adj, vsets):
    """Compute various properties of the graph."""
    k = len(adj)
    if k == 0:
        return {'vertices': 0, 'edges': 0, 'max_degree': 0, 'max_clique': 0,
                'max_indep': 0, 'is_chordal': True, 'is_perfect': True}

    # Basic
    edges = sum(1 for i in range(k) for j in range(i+1, k) if adj[i][j])
    degrees = [sum(1 for j in range(k) if adj[i][j]) for i in range(k)]
    max_degree = max(degrees) if degrees else 0

    # Clique number (max clique)
    max_clique = 0
    for mask in range(2**k):
        verts = [i for i in range(k) if (mask >> i) & 1]
        if len(verts) <= max_clique:
            continue
        is_clique = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if not adj[verts[a]][verts[b]]:
                    is_clique = False
                    break
            if not is_clique:
                break
        if is_clique:
            max_clique = len(verts)

    # Independence number (max independent set)
    max_indep = 0
    for mask in range(2**k):
        verts = [i for i in range(k) if (mask >> i) & 1]
        if len(verts) <= max_indep:
            continue
        ok = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if adj[verts[a]][verts[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            max_indep = len(verts)

    # Chordal check: a graph is chordal iff it has no induced cycle of length ≥ 4
    # Quick check: look for induced C4
    has_c4 = False
    if k >= 4:
        for quad in combinations(range(k), 4):
            a, b, c, d = quad
            # Check all 3 possible cycle orderings of 4 vertices
            for perm in [(a,b,c,d), (a,b,d,c), (a,c,b,d)]:
                # Check if this is an induced 4-cycle
                p = perm
                edges_present = [adj[p[i]][p[(i+1)%4]] for i in range(4)]
                diag1 = adj[p[0]][p[2]]
                diag2 = adj[p[1]][p[3]]
                if all(edges_present) and not diag1 and not diag2:
                    has_c4 = True
                    break
            if has_c4:
                break

    # Check for odd holes (induced C5, C7, ...)
    has_odd_hole = False
    if k >= 5:
        for size in range(5, min(k+1, 8), 2):  # check C5, C7 only
            for verts in combinations(range(k), size):
                for perm in permutations(verts[1:]):
                    cycle = (verts[0],) + perm
                    # Check induced cycle
                    cycle_edges = all(adj[cycle[i]][cycle[(i+1)%size]] for i in range(size))
                    if not cycle_edges:
                        continue
                    # Check no chords
                    has_chord = False
                    for i in range(size):
                        for j in range(i+2, size):
                            if (j - i) % size != 1 and (i - j) % size != 1:
                                if adj[cycle[i]][cycle[j]]:
                                    has_chord = True
                                    break
                        if has_chord:
                            break
                    if not has_chord:
                        has_odd_hole = True
                        break
                if has_odd_hole:
                    break
            if has_odd_hole:
                break

    is_chordal = not has_c4  # simplified
    is_perfect = not has_odd_hole  # by strong perfect graph theorem (no odd hole or odd antihole)

    return {
        'vertices': k,
        'edges': edges,
        'max_degree': max_degree,
        'max_clique': max_clique,
        'max_indep': max_indep,
        'is_chordal': is_chordal,
        'has_c4': has_c4,
        'has_odd_hole': has_odd_hole,
    }

# ═══════════════════════════════════════════
# n=5: EXHAUSTIVE
# ═══════════════════════════════════════════
print("="*60)
print("Ω(T) properties at n=5")
print("="*60)

iso_classes_5 = {}
for T in all_tournaments(5):
    canon = tournament_to_canon(T)
    if canon not in iso_classes_5:
        iso_classes_5[canon] = T

for canon, T in iso_classes_5.items():
    cycles = find_all_directed_odd_cycles(T)
    adj, vsets = build_conflict_graph(cycles)
    props = graph_properties(adj, vsets)
    scores = tuple(sorted(sum(T[i]) for i in range(5)))
    print(f"  scores={scores}: |Ω|={props['vertices']}, edges={props['edges']}, "
          f"ω={props['max_clique']}, α={props['max_indep']}, "
          f"Δ={props['max_degree']}, chordal={props['is_chordal']}")

# ═════════���═════════════════════════════════
# n=6: EXHAUSTIVE
# ═══════════���═══════════════════════════════
print(f"\n{'='*60}")
print("Ω(T) properties at n=6")
print("="*60)

iso_classes_6 = {}
for T in all_tournaments(6):
    canon = tournament_to_canon(T)
    if canon not in iso_classes_6:
        iso_classes_6[canon] = T

prop_dist = Counter()
chordal_count = 0
not_chordal_examples = []

for canon, T in iso_classes_6.items():
    cycles = find_all_directed_odd_cycles(T)
    adj, vsets = build_conflict_graph(cycles)
    props = graph_properties(adj, vsets)
    scores = tuple(sorted(sum(T[i]) for i in range(6)))

    key = (props['vertices'], props['edges'], props['max_clique'],
           props['max_indep'], props['is_chordal'])
    prop_dist[key] += 1

    if props['is_chordal']:
        chordal_count += 1
    elif len(not_chordal_examples) < 3:
        not_chordal_examples.append((scores, props))

print(f"Total iso classes: {len(iso_classes_6)}")
print(f"Chordal: {chordal_count}/{len(iso_classes_6)}")

if not_chordal_examples:
    print(f"\nNon-chordal examples:")
    for scores, props in not_chordal_examples:
        print(f"  scores={scores}: |Ω|={props['vertices']}, ω={props['max_clique']}, "
              f"α={props['max_indep']}, has_C4={props['has_c4']}")

# Summary stats
print(f"\n(|V|, |E|, ω, α, chordal) distribution:")
for key, count in sorted(prop_dist.items()):
    v, e, omega, alpha, chordal = key
    print(f"  |Ω|={v:2d}, |E|={e:3d}, ω={omega}, α={alpha}, chordal={chordal}: {count} classes")

# KEY: check if ω(Ω) = χ(Ω) always (perfect graph)
# For perfect graphs, clique number = chromatic number
# At n≤6 with small Ω, this is easy to check
print(f"\nClique number ω = independence number α relationship:")
omega_alpha = Counter()
for canon, T in iso_classes_6.items():
    cycles = find_all_directed_odd_cycles(T)
    adj, vsets = build_conflict_graph(cycles)
    props = graph_properties(adj, vsets)
    omega_alpha[(props['max_clique'], props['max_indep'])] += 1

for (w, a), count in sorted(omega_alpha.items()):
    print(f"  (ω, α) = ({w}, {a}): {count} classes")

# Density of Ω
print(f"\nΩ(T) density:")
for canon, T in list(iso_classes_6.items())[:5]:
    cycles = find_all_directed_odd_cycles(T)
    if len(cycles) > 1:
        adj, vsets = build_conflict_graph(cycles)
        k = len(cycles)
        edges = sum(1 for i in range(k) for j in range(i+1, k) if adj[i][j])
        density = 2*edges / (k*(k-1)) if k > 1 else 0
        scores = tuple(sorted(sum(T[i]) for i in range(6)))
        print(f"  scores={scores}: |Ω|={k}, density={density:.3f}")

print("\n\nDONE.")
