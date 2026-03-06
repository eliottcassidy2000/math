"""
Class transition graph: which classes are connected by single tile flips?

A single tile flip = flipping one bit in {0,1}^m = reversing one arc.
This changes the tournament and possibly its isomorphism class.

Key questions:
1. What's the structure of the class transition graph?
2. Are there "islands" (disconnected components)?
3. How does H(T) change under single flips? (connects to ΔH formula)
4. What fraction of flips STAY in the same class vs leave?
"""
from itertools import permutations
from collections import defaultdict, Counter

def build_full(n):
    arcs = [(i,j) for i in range(n) for j in range(i+2, n)]
    m = len(arcs)

    results = {}
    for bits in range(1 << m):
        adj = [[False]*n for _ in range(n)]
        for i in range(n-1):
            adj[i][i+1] = True
        for k, (i,j) in enumerate(arcs):
            if bits & (1 << k):
                adj[j][i] = True
            else:
                adj[i][j] = True
        canon = None
        for perm in permutations(range(n)):
            mat = tuple(tuple(adj[perm[a]][perm[b]] for b in range(n)) for a in range(n))
            if canon is None or mat < canon:
                canon = mat
        results[bits] = canon

    return results, arcs, m

def analyze_transitions(n):
    print(f"\n{'='*70}")
    print(f"CLASS TRANSITION GRAPH: n = {n}, m = {(n-1)*(n-2)//2}")
    print(f"{'='*70}")

    results, arcs, m = build_full(n)

    classes = defaultdict(list)
    for bits, canon in results.items():
        classes[canon].append(bits)

    clist = list(classes.keys())
    cidx = {c: i for i, c in enumerate(clist)}
    nc = len(clist)

    # Compute H(T) for each class
    class_h = {}
    for canon in clist:
        adj = [list(row) for row in canon]
        nv = len(adj)
        h = sum(1 for p in permutations(range(nv))
                if all(adj[p[i]][p[i+1]] for i in range(nv-1)))
        class_h[canon] = h

    # Build transition graph: for each tiling, flip each bit, see what class we land in
    transitions = defaultdict(Counter)  # transitions[class_i][class_j] = count
    self_loops = defaultdict(int)

    for bits, canon in results.items():
        ci = cidx[canon]
        for k in range(m):
            flipped = bits ^ (1 << k)
            fcanon = results[flipped]
            cj = cidx[fcanon]
            if ci == cj:
                self_loops[ci] += 1
            else:
                transitions[ci][cj] += 1

    # Transition graph edges (undirected: merge i->j and j->i)
    edges = set()
    edge_weight = {}
    for ci, targets in transitions.items():
        for cj, cnt in targets.items():
            e = (min(ci, cj), max(ci, cj))
            edges.add(e)
            if e not in edge_weight:
                edge_weight[e] = 0
            edge_weight[e] += cnt

    print(f"\nClasses: {nc}")
    print(f"Transition edges: {len(edges)}")
    print(f"Graph density: {2*len(edges)/(nc*(nc-1)):.3f}" if nc > 1 else "")

    # Is the graph connected?
    # Simple BFS
    if nc > 1:
        visited = {0}
        queue = [0]
        while queue:
            node = queue.pop(0)
            for ci, targets in transitions.items():
                if ci == node:
                    for cj in targets:
                        if cj not in visited:
                            visited.add(cj)
                            queue.append(cj)
            for ci, targets in transitions.items():
                for cj in targets:
                    if cj == node and ci not in visited:
                        visited.add(ci)
                        queue.append(ci)
        print(f"Connected: {len(visited) == nc} ({len(visited)}/{nc} reachable)")

    # Self-loop fraction: what fraction of flips stay in same class?
    total_flips = sum(len(classes[c]) * m for c in clist)
    total_self = sum(self_loops.values())
    print(f"\nSelf-loop flips: {total_self}/{total_flips} = {total_self/total_flips:.4f}")

    # Per class: self-loop fraction and degree
    print(f"\nPer-class statistics (sorted by H):")
    for canon in sorted(clist, key=lambda c: class_h[c]):
        ci = cidx[canon]
        size = len(classes[canon])
        h = class_h[canon]
        sl = self_loops.get(ci, 0)
        total = size * m
        sl_frac = sl / total if total > 0 else 0
        neighbors = len(transitions[ci])
        # ΔH for each neighbor
        dh_values = []
        for cj, cnt in transitions[ci].items():
            h2 = class_h[clist[cj]]
            dh_values.append((h2 - h, cnt))
        dh_str = ', '.join(f"Δ{d:+d}×{c}" for d, c in sorted(dh_values)[:6])
        if len(dh_values) > 6:
            dh_str += f", ... ({len(dh_values)} total)"
        print(f"  H={h:3d}, size={size:4d}, self={sl_frac:.3f}, "
              f"nbrs={neighbors:2d}: {dh_str}")

    # ΔH distribution across ALL single-arc flips
    print(f"\nGlobal ΔH distribution under single flips:")
    dh_global = Counter()
    for bits, canon in results.items():
        h1 = class_h[canon]
        for k in range(m):
            flipped = bits ^ (1 << k)
            h2 = class_h[results[flipped]]
            dh_global[h2 - h1] += 1
    for dh in sorted(dh_global.keys()):
        cnt = dh_global[dh]
        print(f"  ΔH = {dh:+3d}: {cnt:6d} ({cnt/total_flips:.4f})")

    # Which arcs cause which ΔH patterns?
    print(f"\nΔH by arc position:")
    for k, (i, j) in enumerate(arcs):
        dh_arc = Counter()
        for bits, canon in results.items():
            h1 = class_h[canon]
            flipped = bits ^ (1 << k)
            h2 = class_h[results[flipped]]
            dh_arc[h2 - h1] += 1
        total_arc = 1 << m
        avg_dh = sum(d * c for d, c in dh_arc.items()) / total_arc
        max_dh = max(dh_arc.keys())
        min_dh = min(dh_arc.keys())
        print(f"  arc ({i},{j}), gap={j-i}: avg_ΔH={avg_dh:+.3f}, "
              f"range=[{min_dh:+d}, {max_dh:+d}], "
              f"ΔH=0: {dh_arc.get(0,0)/total_arc:.3f}")

if __name__ == '__main__':
    for n in range(3, 7):
        analyze_transitions(n)
