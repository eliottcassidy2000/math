"""
asymmetry_of_flips.py — opus-2026-04-04-S15

THE FUNDAMENTAL ASYMMETRY: at n=3, the cyclic tournament is fragile
(every flip → transitive) while the transitive is robust (only 1/3 flips → cyclic).

EXTEND THIS TO ALL n:
1. For each iso class, count: how many arc flips LEAVE this class vs ESCAPE it?
2. "Fragility" = fraction of flips that change the iso class
3. "Attractiveness" = how many flips from OTHER classes lead INTO this class
4. The transitive should be maximally attractive and minimally fragile
5. The H-maximizer should be maximally fragile (any perturbation reduces H)

ABSTRACT QUESTIONS:
- Is the metagraph directed? (flip from A→B doesn't imply flip from B→A)
- Is there a "gradient" structure? (most flips go downhill in H?)
- Which classes are "sinks" (attractors) and "sources" (repellers)?
- Connection to the Fixed/Boundary/Free framework
"""
import numpy as np
from itertools import permutations
from collections import defaultdict, Counter
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import make_tiles, tiling_to_tournament, count_hamiltonian_paths

def tournament_canonical(T, n):
    best = None
    for perm in permutations(range(n)):
        key = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(i+1, n))
        if best is None or key < best:
            best = key
    return best

def compute_flip_asymmetry(n):
    """For each iso class, compute flip-in and flip-out counts."""
    print(f"\n{'='*60}")
    print(f"FLIP ASYMMETRY AT n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    N = 1 << m

    # Compute iso class for each tiling
    tiling_class = {}
    class_H = {}
    class_tilings = defaultdict(list)

    for t in range(N):
        T = tiling_to_tournament(n, tiles, t)
        canon = tournament_canonical(T, n)
        tiling_class[t] = canon
        class_tilings[canon].append(t)
        if canon not in class_H:
            class_H[canon] = count_hamiltonian_paths(T, n)

    V = len(class_tilings)
    print(f"  {V} iso classes, {m} tiles, {N} tilings")

    # For each tiling, flip each tile and record the class transition
    # This gives the DIRECTED metagraph with multiplicities

    # Per-class statistics
    class_self_loops = defaultdict(int)  # flips that stay in class
    class_exits = defaultdict(int)  # flips that leave class
    class_entries = defaultdict(int)  # flips from other classes INTO this class
    class_exits_to = defaultdict(lambda: defaultdict(int))  # class_A → class_B count

    for t in range(N):
        src_class = tiling_class[t]
        for k in range(m):
            t_flipped = t ^ (1 << k)
            dst_class = tiling_class[t_flipped]
            if src_class == dst_class:
                class_self_loops[src_class] += 1
            else:
                class_exits[src_class] += 1
                class_entries[dst_class] += 1
                class_exits_to[src_class][dst_class] += 1

    # Compute fragility and attractiveness
    print(f"\n  {'H':>4s} {'size':>5s} {'self':>6s} {'exit':>6s} {'entry':>6s} "
          f"{'fragility':>10s} {'attract':>10s} {'net_flow':>10s}")

    results = []
    for canon in sorted(class_tilings.keys(), key=lambda c: class_H[c]):
        h = class_H[canon]
        size = len(class_tilings[canon])
        self_loops = class_self_loops[canon]
        exits = class_exits[canon]
        entries = class_entries[canon]
        total_flips = self_loops + exits  # = size * m

        fragility = exits / total_flips if total_flips > 0 else 0
        attractiveness = entries / (exits + entries) if exits + entries > 0 else 0
        net_flow = entries - exits  # positive = net attractor

        results.append({
            'canon': canon, 'H': h, 'size': size,
            'self': self_loops, 'exits': exits, 'entries': entries,
            'fragility': fragility, 'attract': attractiveness,
            'net_flow': net_flow
        })

        if V <= 60:
            print(f"  {h:4d} {size:5d} {self_loops:6d} {exits:6d} {entries:6d} "
                  f"{fragility:10.4f} {attractiveness:10.4f} {net_flow:10d}")

    # KEY ANALYSIS
    print(f"\n  KEY OBSERVATIONS:")

    # 1. Is the transitive (H=1) the most attractive?
    transitive = [r for r in results if r['H'] == 1][0]
    max_attract = max(r['attract'] for r in results)
    print(f"\n  Transitive (H=1): fragility={transitive['fragility']:.4f}, "
          f"attract={transitive['attract']:.4f}")
    print(f"  Most attractive class: H={max(results, key=lambda r: r['attract'])['H']}")

    # 2. Is the H-maximizer the most fragile?
    max_H_class = max(results, key=lambda r: r['H'])
    most_fragile = max(results, key=lambda r: r['fragility'])
    print(f"\n  H-maximizer (H={max_H_class['H']}): fragility={max_H_class['fragility']:.4f}")
    print(f"  Most fragile class: H={most_fragile['H']}, "
          f"fragility={most_fragile['fragility']:.4f}")

    # 3. Net flow analysis: which classes are net attractors vs repellers?
    attractors = [r for r in results if r['net_flow'] > 0]
    repellers = [r for r in results if r['net_flow'] < 0]
    neutral = [r for r in results if r['net_flow'] == 0]
    print(f"\n  Net attractors: {len(attractors)} classes "
          f"(H values: {sorted(r['H'] for r in attractors)[:10]})")
    print(f"  Net repellers: {len(repellers)} classes "
          f"(H values: {sorted(r['H'] for r in repellers)[:10]})")
    print(f"  Neutral: {len(neutral)} classes")

    # 4. Correlation between H and fragility
    h_vals = [r['H'] for r in results]
    frag_vals = [r['fragility'] for r in results]
    if len(set(h_vals)) > 1:
        corr = np.corrcoef(h_vals, frag_vals)[0, 1]
        print(f"\n  Corr(H, fragility) = {corr:.4f}")

    # 5. Correlation between H and net flow
    nf_vals = [r['net_flow'] for r in results]
    corr_nf = np.corrcoef(h_vals, nf_vals)[0, 1]
    print(f"  Corr(H, net_flow) = {corr_nf:.4f}")

    # 6. The n=3 observation extended: for the transitive, which arc flips
    #    give which destination?
    if n <= 6:
        print(f"\n  TRANSITIVE FLIP DESTINATIONS:")
        trans_tilings = class_tilings[[c for c in class_tilings if class_H[c] == 1][0]]
        for t in trans_tilings[:1]:  # just one representative
            T = tiling_to_tournament(n, tiles, t)
            for k in range(m):
                tile = tiles[k]
                skip = tile[0] - tile[1]
                t_flip = t ^ (1 << k)
                dst_h = class_H[tiling_class[t_flip]]
                print(f"    Flip tile {tile} (skip={skip}): → H={dst_h}")

    # 7. The DIRECTED flow: which H-transitions dominate?
    print(f"\n  DIRECTED FLOW (class transitions by H change):")
    h_transitions = Counter()
    for src_canon, destinations in class_exits_to.items():
        src_h = class_H[src_canon]
        for dst_canon, count in destinations.items():
            dst_h = class_H[dst_canon]
            delta_h = dst_h - src_h
            h_transitions[delta_h] += count

    for dh in sorted(h_transitions.keys()):
        count = h_transitions[dh]
        direction = "↑" if dh > 0 else ("↓" if dh < 0 else "=")
        frac = count / sum(h_transitions.values())
        print(f"    ΔH={dh:+4d}: {count:8d} transitions ({frac:.3f}) {direction}")

    return results

# ==============================================================
# Run for small n
# ==============================================================
if __name__ == "__main__":
    for n in [3, 4, 5, 6]:
        compute_flip_asymmetry(n)
