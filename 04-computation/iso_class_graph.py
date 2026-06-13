#!/usr/bin/env python3
"""
iso_class_graph.py — opus-2026-03-13-S67k
Full investigation of tournament isomorphism class structure.

For each n=3..7:
- Enumerate ALL isomorphism classes
- Compute: H, score_seq, c3 (3-cycles), c5, alpha1, alpha2
- Determine: blueself, blackself, self-complement, grid-symmetric
- Build the FLIP GRAPH: iso classes connected by single arc reversal
- Study how structure changes as n increases
- Look for self-similar/fractal patterns
"""

from itertools import combinations, permutations
from collections import defaultdict, Counter
import sys

def tournament_from_bits(n, bits):
    """Create adjacency matrix from bit encoding of upper triangle."""
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def count_ham_paths(A, n):
    """Count Hamiltonian paths via DP with bitmask."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_3cycles(A, n):
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] + A[j][i] == 1:  # always true
                    # Check if i->j->k->i or i->k->j->i
                    if A[i][j] and A[j][k] and A[k][i]:
                        c += 1
                    if A[i][k] and A[k][j] and A[j][i]:
                        c += 1
    return c

def count_odd_cycles(A, n):
    """Count all directed odd cycles (lengths 3,5,...,n or n-1)."""
    total = 0
    for length in range(3, n+1, 2):
        for combo in combinations(range(n), length):
            for perm in permutations(combo):
                is_cycle = True
                for idx in range(length):
                    if not A[perm[idx]][perm[(idx+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    total += 1
    # Each directed cycle of length L is counted L times (cyclic rotations)
    # Actually permutations gives all orderings, cycles of length L have L rotations
    # so we divide by L for each length... but we want total directed cycles
    # A directed cycle i1->i2->...->iL->i1 is found L times in permutations
    result = 0
    for length in range(3, n+1, 2):
        count_L = 0
        for combo in combinations(range(n), length):
            for perm in permutations(combo):
                is_cycle = True
                for idx in range(length):
                    if not A[perm[idx]][perm[(idx+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    count_L += 1
        result += count_L // length
    return result

def count_disjoint_pairs(A, n):
    """Count vertex-disjoint pairs of odd cycles (alpha_2)."""
    # For small n, enumerate all odd cycles then check disjointness
    cycles = []
    for length in range(3, n+1, 2):
        for combo in combinations(range(n), length):
            vset = frozenset(combo)
            for perm in permutations(combo):
                is_cycle = True
                for idx in range(length):
                    if not A[perm[idx]][perm[(idx+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.append(vset)
                    break  # one representative per vertex set per combo
    # Deduplicate by vertex set (a vertex set can support multiple directed cycles)
    # Actually we want: for each vertex set, how many distinct directed cycles exist?
    # For alpha_2 we count pairs of vertex-disjoint cycles
    # Simplification: unique vertex sets with at least one directed cycle
    unique_vsets = list(set(cycles))
    count = 0
    for i in range(len(unique_vsets)):
        for j in range(i+1, len(unique_vsets)):
            if not unique_vsets[i] & unique_vsets[j]:
                count += 1
    return count

def canonical_form(A, n):
    """Canonical form of tournament under vertex permutation."""
    best = None
    for perm in permutations(range(n)):
        form = []
        for i in range(n):
            for j in range(i+1, n):
                form.append(A[perm[i]][perm[j]])
        form = tuple(form)
        if best is None or form < best:
            best = form
        # Also try with all edges reversed (complement)
    return best

def complement(A, n):
    """Reverse all arcs."""
    B = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                B[i][j] = 1 - A[i][j]
    return B

def flip_tournament(A, n):
    """Flip = reverse all arcs (complement/T^op)."""
    return complement(A, n)

def is_grid_symmetric(A, n):
    """Grid symmetric: A[i][j] = A[n-1-j][n-1-i] for all i,j."""
    for i in range(n):
        for j in range(n):
            if i != j:
                ii, jj = n-1-j, n-1-i
                if A[i][j] != A[ii][jj]:
                    return False
    return True

def automorphism_count(A, n):
    """Count automorphisms of tournament."""
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(i+1, n):
                if A[i][j] != A[perm[i]][perm[j]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            count += 1
    return count

def pos_vector(A, n):
    """Compute M[v,v] = sum_P (-1)^{pos(v,P)} for each vertex v."""
    # DP to track position assignments
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = {v: 1}  # pos(v) = 0 in this partial path

    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                contributions = defaultdict(int)
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        if (prev_mask, u) in dp:
                            for vertex, coeff in dp[(prev_mask, u)].items():
                                contributions[vertex] += coeff
                            # vertex v is at position mask_size-1
                            sign = (-1) ** (mask_size - 1)
                            contributions[v] += dp.get((prev_mask, u), {}).get('_count', 0) if False else 0
                # This approach is getting complicated. Use simpler method.
                pass

    # Simpler: enumerate all Ham paths, compute pos
    # For n<=7, this is feasible with DP that tracks the actual path
    # Actually let's just use the DP to count, and separately compute pos_vector
    # via a different method

    # For n<=6, we can enumerate paths via backtracking
    M = [0] * n

    def backtrack(path, visited):
        if len(path) == n:
            for pos_idx, v in enumerate(path):
                M[v] += (-1) ** pos_idx
            return
        last = path[-1]
        for v in range(n):
            if not visited[v] and A[last][v]:
                visited[v] = True
                path.append(v)
                backtrack(path, visited)
                path.pop()
                visited[v] = False

    for start in range(n):
        visited = [False] * n
        visited[start] = True
        backtrack([start], visited)

    return tuple(M)

def single_flip_neighbors(A, n):
    """Generate all tournaments reachable by flipping one arc."""
    neighbors = []
    for i in range(n):
        for j in range(i+1, n):
            B = [row[:] for row in A]
            B[i][j], B[j][i] = B[j][i], B[i][j]
            neighbors.append(B)
    return neighbors

print("=" * 70)
print("TOURNAMENT ISOMORPHISM CLASS GRAPH")
print("=" * 70)

for n in range(3, 8):
    m = n * (n-1) // 2
    num_tournaments = 1 << m

    if n >= 8:
        print(f"\nn={n}: skipping (too large)")
        continue

    print(f"\n{'='*70}")
    print(f"n = {n}, m = {m}, total tournaments = {num_tournaments}")
    print(f"{'='*70}")

    # Step 1: Enumerate all tournaments and group by isomorphism class
    class_of = {}  # bits -> canonical form
    classes = defaultdict(list)  # canonical -> list of bit encodings

    for bits in range(num_tournaments):
        A = tournament_from_bits(n, bits)
        cf = canonical_form(A, n)
        class_of[bits] = cf
        classes[cf].append(bits)

    iso_classes = sorted(classes.keys())
    num_classes = len(iso_classes)
    class_index = {cf: i for i, cf in enumerate(iso_classes)}

    print(f"Number of isomorphism classes: {num_classes}")

    # Step 2: Compute properties of each class
    class_props = []
    for idx, cf in enumerate(iso_classes):
        bits = classes[cf][0]  # representative
        A = tournament_from_bits(n, bits)

        H = count_ham_paths(A, n)
        ss = score_seq(A, n)
        c3 = count_3cycles(A, n)
        aut = automorphism_count(A, n)
        gs = is_grid_symmetric(A, n)
        class_size = len(classes[cf])

        # Check self-complement
        Ac = complement(A, n)
        cf_comp = canonical_form(Ac, n)
        is_sc = (cf_comp == cf)

        # Check blueself/blackself
        # Self-flip: flip(T) isomorphic to T
        is_self_flip = (cf_comp == cf)  # flip = complement for tournaments
        is_blueself = is_self_flip and gs
        is_blackself = is_self_flip and not gs

        # Pos vector (skip for n=7, too slow with backtracking)
        if n <= 6:
            pv = pos_vector(A, n)
            pv_sorted = tuple(sorted(pv))
        else:
            pv = None
            pv_sorted = None

        # alpha_1 (total odd cycles) — skip for n>=7 (too slow)
        if n <= 6:
            alpha1 = count_odd_cycles(A, n)
        else:
            alpha1 = None

        # alpha_2 (disjoint pairs) — skip for n>=7
        if n <= 6:
            alpha2 = count_disjoint_pairs(A, n)
        else:
            alpha2 = None

        props = {
            'idx': idx,
            'cf': cf,
            'H': H,
            'score': ss,
            'c3': c3,
            'aut': aut,
            'size': class_size,
            'gs': gs,
            'sc': is_sc,
            'self_flip': is_self_flip,
            'blueself': is_blueself,
            'blackself': is_blackself,
            'pos': pv_sorted,
            'alpha1': alpha1,
            'alpha2': alpha2,
        }
        class_props.append(props)

    # Sort by H for display
    class_props.sort(key=lambda p: (p['H'], p['score']))

    print(f"\nIsomorphism classes (sorted by H):")
    print(f"{'Cls':>4} {'H':>6} {'Score':>20} {'c3':>4} {'α1':>5} {'α2':>4} {'|Aut|':>5} {'Size':>5} {'GS':>3} {'SC':>3} {'BS':>3} {'BkS':>3} {'Pos':>20}")
    for p in class_props:
        pos_str = str(p['pos']) if p['pos'] else '—'
        a1_str = str(p['alpha1']) if p['alpha1'] is not None else '—'
        a2_str = str(p['alpha2']) if p['alpha2'] is not None else '—'
        bs_str = 'Y' if p['blueself'] else ''
        bks_str = 'Y' if p['blackself'] else ''
        print(f"{p['idx']:4d} {p['H']:6d} {str(p['score']):>20} {p['c3']:4d} {a1_str:>5} {a2_str:>4} {p['aut']:5d} {p['size']:5d} {'Y' if p['gs'] else '':>3} {'Y' if p['sc'] else '':>3} {bs_str:>3} {bks_str:>3} {pos_str:>20}")

    # Step 3: Build flip graph (iso classes connected by single arc reversal)
    print(f"\nFlip graph (adjacency between iso classes):")
    adj = defaultdict(set)
    for cf in iso_classes:
        bits = classes[cf][0]
        A = tournament_from_bits(n, bits)
        neighbors = single_flip_neighbors(A, n)
        for B in neighbors:
            cf_B = canonical_form(B, n)
            if cf_B != cf:
                i = class_index[cf]
                j = class_index[cf_B]
                adj[i].add(j)
                adj[j].add(i)

    for i, p in enumerate(class_props):
        ci = class_index[p['cf']]
        nbrs = sorted(adj[ci])
        nbr_H = [class_props[class_index[iso_classes[j]]]['H'] if j < len(iso_classes) else '?' for j in nbrs]
        # Map to the sorted index
        nbr_info = []
        for j in nbrs:
            for p2 in class_props:
                if class_index[p2['cf']] == j:
                    nbr_info.append(f"H={p2['H']}")
                    break
        print(f"  Class {p['idx']} (H={p['H']}): neighbors = [{', '.join(nbr_info)}]")

    # Step 4: Score class decomposition
    score_groups = defaultdict(list)
    for p in class_props:
        score_groups[p['score']].append(p)

    print(f"\nScore class decomposition:")
    for ss in sorted(score_groups.keys()):
        group = score_groups[ss]
        H_vals = [p['H'] for p in group]
        sc_count = sum(1 for p in group if p['sc'])
        bs_count = sum(1 for p in group if p['blueself'])
        print(f"  {ss}: {len(group)} classes, H = {H_vals}, SC={sc_count}, blueself={bs_count}")

    # Step 5: H distribution
    H_vals = [p['H'] for p in class_props]
    print(f"\nH values: min={min(H_vals)}, max={max(H_vals)}, unique={len(set(H_vals))}")
    print(f"H distribution: {sorted(Counter(H_vals).items())}")

    # Step 6: Blueself/blackself summary
    bs_classes = [p for p in class_props if p['blueself']]
    bks_classes = [p for p in class_props if p['blackself']]
    sf_classes = [p for p in class_props if p['self_flip']]
    print(f"\nSelf-flip: {len(sf_classes)}, Blueself: {len(bs_classes)}, Blackself: {len(bks_classes)}")
    if bs_classes:
        print(f"  Blueself classes: H = {[p['H'] for p in bs_classes]}, scores = {[p['score'] for p in bs_classes]}")
    if bks_classes:
        print(f"  Blackself classes: H = {[p['H'] for p in bks_classes]}, scores = {[p['score'] for p in bks_classes]}")

    # Step 7: Pos analysis
    if n <= 6:
        print(f"\nPos vector analysis:")
        pos_uniform = [p for p in class_props if p['pos'] and len(set(p['pos'])) == 1]
        print(f"  Pos-uniform classes: {len(pos_uniform)} / {num_classes}")
        for p in class_props:
            if p['pos']:
                uniform = 'UNIFORM' if len(set(p['pos'])) == 1 else f'range={max(p["pos"])-min(p["pos"])}'
                print(f"    H={p['H']:4d}, pos={p['pos']}, {uniform}")

print("\n" + "=" * 70)
print("CROSS-n STRUCTURE ANALYSIS")
print("=" * 70)

# Collect all class counts
for n in range(3, 8):
    m = n*(n-1)//2
    num_t = 1 << m
    # Count iso classes (already computed above, but let's summarize)

print("\nIsomorphism class counts by n:")
print("  n=3: 2 classes (transitive, 3-cycle)")
print("  n=4: 4 classes")
print("  n=5: 12 classes")
print("  n=6: 56 classes")
print("  n=7: 456 classes")
print("  Ratios: 4/2=2.0, 12/4=3.0, 56/12=4.67, 456/56=8.14")
print("  OEIS A000568: 1, 1, 1, 2, 4, 12, 56, 456, 6880, ...")

print("\n" + "=" * 70)
print("SELF-SIMILAR STRUCTURE SEARCH")
print("=" * 70)
print("""
Looking for groups of iso classes at n+1 that collectively behave like
a single iso class at n. Key metrics to compare:
- H value relative to maximum
- Score sequence pattern
- Flip graph connectivity
- Blueself/blackself status
- Pos-uniformity
""")
