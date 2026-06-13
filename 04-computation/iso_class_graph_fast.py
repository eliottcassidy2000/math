#!/usr/bin/env python3
"""
iso_class_graph_fast.py — opus-2026-03-13-S67k
Fast version: n=3..6 only (n=7 needs separate optimized script).
Full investigation of tournament isomorphism class structure.
"""

from itertools import combinations, permutations
from collections import defaultdict, Counter
import time

def tournament_from_bits(n, bits):
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

def bits_from_matrix(A, n):
    bits = 0
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j]:
                bits |= (1 << idx)
            idx += 1
    return bits

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def score_vec(A, n):
    """Unsorted score vector."""
    return tuple(sum(A[i]) for i in range(n))

def canonical_form(A, n):
    best = None
    for perm in permutations(range(n)):
        form = []
        for i in range(n):
            for j in range(i+1, n):
                form.append(A[perm[i]][perm[j]])
        form = tuple(form)
        if best is None or form < best:
            best = form
    return best

def count_ham_paths(A, n):
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
                if A[i][j] and A[j][k] and A[k][i]:
                    c += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    c += 1
    return c

def count_5cycles(A, n):
    """Count directed 5-cycles."""
    if n < 5:
        return 0
    c = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            ok = True
            for idx in range(5):
                if not A[perm[idx]][perm[(idx+1) % 5]]:
                    ok = False
                    break
            if ok:
                c += 1
    return c // 5

def count_all_odd_cycles(A, n):
    """Count total directed odd cycles of all lengths."""
    total = 0
    for length in range(3, n+1, 2):
        for combo in combinations(range(n), length):
            for perm in permutations(combo):
                ok = True
                for idx in range(length):
                    if not A[perm[idx]][perm[(idx+1) % length]]:
                        ok = False
                        break
                if ok:
                    total += 1
        # Divide out cyclic rotations
    # Redo properly
    total = 0
    for length in range(3, n+1, 2):
        count_L = 0
        for combo in combinations(range(n), length):
            for perm in permutations(combo):
                ok = True
                for idx in range(length):
                    if not A[perm[idx]][perm[(idx+1) % length]]:
                        ok = False
                        break
                if ok:
                    count_L += 1
        total += count_L // length
    return total

def count_disjoint_pairs(A, n):
    """Count vertex-disjoint pairs of odd cycles."""
    cycle_vsets = []
    for length in range(3, n+1, 2):
        for combo in combinations(range(n), length):
            vset = frozenset(combo)
            has_cycle = False
            for perm in permutations(combo):
                ok = True
                for idx in range(length):
                    if not A[perm[idx]][perm[(idx+1) % length]]:
                        ok = False
                        break
                if ok:
                    has_cycle = True
                    break
            if has_cycle:
                cycle_vsets.append(vset)

    count = 0
    for i in range(len(cycle_vsets)):
        for j in range(i+1, len(cycle_vsets)):
            if not cycle_vsets[i] & cycle_vsets[j]:
                count += 1
    return count

def complement_bits(n, bits):
    m = n*(n-1)//2
    return bits ^ ((1 << m) - 1)

def is_grid_symmetric(A, n):
    for i in range(n):
        for j in range(n):
            if i != j:
                ii, jj = n-1-j, n-1-i
                if A[i][j] != A[ii][jj]:
                    return False
    return True

def automorphism_count(A, n):
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
    """Compute M[v,v] = sum_P (-1)^{pos(v,P)} for each vertex."""
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

print("=" * 70)
print("TOURNAMENT ISOMORPHISM CLASS GRAPH — FAST VERSION (n≤6)")
print("=" * 70)

all_results = {}

for n in range(3, 7):
    m = n * (n-1) // 2
    num_tournaments = 1 << m
    t0 = time.time()

    print(f"\n{'='*70}")
    print(f"n = {n}, m = {m}, total tournaments = {num_tournaments}")
    print(f"{'='*70}")

    # Group by isomorphism class
    class_of = {}
    classes = defaultdict(list)

    for bits in range(num_tournaments):
        A = tournament_from_bits(n, bits)
        cf = canonical_form(A, n)
        class_of[bits] = cf
        classes[cf].append(bits)

    iso_classes = sorted(classes.keys())
    num_classes = len(iso_classes)
    class_index = {cf: i for i, cf in enumerate(iso_classes)}

    print(f"Number of isomorphism classes: {num_classes}")

    # Compute properties
    class_props = []
    for idx, cf in enumerate(iso_classes):
        bits = classes[cf][0]
        A = tournament_from_bits(n, bits)

        H = count_ham_paths(A, n)
        ss = score_seq(A, n)
        c3 = count_3cycles(A, n)
        c5 = count_5cycles(A, n)
        aut = automorphism_count(A, n)
        gs = is_grid_symmetric(A, n)
        alpha1 = count_all_odd_cycles(A, n)
        alpha2 = count_disjoint_pairs(A, n)
        pv = pos_vector(A, n)
        pv_sorted = tuple(sorted(pv))
        class_size = len(classes[cf])

        # Self-complement check
        comp_bits = complement_bits(n, bits)
        cf_comp = class_of.get(comp_bits, canonical_form(tournament_from_bits(n, comp_bits), n))
        is_sc = (cf_comp == cf)
        is_self_flip = is_sc
        is_blueself = is_self_flip and gs
        is_blackself = is_self_flip and not gs

        props = {
            'idx': idx, 'cf': cf, 'H': H, 'score': ss, 'c3': c3, 'c5': c5,
            'aut': aut, 'size': class_size, 'gs': gs, 'sc': is_sc,
            'self_flip': is_self_flip, 'blueself': is_blueself, 'blackself': is_blackself,
            'pos': pv, 'pos_sorted': pv_sorted,
            'alpha1': alpha1, 'alpha2': alpha2,
        }
        class_props.append(props)

    class_props.sort(key=lambda p: (p['H'], p['score']))

    # Reassign display indices after sorting
    for i, p in enumerate(class_props):
        p['display_idx'] = i

    # Build flip graph
    adj = defaultdict(set)
    for cf in iso_classes:
        bits = classes[cf][0]
        A = tournament_from_bits(n, bits)
        ci = class_index[cf]
        for i in range(n):
            for j in range(i+1, n):
                B = [row[:] for row in A]
                B[i][j], B[j][i] = B[j][i], B[i][j]
                cf_B = canonical_form(B, n)
                cj = class_index[cf_B]
                if ci != cj:
                    adj[ci].add(cj)
                    adj[cj].add(ci)

    # Display
    print(f"\n{'Cls':>3} {'H':>5} {'Score':>18} {'c3':>3} {'c5':>4} {'α1':>5} {'α2':>3} {'|A|':>3} {'#T':>4} {'GS':>2} {'SC':>2} {'BS':>2} {'BkS':>3} Pos")
    for p in class_props:
        bs = 'Y' if p['blueself'] else ''
        bks = 'Y' if p['blackself'] else ''
        print(f"{p['display_idx']:3d} {p['H']:5d} {str(p['score']):>18} {p['c3']:3d} {p['c5']:4d} {p['alpha1']:5d} {p['alpha2']:3d} {p['aut']:3d} {p['size']:4d} {'Y' if p['gs'] else '':>2} {'Y' if p['sc'] else '':>2} {bs:>2} {bks:>3} {p['pos']}")

    # Flip graph
    print(f"\nFlip graph (single arc reversal connects classes):")
    cf_to_display = {}
    for p in class_props:
        cf_to_display[class_index[p['cf']]] = p['display_idx']

    for p in class_props:
        ci = class_index[p['cf']]
        nbrs = sorted(adj[ci])
        nbr_strs = []
        for j in nbrs:
            for p2 in class_props:
                if class_index[p2['cf']] == j:
                    nbr_strs.append(str(p2['display_idx']))
                    break
        print(f"  {p['display_idx']} (H={p['H']}, {p['score']}): → [{', '.join(nbr_strs)}]")

    # Self-flip structure
    print(f"\nSelf-flip structure:")
    sf = [p for p in class_props if p['self_flip']]
    print(f"  Self-complement classes: {len(sf)} / {num_classes}")
    for p in sf:
        print(f"    H={p['H']}, score={p['score']}, {'blueself' if p['blueself'] else 'blackself'}")

    # OCF decomposition: H = 1 + 2*alpha1 + 4*alpha2
    print(f"\nOCF verification: H = 1 + 2·α₁ + 4·α₂")
    for p in class_props:
        predicted = 1 + 2*p['alpha1'] + 4*p['alpha2']
        match = '✓' if predicted == p['H'] else '✗'
        print(f"  H={p['H']}: 1 + 2·{p['alpha1']} + 4·{p['alpha2']} = {predicted} {match}")

    # Score class decomposition
    score_groups = defaultdict(list)
    for p in class_props:
        score_groups[p['score']].append(p)
    print(f"\nScore class decomposition:")
    for ss in sorted(score_groups.keys()):
        group = score_groups[ss]
        Hs = [p['H'] for p in group]
        print(f"  {ss}: {len(group)} class(es), H = {Hs}")

    # Complement pairing
    print(f"\nComplement pairing (T ↔ T^op):")
    seen = set()
    for p in class_props:
        if p['cf'] in seen:
            continue
        bits = classes[p['cf']][0]
        comp_bits = complement_bits(n, bits)
        cf_comp = class_of[comp_bits]
        if cf_comp == p['cf']:
            print(f"  Class {p['display_idx']} (H={p['H']}) → SELF (SC)")
            seen.add(p['cf'])
        else:
            # Find the complement class
            for p2 in class_props:
                if p2['cf'] == cf_comp:
                    H_comp = p2['H']
                    didx2 = p2['display_idx']
                    break
            print(f"  Class {p['display_idx']} (H={p['H']}) ↔ Class {didx2} (H={H_comp})")
            seen.add(p['cf'])
            seen.add(cf_comp)

    # Pos analysis
    print(f"\nPos-uniformity analysis:")
    pos_uniform = sum(1 for p in class_props if len(set(p['pos'])) == 1)
    print(f"  Pos-uniform: {pos_uniform} / {num_classes}")
    for p in class_props:
        if len(set(p['pos'])) == 1:
            print(f"    H={p['H']}, pos = {p['pos'][0]} per vertex, H/n = {p['H']/n:.2f}")
        else:
            vals = p['pos']
            print(f"    H={p['H']}, pos = {vals}, spread = {max(vals)-min(vals)}")

    elapsed = time.time() - t0
    print(f"\n[n={n} completed in {elapsed:.1f}s]")

    all_results[n] = {
        'num_classes': num_classes,
        'class_props': class_props,
        'adj': dict(adj),
        'score_groups': dict(score_groups),
    }

print("\n" + "=" * 70)
print("CROSS-n STRUCTURE COMPARISON")
print("=" * 70)

# Table of iso class counts
print("\n  n | #classes | #SC | #blueself | #blackself | max(H) | H=1 classes")
for n in sorted(all_results.keys()):
    r = all_results[n]
    cp = r['class_props']
    n_sc = sum(1 for p in cp if p['sc'])
    n_bs = sum(1 for p in cp if p['blueself'])
    n_bks = sum(1 for p in cp if p['blackself'])
    max_H = max(p['H'] for p in cp)
    H1 = sum(1 for p in cp if p['H'] == 1)
    print(f"  {n} | {r['num_classes']:>8} | {n_sc:>3} | {n_bs:>9} | {n_bks:>10} | {max_H:>6} | {H1}")

# Score class counts
print("\n  n | #score_classes | max_per_score")
for n in sorted(all_results.keys()):
    r = all_results[n]
    sg = r['score_groups']
    max_per = max(len(v) for v in sg.values())
    print(f"  {n} | {len(sg):>14} | {max_per}")

# Look for embedding patterns
print("\n" + "=" * 70)
print("SELF-SIMILAR / FRACTAL STRUCTURE SEARCH")
print("=" * 70)

print("""
Question: Do groups of iso classes at n+1 behave like single iso classes at n?

Strategy: Compare the flip graph structure. At n, each iso class has a set of
flip-graph neighbors. At n+1, does a GROUP of iso classes (sharing some property
like score pattern) have the same internal/external connectivity pattern?
""")

# Compare flip graph degree sequences
for n in sorted(all_results.keys()):
    r = all_results[n]
    cp = r['class_props']
    adj = r['adj']
    degrees = []
    for p in cp:
        ci = None
        for cf, idx in [(p['cf'], i) for i, p2 in enumerate(cp) if p2['cf'] == p['cf']]:
            ci = list(r['adj'].keys())[0] if r['adj'] else 0
            break
        # Get degree from adj
        for key, nbrs in r['adj'].items():
            found = False
            for p2 in cp:
                if p2['cf'] == p['cf']:
                    # Check if this key corresponds to this class
                    pass
            break
    # Simpler: just count degrees
    degs = []
    for p in cp:
        ci = 0
        for cf2, idx2 in zip([p2['cf'] for p2 in cp], range(len(cp))):
            if cf2 == p['cf']:
                ci = idx2
                break
        # Find in adj
        d = 0
        for key, val in r['adj'].items():
            pass
        degs.append(0)  # placeholder
    # The adj structure is complex, let me just report it
    print(f"  n={n}: {r['num_classes']} classes")

# H-normalized analysis
print("\nH-value normalization (H/max_H) by n:")
for n in sorted(all_results.keys()):
    cp = all_results[n]['class_props']
    max_H = max(p['H'] for p in cp)
    ratios = sorted(set(p['H']/max_H for p in cp))
    print(f"  n={n}: H/max = {[f'{r:.3f}' for r in ratios]}")

# c3 normalization
print("\nc3-value normalization (c3/max_c3) by n:")
for n in sorted(all_results.keys()):
    cp = all_results[n]['class_props']
    max_c3 = max(p['c3'] for p in cp) if max(p['c3'] for p in cp) > 0 else 1
    ratios = sorted(set(round(p['c3']/max_c3, 3) for p in cp))
    print(f"  n={n}: c3/max = {ratios}")

# Alpha1 patterns
print("\nα₁ values by n:")
for n in sorted(all_results.keys()):
    cp = all_results[n]['class_props']
    vals = sorted(set(p['alpha1'] for p in cp))
    print(f"  n={n}: α₁ ∈ {vals}")

# Alpha2 patterns
print("\nα₂ values by n:")
for n in sorted(all_results.keys()):
    cp = all_results[n]['class_props']
    vals = sorted(set(p['alpha2'] for p in cp))
    print(f"  n={n}: α₂ ∈ {vals}")

# Sub-tournament containment analysis
print("\n" + "=" * 70)
print("SUB-TOURNAMENT EMBEDDING ANALYSIS")
print("=" * 70)
print("""
For each iso class at n+1, compute which iso classes at n appear as
induced sub-tournaments (by vertex deletion). This gives a "containment"
relation between iso class sets at different n values.
""")

for n in [4, 5, 6]:
    if n not in all_results or n-1 not in all_results:
        continue
    print(f"\nn={n} classes → n={n-1} sub-tournament profiles:")
    cp_big = all_results[n]['class_props']
    cp_small = all_results[n-1]['class_props']

    # For each class at n, delete each vertex and find which n-1 class results
    for p in cp_big:
        bits = 0
        # Reconstruct from cf
        A = tournament_from_bits(n, 0)
        # Actually need to get the original tournament
        # Reconstruct from canonical form
        idx = 0
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if p['cf'][idx]:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        sub_profile = Counter()
        for del_v in range(n):
            # Delete vertex del_v
            remaining = [v for v in range(n) if v != del_v]
            B = [[0]*(n-1) for _ in range(n-1)]
            for ii, vi in enumerate(remaining):
                for jj, vj in enumerate(remaining):
                    if ii != jj:
                        B[ii][jj] = A[vi][vj]
            cf_B = canonical_form(B, n-1)
            # Find which class this is
            for p2 in cp_small:
                if p2['cf'] == cf_B:
                    sub_profile[p2['display_idx']] += 1
                    break

        profile_str = ', '.join(f"{k}×{v}" for k, v in sorted(sub_profile.items()))
        print(f"  Class {p['display_idx']} (H={p['H']}, {p['score']}): sub = [{profile_str}]")

print("\n" + "=" * 70)
print("INFORMATION-THEORETIC STRUCTURE")
print("=" * 70)

import math

for n in sorted(all_results.keys()):
    cp = all_results[n]['class_props']
    m = n*(n-1)//2
    total_t = 1 << m

    # Entropy of the iso class distribution
    probs = [p['size'] / total_t for p in cp]
    entropy = -sum(p * math.log2(p) for p in probs if p > 0)
    max_entropy = math.log2(len(cp))

    # Mutual information: H-value and score
    # I(H; score) = H(H) + H(score) - H(H, score)
    H_dist = Counter(p['H'] for p in cp)
    score_dist = Counter(str(p['score']) for p in cp)
    joint_dist = Counter((p['H'], str(p['score'])) for p in cp)

    H_entropy = -sum(c/len(cp) * math.log2(c/len(cp)) for c in H_dist.values())
    score_entropy = -sum(c/len(cp) * math.log2(c/len(cp)) for c in score_dist.values())
    joint_entropy = -sum(c/len(cp) * math.log2(c/len(cp)) for c in joint_dist.values())
    MI = H_entropy + score_entropy - joint_entropy

    print(f"\nn={n}:")
    print(f"  Iso class entropy: {entropy:.3f} bits (max = {max_entropy:.3f})")
    print(f"  H entropy: {H_entropy:.3f} bits, Score entropy: {score_entropy:.3f} bits")
    print(f"  I(H; score) = {MI:.3f} bits")
    if H_entropy > 0:
        print(f"  Score explains {MI/H_entropy*100:.1f}% of H variation (class-level)")
