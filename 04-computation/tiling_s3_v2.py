"""
S3 symmetry on the triangular tiling grid — corrected sigma.

sigma = converse + relabeling = reverse tournament + reverse vertex order.

For arc (i,j) with bit b:
- b=0 means i->j; T^op has j->i; relabel k->n-1-k: (n-1-j)->(n-1-i), forward. Bit=0.
- b=1 means j->i; T^op has i->j; relabel: (n-1-i)->(n-1-j), backward. Bit=1.

So sigma PRESERVES bits and PERMUTES arc positions:
  arc (i,j) -> arc (n-1-j, n-1-i), same bit.

This is the key: sigma is a pure permutation on {0,1}^m, no complement!

The complement (all bits flip) corresponds to: reverse all non-path arcs,
which is NOT the converse tournament.

What IS the converse in tiling language?
T^op relabeled by vertex reversal just permutes the tiling bits.
"""
from itertools import permutations
from collections import defaultdict, Counter

def build_full_data(n):
    arcs = []
    for i in range(n):
        for j in range(i+2, n):
            arcs.append((i, j))
    m = len(arcs)
    arc_idx = {a: k for k, a in enumerate(arcs)}

    # sigma permutation: arc (i,j) -> arc (n-1-j, n-1-i)
    sigma_perm = []
    for k, (i, j) in enumerate(arcs):
        new_arc = (n-1-j, n-1-i)
        sigma_perm.append(arc_idx[new_arc])

    def apply_sigma(bits):
        new_bits = 0
        for k in range(m):
            if bits & (1 << k):
                new_bits |= (1 << sigma_perm[k])
        return new_bits

    # Verify involution
    for bits in range(min(32, 1<<m)):
        assert apply_sigma(apply_sigma(bits)) == bits, f"Not involution at {bits}"

    # Also define complement (all bits flip) — this is a DIFFERENT operation
    def apply_complement(bits):
        return ((1 << m) - 1) ^ bits

    # And the "full reversal" = complement + sigma (reverse all arcs AND relabel)
    def apply_full_rev(bits):
        return apply_sigma(apply_complement(bits))

    # Build all tournaments
    results = []
    for bits in range(1 << m):
        adj = [[False]*n for _ in range(n)]
        for i in range(n-1):
            adj[i][i+1] = True
        for k, (i, j) in enumerate(arcs):
            if bits & (1 << k):
                adj[j][i] = True
            else:
                adj[i][j] = True

        canon = None
        for perm in permutations(range(n)):
            mat = tuple(tuple(adj[perm[a]][perm[b]] for b in range(n)) for a in range(n))
            if canon is None or mat < canon:
                canon = mat

        results.append({
            'bits': bits,
            'canon': canon,
            'sigma_bits': apply_sigma(bits),
            'comp_bits': apply_complement(bits),
            'fullrev_bits': apply_full_rev(bits),
        })

    return results, arcs, sigma_perm, apply_sigma

def analyze(n):
    print(f"\n{'='*70}")
    print(f"n = {n}, m = {(n-1)*(n-2)//2}")
    print(f"{'='*70}")

    results, arcs, sigma_perm, apply_sigma = build_full_data(n)
    m = len(arcs)

    print(f"Sigma permutation: {sigma_perm}")
    # Fixed points of sigma
    sigma_fixed_arcs = [k for k in range(m) if sigma_perm[k] == k]
    print(f"Sigma-fixed arcs: {[arcs[k] for k in sigma_fixed_arcs]}")

    classes = defaultdict(list)
    bits_to_canon = {}
    for r in results:
        classes[r['canon']].append(r)
        bits_to_canon[r['bits']] = r['canon']

    canon_list = list(classes.keys())

    # Check: does sigma map classes to classes?
    class_sigma = {}
    for canon in canon_list:
        target_canons = set()
        for r in classes[canon]:
            target_canons.add(bits_to_canon[r['sigma_bits']])
        if len(target_canons) != 1:
            print(f"  BUG: sigma maps class (size={len(classes[canon])}) to {len(target_canons)} classes!")
            # Debug: show the targets
            tc = Counter(bits_to_canon[r['sigma_bits']] for r in classes[canon])
            for c, cnt in tc.most_common(3):
                print(f"    -> class (size={len(classes[c])}): {cnt} tilings")
        else:
            class_sigma[canon] = target_canons.pop()

    if len(class_sigma) == len(canon_list):
        self_conv = sum(1 for c in canon_list if class_sigma[c] == c)
        paired = (len(canon_list) - self_conv) // 2
        print(f"\nSigma maps classes to classes cleanly!")
        print(f"Self-converse: {self_conv}, Paired: {paired}, Total: {len(canon_list)}")

        # Show self-converse classes
        print(f"\nSelf-converse classes:")
        for canon in sorted(canon_list, key=lambda c: len(classes[c]), reverse=True):
            if class_sigma[canon] != canon:
                continue
            members = classes[canon]
            # Count sigma-fixed tilings
            fixed = sum(1 for r in members if r['sigma_bits'] == r['bits'])
            # H(T)
            adj = [list(row) for row in canon]
            nv = len(adj)
            h = sum(1 for perm in permutations(range(nv))
                    if all(adj[perm[i]][perm[i+1]] for i in range(nv-1)))
            aut = h // len(members) if len(members) > 0 else 0
            print(f"  size={len(members):4d}, H={h:5d}, |Aut|={aut:3d}, "
                  f"sigma-fixed={fixed}")

        # Show paired classes
        print(f"\nPaired classes:")
        seen = set()
        for canon in sorted(canon_list, key=lambda c: len(classes[c]), reverse=True):
            if canon in seen or class_sigma[canon] == canon:
                continue
            partner = class_sigma[canon]
            seen.add(partner)
            m1 = classes[canon]
            m2 = classes[partner]
            h1 = sum(1 for perm in permutations(range(n))
                     if all(list(canon[perm[i]])[perm[i+1]] for i in range(n-1)))
            h2 = sum(1 for perm in permutations(range(n))
                     if all(list(partner[perm[i]])[perm[i+1]] for i in range(n-1)))
            print(f"  size={len(m1):4d} (H={h1}) <-> size={len(m2):4d} (H={h2})")

    # Check complement and full reversal too
    print(f"\nComplement (flip all non-path arcs):")
    class_comp = {}
    for canon in canon_list:
        target = set(bits_to_canon[r['comp_bits']] for r in classes[canon])
        if len(target) == 1:
            class_comp[canon] = target.pop()
    if len(class_comp) == len(canon_list):
        self_comp = sum(1 for c in canon_list if class_comp[c] == c)
        print(f"  Complement maps classes to classes: self-complement={self_comp}")
    else:
        print(f"  Complement does NOT map classes to classes cleanly!")

    print(f"\nFull reversal (complement + sigma):")
    class_fr = {}
    for canon in canon_list:
        target = set(bits_to_canon[r['fullrev_bits']] for r in classes[canon])
        if len(target) == 1:
            class_fr[canon] = target.pop()
    if len(class_fr) == len(canon_list):
        self_fr = sum(1 for c in canon_list if class_fr[c] == c)
        print(f"  Full reversal maps classes to classes: self-fixed={self_fr}")
    else:
        print(f"  Full reversal does NOT map classes to classes cleanly!")

    # Sigma-fixed tilings: these encode self-converse tournaments
    # The number should be 2^floor((n-1)^2/4) per the tex
    sigma_fixed = sum(1 for r in results if r['sigma_bits'] == r['bits'])
    expected = 2**((n-1)**2 // 4)
    print(f"\nSigma-fixed tilings: {sigma_fixed} (expected 2^floor((n-1)^2/4) = {expected})")

    # Weight of sigma(t): sigma preserves weight (just permutes bits)
    weight_preserved = all(
        bin(r['bits']).count('1') == bin(r['sigma_bits']).count('1')
        for r in results
    )
    print(f"Sigma preserves weight: {weight_preserved}")

    # Complement changes weight: w -> m-w
    comp_weight = all(
        bin(r['bits']).count('1') + bin(r['comp_bits']).count('1') == m
        for r in results
    )
    print(f"Complement: weight(t) + weight(comp(t)) = m: {comp_weight}")

    return classes, class_sigma

def deep_pattern_analysis(n):
    """For n=5,6: what distinguishes classes that share score and c3?"""
    print(f"\n{'='*70}")
    print(f"DEEP PATTERN ANALYSIS: n = {n}")
    print(f"{'='*70}")

    results, arcs, sigma_perm, apply_sigma = build_full_data(n)
    m = len(arcs)

    classes = defaultdict(list)
    bits_to_canon = {}
    for r in results:
        classes[r['canon']].append(r)
        bits_to_canon[r['bits']] = r['canon']

    # For each class, compute:
    # 1. Score sequence
    # 2. c3 (3-cycle count)
    # 3. c5 (5-cycle count if applicable)
    # 4. Degree sequence of the "3-cycle graph" Omega_3
    class_info = {}
    for canon, members in classes.items():
        adj = [list(row) for row in canon]
        nv = len(adj)
        scores = tuple(sorted(sum(adj[v]) for v in range(nv)))

        # 3-cycles
        c3 = 0
        for a in range(nv):
            for b in range(a+1, nv):
                for c in range(b+1, nv):
                    if (adj[a][b] and adj[b][c] and adj[c][a]) or \
                       (adj[a][c] and adj[c][b] and adj[b][a]):
                        c3 += 1

        # 5-cycles (directed)
        c5 = 0
        if nv >= 5:
            from itertools import combinations
            for verts in combinations(range(nv), 5):
                # Check all 12 directed 5-cycles on these 5 vertices
                for perm in permutations(verts):
                    if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
                        c5 += 1
            c5 //= 5  # each cycle counted 5 times

        # H(T)
        h = sum(1 for perm in permutations(range(nv))
                if all(adj[perm[i]][perm[i+1]] for i in range(nv-1)))

        # Omega_3 degree sequence
        cycles_3 = []
        for a in range(nv):
            for b in range(a+1, nv):
                for c in range(b+1, nv):
                    if (adj[a][b] and adj[b][c] and adj[c][a]) or \
                       (adj[a][c] and adj[c][b] and adj[b][a]):
                        cycles_3.append(frozenset([a,b,c]))

        # Build Omega_3: vertices are 3-cycles, edges when they share a vertex
        omega_adj = [[False]*len(cycles_3) for _ in range(len(cycles_3))]
        for i in range(len(cycles_3)):
            for j in range(i+1, len(cycles_3)):
                if cycles_3[i] & cycles_3[j]:
                    omega_adj[i][j] = omega_adj[j][i] = True
        omega_deg = tuple(sorted(sum(omega_adj[i]) for i in range(len(cycles_3))))

        class_info[canon] = {
            'scores': scores,
            'c3': c3,
            'c5': c5,
            'h': h,
            'size': len(members),
            'omega_deg': omega_deg,
        }

    # Group by (scores, c3)
    groups = defaultdict(list)
    for canon, info in class_info.items():
        key = (info['scores'], info['c3'])
        groups[key].append((canon, info))

    print(f"\nGroups with same (score, c3) but multiple classes:")
    for key, class_list in sorted(groups.items()):
        if len(class_list) < 2:
            continue
        scores, c3 = key
        print(f"\n  Score={scores}, c3={c3}: {len(class_list)} classes")
        for canon, info in class_list:
            print(f"    H={info['h']:5d}, size={info['size']:4d}, c5={info['c5']}, "
                  f"omega_deg={info['omega_deg']}")

    # Can (scores, c3, c5) distinguish?
    groups2 = defaultdict(list)
    for canon, info in class_info.items():
        key = (info['scores'], info['c3'], info['c5'])
        groups2[key].append((canon, info))
    multi2 = sum(1 for g in groups2.values() if len(g) > 1)
    print(f"\nGroups with same (score, c3, c5) and multiple classes: {multi2}")
    for key, class_list in sorted(groups2.items()):
        if len(class_list) < 2:
            continue
        print(f"  {key}: {len(class_list)} classes")
        for canon, info in class_list:
            print(f"    H={info['h']:5d}, omega_deg={info['omega_deg']}")

    # Can omega_deg distinguish remaining ones?
    groups3 = defaultdict(list)
    for canon, info in class_info.items():
        key = (info['scores'], info['c3'], info['c5'], info['omega_deg'])
        groups3[key].append((canon, info))
    multi3 = sum(1 for g in groups3.values() if len(g) > 1)
    print(f"\nGroups with same (score, c3, c5, omega_deg) and multiple classes: {multi3}")
    for key, class_list in sorted(groups3.items()):
        if len(class_list) < 2:
            continue
        print(f"  {key}: {len(class_list)} classes")
        for canon, info in class_list:
            print(f"    H={info['h']:5d}, size={info['size']:4d}")

if __name__ == '__main__':
    for n in range(3, 7):
        analyze(n)

    for n in range(5, 7):
        deep_pattern_analysis(n)
