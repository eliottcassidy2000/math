"""
Verify: are indistinguishable classes exactly sigma (converse) pairs?

Key hypothesis: All classes sharing (score, c3, c5, omega_deg) are either:
1. A sigma pair {C, C^op} — converse-paired, NOT self-converse
2. Multiple self-converse classes that happen to share invariants
3. A mixture

Also: sigma preserves weight, so sigma-paired classes have identical
weight distributions. But we saw at n=5 that some same-invariant
classes have DIFFERENT weight distributions — those must both be
self-converse, not a pair!
"""
from itertools import permutations, combinations
from collections import defaultdict, Counter

def build_full(n):
    arcs = [(i,j) for i in range(n) for j in range(i+2, n)]
    m = len(arcs)
    arc_idx = {a: k for k, a in enumerate(arcs)}
    sigma_perm = [arc_idx[(n-1-j, n-1-i)] for (i,j) in arcs]

    def apply_sigma(bits):
        nb = 0
        for k in range(m):
            if bits & (1 << k):
                nb |= (1 << sigma_perm[k])
        return nb

    results = []
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
        results.append({'bits': bits, 'canon': canon, 'sigma_bits': apply_sigma(bits)})

    return results, arcs, m

def analyze(n):
    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")

    results, arcs, m = build_full(n)

    classes = defaultdict(list)
    b2c = {}
    for r in results:
        classes[r['canon']].append(r)
        b2c[r['bits']] = r['canon']

    # sigma map on classes
    clist = list(classes.keys())
    csigma = {}
    for c in clist:
        csigma[c] = b2c[classes[c][0]['sigma_bits']]

    # Compute invariants for each class
    def class_invariants(canon):
        adj = [list(row) for row in canon]
        nv = len(adj)
        scores = tuple(sorted(sum(adj[v]) for v in range(nv)))

        c3 = 0
        triples = []
        for a in range(nv):
            for b in range(a+1, nv):
                for cc in range(b+1, nv):
                    if (adj[a][b] and adj[b][cc] and adj[cc][a]) or \
                       (adj[a][cc] and adj[cc][b] and adj[b][a]):
                        c3 += 1
                        triples.append(frozenset([a,b,cc]))

        c5 = 0
        if nv >= 5:
            for vs in combinations(range(nv), 5):
                for p in permutations(vs):
                    if all(adj[p[i]][p[(i+1)%5]] for i in range(5)):
                        c5 += 1
            c5 //= 5

        h = sum(1 for p in permutations(range(nv))
                if all(adj[p[i]][p[i+1]] for i in range(nv-1)))

        # Omega_3 degree seq
        omega_adj = [[False]*len(triples) for _ in range(len(triples))]
        for i in range(len(triples)):
            for j in range(i+1, len(triples)):
                if triples[i] & triples[j]:
                    omega_adj[i][j] = omega_adj[j][i] = True
        omega_deg = tuple(sorted(sum(omega_adj[i]) for i in range(len(triples))))

        return (scores, c3, c5, omega_deg, h)

    cinv = {c: class_invariants(c) for c in clist}

    # Group by invariants
    inv_groups = defaultdict(list)
    for c in clist:
        inv_groups[cinv[c]].append(c)

    multi = {k: v for k, v in inv_groups.items() if len(v) > 1}
    print(f"\nInvariant groups with >1 class: {len(multi)}")

    for inv, cs in sorted(multi.items(), key=lambda x: -len(x[1])):
        print(f"\n  Invariants: score={inv[0]}, c3={inv[1]}, c5={inv[2]}, "
              f"omega={inv[3]}, H={inv[4]}")
        print(f"  {len(cs)} classes:")
        for c in cs:
            is_sc = csigma[c] == c
            partner = csigma[c]
            members = classes[c]
            wdist = Counter(bin(r['bits']).count('1') for r in members)
            sf = sum(1 for r in members if r['sigma_bits'] == r['bits'])
            # Is partner in same group?
            partner_in_group = partner in cs
            print(f"    self-conv={is_sc}, partner_in_group={partner_in_group}, "
                  f"size={len(members)}, sigma-fixed={sf}")
            print(f"    weights: {dict(sorted(wdist.items()))}")

    # The KEY question: for each multi-group, classify the structure
    print(f"\n--- Classification of multi-invariant groups ---")
    all_sigma_pairs = 0
    all_self_conv = 0
    mixed = 0
    for inv, cs in multi.items():
        sc_count = sum(1 for c in cs if csigma[c] == c)
        pair_count = 0
        seen = set()
        for c in cs:
            if c in seen or csigma[c] == c:
                continue
            if csigma[c] in cs:
                pair_count += 1
                seen.add(csigma[c])
        if sc_count == 0 and pair_count * 2 == len(cs):
            all_sigma_pairs += 1
            t = "ALL SIGMA PAIRS"
        elif sc_count == len(cs):
            all_self_conv += 1
            t = "ALL SELF-CONVERSE"
        else:
            mixed += 1
            t = "MIXED"
        print(f"  {t}: {inv[0]}, c3={inv[1]}, c5={inv[2]}, "
              f"H={inv[4]}: {sc_count} SC + {pair_count} pairs = {len(cs)} classes")

    print(f"\nSummary: {all_sigma_pairs} all-pair groups, "
          f"{all_self_conv} all-SC groups, {mixed} mixed groups")

if __name__ == '__main__':
    for n in range(3, 7):
        analyze(n)
