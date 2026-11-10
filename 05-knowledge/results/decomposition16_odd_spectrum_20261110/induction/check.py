#!/usr/bin/env python3
"""Joint deletion cores: scoped exact small census and rational-series controls."""
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import permutations, product
from math import gcd
from pathlib import Path
import json
import gzip

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = ROOT / '05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/n7_representatives.txt'
WINDOW = ROOT / '05-knowledge/results/tilted_ensemble_local_limit_20261110/out/pA_le10.bin.gz'

def require(q, msg):
    if not q:
        raise RuntimeError(msg)

def rows_from_code(n, code):
    rows = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code >> k & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
            k += 1
    return tuple(rows)

def induced(rows, vertices):
    return tuple(sum(((rows[v] >> w) & 1) << j for j, w in enumerate(vertices)) for v in vertices)

def hp_words(rows):
    n = len(rows)
    if not n:
        return [()]
    ans = []
    def walk(word, unused):
        if not unused:
            ans.append(word)
            return
        m = unused if not word else unused & rows[word[-1]]
        while m:
            b = m & -m
            m -= b
            walk(word + (b.bit_length() - 1,), unused ^ b)
    walk((), (1 << n) - 1)
    return ans

@lru_cache(None)
def hdp(rows):
    n = len(rows)
    if not n:
        return 1
    d = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        d[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for i, x in enumerate(d[mask]):
            if x:
                out = rows[i] & ~mask
                while out:
                    b = out & -out
                    out -= b
                    d[mask | b][b.bit_length() - 1] += x
    return sum(d[-1])

def code_at(rows, p):
    return sum(((rows[p[i]] >> p[j]) & 1) << k for k, (i, j) in enumerate((i, j) for i in range(len(p)) for j in range(i + 1, len(p))))

@lru_cache(None)
def canon(rows):
    return min(code_at(rows, p) for p in hp_words(rows))

def sccs(rows):
    n = len(rows)
    reach = [r | (1 << i) for i, r in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    left = set(range(n)); ans = []
    while left:
        v = min(left)
        part = tuple(w for w in sorted(left) if reach[v] >> w & 1 and reach[w] >> v & 1)
        ans.append(part)
        left.difference_update(part)
    return ans

def root_orbits(rows):
    n = len(rows)
    groups = defaultdict(list)
    for i, row in enumerate(rows):
        groups[row.bit_count()].append(i)
    groups = list(groups.values())
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for image_blocks in product(*(list(permutations(g)) for g in groups)):
        p = [0] * n
        for block, image in zip(groups, image_blocks):
            for v, w in zip(block, image):
                p[v] = w
        if all(((rows[i] >> j) & 1) == ((rows[p[i]] >> p[j]) & 1) for i in range(n) for j in range(i + 1, n)):
            for v, w in enumerate(p):
                parent[find(v)] = find(w)
    return sorted({find(i) for i in range(n)})

def convolution(a, b, degree=7):
    out = Counter()
    for i, x in a.items():
        for j, y in b.items():
            if i + j <= degree:
                out[i + j] += x * y
    return out

def main():
    counts = [1, 1, 1, 2, 4, 12, 56, 456]
    reps = {0: [()]}
    presentation_count = 0
    # Every tournament has a Hamilton path; force 0->1->... and quotient by all actual Hamilton words.
    for n in range(1, 7):
        free = [(i, j) for i in range(n) for j in range(i + 2, n)]
        classes = {}
        for mask in range(1 << len(free)):
            rows = [0] * n
            for i in range(n - 1):
                rows[i] |= 1 << (i + 1)
            for k, (i, j) in enumerate(free):
                a, b = (i, j) if mask >> k & 1 else (j, i)
                rows[a] |= 1 << b
            rows = tuple(rows)
            classes[canon(rows)] = rows
            presentation_count += 1
        reps[n] = list(classes.values())
        require(len(reps[n]) == counts[n], ('class count', n))
    reps[7] = [tuple(map(int, line.split())) for line in SOURCE.read_text().splitlines() if line.strip()]
    require(len(reps[7]) == counts[7], 'order-seven source length')
    keys7 = {canon(r) for r in reps[7]}
    require(len(keys7) == counts[7], 'order-seven pairwise inequivalence')
    # Independent exhaustive augmentation verifies the inherited representative file's coverage.
    augmented = set()
    for rows in reps[6]:
        for cut in range(64):
            child = tuple(rows[i] | (0 if cut >> i & 1 else 1 << 6) for i in range(6)) + (cut,)
            augmented.add(canon(child))
    require(augmented == keys7, 'order-seven complete augmentation')
    F = defaultdict(Counter); Q = defaultdict(Counter); G = defaultdict(Counter)
    F[1][0] = 1
    core_pairs = defaultdict(list)
    factor_checks = marked_checks = path_checks = 0
    for n in range(1, 8):
        for rows in reps[n]:
            N = hdp(rows)
            require(N == len(hp_words(rows)), 'independent path count')
            path_checks += 1
            F[N][n] += 1
            components = sccs(rows)
            component_count = {c: hdp(induced(rows, c)) for c in components}
            for v in range(n):
                h = hdp(induced(rows, tuple(i for i in range(n) if i != v)))
                c = next(c for c in components if v in c)
                hc = component_count[c]
                dc = hdp(induced(rows, tuple(i for i in c if i != v)))
                outside = 1
                for part in components:
                    if part != c:
                        outside *= component_count[part]
                require((N, h) == (outside * hc, outside * dc), 'SCC pair factorization')
                require(gcd(N, h) % outside == 0, 'outside divisor')
                if N > h and gcd(N, h) == 1:
                    require(len(c) >= 3 and all(len(part) == 1 for part in components if part != c), 'coprime core')
                factor_checks += 1
            for v in root_orbits(rows):
                h = hdp(induced(rows, tuple(i for i in range(n) if i != v)))
                G[h, N][n] += 1
                marked_checks += 1
                if len(components) == 1:
                    Q[h, N][n] += 1
                    if N < 45:
                        core_pairs[h, N].append({'order': n, 'rows': list(rows), 'root': v})
    # G_(h,N) = sum_ab r, r|gcd F_a F_b Q_(h/r,N/r), with C singleton allowed.
    predicted = defaultdict(Counter)
    for (h, N), poly in Q.items():
        for a, left in F.items():
            for b, right in F.items():
                term = convolution(poly, convolution(left, right))
                if term:
                    r = a * b
                    predicted[h * r, N * r].update(term)
    require(dict(predicted) == dict(G), 'full marked generating series through order seven')
    # Complete in all orders for coprime N<45, since strong order eight has minimum45.
    small_branches = []
    for h in [1, 3, 5, 9, 11, 13]:
        for offset in [2, 4]:
            N = 3 * h + offset
            small_branches.append({'h': h, 'N': N, 'exists': (h, N) in core_pairs, 'core_polynomial': dict(sorted(Q[h, N].items()))})
    require((13, 27) not in core_pairs, '13 to27 hostile')
    require((5, 11) in core_pairs and (5, 13) in core_pairs, 'positive two-branch controls')
    sample = {(h,N):dict(sorted(Q[h,N].items())) for h,N in [(1,3),(1,5),(3,5),(3,11),(3,13),(5,11),(5,13),(5,17),(5,19),(9,29),(9,31),(11,35),(11,37),(13,41),(13,43),(13,27)]}
    # Canonical Qs are marked isomorphism classes, not n raw root choices.
    # Read the separately audited inherited all-base-order<=10 census; do not claim to rerun it.
    raw = gzip.decompress(WINDOW.read_bytes())
    require(len(raw) == 512 * 512, 'inherited pair table dimensions')
    inherited = {h: int.from_bytes(raw[(h//2)*512:(h//2+1)*512], 'little') for h in range(1,1024,2)}
    require([h for h, row in inherited.items() if not row] == [7,21], 'inherited table attained rows')
    misses = {offset: [(h, 3*h+offset) for h, row in inherited.items() if row and not (row >> (h+offset//2)) & 1] for offset in [2,4]}
    require(misses == {2:[(3,11),(15,47),(35,107),(93,281)],4:[(1,7),(17,55),(35,109),(41,127),(87,265)]}, 'inherited ternary failure list')
    require(not (inherited[13] >> 7) & 1, 'inherited13to27 agrees with independent small census')
    def busch(n):
        q, r = divmod(n,3)
        return [3*5**(q-1),5**q,9*5**(q-1)][r]
    require(busch(12) == 375 and all(N < busch(12) and gcd(h,N)==1 for pairs in misses.values() for h,N in pairs), 'all-order cutoff for every inherited missing branch')
    result = {
        'status':'PASS',
        'inherited_not_rerun_window_census': {'base_order_max':10,'h_max':1023,'missing_ternary_pairs':misses,'both_branches_h_interval':[95,1023],'both_branch_rows':465,'positive_branch_pairs':1011,'global_absence_cutoff':{'f12':375,'largest_missing_N':281}},
        'universe':'all unlabelled tournament classes of orders0..7; marked classes use automorphism root orbits',
        'class_counts': counts,
        'hamiltonian_frame_presentations_through6': presentation_count,
        'order7_augmentations': len(reps[6])*64,
        'independent_path_count_checks':path_checks,
        'SCC_factorization_checks':factor_checks,
        'marked_classes':marked_checks,
        'joint_series_pairs_checked':len(G),
        'small_ternary_branches_complete_all_orders':small_branches,
        'selected_strong_marked_polynomials': {f'{h},{N}':p for (h,N),p in sample.items()},
        'hostile_13_to27': {'possible_strong_core_orders':[3,4,5,6,7], 'matches':0},
        'source_sha256': {str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [SOURCE,WINDOW]},
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()
    }
    print(json.dumps(result, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
