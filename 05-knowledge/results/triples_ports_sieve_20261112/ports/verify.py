"""Exact controls for short quiver returns and the inherited paired-port code.

No randomness, no assert-dependent checks, no writes outside this directory.
Run from any directory: python verify.py [--output receipt.json]
"""
from collections import Counter, deque
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
import argparse
import gzip
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = ROOT / '05-knowledge/results/session290_core_descent_20261113/descent/profiles.json.gz'


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def validate(B):
    n = len(B)
    need(n >= 1 and all(len(r) == n for r in B), 'nonempty square matrix')
    need(all(type(x) is int for r in B for x in r), 'integer matrix')
    need(all(B[i][j] == -B[j][i] for i in range(n) for j in range(n)), 'skew matrix')


def mutate(B, k):
    n = len(B)
    return tuple(tuple(-B[i][j] if i == k or j == k else B[i][j] +
                       (abs(B[i][k])*B[k][j] + B[i][k]*abs(B[k][j])) // 2
                       for j in range(n)) for i in range(n))


def mutate_paths(B, k):
    """Independent directed-path addition, cancellation, and arrow reversal."""
    C = [list(r) for r in B]
    n = len(B)
    for i in range(n):
        for j in range(n):
            if i != j and B[i][k] > 0 and B[k][j] > 0:
                z = B[i][k]*B[k][j]
                C[i][j] += z
                C[j][i] -= z
    for i in range(n):
        C[i][k] = -B[i][k]
        C[k][i] = -B[k][i]
    return tuple(map(tuple, C))


def tournament(B):
    return all(abs(B[i][j]) == 1 for i in range(len(B)) for j in range(i))


def matrix(n, mask):
    B = [[0]*n for _ in range(n)]
    for q, (i, j) in enumerate(combinations(range(n), 2)):
        B[i][j] = 1 if mask >> q & 1 else -1
        B[j][i] = -B[i][j]
    return tuple(map(tuple, B))


def maskof(B):
    return sum((B[i][j] > 0) << q for q, (i, j) in
               enumerate(combinations(range(len(B)), 2)))


def hamilton(B):
    n = len(B)
    d = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        d[1 << i][i] = 1
    for m in range(1, 1 << n):
        for i in range(n):
            if not (m >> i & 1):
                continue
            for j in range(n):
                if not (m >> j & 1) and B[i][j] > 0:
                    d[m | 1 << j][j] += d[m][i]
    return sum(d[-1])


def extreme(B, k):
    return len({B[k][i] for i in range(len(B)) if i != k}) <= 1


def extreme_pair(B, k, l):
    """R -> k -> l or l -> k -> R; k is the inner member."""
    a = B[k][l]
    return all(B[k][v] == -a and B[l][v] == -a
               for v in range(len(B)) if v not in (k, l))


def sccs(B):
    n = len(B)
    reach = [[i == j or B[i][j] > 0 for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] |= reach[i][k] and reach[k][j]
    unused = set(range(n))
    parts = []
    while unused:
        i = min(unused)
        c = tuple(j for j in unused if reach[i][j] and reach[j][i])
        unused.difference_update(c)
        parts.append(c)
    unsorted_parts = tuple(parts)
    parts.sort(key=lambda c: -sum(B[c[0]][d[0]] > 0 for d in unsorted_parts if c != d))
    return parts


def mutation_checks():
    previous_a = 1
    rows = []
    for n in range(1, 6):
        universe = {m: matrix(n, m) for m in range(1 << (n*(n-1)//2))}
        hs = {m: hamilton(B) for m, B in universe.items()}
        iso = {}
        for m, B in universe.items():
            if m in iso:
                continue
            orbit = {maskof(tuple(tuple(B[i][j] for j in p) for i in p))
                     for p in permutations(range(n))}
            rep = min(orbit)
            iso.update({q: rep for q in orbit})
        a = len(set(iso.values()))
        adjacency = {}
        exits = three_returns = mutation_count = 0
        for mask, B in universe.items():
            brute_h = sum(all(B[p[i]][p[i+1]] > 0 for i in range(n-1))
                          for p in permutations(range(n)))
            need(brute_h == hs[mask], 'subset DP versus permutation count')
            adjacency[mask] = set()
            for k in range(n):
                C = mutate(B, k)
                need(C == mutate_paths(B, k), 'matrix/path mutation cross-check')
                need(mutate(C, k) == B, 'involution')
                need(tournament(C) == extreme(B, k), 'one-step iff source/sink')
                mutation_count += 1
                if tournament(C):
                    q = maskof(C)
                    adjacency[mask].add(q)
                    need(hs[q] == hs[mask], 'native H invariance')
                    continue
                exits += 1
                for l in range(n):
                    if l == k:
                        continue
                    D = mutate(C, l)
                    need(D == mutate_paths(C, l), 'weighted path mutation cross-check')
                    need(not tournament(D), 'no non-backtracking two-step return')
                    for m in range(n):
                        if m == l:
                            continue
                        E = mutate(D, m)
                        expected = m == k and extreme_pair(B, k, l)
                        need(tournament(E) == expected, 'complete three-step classification')
                        if expected:
                            three_returns += 1
                            need(hamilton(E) == hs[mask], 'three-step H invariance')
        unseen = set(universe)
        components = []
        while unseen:
            seed = min(unseen)
            todo = [seed]
            comp = set()
            while todo:
                m = todo.pop()
                if m in comp:
                    continue
                comp.add(m)
                todo.extend(adjacency[m] - comp)
            unseen.difference_update(comp)
            parts = sccs(universe[seed])
            nontrivial = [i for i, c in enumerate(parts) if len(c) > 1]
            if nontrivial:
                t = nontrivial[0] + len(parts) - 1 - nontrivial[-1]
                need(len(comp) == t+1, 'fixed-hull path size')
                degrees = sorted(len(adjacency[m] - {m}) for m in comp)
                need(degrees == ([0] if t == 0 else [1, 1] + [2]*(t-1)), 'path degrees')
            else:
                need(len(comp) == n, 'transitive labelled cyclic rotations')
            need(len({hs[m] for m in comp}) == 1, 'component H')
            components.append(comp)
        ordinary_components = {tuple(sorted({iso[m] for m in c})) for c in components}
        ordinary_count = len(ordinary_components)
        need(ordinary_count == a-previous_a+1, 'A000568 first-difference formula')
        if n >= 3:
            need(three_returns == {3: 12, 4: 48, 5: 320}[n], 'three-step census')
        rows.append(dict(n=n, labelled=len(universe), isotypes=a,
                         labelled_native_components=len(components),
                         ordinary_native_components=ordinary_count,
                         one_step_checks=mutation_count, nonnative_first_steps=exits,
                         three_step_excursions=three_returns))
        previous_a = a
    # Minimal hostile: a cyclic triangle loses an edge at every mutation.
    cyclic = ((0, 1, -1), (-1, 0, 1), (1, -1, 0))
    need(hamilton(cyclic) == 3 and all(not tournament(mutate(cyclic, k)) for k in range(3)),
         'cyclic hostile')
    rejects = 0
    for bad in [(), ((0, 1),), ((0, 1), (1, 0)), ((0, 1.0), (-1.0, 0))]:
        try:
            validate(bad)
        except RuntimeError:
            rejects += 1
    need(rejects == 4, 'malformed-input controls remain active under -O')
    return dict(census=rows, malformed_rejections=rejects)


def port_checks():
    raw = SOURCE.read_bytes()
    r = json.loads(gzip.decompress(raw))['even_records']
    need(len(r) == 24, '24 inherited even charts')
    bases = [1, 2, 3, 4, 6, 8, 9, 12]
    D = [[dict(r[2*j+1]['U']).get(b, 0)-dict(r[2*j]['U']).get(b, 0)
          for j in range(12)] for b in bases]
    need(len({(tuple(q['rows']), q['selected']) for q in r}) == 12, '12 marked bases')
    for j in range(12):
        a, b = r[2*j:2*j+2]
        need(a['rows'] == b['rows'] and a['selected'] == b['selected'], 'paired terminal graph')
        need(a['ports'] != b['ports'], 'distinct port marking')
    # Physical lengths 2 through 12 are all nonconstant phases modulo 12.
    times = tuple(range(2, 13))
    columns = [tuple(sum(D[i][j]*pow(b, m, 13) for i, b in enumerate(bases)) % 13
                     for m in times) for j in range(12)]
    need(all(sum(D[i][j]*b for i, b in enumerate(bases)) % 13 == 0
             for j in range(12)), 'omitted length-one phase is constant')
    words = [tuple(sum(columns[j][i] for j in range(12) if e >> j & 1) % 13
                   for i in range(11)) for e in range(4096)]
    need(len(set(words)) == 4096, 'full observation code injectivity')
    # Independent evaluation from entire U-polynomials, after subtracting the zero-choice baseline.
    for e in (0, 1, 37, 2048, 4095):
        direct = tuple(sum(sum(c*pow(b, m, 13) for b, c in r[2*j+((e >> j) & 1)]['U'])-
                           sum(c*pow(b, m, 13) for b, c in r[2*j]['U'])
                           for j in range(12)) % 13 for m in times)
        need(direct == words[e], 'direct polynomial versus difference-column observer')
    best_four = 0
    best_four_times = None
    for slots in combinations(range(11), 4):
        count = len({tuple(w[i] for i in slots) for w in words})
        if count > best_four:
            best_four, best_four_times = count, tuple(times[i] for i in slots)
    need(best_four == 4064, 'four observations never distinguish all choices')
    good_five = []
    for slots in combinations(range(11), 5):
        if len({tuple(w[i] for i in slots) for w in words}) == 4096:
            good_five.append(tuple(times[i] for i in slots))
    need(len(good_five) == 225 and (2, 3, 4, 6, 12) in good_five, 'optimal five observations')
    # Every binary difference is ternary, with exactly 2^(number of zero digits)
    # ordered preimage pairs. Enumerate ternary differences incrementally.
    states = [(tuple([0]*11), 0, 1)]
    power = 1
    for col in columns:
        nxt = []
        for v, label, multiplicity in states:
            nxt.append((tuple((a-b) % 13 for a, b in zip(v, col)), label, multiplicity))
            nxt.append((v, label+power, 2*multiplicity))
            nxt.append((tuple((a+b) % 13 for a, b in zip(v, col)), label+2*power, multiplicity))
        states = nxt
        power *= 3
    zero_label = (3**12-1)//2
    ternary_hist = Counter()
    pair_hist = Counter()
    witness = None
    for v, label, multiplicity in states:
        if label == zero_label:
            continue
        weight = sum(a != 0 for a in v)
        ternary_hist[weight] += 1
        pair_hist[weight] += multiplicity
        if weight == 5 and witness is None:
            q = label
            digits = []
            for _ in range(12):
                digits.append(q % 3-1)
                q //= 3
            positive = sum(1 << i for i, a in enumerate(digits) if a == 1)
            negative = sum(1 << i for i, a in enumerate(digits) if a == -1)
            need(tuple((a-b) % 13 for a, b in zip(words[positive], words[negative])) == v,
                 'distance witness has actual binary preimages')
            witness = dict(positive_choice=positive, negative_choice=negative,
                           ternary_difference=digits, observation_difference=list(v))
    need(min(ternary_hist) == 5 and sum(ternary_hist.values()) == 3**12-1, 'exact distance five')
    need(sum(pair_hist.values()) == 4096*4095, 'ordered-pair histogram mass')
    need(ternary_hist == Counter({5: 28, 6: 468, 7: 3416, 8: 21446,
                                9: 84972, 10: 202604, 11: 218506}), 'distance census')
    # Independent all-pairs route: pack eleven 13-ary digits into four-bit slots.
    # ORing the four shifted XORs marks precisely the differing slots.
    packed = [sum(x << (4*i) for i, x in enumerate(w)) for w in words]
    slot_mask = sum(1 << (4*i) for i in range(11))
    direct_hist = Counter()
    for i, a in enumerate(packed):
        for b in packed[:i]:
            x = a ^ b
            weight = ((x | x >> 1 | x >> 2 | x >> 3) & slot_mask).bit_count()
            direct_hist[weight] += 2
    need(direct_hist == pair_hist, 'independent complete all-pairs distance census')
    # The canonical paired binary code independently has weight-two words.
    pair_weights = Counter(2*e.bit_count() for e in range(4096))
    need(min(w for w in pair_weights if w) == 2, 'paired binary code is not Golay')
    return dict(source=str(SOURCE.relative_to(ROOT)).replace('\\', '/'),
                source_sha256=sha256(raw).hexdigest(), difference_matrix=D,
                times=list(times), modulus=13, words=len(words), minimum_distance=5,
                best_four_images=best_four, best_four_times=best_four_times,
                working_five_sets=good_five, working_five_count=len(good_five),
                ternary_distance_histogram=dict(sorted(ternary_hist.items())),
                ordered_pair_distance_histogram=dict(sorted(pair_hist.items())),
                independent_all_pairs_checked=4096*4095//2,
                distance_witness=witness,
                canonical_paired_binary_weight_histogram=dict(sorted(pair_weights.items())))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='receipt.json')
    args = parser.parse_args()
    output = (HERE / args.output).resolve()
    need(output.parent == HERE, 'receipt stays inside owned directory')
    result = dict(status='FINITE-EXACT; theorem proofs in PROOF.md',
                  mutation=mutation_checks(), ports=port_checks())
    output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    print(json.dumps(dict(status='PASS', mutation_census=result['mutation']['census'],
                          port_code=[11, 4096, 5, 13], minimum_observations=5,
                          working_five_sets=225), sort_keys=True))


if __name__ == '__main__':
    main()
