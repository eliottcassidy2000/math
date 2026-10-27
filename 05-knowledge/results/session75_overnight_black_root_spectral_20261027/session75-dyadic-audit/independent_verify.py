#!/usr/bin/env python3
"""Independent Session75 audit, Python >=3.9 standard library only.

No producer code is imported or executed. Source formulas are implemented
directly; producer JSON is comparison data, not construction input.
Run with --producer /path/to/producer to compare every claimed source record.
"""
import argparse
import hashlib
import itertools
import json
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PIN = '96e35642110104cc547516a63f37eebcf930f707'
FRONTIER_SHA = 'df88f03b576b13c98c0c5378fc1cf61685eeb6817073681266d2a781dc329b26'
PATTERNS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))


def require(condition, detail):
    if not condition:
        raise ValueError(detail)


def plus(*vectors):
    return tuple(sum(column) for column in zip(*vectors)) if vectors else (0, 0)


def norm(x):
    require(len(x) % 2 == 0, 'terminal direct-sum coordinate pairs')
    return sum(x[j]**2 - x[j]*x[j+1] + x[j+1]**2 for j in range(0, len(x), 2))


def cross(x, y):
    return norm(plus(x, y)) - norm(x) - norm(y)


def esum(vectors):
    return sum(map(norm, vectors))


def debt(path):
    return sum(max(0, left-right) for left, right in zip(path, path[1:]))


def excess(mid, left, right):
    return max(min(left, right)-mid, mid-max(left, right), 0)


def matchings(indices):
    if not indices:
        yield ()
    else:
        first, *tail = indices
        for second in tail:
            for rest in matchings(tuple(i for i in tail if i != second)):
                yield ((first, second),) + rest


def merge(vectors, matching):
    return [plus(vectors[i], vectors[j]) for i, j in matching]


def four_test(vectors):
    p, c = norm(plus(*vectors)), esum(vectors)
    mids = [esum(merge(vectors, pairing)) for pairing in PATTERNS]
    require(sum(mids) == p+2*c, 'three-pair exact identity')
    require(all(F(p, 2) <= m <= 2*c for m in mids), 'middle energy bounds')
    d = [debt([p, m, c]) for m in mids]
    ex = [excess(m, p, c) for m in mids]
    require(d == [max(p-c, 0)+e for e in ex], 'exact insertion formula')
    bound = min(F(p, 2), max(c-F(p, 2), 0))
    mean_bound = (excess(2*c, p, c)+2*excess(F(p, 2), p, c))/3
    require(min(ex) <= bound, 'deterministic excess bound')
    require(F(sum(ex), 3) <= mean_bound, 'mean excess bound')
    return {'P': p, 'C': c, 'middles': mids, 'debts': d}


def selected_tree(vectors):
    n = len(vectors)
    require(n > 0 and not n & (n-1), 'positive power of two')
    p = norm(plus(*vectors))
    bins = [(v, (j,)) for j, v in enumerate(vectors)]
    energies = [esum(vectors)]
    levels = []
    while len(bins) > 1:
        m = len(bins)
        candidates = list(matchings(tuple(range(m))))
        pairing = min(candidates, key=lambda t: esum(merge([b[0] for b in bins], t)))
        new = [(plus(bins[i][0], bins[j][0]), bins[i][1]+bins[j][1]) for i, j in pairing]
        e, ep = energies[-1], esum([b[0] for b in new])
        require(ep <= F(p+(m-2)*e, m-1), 'minimum versus perfect matching mean')
        require(max(ep-e, 0) <= F(max(p-e, 0), m-1) <= F(p, m), 'coarsening gain direction')
        flat = [label for _, labels in new for label in labels]
        require(sorted(flat) == list(range(n)), 'all labels exactly once after every merge')
        levels.append([list(labels) for _, labels in new])
        bins = new
        energies.append(ep)
    path = list(reversed(energies))
    require(debt(path) <= F((n-1)*p, n), 'selected compatible tree theorem')
    return {'P': p, 'energies_root_to_leaf': path, 'D': debt(path), 'blocks_bottom_up': levels}


def algebra():
    grid = tuple(itertools.product(range(-1, 2), repeat=2))
    count = 0
    for xs in itertools.product(grid, repeat=4):
        four_test(xs)
        count += 1
    hostile = four_test(((2, -1), (2, -1), (0, 1), (0, 1)))
    require(hostile == {'P': 16, 'C': 16, 'middles': [32, 8, 8], 'debts': [16, 8, 8]}, 'sharp P=C16 witness')
    # All feasible rational (P,C) regimes via x1=x2=(p+a)/4, x3=x4=(p-a)/4.
    sharp_samples = 0
    for scalar_p, scalar_a in itertools.product(range(7), repeat=2):
        vectors = ((F(scalar_p+scalar_a, 4), 0),)*2 + ((F(scalar_p-scalar_a, 4), 0),)*2
        result = four_test(vectors)
        p, c = result['P'], result['C']
        ds = [excess(m, p, c) for m in result['middles']]
        require(min(ds) == min(F(p, 2), max(c-F(p, 2), 0)), 'deterministic sharp family')
        require(F(sum(ds), 3) == (excess(2*c, p, c)+2*excess(F(p, 2), p, c))/3, 'expected sharp family')
        sharp_samples += 1
    high = []
    for n in (2, 4, 6, 8, 10):
        xs = [(j-3, (j*j+1) % 7-3, (-1)**j, j % 3-1) for j in range(n)]
        p, c = norm(plus(*xs)), esum(xs)
        energies = [esum(merge(xs, m)) for m in matchings(tuple(range(n)))]
        mean = F(sum(energies), len(energies))
        require(mean == F(p+(n-2)*c, n-1), 'labelled matching mean in terminal direct sum')
        high.append({'n': n, 'matchings': len(energies), 'P': p, 'C': c, 'mean': str(mean)})
    # One shared permutation gives nested blocks; compute all 8! labelled orders.
    xs = [(1, 0)]*4 + [(-1, 0)]*4
    sums = [0]*4
    ds = 0
    permutations = 0
    for order in itertools.permutations(range(8)):
        path = [sum(norm(plus(*(xs[j] for j in order[start:start+b]))) for start in range(0, 8, b)) for b in (8, 4, 2, 1)]
        sums = [a+b for a, b in zip(sums, path)]
        ds += debt(path)
        permutations += 1
    means = [F(s, permutations) for s in sums]
    require(means == [F(0), F(32, 7), F(48, 7), F(8)], 'uniform equal-block means')
    require(debt(means) == 0 and F(ds, permutations) == F(8, 7), 'actual Jensen debt obstruction')
    sign_tests = 0
    for signs in itertools.product((-1, 1), repeat=8):
        selected_tree([(x, 0) for x in signs])
        sign_tests += 1
    terminal_controls = [
        selected_tree([(1, 0, 0, 0), (0, 0, -1, 0)]),
        selected_tree([(1, 0, 0, 0), (-1, 0, 0, 0), (0, 0, 1, 0), (0, 0, -1, 0)]),
        selected_tree([(0, 0)]*8),
        selected_tree([(1, 0), (-1, 0), (0, 0), (0, 0)]),
        selected_tree([(1, 0)]),
    ]
    require(terminal_controls[0]['P'] == 2, 'different terminals do not coherently cancel')
    coherent = [selected_tree([(2, -1)]*n) for n in (1, 2, 4, 8)]
    require(all(r['D'] == F(len(r['blocks_bottom_up'][-1][0])-1, len(r['blocks_bottom_up'][-1][0]))*r['P'] if r['blocks_bottom_up'] else r['D'] == 0 for r in coherent), 'tree sharpness')
    fixed = []
    for h in range(1, 10):
        n = 2**h
        values = [(-1)**j for j in range(n)]
        path = [sum(sum(values[r::q])**2 for r in range(q)) for q in [2**k for k in range(h+1)]]
        require(path == [0]+[n*n//(2**k) for k in range(1, h+1)], 'independent fixed low-bit bins')
        require(debt(path) == n*n//2-n, 'unbounded fixed-tree family')
        fixed.append({'n': n, 'S': n, 'P': path[0], 'D': debt(path)})
    # The insertion and endpoint-fixed stability identities include all signs/orderings.
    for p, m, c in itertools.product(range(11), repeat=3):
        require(debt([p, m, c]) == max(p-c, 0)+excess(m, p, c), 'three energy nonnegative exhaustive identity')
    return {'four_tuple_count': count, 'sharp_family_samples': sharp_samples, 'hostile': hostile,
            'matching_controls': high, 'permutations': permutations, 'mean_energies': list(map(str, means)),
            'expected_actual_debt': str(F(ds, permutations)), 'sign_tree_controls': sign_tests,
            'terminal_padding_empty_controls': terminal_controls, 'coherent_controls': coherent,
            'fixed_low_bit_family': fixed}


def multiply(x, y):
    # Quotient ring Z[T]/(T^2+T+1), not producer's lookup-table source implementation.
    a, b = x
    c, d = y
    return (a*c-b*d, a*d+b*c-b*d)


def omega_power(k):
    answer = (1, 0)
    for _ in range(k % 3):
        answer = multiply(answer, (0, 1))
    return answer


def source_charge(seed, scale):
    chain = []
    charge = (0, 0)
    rows = [0, 0, 0]
    x = seed
    while x < 4*scale:
        factor = (1, 2) if x < 2*scale else (0, 1)  # omega-omega^2 or omega.
        term = multiply(omega_power(x), factor)
        chain.append({'x': x, 'lower': x < 2*scale, 'term': list(term)})
        charge = plus(charge, term)
        weights = (0, -3, 3) if x < 2*scale else (-1, -1, 2)
        rows = [rows[r]+weights[(x+r) % 3] for r in range(3)]
        x = 4*x+1
    a, b = charge
    require(rows == [2*a-b, -a-b, -a+2*b], 'literal B-chain versus all three real rows')
    require(norm(charge) in (0, 1, 3), 'native charge norm palette')
    return charge, chain


def edge(seed):
    require(seed > 0 and seed % 2 == 1 and seed % 8 != 5, 'native primitive seed guard')
    a = 1 + int(seed % 8 == 1)
    y, remainder = divmod(3*seed+1, 2**a)
    require(remainder == 0, 'exact native initial division')
    braid = 0
    while y % 8 == 5:
        y, remainder = divmod(y-1, 4)
        require(remainder == 0, 'exact inverse braid')
        braid += 1
    t = a+2*braid
    require(2**t*(3*y+1) == 3*(3*seed+1)+2**a, 'corrected native integer-clock identity')
    return y, t, a, braid


def load_frontier():
    data = (ROOT/'frontier.txt').read_bytes()
    require(hashlib.sha256(data).hexdigest() == FRONTIER_SHA, 'frontier content pin')
    rows = [tuple(map(int, line.split())) for line in data.decode().splitlines() if line and not line.startswith('#')]
    ranks = dict(rows)
    require(len(rows) == len(ranks) == 627, 'unique complete frontier')
    for s, rank in rows:
        y, _, _, _ = edge(s)
        require(y in ranks, 'frontier forward closed')
        require(s == 1 and y == 1 or ranks[y] < rank, 'frontier descending rank')
    return ranks


def reconstruct(scale, frontier):
    records, detailed = [], []
    for s in range(1, 4*scale, 2):
        if s % 8 == 5:
            continue
        x = s
        elapsed = 0
        visited = set()
        route = []
        while x not in frontier:
            require(x not in visited and len(route) < 10000, 'finite route unresolved, never silently omitted')
            visited.add(x)
            y, t, a, braid = edge(x)
            route.append({'from': x, 'to': y, 'clock': t, 'a': a, 'braid': braid})
            elapsed += t
            x = y
        charge, chain = source_charge(s, scale)
        sign = (-1)**elapsed
        record = {'s': s, 'u': x, 'p': route[-1]['from'] if route else 0,
                  't': elapsed, 'd': len(route), 'A': list(charge),
                  'v': [sign*c for c in charge], 'm': (s-1)//2}
        records.append(record)
        detailed.append(dict(record, source_B_chain=chain, first_entry_route=route))
    return records, detailed


def ledger(records, scale, cutoff, full_sources, include_pairs=True):
    arrived = [r for r in records if r['d'] <= cutoff]
    qs = [2**k for k in range((2*scale-1).bit_length()+1)]
    levels, energies, bin_counts = [], [], []
    terminals = {r['u'] for r in arrived}
    for q in qs:
        bins = defaultdict(list)
        for r in arrived:
            bins[(r['u'], r['m'] % q)].append(r)
        sums = {key: plus(*(r['v'] for r in group)) for key, group in bins.items()}
        levels.append(sums)
        energies.append(esum(sums.values()))
        bin_counts.append({'q': q, 'conceptual': len(terminals)*q, 'label_occupied': len(bins),
                           'nonempty_zero_sum': sum(norm(v) == 0 for v in sums.values())})
    deltas = [a-b for a, b in zip(energies, energies[1:])]
    sarr = sum(norm(r['A']) for r in arrived)
    require(energies[-1] == sarr, 'separated-source endpoint')
    require(energies[0]-sarr == sum(deltas), 'signed telescope')
    dplus = sum(max(d, 0) for d in deltas)
    require(energies[0] <= full_sources+dplus, 'full source denominator, including moving sources')
    result = {'N': scale, 'cutoff': cutoff, 'arrived': len(arrived), 'S_full': full_sources,
              'S_arrived': sarr, 'P': energies[0], 'qs': qs, 'energies': energies, 'deltas': deltas,
              'D': dplus, 'zero_charge_labels': sum(norm(r['A']) == 0 for r in arrived), 'bin_counts': bin_counts}
    if not include_pairs:
        return result
    owned = dict.fromkeys(qs[:-1], 0)
    pair_positive = 0
    pair_count = 0
    grouped = defaultdict(list)
    for r in arrived:
        grouped[r['u']].append(r)
    for group in grouped.values():
        for left, right in itertools.combinations(group, 2):
            difference = abs(left['m']-right['m'])
            require(difference > 0, 'labelled distinct addresses')
            owner = 1
            while difference % 2 == 0:
                owner *= 2
                difference //= 2
            b = cross(left['v'], right['v'])
            owned[owner] += b
            pair_positive += max(b, 0)
            pair_count += 1
    require([owned[q] for q in qs[:-1]] == deltas, 'literal v2 pair owner at every level')
    parent_positive = 0
    three_patterns = []
    for j, q in enumerate(qs[:-1]):
        children = levels[j+1]
        parents = {(u, r % q) for u, r in children}
        bs = [cross(children.get((u, r), (0, 0)), children.get((u, r+q), (0, 0))) for u, r in parents]
        require(sum(bs) == deltas[j], 'parent-bin polarized difference')
        parent_positive += sum(max(b, 0) for b in bs)
        if j+2 < len(levels):
            fine = levels[j+2]
            mids = [0, 0, 0]
            for u, r in {(u, r % q) for u, r in fine}:
                vectors = [fine.get((u, r+k*q), (0, 0)) for k in range(4)]
                for k, pairing in enumerate(PATTERNS):
                    mids[k] += esum(merge(vectors, pairing))
            require(mids[1] == energies[j+1], 'native pairing is 02|13, not 01|23')
            require(sum(mids) == energies[j]+2*energies[j+2], 'native global three-pattern mean')
            three_patterns.append({'q': q, 'P': energies[j], 'C': energies[j+2], 'middles': mids})
    require(dplus <= parent_positive <= pair_positive, 'three distinct positive-part placements')
    result.update(pair_positive=pair_positive, parent_positive=parent_positive, pairs=pair_count,
                  three_patterns=three_patterns)
    return result


def json_bytes(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'))+'\n').encode()


def native(producer, write_dir):
    frontier = load_frontier()
    sets, controls, cutoffs, files = {}, [], [], {}
    for scale in (512, 1024, 4096):
        records, detailed = reconstruct(scale, frontier)
        sets[scale] = records
        require(len(records) == 3*scale//2, 'all primitive original source labels present')
        payload = json_bytes(records)
        files[f'native_sources_N{scale}.json'] = hashlib.sha256(payload).hexdigest()
        if producer:
            given = json.loads((producer/f'native_sources_N{scale}.json').read_text())
            require(given == records, f'every producer source, charge, port, clock, depth and terminal N{scale}')
            require((producer/f'native_sources_N{scale}.json').read_bytes() == payload, 'canonical ledger bytes')
        if write_dir:
            (write_dir/f'independent_sources_N{scale}.json').write_bytes(json_bytes(detailed))
        full = sum(norm(r['A']) for r in records)
        maximum = max(r['d'] for r in records)
        declared = {maximum} | ({1, 3, 10} if scale == 4096 else set())
        for d in range(maximum+1):
            r = ledger(records, scale, d, full, d in declared)
            cutoffs.append({k: r[k] for k in ('N', 'cutoff', 'arrived', 'S_full', 'S_arrived', 'P', 'D')})
            if d in declared:
                controls.append(r)
    small = next(r for r in controls if r['N'] == 512)
    require([small[k] for k in ('arrived', 'S_full', 'P', 'pair_positive', 'parent_positive', 'D')] == [768, 1244, 1313, 433, 293, 79], 'recovered 433/293/79')
    require(small['energies'] == [1313, 1274, 1268, 1261, 1255, 1242, 1252, 1247, 1244, 1244, 1244], 'complete N512 level energies')
    witnesses = []
    for scale, q, terminal, residue, expected_labels in ((512, 1, 435, 0, [[1833], [435], [], [1375]]), (4096, 4, 443, 0, [[12609], [13993], [9457], []])):
        children = [[r for r in sets[scale] if r['u'] == terminal and r['m'] % (4*q) == residue+k*q] for k in range(4)]
        require([[r['s'] for r in group] for group in children] == expected_labels, 'complete parent native membership')
        vectors = [plus(*(r['v'] for r in group)) for group in children]
        result = four_test(vectors)
        require(result['P'] == result['C'] == 3 and sorted(result['middles']) == [2, 2, 5], 'native three-source obstruction')
        witnesses.append(dict(result, N=scale, q=q, terminal=terminal, residue=residue, children=children, vectors=vectors))
    by_label = {r['s']: r for r in sets[512]}
    observer_records = [by_label[s] for s in (799, 1065)]
    first, second = observer_records
    qcross = cross(first['v'], second['v'])
    marked = 2*(2*first['v'][0]-first['v'][1])*(2*second['v'][0]-second['v'][1])
    require(first['u'] == second['u'] == 799 and qcross == 0 and marked == 6, 'marked observer pair fully authenticated')
    observer = {'N': 512, 'records': observer_records, 'beta_Q': qcross, 'marked_cross': marked}
    receipt_status = 'not supplied'
    if producer:
        receipt = json.loads((producer/'receipt.json').read_text())['native']
        require(receipt['native_marked_cross_hostile'] == observer, 'producer receipt must serialize BOTH authenticated source records')
        require(len(controls) == len(receipt['controls']), 'complete declared native control list')
        for r, given in zip(controls, receipt['controls']):
            mapping = {'N': 'N', 'cutoff': 'depth_cutoff', 'arrived': 'arrived_sources', 'S_full': 'S_full', 'S_arrived': 'S_arrived', 'P': 'terminal_energy', 'D': 'D_plus', 'energies': 'energies', 'deltas': 'level_deltas', 'pair_positive': 'sum_positive_labelled_pair_interactions', 'parent_positive': 'sum_positive_parent_bin_interactions', 'pairs': 'same_terminal_labelled_pairs'}
            require(all(r[k] == given[v] for k, v in mapping.items()), 'all declared producer native controls')
        receipt_status = 'all declared native controls and both marked records agree'
    return {'frontier_count': len(frontier), 'frontier_sha256': FRONTIER_SHA, 'source_hashes': files,
            'declared_controls': controls, 'all_depth_cutoffs': cutoffs, 'local_witnesses': witnesses,
            'marked_observer': observer, 'producer_receipt': receipt_status}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--producer', type=Path)
    parser.add_argument('--write-sources', type=Path)
    args = parser.parse_args()
    if args.write_sources:
        args.write_sources.mkdir(parents=True, exist_ok=True)
    if args.producer:
        manifest = json.loads((args.producer/'MANIFEST.json').read_text())
        require(manifest.get('revision') == 2, 'audited producer revision2')
        for name, expected in manifest['files'].items():
            require(hashlib.sha256((args.producer/name).read_bytes()).hexdigest() == expected, 'producer frozen manifest: '+name)
        pins = json.loads((args.producer/'SOURCE_PINS.json').read_text())
        require(pins['commit'] == PIN, 'repository source commit')
        for source in pins['sources']:
            require(hashlib.sha256((args.producer/source['local_path']).read_bytes()).hexdigest() == source['sha256'], 'pinned canonical source: '+source['local_path'])
    result = {'status': 'PASS', 'source_commit': PIN, 'algebra': algebra(), 'native': native(args.producer, args.write_sources)}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
