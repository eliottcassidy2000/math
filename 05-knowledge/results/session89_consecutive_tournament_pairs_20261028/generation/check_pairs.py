#!/usr/bin/env python3
"""Aligned consecutive tournaments: full representatives and full stabilizers.

Production reuses the byte-pinned Session88 block/root certificates. Audit
independently exhausts the literal 2^(n*n) labelled universe through n=4.
No third-party dependencies; no baseline writes. Run with Python -B.
"""
from collections import Counter
from hashlib import sha256
from itertools import permutations
from math import factorial
from pathlib import Path
from statistics import median
from time import perf_counter
import argparse
import json
import sys
import types

PIN = '03408f84ac43fcdcfb0887e5d6bc1fc8792415cf'
SOURCES = {
 'generation': ('05-knowledge/results/session87_four_vertices_enumeration_20261028/check_generation.py',
  'f5d3da961f34dc1a32eaa0c1b69fad9bda983154578ea0187c5fd6b2e2d10fae'),
 'blockfill': ('05-knowledge/results/session88_symmetry_certificates_20261028/blockfill/check_blockfill.py',
  'cb246f280786f6a7f8b1f829d9d9b9e0993cb2298eb399783809bedb27780bdb')}

def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_baseline(root):
    if root is None:
        candidates = list(Path(__file__).resolve().parents)
    else:
        candidates = [Path(root).resolve()]
    root = next((p for p in candidates if (p / SOURCES['generation'][0]).is_file()), None)
    need(root is not None, 'Use --baseline-root PATH to the source repository')
    modules = {}
    for name, (relative, expected) in SOURCES.items():
        path = root / relative
        data = path.read_bytes()
        need(sha256(data).hexdigest() == expected, 'Changed baseline: ' + str(path))
        module = types.ModuleType('session89_' + name)
        module.__file__ = str(path)
        if name == 'blockfill':
            bridge = types.ModuleType('baseline_sources')
            bridge.generation = modules['generation']
            saved = sys.modules.get('baseline_sources')
            sys.modules['baseline_sources'] = bridge
        try:
            exec(compile(data, str(path), 'exec'), module.__dict__)
        finally:
            if name == 'blockfill':
                if saved is None:
                    del sys.modules['baseline_sources']
                else:
                    sys.modules['baseline_sources'] = saved
        modules[name] = module
    return modules['generation'], modules['blockfill']


def inverse(p):
    result = [0] * len(p)
    for i, j in enumerate(p):
        result[j] = i
    return tuple(result)


def prepare(n, old, block):
    """Identical charged setup for all methods, including S88 through n+1.

    Each p is a destination-to-source ordering: action sends edge A_ij to
    A_{p(i),p(j)}. Automorphism tuples can use either direction, as full groups.
    """
    start = perf_counter()
    parents = [()]
    levels = []
    core_result = None
    for k in range(n + 1):
        result = block.build(parents, k, 'blockfill')
        levels.append({key: value for key, value in result['counts'].items()
                       if key not in {'seconds', 'preparation_seconds'}})
        parents = list(result['bank'])
        if k + 1 == n:
            core_result = result
    pointed_result = result
    core = []
    engine = old.ParentCanonicalizer()
    if n == 0:
        core.append((0, (), ((),)))
    else:
        states_reverse = {child: key for key, child in core_result['states'].items()}
        for child in core_result['bank']:
            code, maps = engine.maps(child)
            p = maps[0]
            ip = inverse(p)
            raw_group = core_result['automorphisms'][states_reverse[child]]
            group = tuple(sorted(tuple(ip[g[p[i]]] for i in range(n)) for g in raw_group))
            rows = old.matrix(n, code)
            need(all(old.encode(rows, g) == code for g in group), 'bad transported S88 group')
            need(len(set(group)) == len(group), 'duplicated S88 group')
            core.append((code, rows, group))
    core.sort()
    perms = tuple(permutations(range(n)))
    all_group_elements = sorted({p for _, _, group in core for p in group})
    mask_action = {p: tuple(old.move_mask(mask, p) for mask in range(1 << n))
                   for p in all_group_elements}
    # Reuse every pointed state once. Record its certified root stabilizer.
    pointed = []
    for key, child in sorted(pointed_result['states'].items()):
        root_maps = pointed_result['root_transports'][key][n][1]
        H = tuple(sorted(p[:n] for p in root_maps))
        need(all(p[-1] == n for p in root_maps), 'root stabilizer moved root')
        pointed.append((key, H))
    # One named lookup, shared by all outputs, links to the S88 rooted dictionary.
    pointed_lookup = {}
    for (scode, mask), H in pointed:
        rows = old.matrix(n, scode)
        images = {(old.encode(rows, p), old.move_mask(mask, p)) for p in perms}
        need(len(images) * len(H) == factorial(n), 'S88 rooted orbit-stabilizer failed')
        for image in images:
            need(image not in pointed_lookup, 'S88 pointed classes overlap')
            pointed_lookup[image] = (scode, mask)
    expected = 1 << (n * (n - 1) // 2 + n)
    need(len(pointed_lookup) == expected, 'incomplete S88 rooted-labelled universe')
    return {'n': n, 'core': tuple(core), 'perms': perms, 'mask_action': mask_action,
            'pointed': tuple(pointed), 'pointed_lookup': pointed_lookup,
            'levels': levels, 'seconds': perf_counter() - start,
            'pointed_dictionary_sha256': digest(tuple(pointed)),
            'normalization_requests': engine.requests}


def digest(value):
    return sha256(repr(value).encode()).hexdigest()


def reference(prep, old):
    """Fair simple baseline: fix T's class, extend every labelled S by every mask.

    Reuse group actions on each S and precomputed mask permutations; reject a
    candidate immediately after finding a smaller image. All tests are joint.
    """
    start = perf_counter()
    n = prep['n']
    bank = {}
    count = Counter()
    for tcode, _, G in prep['core']:
        for scode in range(1 << (n * (n - 1) // 2)):
            rows = old.matrix(n, scode)
            actions = [(p, old.encode(rows, p)) for p in G]
            count['core_action_evaluations'] += len(G)
            for mask in range(1 << n):
                count['joint_candidates'] += 1
                stabilizer = []
                for p, image in actions:
                    count['joint_comparisons'] += 1
                    transformed = (image, prep['mask_action'][p][mask])
                    if transformed < (scode, mask):
                        break
                    if transformed == (scode, mask):
                        stabilizer.append(p)
                else:
                    bank[(tcode, scode, mask)] = (prep['pointed_lookup'][(scode, mask)], tuple(stabilizer))
    return bank, {**dict(count), 'seconds': perf_counter() - start}


def core_first(prep, old):
    """G\\S_n/H core alignments, then masks under K=G intersect gHg^-1.

    Materialize S_n/H once per ordinary S using its distinct labelled images.
    Traversal of each G-orbit uses all G, hence obtains the complete K exactly.
    Reuse the mask-orbit/stabilizer recipe for identical concrete subgroup K.
    """
    start = perf_counter()
    n = prep['n']
    bank = {}
    count = Counter()
    labelled = []
    for _, rows, H in prep['core']:
        images = tuple(sorted({old.encode(rows, p) for p in prep['perms']}))
        count['right_coset_permutation_visits'] += len(prep['perms'])
        need(len(images) * len(H) == factorial(n), 'ordinary right-coset size mismatch')
        labelled.append(images)
    recipes = {}
    for tcode, _, G in prep['core']:
        for images in labelled:
            seen = set()
            for scode in images:
                count['core_image_visits'] += 1
                if scode in seen:
                    continue
                rows = old.matrix(n, scode)
                actions = tuple((p, old.encode(rows, p)) for p in G)
                count['core_action_evaluations'] += len(G)
                orbit = {image for p, image in actions}
                need(min(orbit) == scode, 'bad core double-coset ordering')
                seen.update(orbit)
                K = tuple(p for p, image in actions if image == scode)
                count['core_double_cosets'] += 1
                if K not in recipes:
                    recipes[K] = []
                    seen_masks = set()
                    for mask in range(1 << n):
                        if mask in seen_masks:
                            continue
                        orbit_m = {prep['mask_action'][p][mask] for p in K}
                        count['mask_action_evaluations'] += len(K)
                        need(min(orbit_m) == mask, 'bad mask orbit ordering')
                        seen_masks.update(orbit_m)
                        stabilizer = tuple(p for p in K if prep['mask_action'][p][mask] == mask)
                        recipes[K].append((mask, stabilizer))
                    need(len(seen_masks) == 1 << n, 'mask orbit coverage failed')
                for mask, stabilizer in recipes[K]:
                    key = (tcode, scode, mask)
                    need(key not in bank, 'duplicate production pair')
                    bank[key] = (prep['pointed_lookup'][(scode, mask)], stabilizer)
            need(seen == set(images), 'double-coset coverage failed')
    count['distinct_residual_group_recipes'] = len(recipes)
    return bank, {**dict(count), 'seconds': perf_counter() - start}


def rooted_first(prep, old):
    """G\\S_n/H with H the certified rooted-child stabilizer.

    This is the other exact decomposition; it is not generally the fastest.
    """
    start = perf_counter()
    n = prep['n']
    bank = {}
    count = Counter()
    catalogues = []
    for (scode, mask), H in prep['pointed']:
        rows = old.matrix(n, scode)
        images = tuple(sorted({(old.encode(rows, p), old.move_mask(mask, p))
                               for p in prep['perms']}))
        count['right_coset_permutation_visits'] += len(prep['perms'])
        need(len(images) * len(H) == factorial(n), 'rooted right-coset mismatch')
        catalogues.append(((scode, mask), images))
    for tcode, _, G in prep['core']:
        for rooted_key, images in catalogues:
            seen = set()
            for scode, mask in images:
                count['rooted_image_visits'] += 1
                if (scode, mask) in seen:
                    continue
                rows = old.matrix(n, scode)
                actions = [(p, (old.encode(rows, p), prep['mask_action'][p][mask])) for p in G]
                count['joint_action_evaluations'] += len(G)
                orbit = {image for p, image in actions}
                need(min(orbit) == (scode, mask), 'bad rooted double-coset ordering')
                seen.update(orbit)
                stabilizer = tuple(p for p, image in actions if image == (scode, mask))
                key = (tcode, scode, mask)
                need(key not in bank, 'duplicate rooted-first pair')
                bank[key] = (rooted_key, stabilizer)
            need(seen == set(images), 'rooted double-coset coverage failed')
    count['rooted_double_cosets'] = len(bank)
    return bank, {**dict(count), 'seconds': perf_counter() - start}


# Audit-only implementation: no S87/S88 action or canonicalization functions.
def direct_edge_action(n, code, p):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    bits = {edge: (code >> k) & 1 for k, edge in enumerate(edges)}
    answer = 0
    for k, (i, j) in enumerate(edges):
        a, b = p[i], p[j]
        value = bits[(a, b)] if a < b else 1 - bits[(b, a)]
        answer |= value << k
    return answer


def direct_mask_action(n, mask, p):
    return sum(((mask >> p[i]) & 1) << i for i in range(n))


def literal_audit(n, bank):
    """Mark all literal pair orbits; independently recover every key/group.

    The arbitrary ordinary representative convention is imported only as target
    labels; its coverage and uniqueness are independently verified first.
    """
    need(n <= 4, 'literal independent scope limited to n<=4')
    start = perf_counter()
    perms = tuple(permutations(range(n)))
    k = n * (n - 1) // 2
    L = 1 << k
    edge = {(p, code): direct_edge_action(n, code, p) for p in perms for code in range(L)}
    masks = {(p, mask): direct_mask_action(n, mask, p) for p in perms for mask in range(1 << n)}
    ttargets = {key[0] for key in bank}
    normalize = {}
    for raw in range(L):
        targets = {edge[(p, raw)] for p in perms} & ttargets
        need(len(targets) == 1, 'ordinary bank not a full independent transversal')
        target = next(iter(targets))
        normalize[raw] = (target, tuple(p for p in perms if edge[(p, raw)] == target))
    universe = 1 << (n * n)
    visited = bytearray(universe)
    output = {}
    hist = Counter()
    for packed in range(universe):
        if visited[packed]:
            continue
        t, s, m = packed & (L - 1), (packed >> k) & (L - 1), packed >> (2 * k)
        orbit = {edge[(p, t)] | (edge[(p, s)] << k) | (masks[(p, m)] << (2 * k)) for p in perms}
        need(all(not visited[z] for z in orbit), 'independent orbits overlap')
        for z in orbit:
            visited[z] = 1
        target, maps = normalize[t]
        s0, m0 = min((edge[(p, s)], masks[(p, m)]) for p in maps)
        key = (target, s0, m0)
        stabilizer = tuple(p for p in perms if edge[(p, target)] == target
                           and edge[(p, s0)] == s0 and masks[(p, m0)] == m0)
        score, smaps = normalize[s0]
        rooted_key = (score, min(masks[(p, m0)] for p in smaps))
        need(key not in output, 'duplicate independent normalized output')
        output[key] = (rooted_key, stabilizer)
        hist[len(orbit)] += 1
    need(all(visited), 'literal universe incomplete')
    need(output == bank, 'independent FULL key/group/root dictionary mismatch')
    return {'universe_size': universe, 'class_count': len(output),
            'orbit_size_histogram': dict(sorted(hist.items())),
            'every_labelled_input_covered_once': True,
            'all_keys_root_classes_and_complete_stabilizers_equal': True,
            'seconds': perf_counter() - start}


def semantic(report):
    if isinstance(report, dict):
        return {k: semantic(v) for k, v in report.items()
                if k not in {'seconds', 'setup_seconds', 'total_seconds', 'benchmarks'}}
    if isinstance(report, (list, tuple)):
        return [semantic(x) for x in report]
    return report


def semantic_hash(report):
    return sha256(json.dumps(semantic(report), sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def run(n, old, block, repeats=0, include_rooted=True):
    prep = prepare(n, old, block)
    methods = {'reference': reference, 'core_first': core_first}
    if include_rooted:
        methods['rooted_first'] = rooted_first
    records = {}
    expected = None
    for name, method in methods.items():
        bank, counts = method(prep, old)
        if expected is None:
            expected = bank
        else:
            need(bank == expected, 'method full-output mismatch: ' + name)
        records[name] = {**counts, 'setup_seconds': prep['seconds'],
                         'total_seconds': counts['seconds'] + prep['seconds']}
    bank = expected
    nested_count = sum(tcode == scode for tcode, scode, mask in bank)
    need(nested_count == len(prep['pointed']), 'nested child/core extraction mismatch')
    histogram = Counter(len(group) for root, group in bank.values())
    need(sum(factorial(n) // order * count for order, count in histogram.items()) == 1 << (n*n),
         'labelled pair mass mismatch')
    report = {'n': n, 'ordinary_core_classes': len(prep['core']),
              'pointed_child_classes': len(prep['pointed']), 'pair_classes': len(bank),
              'nested_U_minus_root_equals_T_classes': nested_count,
              'full_output_sha256': digest(tuple(sorted(bank.items()))),
              'stabilizer_order_histogram': dict(sorted(histogram.items())),
              'labelled_mass': 1 << (n*n),
              'same_full_keys_named_root_classes_and_full_stabilizers': True,
              'setup': {'s88_levels': prep['levels'],
                        'normalization_requests': prep['normalization_requests'],
                        'pointed_dictionary_sha256': prep['pointed_dictionary_sha256']},
              'methods': records}
    if n <= 4:
        report['independent_literal_audit'] = literal_audit(n, bank)
    if repeats:
        timings = {name: [] for name in methods}
        setups = {name: [] for name in methods}
        for i in range(repeats):
            names = list(methods)
            if i % 2:
                names.reverse()
            for name in names:
                # Charge every fresh S88 setup to every individual trial.
                trial_prep = prepare(n, old, block)
                trial_bank, trial_counts = methods[name](trial_prep, old)
                need(trial_bank == bank, 'benchmark full-output changed')
                timings[name].append(trial_prep['seconds'] + trial_counts['seconds'])
                setups[name].append(trial_prep['seconds'])
                del trial_bank, trial_prep
        medians = {name: median(values) for name, values in timings.items()}
        report['benchmarks'] = {'repeats': repeats, 'fresh_setup_each_trial': True,
                                'total_seconds': timings, 'setup_seconds': setups,
                                'median_total_seconds': medians,
                                'reference_over_core_first': medians['reference'] / medians['core_first']}
    return report


def hostile_controls(old, block):
    # C3 has three automorphisms; a transitive aligned S has only identity.
    prep = prepare(3, old, block)
    bank, _ = core_first(prep, old)
    G = next(group for code, rows, group in prep['core'] if code == 2)
    K = tuple(p for p in G if direct_edge_action(3, 0, p) == 0)
    need(len(G) == 3 and len(K) == 1, 'residual-intersection hostile changed')
    need(any(direct_mask_action(3, 1, p) == 2 for p in G), 'G must identify the hostile masks')
    need((2, 0, 1) in bank and (2, 0, 2) in bank,
         'naively quotienting by G lost the distinguishing hostile pairs')
    prep4 = prepare(4, old, block)
    groups3 = [G for code, rows, G in prep4['core'] if len(G) == 3]
    need(len(groups3) >= 2, 'missing equal-order concrete groups')
    A, B = groups3[:2]
    witness = next(mask for mask in range(16)
                   if {direct_mask_action(4, mask, p) for p in A}
                   != {direct_mask_action(4, mask, p) for p in B})
    edges7 = [(i, j) for i in range(7) for j in range(i + 1, 7)]
    code7 = sum(int((j - i) % 7 in {1, 2, 4}) << k for k, (i, j) in enumerate(edges7))
    group7 = tuple(p for p in permutations(range(7)) if direct_edge_action(7, code7, p) == code7)
    root7 = tuple(p for p in group7 if p[0] == 0)
    need(len(group7) == 21 and len(root7) == 3, 'P7 root-stabilizer hostile changed')
    return {'class_product_is_not_pair_count':
                {'n': 2, 'ordinary_T_times_pointed_U': 4, 'actual_aligned_pairs': 8},
            'must_use_residual_intersection_not_first_group':
                {'n': 3, 'T_code': 2, 'S_code': 0, 'G_order': len(G), 'K_order': len(K),
                 'masks_identified_by_G_but_distinct_as_pairs': [1, 2]},
            'pointed_type_does_not_fix_a_labelled_orientation':
                {'n': 7, 'arcs': 'i -> j iff (j-i) mod 7 is in {1,2,4}',
                 'root': 0, 'full_automorphism_order': len(group7),
                 'full_root_stabilizer': root7,
                 'labelled_realizations_with_root_fixed': factorial(6) // len(root7)},
            'equal_group_orders_do_not_identify_mask_recipes':
                {'n': 4, 'group_order': 3, 'groups': [A, B], 'mask': witness,
                 'mask_orbits': [sorted({direct_mask_action(4, witness, p) for p in G})
                                 for G in [A, B]]}}



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline-root')
    parser.add_argument('--max-n', type=int, default=4)
    parser.add_argument('--benchmark', type=int, default=0)
    parser.add_argument('--skip-rooted-first', action='store_true')
    args = parser.parse_args()
    need(0 <= args.max_n <= 5, 'finite tested materialization limit is n=5')
    old, block = load_baseline(args.baseline_root)
    results = [run(n, old, block, args.benchmark, not args.skip_rooted_first)
               for n in range(args.max_n + 1)]
    need([r['pair_classes'] for r in results] == [1, 2, 8, 88, 2752, 279968][:args.max_n + 1],
         'consecutive pair count regression')
    return {'status': 'PASS', 'scope': 'ordinary tournaments, root fixed, diagonal shared-core S_n; all masks',
            'not_standard_binary_matrix_conjugation': True,
            'baseline_commit': PIN, 'baseline_source_sha256': {k: v[1] for k,v in SOURCES.items()},
            'results': results, 'hostile_controls': hostile_controls(old, block),
            'semantic_sha256': semantic_hash(results)}

if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
