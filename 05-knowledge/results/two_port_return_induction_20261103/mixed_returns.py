"""Small exact certificate for the four-macro mixed-return obstruction.

No support census and no producer imports. Read only the pinned atom table
and the published 30-node repair receipt. Words below use processing order.
Run: python3 mixed_returns.py > mixed_returns.json
"""
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
ATLAS_PATH = RESULTS / 'golden_seam_recurrence_20261103/atom_atlas.json'
REPAIR_PATH = RESULTS / 'ternary_repair_clock_20261103/receipt.json'
ATLAS_SHA256 = 'b7f7715f5e6577c056da67ea2752f4643bd4ae16a8ef06707c2286032b434eba'
TAU = (2, -1, 0)
B = ((0, 1), (0, 2), (1, -1))
ROOT = (1, 0, 0)


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def matmul(a, b):
    return tuple(tuple(sum(x * b[k][j] for k, x in enumerate(row))
                       for j in range(len(b[0]))) for row in a)


def step(v, p):
    return matmul((v,), p)[0]


def identity(n):
    return tuple(tuple(int(i == j) for j in range(n)) for i in range(n))


def det2(q):
    return q[0][0] * q[1][1] - q[0][1] * q[1][0]


def inv(q):
    n = len(q)
    a = [[Fraction(x) for x in row] + [Fraction(i == j) for j in range(n)]
         for i, row in enumerate(q)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        z = a[j][j]
        a[j] = [x / z for x in a[j]]
        for i in range(n):
            if i != j:
                z = a[i][j]
                a[i] = [x - z * y for x, y in zip(a[i], a[j])]
    return tuple(tuple(row[n:]) for row in a)


def val3(z):
    z = abs(z)
    need(z != 0, 'valuation called at zero')
    k = 0
    while z % 3 == 0:
        z //= 3
        k += 1
    return k


def level(beta):
    return min(val3(z) for z in beta if z)


def height(beta):
    return beta[0] // 3 ** level(beta)


def cone(beta):
    y, x = beta
    return y > 0 and 0 <= x <= 3 * y


def physical(v):
    h, x, y = v
    return h >= 1 and h % 2 == 1 and 0 <= x <= y <= h + x and h <= x + y


def ratio(beta):
    return Fraction(beta[1], beta[0])


def frac(z):
    return str(z)


def matrix_word(word, matrices):
    p = identity(3)
    for label in word:
        p = matmul(p, matrices[label])
    return p


def value(word, matrices):
    return step(ROOT, matrix_word(word, matrices))


def actual_rows(word, atlas):
    # Construct the ordered source SCC word, retain its adjacent spine arcs,
    # complement all other arcs, then reverse labels into processing order.
    rows = []
    for label in reversed(word):
        block = atlas[label]['source_rows']
        n, m = len(rows), len(block)
        rows = [row | (((1 << m) - 1) << n) for row in rows]
        rows += [row << n for row in block]
    n = len(rows)
    for i in range(n):
        for j in range(i + 2, n):
            rows[i] ^= 1 << j
            rows[j] ^= 1 << i
    return tuple(sum(1 << (n - 1 - j) for j in range(n)
                     if rows[n - 1 - i] >> j & 1) for i in range(n))


def count_paths(rows):
    @lru_cache(None)
    def paths(last, left):
        if not left:
            return 1
        choices = rows[last] & left
        total = 0
        while choices:
            bit = choices & -choices
            choices -= bit
            total += paths(bit.bit_length() - 1, left ^ bit)
        return total
    allv = (1 << len(rows)) - 1
    return sum(paths(i, allv ^ (1 << i)) for i in range(len(rows)))


def interval(q):
    endpoints = [ratio(step(beta, q)) for beta in ((1, 0), (1, 3))]
    return min(endpoints), max(endpoints)


def charge_interval(beta):
    y, x = beta
    return {'g_min': max(x - 2 * y, -x, 1 - y), 'g_max': x // 3,
            'g_mod_2': x % 2, 'odd_h_possible': (x + y) % 2 == 1}


def decode_bulk(target, qs, ps, intervals):
    """All possible nonempty suffixes; only physical integer bulk guards.

    Every real root is represented, but acceptance is not actual membership.
    Disjoint last-letter image intervals make this a single reverse chain.
    """
    current = tuple(3 * z for z in step(target, B))
    suffix = ()
    p = identity(3)
    q_product = identity(2)
    stages = []
    while True:
        candidates = [label for label, (lo, hi) in intervals.items()
                      if lo <= ratio(current) <= hi]
        need(len(candidates) <= 1, 'last-letter intervals overlap')
        if not candidates:
            return stages, 'OUTSIDE_ALL_LAST_MACRO_CONES'
        label = candidates[0]
        parent = step(current, inv(qs[label]))
        if any(z.denominator != 1 for z in parent):
            return stages, 'NONINTEGRAL_LAST_MACRO_INVERSE'
        parent = tuple(int(z) for z in parent)
        if not cone(parent):
            return stages, 'OUTSIDE_PHYSICAL_BULK_CONE'
        need(parent[0] < current[0], 'inverse descent')
        suffix = (label,) + suffix
        p = matmul(ps[label], p)
        q_product = matmul(qs[label], q_product)
        need(step(parent, q_product) == tuple(3 * z for z in step(target, B)),
             'decoded bulk intertwining')
        fixed_offsets = []
        for m in (-2, -1, 1, 2):
            child = tuple(3 * z + m * t for z, t in zip(target, TAU))
            w = step(child, inv(p))
            if all(z.denominator == 1 for z in w):
                w = tuple(int(z) for z in w)
                if physical(w):
                    need(step(w, p) == child, 'full-state inverse')
                    fixed_offsets.append({'m': m, 'physical_integer_parent': w})
        stages.append({'return_word': suffix, 'parent_bulk': parent,
                       'charge_fibre': charge_interval(parent),
                       'offset_candidates_NOT_actual_certificates': fixed_offsets})
        current = parent


def source_inverse(atlas, matrices):
    """Recheck the inherited complete inverse guards for exactly 66 letters.

    A negative answer is relative to this alphabet, never to all tournaments.
    No scalar saturation is performed.
    """
    six = ('s14', 's268', 's330', 's248', 's261', 's85', 's262', 's453', 's78')
    labels = tuple(k for k, data in atlas.items() if len(data['source_rows']) <= 5) + six
    need(len(labels) == 66, 'declared exact inverse alphabet')
    weights = (6, 9, 5)
    norm = lambda v: sum(x * y for x, y in zip(v, weights))
    rays16 = ((1, 0, 0), (1, 16, 0), (1, 0, 16), (1, 16, 16))
    physical_rays = ((2, 1, 1), (1, 0, 1), (0, 1, 1))
    for label in labels:
        p = matrices[label]
        need(all(z >= 0 for row in p for z in row), ('nonnegative matrix', label))
        need(physical(p[0]), ('actual first state', label))
        need(p[0][0] >= 1, ('H monotone', label))
        for ray in rays16:
            h, x, y = step(ray, p)
            need(0 <= x <= 16 * h and 0 <= y <= 16 * h, ('C16 invariant', label))
        for ray in physical_rays:
            h, x, y = step(ray, p)
            need(0 <= x <= y <= h + x and h <= x + y, ('physical cone invariant', label))
        need(all(5 * norm(row) >= 9 * weights[i] for i, row in enumerate(p)),
             ('strict inverse norm', label))
    inverses = {label: inv(matrices[label]) for label in labels}

    @lru_cache(None)
    def decode(v):
        if v == ROOT:
            return ()
        if not physical(v) or v[1] > 16 * v[0] or v[2] > 16 * v[0]:
            return None
        options = []
        for label in labels:
            parent = step(v, inverses[label])
            if any(z.denominator != 1 for z in parent):
                continue
            parent = tuple(int(z) for z in parent)
            if parent != ROOT and (not physical(parent) or parent[1] > 16 * parent[0]
                                   or parent[2] > 16 * parent[0]):
                continue
            if parent[0] > v[0] or 9 * norm(parent) > 5 * norm(v):
                continue
            options.append((norm(parent), len(atlas[label]['source_rows']), label, parent))
        for _, _, label, parent in sorted(options):
            word = decode(parent)
            if word is not None:
                result = word + (label,)
                need(value(result, matrices) == v, ('source inverse replay', result))
                return result
        return None
    return decode, labels


def main():
    need(sha(ATLAS_PATH) == ATLAS_SHA256, 'pinned atlas changed')
    atlas = json.loads(ATLAS_PATH.read_text())
    repair = json.loads(REPAIR_PATH.read_text())
    matrices = {k: tuple(map(tuple, v['matrix'])) for k, v in atlas.items()}
    macro_words = {'b': ('b',), 'q2': ('q2',), 's78': ('s78',),
                   'D': ('r24', 'a', 'a', 'a', 'b', 'r24', 'a', 'a', 'a')}
    expected_q = {'b': ((5, 8), (2, 3)), 'q2': ((11, 18), (4, 7)),
                  's78': ((77, 152), (28, 55)),
                  'D': ((715050, 1009641), (268050, 378483))}
    expected_gain = {'b': -1, 'q2': 1, 's78': 1, 'D': -4}
    ps = {k: matrix_word(w, matrices) for k, w in macro_words.items()}
    qs = expected_q
    intervals = {k: interval(q) for k, q in qs.items()}
    ordered = sorted(intervals, key=lambda k: intervals[k][0])
    need(all(intervals[a][1] < intervals[b][0] for a, b in zip(ordered, ordered[1:])),
         'pairwise disjoint closed image intervals')
    macro_data = {}
    for label, p in ps.items():
        q = qs[label]
        need(matmul(p, B) == matmul(B, q), ('saturated quotient', label))
        need(step(TAU, p) == tuple(expected_gain[label] * z for z in TAU),
             ('error line', label))
        need(all(z > 0 for row in q for z in row), ('strict positivity', label))
        determinant, c = det2(q), q[0][0]
        need(c > 3 ** val3(determinant), ('growth defeats determinant', label))
        need(all(cone(step(beta, q)) for beta in ((1, 0), (1, 3))),
             ('physical bulk cone invariant', label))
        linear = tuple(25 * q[i][1] - 34 * q[i][0] for i in range(2))
        need(min(linear) > 0, ('original seed ray separator', label))
        macro_data[label] = {'processing_word': macro_words[label], 'P': p, 'Q': q,
                             'tau_gain': expected_gain[label], 'det_Q': determinant,
                             'v3_det': val3(determinant), 'Y_growth_lower': c,
                             'normalized_height_growth_lower': frac(Fraction(c, 3 ** val3(determinant))),
                             'image_ratio_interval': list(map(frac, intervals[label])),
                             '25X_minus_34Y_after_Q_coefficients': linear}

    # A separate seven-macro obstruction certificate. These three additions
    # are NOT included in the disjoint-interval decoder above or below.
    extra_words = {'L6': ('r24', 'q4'), 'L45': ('r19', 'r32', 'r48', 'q4'),
                   'L13': ('r19', 'r19', 'a', 'a', 'r18')}
    extra_q = {'L6': ((720, 1101), (270, 413)),
               'L45': ((1020432, 1560769), (366879, 561148)),
               'L13': ((211965, 378780), (76220, 136205))}
    extra_gain = {'L6': 6, 'L45': 45, 'L13': 13}
    extra_data = {}
    for label, word in extra_words.items():
        p, q = matrix_word(word, matrices), extra_q[label]
        need(matmul(p, B) == matmul(B, q), ('extra loop saturated quotient', label))
        need(step(TAU, p) == tuple(extra_gain[label] * z for z in TAU),
             ('extra loop error line', label))
        need(all(z > 0 for row in q for z in row), ('extra loop positivity', label))
        need(all(cone(step(beta, q)) for beta in ((1, 0), (1, 3))),
             ('extra loop physical bulk cone', label))
        determinant, c = det2(q), q[0][0]
        need(c > 3 ** val3(determinant), ('extra loop height criterion', label))
        extra_data[label] = {'processing_word': word, 'P': p, 'Q': q,
                             'tau_gain': extra_gain[label], 'det_Q': determinant,
                             'v3_det': val3(determinant), 'Y_growth_lower': c,
                             'normalized_height_growth_lower': frac(Fraction(c, 3 ** val3(determinant))),
                             'vertices_added': sum(len(atlas[k]['source_rows']) for k in word)}

    seeds = {'U': ('r5', 'a', 'a'), 'U2': ('a', 'a', 'q1', 'a', 'a'),
             'U4': ('q0', 'a', 'a', 'a', 'a', 'a')}
    seed_data = {}
    for label, word in seeds.items():
        v = value(word, matrices)
        counts = [count_paths(actual_rows(word + ('a',) * j, atlas)) for j in range(3)]
        need(tuple((counts[0], (counts[1] - counts[0]) // 2,
                    (counts[2] - counts[1]) // 2)) == v, ('independent seed state', label))
        seed_data[label] = {'processing_word': word, 'state': v, 'bulk': step(v, B),
                            'singleton_queries': counts, 'normalized_Y': height(step(v, B))}
    rows = actual_rows(seeds['U'], atlas)
    literal = sum(all(rows[a] >> b & 1 for a, b in zip(p, p[1:]))
                  for p in permutations(range(len(rows))))
    need(literal == 31, '7! literal positive seed')

    targets = [('original_' + k, value(word, matrices)) for k, word in seeds.items()]
    targets += [('compact_' + str(i), tuple(z['target_state']))
                for i, z in enumerate(repair['compact_repairs'])]
    targets += [('clock_twin_' + str(i), tuple(v))
                for i, v in enumerate(repair['clock_twins']['states'])]
    targets += [('tree_' + str(i), tuple(z['state']))
                for i, z in enumerate(repair['four_level_tree'])]
    target_data = []
    for name, v in targets:
        need(physical(v), ('input physical guard', name))
        stages, stop = decode_bulk(v, qs, ps, intervals)
        target_data.append({'name': name, 'state': v, 'bulk': step(v, B),
                            'ratio': frac(ratio(step(v, B))), 'stages': stages, 'stop': stop})
    need(all(not z['stages'] for z in target_data if z['name'].startswith(('original_', 'compact_', 'clock_twin_'))),
         'named hostile endpoints have no physical nonempty inverse')

    decode_source, inverse_labels = source_inverse(atlas, matrices)
    for label, word in seeds.items():
        decoded = decode_source(value(word, matrices))
        need(decoded is not None and value(decoded, matrices) == value(word, matrices),
             ('positive complete-inverse control', label))
    inverse_queries = []
    for target in target_data:
        if not target['name'].startswith('tree_'):
            continue
        for stage in target['stages']:
            for candidate in stage['offset_candidates_NOT_actual_certificates']:
                state = tuple(candidate['physical_integer_parent'])
                word = decode_source(state)
                result = {'target': target['name'], 'return_word': stage['return_word'],
                          'm': candidate['m'], 'candidate_state': state,
                          'status': 'ACTUAL_66_BANK_SOURCE' if word is not None else 'ABSENT_IN_66_BANK'}
                if word is not None:
                    rows = actual_rows(word, atlas)
                    need(count_paths(rows) == state[0], ('literal candidate source H', state))
                    result.update({'processing_word': word, 'rows': rows})
                inverse_queries.append(result)

    positive_data = []
    for macro_word in (('b',), ('q2',), ('s78',), ('b', 'q2'), ('s78', 'b'), ('D',)):
        atoms = sum((macro_words[k] for k in macro_word), ())
        v = value(seeds['U'] + atoms, matrices)
        parent = value(seeds['U2'], matrices)
        p = matrix_word(atoms, matrices)
        m = 1
        for label in macro_word:
            m *= expected_gain[label]
        need(step(parent, p) == tuple(3 * z + m * t for z, t in zip(v, TAU)),
             ('inherited transported positive return', macro_word))
        stages, stop = decode_bulk(v, qs, ps, intervals)
        need(any(tuple(stage['return_word']) == macro_word and
                 tuple(stage['parent_bulk']) == step(parent, B) for stage in stages),
             ('decoder recovers known actual mixed return', macro_word))
        positive_data.append({'target_processing_word': seeds['U'] + atoms,
                              'target_state': v, 'actual_parent_word': seeds['U2'],
                              'actual_parent_state': parent, 'return_word': macro_word,
                              'm': m, 'stages': stages, 'stop': stop})

    print(json.dumps({'status': 'PASS', 'scope': 'seven-macro height obstruction; four-macro physical bulk decoder; no coverage census',
                      'source_sha256_lf': {ATLAS_PATH.relative_to(RESULTS).as_posix(): sha(ATLAS_PATH),
                                           REPAIR_PATH.relative_to(RESULTS).as_posix(): sha(REPAIR_PATH)},
                      'macros': macro_data, 'interval_order': ordered,
                      'separate_seven_macro_obstruction_additions': extra_data,
                      'actual_seed_controls': seed_data, 'literal_U_permutations': 5040,
                      'literal_U_count': literal, 'hostile_targets': target_data,
                      'source_inverse_alphabet': inverse_labels,
                      'source_inverse_queries': inverse_queries,
                      'inherited_transport_positive_controls': positive_data},
                     sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
