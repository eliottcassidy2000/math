"""Small exact audits of correction-channel returns; no spectrum census.

Run from any directory: python3 /absolute/path/to/return_bank_probe.py
Only the adjacent return_bank.json is written. All finite scope is printed.
"""
from functools import lru_cache
from itertools import product, permutations
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
ATLAS_PATH = RESULTS / 'golden_seam_recurrence_20261103/atom_atlas.json'
ATLAS_SHA256 = 'b7f7715f5e6577c056da67ea2752f4643bd4ae16a8ef06707c2286032b434eba'
CHANNELS = {
    'bq2': (2, -1, 0), 'q1': (3, -2, 1),
    'q4': (3, -2, 0), 'q5': (2, -2, 1), 'r21': (3, -1, -1),
}
BASIS = {
    'bq2': ((0, 1), (0, 2), (1, -1)),
    'r21': ((0, 1), (1, 3), (-1, 0)),
    'q5': ((1, 0), (1, 1), (0, 2)),
}


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def row(v, m):
    return tuple(sum(v[i] * m[i][j] for i in range(len(v))) for j in range(len(m[0])))


def mul(a, b):
    return tuple(row(v, b) for v in a)


def det(m):
    if len(m) == 2:
        return m[0][0] * m[1][1] - m[0][1] * m[1][0]
    return sum((-1) ** j * m[0][j] * det(tuple(tuple(m[i][k] for k in range(3) if k != j)
                                                for i in (1, 2))) for j in range(3))


def shifted(p, sign):
    return tuple(tuple(p[i][j] - (sign if i == j else 0) for j in range(3)) for i in range(3))


def quotient(p, b):
    """Recover Q from PB=BQ using a unimodular two-row minor, then audit."""
    z = mul(p, b)
    for i in range(3):
        for j in range(i + 1, 3):
            a = (b[i], b[j])
            d = det(a)
            if abs(d) == 1:
                inv = ((a[1][1] // d, -a[0][1] // d), (-a[1][0] // d, a[0][0] // d))
                q = mul(inv, (z[i], z[j]))
                need(mul(b, q) == z, 'quotient intertwining')
                return q
    raise RuntimeError('No saturated quotient minor')


def valuation(v, p=3):
    need(any(v), 'nonzero bulk vector')
    power = 0
    while all(x % p == 0 for x in v):
        v = tuple(x // p for x in v)
        power += 1
    return power


def main():
    raw = ATLAS_PATH.read_bytes().replace(b'\r\n', b'\n')
    need(hashlib.sha256(raw).hexdigest() == ATLAS_SHA256, 'pinned source atlas')
    atlas = json.loads(raw)
    labels = tuple(k for k, v in atlas.items() if len(v['source_rows']) <= 5)
    labels += ('s14', 's268', 's330', 's248', 's261', 's85', 's262', 's453')
    need(len(labels) == 65, 'exact inherited presentation bank')
    mats = {k: tuple(map(tuple, atlas[k]['matrix'])) for k in labels}

    def vector(word):
        v = (1, 0, 0)
        for l in word:
            v = row(v, mats[l])
        return v

    def graph(word):
        # Reconstruct the actual reversed source-component word, toggle precisely
        # the nonspine pairs, then reverse labels into processing-port convention.
        sizes = [len(atlas[l]['source_rows']) for l in reversed(word)]
        n = sum(sizes)
        arcs = [[False] * n for _ in range(n)]
        at = 0
        for l, size in zip(reversed(word), sizes):
            src = atlas[l]['source_rows']
            for i in range(size):
                for j in range(size):
                    arcs[at + i][at + j] = bool(src[i] >> j & 1)
                for j in range(at + size, n):
                    arcs[at + i][j] = True
            at += size
        for i in range(n):
            for j in range(i + 2, n):
                arcs[i][j], arcs[j][i] = arcs[j][i], arcs[i][j]
        return tuple(sum(1 << j for j in range(n) if arcs[n - 1 - i][n - 1 - j]) for i in range(n))

    def count(rows):
        n = len(rows)
        @lru_cache(None)
        def complete(last, remaining):
            if not remaining:
                return 1
            choices = rows[last] & remaining
            total = 0
            while choices:
                bit = choices & -choices
                choices -= bit
                total += complete(bit.bit_length() - 1, remaining ^ bit)
            return total
        return sum(complete(i, ((1 << n) - 1) ^ (1 << i)) for i in range(n))

    def strong(rows):
        n = len(rows)
        for source in range(n):
            seen = 1 << source
            todo = [source]
            while todo:
                v = todo.pop()
                for u in range(n):
                    if rows[v] >> u & 1 and not seen >> u & 1:
                        seen |= 1 << u
                        todo.append(u)
            if seen != (1 << n) - 1:
                return False
        return True

    quotient_records = []
    for channel, letters in [('bq2', ('b', 'q2')), ('r21', ('r21', 'r38')), ('q5', ('q5',))]:
        b = BASIS[channel]
        need(row(CHANNELS[channel], b) == (0, 0), 'annihilator basis')
        for l in letters:
            q = quotient(mats[l], b)
            quotient_records.append({'channel': channel, 'letter': l, 'basis': b,
                                     'matrix': q, 'determinant': det(q), 'unit_mod3': det(q) % 3 != 0})

    # Complete one-letter graph on these FIVE inherited primitive error lines.
    # Test rational proportionality with integer cross products, not floats.
    bridges = []
    for source, tau in CHANNELS.items():
        for l in labels:
            image = row(tau, mats[l])
            for target, other in CHANNELS.items():
                k = next(i for i, x in enumerate(other) if x)
                if all(image[i] * other[k] == image[k] * other[i] for i in range(3)):
                    need(image[k] % other[k] == 0, 'integral line multiplier')
                    bridges.append({'source': source, 'letter': l, 'target': target,
                                    'multiplier': image[k] // other[k]})
    need(len(bridges) == 10, 'five-channel one-letter graph cardinality')

    two_letter_units = []
    for a, b in product(labels, repeat=2):
        p = mul(mats[a], mats[b])
        for sign in (1, -1):
            if det(shifted(p, sign)) == 0:
                two_letter_units.append({'word': (a, b), 'sign': sign, 'determinant': det(p)})
    need(len(two_letter_units) == 11, 'two-letter signed-unit audit')

    # Small H-bounded saturation supplies exact discovery examples; every lift
    # requires BOTH source and target states inside the bound. This is not a
    # global classification of lifts from source H<=500.
    cap = 500
    states = [(1, 0, 0)]
    words = {states[0]: ()}
    for v in states:
        for l in labels:
            w = row(v, mats[l])
            if w[0] <= cap and w not in words:
                words[w] = words[v] + (l,)
                states.append(w)
    lifts = {}
    for name, tau in CHANNELS.items():
        multipliers = (1, 2) if tau[0] % 2 == 0 else (2, 4)
        hits = []
        for v in states[1:]:
            outputs = []
            for m in multipliers:
                w = tuple(3 * x + m * y for x, y in zip(v, tau))
                if w in words:
                    outputs.append({'m': m, 'state': w, 'word': words[w]})
            if outputs:
                hits.append({'state': v, 'word': words[v], 'outputs': outputs})
        lifts[name] = hits
    need(len(states) == 3592, 'small inherited-state universe')
    need({k: len(v) for k, v in lifts.items()} == {'bq2': 25, 'q1': 5, 'q4': 4, 'q5': 8, 'r21': 3},
         'discovery counts')

    # Actual q5-to-bq2 bridge: a new transported replacement, checked on literal
    # graphs. Generality follows from the displayed exact matrix identities.
    u, v = ('q2',), ('a',) * 6
    need(vector(v) == tuple(3*x+y for x,y in zip(vector(u), CHANNELS['q5'])), 'q5 actual seed lift')
    need(row(CHANNELS['q5'], mats['q5']) == tuple(-x for x in CHANNELS['q5']), 'q5 phase')
    need(row(CHANNELS['q5'], mats['a']) == tuple(-x for x in CHANNELS['bq2']), 'singleton bridge')
    actual = []
    for k, z in [(0, ()), (0, ('b',)), (1, ())]:
        left, right = u + ('q5',)*k + ('a',) + z, v + ('q5',)*k + ('a',) + z
        left_rows, right_rows = graph(left), graph(right)
        h, target = vector(left)[0], vector(right)[0]
        offset = 2 * (-1)**(k + 1 + z.count('b'))
        need(target == 3*h + offset, 'switched affine law')
        need(count(left_rows) == h and count(right_rows) == target, 'independent literal-graph count')
        need(strong(left_rows) and strong(right_rows), 'actual strong controls')
        if len(right_rows) <= 7:
            literal = sum(all(right_rows[p[i]] >> p[i+1] & 1 for i in range(len(p)-1))
                          for p in permutations(range(len(right_rows))))
            need(literal == target, 'literal permutations control')
        actual.append({'k': k, 'tail': z, 'source_word': left, 'target_word': right,
                       'source_state': vector(left), 'target_state': vector(right),
                       'source_H': h, 'target_H': target, 'offset': offset,
                       'source_order': len(left_rows), 'target_order': len(right_rows),
                       'source_rows': left_rows, 'target_rows': right_rows})

    # q5's determinant is divisible by3: it really can change the bulk level.
    q5_bulk = []
    for k in range(3):
        w = u + ('q5',)*k
        b = row(vector(w), BASIS['q5'])
        q5_bulk.append({'word': w, 'bulk': b, 'valuation3': valuation(b)})
    need([x['valuation3'] for x in q5_bulk] == [0, 1, 1], 'singular-channel positive control')

    seeds = []
    for w in [('r5','a','a'), ('a','a','q1','a','a'), ('q0','a','a','a','a','a')]:
        b = row(vector(w), BASIS['bq2'])
        seeds.append({'word': w, 'state': vector(w), 'bulk': b, 'valuation3': valuation(b)})
    need([s['valuation3'] for s in seeds] == [0, 1, 1], 'inherited seed obstruction')

    result = {
        'status': 'PASS', 'base_commit': 'ff94a6d741', 'atlas_sha256': ATLAS_SHA256,
        'scope': '65 supplied matrices; five inherited error lines; all 4225 length-two words; '
                 'both ends of each discovery lift have H<=500; three literal switched-family controls. '
                 'No all-word signed-charge or global-spectrum classification.',
        'quotient_records': quotient_records, 'one_letter_bridges': bridges,
        'two_letter_signed_units': two_letter_units,
        'small_lift_discovery': {'cap': cap, 'state_count': len(states), 'lifts': lifts},
        'actual_switched_family_controls': actual, 'singular_q5_bulk_control': q5_bulk,
        'inherited_bq2_seed_levels': seeds,
    }
    output = HERE / 'return_bank.json'
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'one_letter_bridges': len(bridges),
                      'two_letter_words': 65**2, 'signed_unit_hits': len(two_letter_units),
                      'small_state_count': len(states), 'literal_pairs': len(actual),
                      'output': output.name}))


if __name__ == '__main__':
    main()
