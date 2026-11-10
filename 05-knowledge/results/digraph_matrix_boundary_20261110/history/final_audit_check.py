"""Independent small cross-packet controls; imports no producer code."""
from collections import Counter
from itertools import permutations, product
from pathlib import Path
from hashlib import sha256
import json

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def rows_from_matrix(a):
    return tuple(sum(int(v) << j for j, v in enumerate(row)) for row in a)


def wordset(rows):
    return {p for p in permutations(range(len(rows)))
            if all(rows[u] >> v & 1 for u, v in zip(p, p[1:]))}


def canon(rows):
    n = len(rows)
    return min(tuple(sum(((rows[p[i]] >> p[j]) & 1) << j for j in range(n))
                     for i in range(n)) for p in permutations(range(n)))


def family_degree(rows):
    counts = Counter((rows[i] >> j & 1) + (rows[j] >> i & 1)
                     for i in range(len(rows)) for j in range(i))
    return counts[0], counts[2], counts[1]


def main():
    enum = json.loads((PACKET / 'enumeration/normal.json').read_text())
    comp = json.loads((PACKET / 'completion/check.json').read_text())
    comparisons = []
    for n in (1, 2, 3):
        arcs = [(i, j) for i in range(n) for j in range(n) if i != j]
        classes = set()
        for bits in product((0, 1), repeat=len(arcs)):
            rows = [0] * n
            for (u, v), present in zip(arcs, bits):
                rows[u] |= present << v
            classes.add(canon(rows))
        actual = Counter(family_degree(rows) for rows in classes)
        expected = {(x['missing'], x['digons'], x['single']): x['types']
                    for x in enum['weighted_orbit_coefficients'][str(n)]}
        need(dict(actual) == expected, 'Independent weighted raw class counts')
        comparisons.append(dict(n=n, labelled_inputs=1 << len(arcs), types=len(classes)))

    c4 = (2, 4, 8, 1)
    loopless, looped = set(), set()
    for a, b in product((0, 1), repeat=2):
        rows = list(c4)
        for u, v, bit in ((0, 2, a), (1, 3, b)):
            rows[u] |= bit << v
            rows[v] |= bit << u
        loopless.add(canon(rows))
        for loops in product((0, 1), repeat=4):
            looped.add(canon([row | (bit << i) for i, (row, bit)
                              in enumerate(zip(rows, loops))]))
    need(len(loopless) == 3 and len(looped) == 22, 'C4 fibre counts')
    need(sorted(len(wordset(r)) for r in loopless) == [4, 4, 8], 'C4 ordinary counts')

    collisions = []
    for record in comp['one_digon_collision']:
        rows = tuple(record['rows'])
        raw_words = wordset(rows)
        resolutions = []
        types = []
        for u, v in ((1, 0), (0, 1)):
            altered = list(rows)
            altered[u] &= ~(1 << v)
            resolutions.append(wordset(altered))
            types.append(canon(altered))
        need(len(raw_words) == record['H'], 'Collision raw H')
        need(resolutions[0] | resolutions[1] == raw_words, 'Collision actual union')
        need([len(s) for s in resolutions] == [5, 5], 'Collision scalar responses')
        collisions.append(dict(rows=rows, raw_H=len(raw_words),
                               common_H=len(resolutions[0] & resolutions[1]), types=types))
    need([c['raw_H'] for c in collisions] == [6, 8], 'Different actual unions')
    need([c['common_H'] for c in collisions] == [4, 2], 'Different actual intersections')
    need(len({tuple(t) for c in collisions for t in c['types']}) == 1,
         'Same unmarked tournament type in every resolution')

    # Independent literal order-eight verification of the selected H7 response.
    target = comp['constrained_target_witnesses'][0]
    tournament = tuple(target['tournament_rows'])
    forced = set(map(tuple, target['matching_arcs']))
    paths = wordset(tournament)
    constrained = [p for p in paths if forced <= set(zip(p, p[1:]))]
    need(len(constrained) == 7, 'Literal H7 constrained response')
    need(len(paths) == target['full_tournament_H'], 'Literal full order-eight H')
    need(len(wordset(tuple(target['input_rows']))) == 7, 'Literal compressed input H')

    files = [
        PACKET / 'enumeration/NOTE.md', PACKET / 'enumeration/check.py',
        PACKET / 'enumeration/normal.json', PACKET / 'enumeration/optimized.json',
        PACKET / 'completion/PROOFS.md', PACKET / 'completion/check.py', PACKET / 'completion/check.json',
        PACKET / 'completion/check.optimized.json',
        PACKET / 'construction/PROOF.md', PACKET / 'construction/check.py',
        Path(__file__).resolve(),
    ]
    need((PACKET / 'enumeration/normal.json').read_bytes() ==
         (PACKET / 'enumeration/optimized.json').read_bytes(), 'Enumeration saved normal/O equality')
    need((PACKET / 'completion/check.json').read_bytes() ==
         (PACKET / 'completion/check.optimized.json').read_bytes(), 'Completion saved normal/O equality')
    print(json.dumps(dict(status='PASS', raw_class_controls=comparisons,
        c4_fibres=dict(labelled_loopless=4, labelled_looped=64,
                      types_loopless=len(loopless), types_looped=len(looped)),
        collisions=collisions,
        selected_tournament=dict(order=len(tournament), full_H=len(paths), constrained_H=len(constrained)),
        scope='Small independent controls and saved normal/O equality, not a fresh full producer rerun',
        source_sha256={str(p.relative_to(PACKET)): sha256(p.read_bytes()).hexdigest() for p in files}),
        indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
