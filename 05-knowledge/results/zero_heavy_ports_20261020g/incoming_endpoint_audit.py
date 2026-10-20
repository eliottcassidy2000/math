"""Independent finite-head generation and literal endpoint-two certificate replay.

Imports the incoming word-template banks, never their evaluators or search mains.
The n6..8 universe is independently generated on all n vertices, not by the
incoming two-port construction. The n9..10 control uses independently generated
blue middles and a local endpoint construction. Every K move, directedness test,
positional anti-involution test and literal target comparison is local here.

No files are written. Explicit exceptions survive python -O. This checks the
finite head and bounded proof controls, not the all-order symbolic proofs.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json
import runpy
import sys


HERE = Path(__file__).resolve().parent
INCOMING = HERE.parent / 'endpoint_two_repairs_20261020g'


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def blue_presentations(n):
    """All full-n tournaments whose increasing path is blue and directed."""
    remaining = set(combinations(range(n), 2))
    orbits = []
    while remaining:
        a, b = min(remaining)
        orbit = {(a, b), (n - 1 - b, n - 1 - a)}
        remaining -= orbit
        orbits.append(orbit)
    fixed = [o for o in orbits if any(b == a + 1 for a, b in o)]
    free = [o for o in orbits if o not in fixed]
    for bits in product((0, 1), repeat=len(free)):
        arcs = set().union(*fixed)
        for orbit, bit in zip(free, bits):
            arcs.update(orbit if bit else {(b, a) for a, b in orbit})
        yield frozenset(arcs)


def move(arcs, word):
    keep = {frozenset((a, b)) for a, b in zip(word, word[1:])}
    return frozenset((a, b) if frozenset((a, b)) in keep else (b, a)
                     for a, b in arcs)


def blue(arcs, word):
    reversal = dict(zip(word, word[::-1]))
    return all((reversal[b], reversal[a]) in arcs for a, b in arcs)


def source_controls(arcs, n):
    need(len(arcs) == n * (n - 1) // 2, 'pair count')
    need(all(((a, b) in arcs) != ((b, a) in arcs)
             for a, b in combinations(range(n), 2)), 'tournament pairs')
    need(all((i, i + 1) in arcs for i in range(n - 1)), 'source path')
    need(blue(arcs, tuple(range(n))), 'source blue')
    need((n - 1, 0) in arcs, 'source closed')
    need(sum(a == 0 for a, b in arcs) == n - 3, 'first endpoint score')
    need(sum(a == n - 1 for a, b in arcs) == 2, 'last endpoint score')


def replay(arcs, words, n):
    source = arcs
    need(len(words) == 3, 'exactly three moves')
    for word in words:
        need(len(word) == n and set(word) == set(range(n)), 'permutation')
        if not all(edge in arcs for edge in zip(word, word[1:])):
            return False
        if blue(arcs, word):
            return False
        arcs = move(arcs, word)
    need(arcs == move(source, tuple(range(n))), 'literal labelled target')
    return True


def endpoint_source(n, middle, exceptional):
    """Construct locally from zero-based middle arcs; no incoming generator."""
    z = n - 1
    arcs = {(a + 1, b + 1) for a, b in middle}
    arcs.update((0, i) for i in range(1, z))
    arcs.update((i, z) for i in range(1, z))
    arcs.add((z, 0))
    reflected = n - 1 - exceptional
    arcs.remove((0, exceptional))
    arcs.add((exceptional, 0))
    arcs.remove((reflected, z))
    arcs.add((z, reflected))
    return frozenset(arcs)


def main():
    # run_path does not execute guarded mains; sibling imports provide only the
    # declared candidate words. No discovery search or receipt writer is called.
    sys.path.insert(0, str(INCOMING))
    coverage = runpy.run_path(str(INCOMING / 'probe_coverage.py'))
    unified = runpy.run_path(str(INCOMING / 'unified_check.py'))
    receipt = json.loads((INCOMING / 'probe_coverage.json').read_text())
    for name, expected in sorted(receipt['source_hashes'].items()):
        need(sha256((INCOMING / name).read_bytes()).hexdigest() == expected,
             'incoming finite receipt source hash: ' + name)
    print('Incoming finite receipt source hashes: PASS')
    finite = {tuple(c['rows']): c['words']
              for c in receipt['finite_residual_certificates']}
    used_finite = set()
    for n in (6, 7, 8):
        count = fallback = 0
        covered = set()
        for arcs in blue_presentations(n):
            if (n - 1, 0) not in arcs:
                continue
            if sum(a == 0 for a, b in arcs) != n - 3:
                continue
            count += 1
            source_controls(arcs, n)
            rows = tuple(sum(1 << b for a, b in arcs if a == i)
                         for i in range(n))
            exceptional, = (a for a, b in arcs if b == 0 and 0 < a < n - 1)
            found = False
            for name, words in coverage['candidates'](n, exceptional):
                if replay(arcs, words, n):
                    found = True
                    break
            if not found:
                need(rows in finite, 'missing finite fallback')
                need(replay(arcs, finite[rows], n), 'invalid finite fallback')
                used_finite.add(rows)
                fallback += 1
            covered.add(rows)
        expected = {6: 12, 7: 64, 8: 320}[n]
        need(count == expected and len(covered) == expected,
             'full-n generated finite universe')
        need(fallback == {6: 3, 7: 1, 8: 1}[n], 'finite fallback count')
        print(f'Full-n independent generator/replay: n={n}, '
              f'sources={count}, finite_fallbacks={fallback}')
    need(used_finite == set(finite), 'every stored finite exception replayed')
    bounded_receipt = json.loads((INCOMING / 'unified_check.json').read_text())
    for n in (9, 10):
        counts = Counter()
        coincident = 0
        for middle in blue_presentations(n - 2):
            for exceptional in range(2, n - 1):
                arcs = endpoint_source(n, middle, exceptional)
                source_controls(arcs, n)
                branch, words = unified['choose_ledger'](n, arcs, exceptional)
                need(replay(arcs, words, n), 'unified certificate')
                counts[branch] += 1
                coincident += 2 * exceptional == n - 1
        expected = {9: 3072, 10: 28672}[n]
        need(sum(counts.values()) == expected, 'unified universe')
        need(dict(counts) == bounded_receipt['exhaustive'][str(n)]['counts'],
             'seven-branch receipt counts')
        need(coincident == {9: 512, 10: 0}[n], 'coincident ports included')
        print(f'Independent unified replay: n={n}, sources={expected}, '
              f'coincident={coincident}, branches={len(counts)}')
    print('PASS: no strongness/mixedness filter; no receipt writes; '
          'no discovery search; no all-order inference from finite controls')


if __name__ == '__main__':
    main()
