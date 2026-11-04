#!/usr/bin/env python3
"""Independent literal-permutation audit of the component inverse supplier.

No graph is quotiented before its seven observables are independently counted.
The supplier is imported only as the implementation under test.  All gates
remain active with -O.  Outputs concern the declared finite universes, not an
unbounded coverage or physical-return theorem.
"""
from argparse import ArgumentParser
from functools import lru_cache
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json

import component_inverse as subject

GATES = 0


def check(value, message):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(message)


def graphs(n):
    edges = tuple(combinations(range(n), 2))
    for mask in range(1 << len(edges)):
        rows = [0] * n
        for k, (u, v) in enumerate(edges):
            if (mask >> k) & 1:
                rows[u] |= 1 << v
            else:
                rows[v] |= 1 << u
        yield tuple(rows)


def strong(rows):
    for source in range(len(rows)):
        todo, seen = [source], {source}
        while todo:
            u = todo.pop()
            for v in range(len(rows)):
                if (rows[u] >> v) & 1 and v not in seen:
                    seen.add(v)
                    todo.append(v)
        if len(seen) != len(rows):
            return False
    return True


def canonical(rows, mask):
    n = len(rows)
    return min((tuple(sum(((rows[p[i]] >> p[j]) & 1) << j
                          for j in range(n)) for i in range(n)),
                sum(((mask >> p[i]) & 1) << i for i in range(n)))
               for p in permutations(range(n)))


@lru_cache(None)
def incidence(rows):
    """One full permutation and a cut uniquely specify ordered path slots."""
    n = len(rows)
    h = p = 0
    ff, gg, ee, ss = ([0] * n for _ in range(4))
    jj = [[0] * n for _ in range(n)]
    for perm in permutations(range(n)):
        valid = [bool((rows[perm[i]] >> perm[i + 1]) & 1)
                 for i in range(n - 1)]
        if all(valid):
            h += 1
            ee[perm[-1]] += 1
            ss[perm[0]] += 1
        for k in range(n + 1):
            # The omitted link is exactly the link across the slot boundary.
            if not all(valid[:max(0, k - 1)]) or not all(valid[k:]):
                continue
            p += 1
            if k:
                ff[perm[k - 1]] += 1
            if k < n:
                gg[perm[k]] += 1
            if 0 < k < n:
                jj[perm[k - 1]][perm[k]] += 1
    return h, p, ff, gg, jj, ee, ss


def signature(rows, ports):
    h, p, ff, gg, jj, ee, ss = incidence(rows)
    vs = [i for i in range(len(rows)) if (ports >> i) & 1]
    return (h, p, sum(ff[i] for i in vs), sum(gg[i] for i in vs),
            sum(jj[i][j] for i in vs for j in vs),
            sum(ee[i] for i in vs), sum(ss[i] for i in vs))


def ordered_graph(a, b):
    n = len(a.rows)
    right = ((1 << len(b.rows)) - 1) << n
    return (tuple(row | right for row in a.rows) +
            tuple(row << n for row in b.rows), a.ports | (b.ports << n))


def closure_census():
    """Separate reconstruction using the analytically proved first-step iff.

    This is not a replay of the producer's four-million extension comparisons.
    The literal path counts and terminal reachability are independently formed.
    """
    from collections import Counter
    first = Counter()
    digits = {side: Counter() for side in ('append', 'prepend')}
    states = guarded = small = small_guarded = 0
    for n in range(1, 6):
        for rows in graphs(n):
            reach = []
            for source in range(n):
                todo, seen = [source], {source}
                while todo:
                    u = todo.pop()
                    for v in range(n):
                        if (rows[u] >> v) & 1 and v not in seen:
                            seen.add(v)
                            todo.append(v)
                reach.append(seen)
            terminal = set.intersection(*reach)
            for ports in range(1 << n):
                states += 1
                small += n <= 3
                if not any((ports >> v) & 1 for v in terminal):
                    continue
                guarded += 1
                small_guarded += n <= 3
                h, p, f, g, j, e, s = signature(rows, ports)
                if p == 2 * f:
                    first['append'] += 1
                    digits['append'][2 * j - 2 * h - e - g] += 1
                if f == h:
                    first['prepend'] += 1
                    digits['prepend'][-3 * h - s] += 1
    check((states, guarded, small, small_guarded) ==
          (33866, 27475, 74, 43), 'closure raw state universes')
    check(dict(first) == {'append': 775, 'prepend': 311},
          'closure first-step classification census')
    check(dict(digits['append']) ==
          {-42: 120, -20: 120, -16: 144, -12: 88, -10: 150, -4: 153},
          'closure append digit distribution')
    return {'actual_states': states, 'terminal_guard_states': guarded,
            'small_blocks': small, 'small_terminal_blocks': small_guarded,
            'first_steps': dict(first),
            'digits': {k: dict(sorted(v.items())) for k, v in digits.items()},
            'raw_extension_comparisons_implied': 2 * guarded * small,
            'second_trials_implied': sum(first.values()) * (small + small_guarded)}


def main():
    parser = ArgumentParser()
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    base, bank = subject.make_bank(), subject.make_bank_five()
    counts, marked, shapes = {}, set(), {}
    five = set()
    raw_marked = 0
    for n in range(1, 6):
        count = 0
        shape_set = set()
        for rows in graphs(n):
            if not strong(rows):
                continue
            count += 1
            if n <= 4:
                shape_set.add(canonical(rows, 0)[0])
                for ports in range(1 << n):
                    marked.add(canonical(rows, ports))
            else:
                for ports in range(32):
                    value = signature(rows, ports)
                    check(value == subject.direct(rows, ports),
                          'raw five-vertex literal signature')
                    five.add(value)
                    raw_marked += 1
        counts[n] = count
        if n <= 4:
            shapes[n] = len(shape_set)
    check(counts == {1: 1, 2: 0, 3: 2, 4: 24, 5: 544},
          'complete strong labelled graph counts')
    check(shapes == {1: 1, 2: 0, 3: 1, 4: 1},
          'unique strong shapes through four')
    check(marked == {canonical(a.rows, a.ports) for a in base},
          'all marked isomorphism classes through four')
    check(len(marked) == len(base) == 22, '22 distinct marked classes')
    for a in bank:
        check(strong(a.rows), 'representative strongness')
        check(signature(a.rows, a.ports) == a.sigma,
              'representative literal signature')
    supplied = {a.sigma for a in bank}
    check(supplied == {a.sigma for a in base} | five,
          'five-vertex signature quotient completeness')
    check(len(five) == 148 and len(supplied) == len(bank) == 170,
          'five-vertex quotient cardinality')
    inv = subject.Inverse(base)
    check(inv.C == Fraction(32, 5), 'four-vertex terminal bound')
    check(subject.Inverse(bank).C == Fraction(28, 3),
          'five-vertex terminal bound')
    transfers = 0
    for a, b in product(base, repeat=2):
        rows, ports = ordered_graph(a, b)
        actual = signature(rows, ports)
        check(actual == subject.product(a.sigma, b.sigma),
              'literal seven-coordinate ordered product')
        if b.ports:
            h, p, f, g, j, e, s = actual
            check(h <= f and p <= inv.C * f and j <= f - e,
                  'actual fixed-f and fixed-target bounds')
        transfers += 1
    for f in range(1, 101):
        check(inv.words(f) == inv.words_right(f), 'inverse directions')
    levels = []
    for level in range(5):
        f = 7 * 3 ** level
        words = inv.words(f)
        check(words == inv.words_right(f), 'target inverse directions')
        intercepts = {inv.sigma(w)[0] - inv.sigma(w)[4] for w in words}
        accepted = [digits for digits, a in subject.targets(level)
                    if a in intercepts]
        levels.append({'level': level, 'f': f, 'words': len(words),
                       'accepted': accepted})
    check([r['words'] for r in levels] == [2, 5, 9, 18, 29],
          'four-vertex all-length counts')
    check([r['accepted'] for r in levels] ==
          [[''], ['2', '4'], ['22', '42'], ['422'], ['4222']],
          'four-vertex clock-uniform target acceptance')
    # The old h-divides-f restriction must not leak into the enlarged bank.
    q = next(a for a in base if a.name == 'Q1')
    check(q.sigma[0] == 5 and q.sigma[2] == 7,
          'nonintegral old-radix hostile')
    # An unmarked prefix has f=0 and is valid only inside the right inverse.
    t0 = next(k for k, a in enumerate(base) if a.name == 'T0')
    t1 = next(k for k, a in enumerate(base) if a.name == 'T1')
    check((t0, t1) in inv.words_right(2), 'zero-prefix right-inverse hostile')
    closure = closure_census()
    supplier = Path(subject.__file__)
    record = {
        'status': 'FINITE-EXACT independent audit controls PASS',
        'scope': 'Declared finite graph universes and implementation controls; analytic unbounded inverse audited separately.',
        'supplier_sha256': hashlib.sha256(supplier.read_bytes()).hexdigest(),
        'gates': GATES, 'strong_labelled_counts': counts,
        'strong_shape_counts_through_four': shapes,
        'marked_classes_through_four': len(marked),
        'raw_marked_five_cases': raw_marked,
        'five_signatures': len(five), 'combined_signatures': len(supplied),
        'signature_set_sha256': hashlib.sha256(json.dumps(
            sorted(supplied), separators=(',', ':')).encode()).hexdigest(),
        'literal_two_atom_transfers': transfers,
        'dense_inverse_queries': 100, 'levels': levels,
        'C_four': str(inv.C), 'C_five': str(subject.Inverse(bank).C),
        'closure_independent_census': closure,
    }
    content = json.dumps(record, indent=2, sort_keys=True) + '\n'
    args.output.write_text(content)
    print(json.dumps({'gates': GATES, 'receipt_sha256': hashlib.sha256(
        content.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
