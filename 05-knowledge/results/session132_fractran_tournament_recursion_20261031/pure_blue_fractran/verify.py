#!/usr/bin/env python3
"""Pure-blue 22-fraction stream: literal FRACTRAN and independent finite tests.
Standard library only. No external solver, tournament data, or macro execution.
"""
from collections import Counter
from itertools import permutations, product
from math import comb, factorial, gcd, prod
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
REG = (2, 3, 5, 7, 11, 13)
TEMP, INDEX, CHECKPOINT = 17, 19, 23
MAIN = (29, 37, 43, 53, 61, 71, 79)
AUX = (31, 41, 47, 59, 67, 73, 83)
STATES = (CHECKPOINT,) + MAIN + AUX
PRIMES = REG + (TEMP, INDEX) + STATES
SEED = CHECKPOINT * prod(p ** a for p, a in zip(REG, (1, 1, 1, 2, 1, 3)))
# (source register, destination registers). Source is never a destination.
TRANSFERS = ((2, (17,)), (3, (2,)), (5, (3,)), (7, (5,)),
             (11, (7, 17)), (13, (11,)), (17, (13,)))


def require(condition, message="verification condition failed"):
    """Hard check retained under python -O; never silently skip validation."""
    if not condition:
        raise AssertionError(message)


def compile_program():
    fractions = []
    for i, (source, dests) in enumerate(TRANSFERS):
        following = MAIN[i + 1] if i + 1 < len(MAIN) else CHECKPOINT
        fractions.extend(((AUX[i] * prod(dests), MAIN[i] * source),
                          (MAIN[i], AUX[i]), (following, MAIN[i])))
    fractions.append((MAIN[0] * INDEX, CHECKPOINT))
    require(all(gcd(a, b) == 1 for a, b in fractions))
    return fractions


def load_program():
    return [tuple(map(int, line.split('#')[0].strip().split('/')))
            for line in (HERE / 'pure_blue.fractran').read_text().splitlines()
            if line.split('#')[0].strip()]


def step(integer, fractions):
    """First fraction yielding an integer; literal rational-step semantics."""
    for i, (numerator, denominator) in enumerate(fractions):
        product_ = integer * numerator
        if product_ % denominator == 0:
            return product_ // denominator, i
    raise AssertionError('unexpected halt')


def valuation(integer, prime):
    exponent = 0
    while integer % prime == 0:
        exponent += 1
        integer //= prime
    return exponent


def factor_known(integer):
    result = {}
    for prime in PRIMES:
        value = valuation(integer, prime)
        if value:
            result[prime] = value
            integer //= prime ** value
    require(integer == 1, ('unexpected prime factor', integer))
    return result


def narayana(n):
    # Composition count by number k of weight-three letters; no recurrence.
    return sum(comb(n - 2 * k, k) for k in range(n // 3 + 1)) if n >= 0 else 0


def formula(n):
    return sum(narayana((n - c) // 2) for c in (0, 1, 3, 5)
               if n >= c and (n - c) % 2 == 0)


def words(weight):
    """All ordered {1,3} words, without using a count recurrence."""
    if weight == 0:
        yield ()
    for size in (1, 3):
        if size <= weight:
            for suffix in words(weight - size):
                yield (size,) + suffix


def normal_forms(n):
    for center in (0, 1, 3, 5):
        if center <= n and (n - center) % 2 == 0:
            for half in words((n - center) // 2):
                yield half + ((center,) if center else ()) + half[::-1]


def expected_budget(values):
    a, b, c, d, e, f = values
    return 8 + 2 * (2 * a + b + c + d + 2 * e + f)


def checkpoint_record(integer, n, steps):
    values = [formula(n + j) for j in range(6)]
    expected = CHECKPOINT * INDEX ** n * prod(p ** a for p, a in zip(REG, values))
    require(integer == expected)
    require(factor_known(integer) == {p: e for p, e in
        [(CHECKPOINT, 1), (INDEX, n)] + list(zip(REG, values)) if e})
    return {'n': n, 'output': valuation(integer, 2), 'window': values,
            'rational_steps_from_seed': steps, 'integer_bits': integer.bit_length()}


def stream_test(fractions, count=25):
    integer, steps = SEED, 0
    rows = [checkpoint_record(integer, 0, steps)]
    trace = [{'step': 0, 'integer': str(integer), 'fraction_index_1based': None}]
    fraction_usage = Counter()
    # This is a literal integer FRACTRAN run; recurrence values are assertions only.
    while len(rows) < count:
        integer, selected = step(integer, fractions)
        steps += 1
        fraction_usage[selected + 1] += 1
        if steps <= 30:
            trace.append({'step': steps, 'integer': str(integer),
                          'fraction_index_1based': selected + 1})
        if integer % CHECKPOINT == 0:
            n = len(rows)
            require(steps - rows[-1]['rational_steps_from_seed'] == expected_budget(rows[-1]['window']))
            rows.append(checkpoint_record(integer, n, steps))
    require(set(fraction_usage) == set(range(1, len(fractions) + 1)))
    return rows, trace, dict(sorted(fraction_usage.items()))


def macro_test(fractions):
    # Exhaust all 3^6 windows, including zero branches and duplicated source.
    max_steps = 0
    for values in product(range(3), repeat=6):
        integer = CHECKPOINT * prod(p ** a for p, a in zip(REG, values))
        steps = 0
        while True:
            integer, _ = step(integer, fractions)
            steps += 1
            if integer % CHECKPOINT == 0:
                break
            require(steps <= expected_budget(values))
        target = values[1:] + (values[0] + values[4],)
        require(integer == CHECKPOINT * INDEX * prod(p ** a for p, a in zip(REG, target)))
        require(steps == expected_budget(values))
        max_steps = max(max_steps, steps)
    return {'windows': 3 ** 6, 'value_range': [0, 2], 'max_steps': max_steps}


def adjacency(n, code):
    rows = [0] * n
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code >> bit & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
            bit += 1
    return rows


def all_paths(rows):
    n = len(rows)
    all_mask = (1 << n) - 1
    def extend(path, unused):
        if not unused:
            yield path
            return
        options = unused if not path else unused & rows[path[-1]]
        while options:
            bit = options & -options
            options -= bit
            vertex = bit.bit_length() - 1
            yield from extend(path + (vertex,), unused ^ bit)
    yield from extend((), all_mask)


def anti_reversal(rows, path):
    # Arc-level positional reversal, deliberately not the rejected affine test.
    reverse = dict(zip(path, path[::-1]))
    for u in range(len(rows)):
        for v in range(len(rows)):
            if rows[u] >> v & 1:
                if not (rows[reverse[v]] >> reverse[u] & 1):
                    return False
    return True


def is_pure(rows):
    found = False
    for path in all_paths(rows):
        found = True
        if not anti_reversal(rows, path):
            return False
    require(found)  # includes the empty word for the formal n=0 object
    return True


def relabel_code(rows, permutation):
    result, bit = 0, 0
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            result |= ((rows[permutation[i]] >> permutation[j]) & 1) << bit
            bit += 1
    return result


def canonical(rows):
    return min(relabel_code(rows, perm) for perm in permutations(range(len(rows))))


def stack(sizes):
    rows, offset = [0] * sum(sizes), 0
    for size in sizes:
        require(size in (1, 3, 5))
        for u in range(size):
            for v in range(size):
                if (v - u) % size in range(1, size // 2 + 1):
                    rows[offset + u] |= 1 << (offset + v)
            for v in range(offset + size, len(rows)):
                rows[offset + u] |= 1 << v
        offset += size
    return rows


def census(max_n=6):
    output = []
    for n in range(max_n + 1):
        keys, labelled = set(), 0
        for code in range(1 << (n * (n - 1) // 2)):
            rows = adjacency(n, code)
            if is_pure(rows):
                labelled += 1
                keys.add(canonical(rows))
        expected_keys = {canonical(stack(sizes)) for sizes in normal_forms(n)}
        require(keys == expected_keys)
        for key in keys:
            rows = adjacency(n, key)
            converse = [sum(1 << u for u in range(n) if rows[u] >> v & 1)
                        for v in range(n)]
            require(canonical(converse) == key)
        output.append({'n': n, 'labelled_universe': 1 << (n * (n - 1) // 2),
                       'labelled_pure_blue': labelled, 'isomorphism_classes': len(keys),
                       'converse_merged_classes': len(keys), 'canonical_codes': sorted(keys)})
        require(len(keys) == formula(n))
    return output


def path_inventory(sizes):
    rows = stack(sizes)
    paths = list(all_paths(rows))
    aut = sum(relabel_code(rows, p) == relabel_code(rows, tuple(range(len(rows))))
              for p in permutations(range(len(rows)))) if len(rows) <= 6 else prod(
                  {1: 1, 3: 3, 5: 5}[s] for s in sizes)
    return {'strong_component_sizes': list(sizes), 'H': len(paths),
            'palindromic_H': sum(anti_reversal(rows, p) for p in paths),
            'Aut': aut, 'tc': len(paths) // aut,
            'blue_tilings': sum(anti_reversal(rows, p) for p in paths) // aut}


def audit_affine_hostile():
    # The inherited odd-proof gap must remain a failing affine shortcut control.
    n, connection = 9, {1, 3, 4, 7}
    rows = [sum(1 << v for v in range(n) if (v-u) % n in connection) for u in range(n)]
    path = (0, 1, 2, 3, 7, 5, 6, 4, 8)
    require(all(rows[u] >> v & 1 for u, v in zip(path, path[1:])))
    sums = sorted({(u+v) % n for u, v in zip(path, path[::-1])})
    require(sums == [5, 8] and anti_reversal(rows, path))
    return {'order': n, 'connection_set': sorted(connection), 'path': list(path),
            'partner_sums': sums, 'actual_arc_reversal_passes': True,
            'affine_only_nonpalindromicity_inference_is_invalid': True}


def main():
    fractions = load_program()
    require(fractions == compile_program())
    require(len(fractions) == 22)
    rows, trace, usage = stream_test(fractions)
    for n in range(25):
        forms = list(normal_forms(n))
        require(len(forms) == len(set(forms)) == formula(n))
        require(formula(2*n) == narayana(n))
        require(formula(2*n + 1) == narayana(n + 2))
    controls = [path_inventory(x) for x in [(1,), (3,), (5,), (1,3,1), (3,3), (5,5)]]
    require(controls[-1]['tc'] == 9 and controls[-1]['blue_tilings'] == 3)
    report = {'status': 'ALL CHECKS PASSED', 'source_commit': 'fae72243028c39ce766cd9c59df763c9afd883e2',
              'program_sha256': hashlib.sha256((HERE/'pure_blue.fractran').read_bytes()).hexdigest(),
              'fraction_count': len(fractions), 'seed': SEED,
              'first_25_outputs_n_0_through_24': [row['output'] for row in rows],
              'checkpoints': rows, 'fraction_usage': usage, 'macro_exhaustive': macro_test(fractions),
              'direct_labelled_census': census(), 'controls': controls,
              'inherited_affine_shortcut_hostile': audit_affine_hostile()}
    (HERE/'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    (HERE/'first_transition_trace.json').write_text(json.dumps(trace, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ('status','fraction_count','seed','first_25_outputs_n_0_through_24','macro_exhaustive','direct_labelled_census','controls')}, indent=2))


if __name__ == '__main__':
    main()
