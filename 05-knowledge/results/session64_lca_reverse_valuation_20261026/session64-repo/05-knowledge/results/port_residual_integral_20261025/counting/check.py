"""Exact integer-lattice controls; all gates survive python -O.

No timing claim. The rational oracle retains f rather than D and has no
moment recurrence. Its intermediate states are compared coordinatewise.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb, factorial, gcd, prod
from pathlib import Path
from random import Random
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def valid_n(n):
    if type(n) is not int or n < 0:
        raise ValueError('target must be a nonnegative integer')


def beta(k, e):
    return (-1) ** e * comb(k + e - 1, e) if k else int(e == 0)


def vfact(k):
    return k - k.bit_count()


@lru_cache(None)
def oddfact(k):
    return factorial(k) >> vfact(k)


def bound(defect):
    need(all(type(x) is int and x >= 0 for x in defect), 'nonnegative integer defects')
    need(not defect or defect[0] == 0, 'D(1)=0')
    return max(((value + 2 * i - 1) // (2 * i)
                for i, value in enumerate(defect) if i), default=0)


def intrinsic(ell, k):
    return k * (ell - 1) // 2 + ell * k * (k - 1) // 2


@lru_cache(None)
def transfer(u, ell, k, e):
    up = u + ell * (k + e)
    divisor = oddfact(u) * ell ** k * oddfact(k)
    quotient, remainder = divmod(oddfact(up), divisor)
    need(remainder == 0, ('odd transfer divisibility', u, ell, k, e))
    return quotient * beta(k, e)


def sparse_vector(s, defect, support):
    if not s:
        return [sum(support.values())]
    b = bound(defect)
    need(all(m >= b for m in support), 'sparse support respects b')
    return [0] + [sum(x << (r * (m - b)) for m, x in support.items())
                  for r in range(1, s + 1)]


def dense_contribution(n, ell, s, defect, data, k, e, counters=None):
    rest = s - ell * (k + e)
    dn = tuple(defect[(d - 1) // 2] + k * (d - gcd(ell, d))
               for d in range(1, min(ell - 2, rest) + 1, 2))
    b, bp = bound(defect), bound(dn)
    selected = defect[(ell - 1) // 2] if ell <= s else 0
    need(bp <= b + k, ('destination lower bound', ell, s, defect, k, e))
    scale = transfer(n - s, ell, k, e)
    base = intrinsic(ell, k) - vfact(k)
    out = [0] * (rest + 1)
    for r in (range(1, rest + 1) if rest else (0,)):
        shift = base + k * (ell * b - selected) + r * (k + b - bp)
        need(shift >= 0, ('nonnegative dense shift', ell, s, defect, k, e, r))
        source_index = ell * k + r
        need(source_index <= s and (source_index > 0 or s == 0), 'affordable source moment')
        need(type(data[source_index]) is int, 'integer source moment')
        out[r] = (scale * data[source_index]) << shift
        if counters is not None:
            counters['dense_channels'] += 1
            counters['negative_ungauged_channels'] += int(base - k * selected + r * k < 0)
    return (rest, dn), out


def rational_stages(n, core):
    """Independent original f-response recurrence, including intermediate stages."""
    if not n:
        return {}, Fraction(1)
    top = n if n % 2 else n - 1
    states = {(n, (0,) * ((top + 1) // 2)): Fraction(1)}
    snapshots = {}
    for ell in range(top, 0, -2):
        snapshots[ell] = states
        nxt = defaultdict(Fraction)
        for (s, response), w in states.items():
            selected = response[(ell - 1) // 2] if ell <= s else 0
            for k in range(s // ell + 1):
                choices = range((s - ell * k) // ell + 1) if core and k else (0,)
                for e in choices:
                    rest = s - ell * (k + e)
                    if ell == 1 and rest:
                        continue
                    fn = tuple(response[(d - 1) // 2] + k * gcd(ell, d)
                               for d in range(1, min(ell - 2, rest) + 1, 2))
                    factor = Fraction(beta(k, e) << (intrinsic(ell, k) + k * selected),
                                      ell ** k * factorial(k))
                    nxt[(rest, fn)] += w * factor
        states = {key: w for key, w in nxt.items() if w}
    snapshots[-1] = states
    return snapshots, states.get((0, ()), Fraction(0))


def audit_state(n, ell, states, reference, core):
    banks = {}
    for (s, response), w in reference.items():
        if s:
            m = response[0]
            defect = tuple(d * m - value for d, value in zip(range(1, min(ell, s) + 1, 2), response))
            bank = banks.setdefault((s, defect), {})
            need(m not in bank, 'f determines m at a D fibre')
            bank[m] = w * oddfact(n - s)
            need(bank[m].denominator == 1, 'integer scaled reference weight')
        else:
            banks[(0, ())] = {0: w * oddfact(n)}
    need(set(states) == set(banks), ('same stagewise defect fibres', n, ell, core))
    for (s, defect), (kind, data) in states.items():
        bank = {m: int(x) for m, x in banks[(s, defect)].items()}
        if kind == 'sparse':
            need(data == bank, ('sparse reference equality', n, ell, s, defect))
        elif kind == 'dense':
            need(data == sparse_vector(s, defect, bank), ('dense reference equality', n, ell, s, defect))
        else:
            expected = sparse_vector(s, defect, bank)[s]
            need(kind == 'ordinary' and data == expected, ('ordinary reference equality', n, ell, s, defect))


def integral_target(n, core=True, reference=None):
    valid_n(n)
    if n == 0:
        return 1, {'keys': 1, 'cells': 1, 'rows': [], 'dense_channels': 0,
                   'negative_ungauged_channels': 0, 'conversions': 0}
    top = n if n % 2 else n - 1
    initial = ('sparse', {0: 1}) if core else ('ordinary', 1)
    states = {(n, (0,) * ((top + 1) // 2)): initial}
    stats = {'dense_channels': 0, 'negative_ungauged_channels': 0, 'conversions': 0}
    rows = []
    for ell in range(top, 0, -2):
        if reference is not None:
            audit_state(n, ell, states, reference[ell], core)
        rows.append({'length': ell, 'keys': len(states),
                     'cells': sum(len(data) if kind == 'sparse' else max(1, key[0]) if kind == 'dense' else 1
                                  for key, (kind, data) in states.items())})
        nxt = {}

        def add_sparse(key, m, x):
            rest, dn = key
            need(type(x) is int and x != 0, 'nonzero integer sparse contribution')
            if rest:
                need((x > 0) == ((n - rest - m) % 2 == 0), 'prefix sign coherence')
            if key not in nxt:
                nxt[key] = ('sparse', {}) if rest else ('dense', [0])
            kind, data = nxt[key]
            if kind == 'sparse':
                if m in data:
                    need((data[m] > 0) == (x > 0), 'same-m terms never cancel')
                data[m] = data.get(m, 0) + x
                if len(data) > rest:
                    nxt[key] = ('dense', sparse_vector(rest, dn, data))
                    stats['conversions'] += 1
            elif rest:
                bp = bound(dn)
                need(m >= bp, 'destination sparse-to-dense bound')
                for r in range(1, rest + 1):
                    data[r] += x << (r * (m - bp))
            else:
                data[0] += x

        def add_dense(key, vector):
            rest, dn = key
            if key not in nxt:
                nxt[key] = ('dense', [0] * (rest + 1))
            kind, data = nxt[key]
            if kind == 'sparse':
                data = sparse_vector(rest, dn, data)
                nxt[key] = ('dense', data)
                stats['conversions'] += 1
            for r in (range(1, rest + 1) if rest else (0,)):
                data[r] += vector[r]

        for (s, defect), (kind, data) in states.items():
            selected = defect[(ell - 1) // 2] if ell <= s else 0
            for k in range(s // ell + 1):
                choices = range((s - ell * k) // ell + 1) if core and k else (0,)
                for e in choices:
                    rest = s - ell * (k + e)
                    if ell == 1 and rest:
                        continue
                    dn = tuple(defect[(d - 1) // 2] + k * (d - gcd(ell, d))
                               for d in range(1, min(ell - 2, rest) + 1, 2))
                    key = (rest, dn)
                    scale = transfer(n - s, ell, k, e)
                    if kind == 'sparse':
                        for m, x in data.items():
                            shift = intrinsic(ell, k) - vfact(k) + k * (ell * m - selected)
                            need(m >= bound(defect) and shift >= 0, 'nonnegative sparse shift')
                            add_sparse(key, m + k, (x * scale) << shift)
                    elif kind == 'dense':
                        checked_key, vector = dense_contribution(n, ell, s, defect, data, k, e, stats)
                        need(checked_key == key, 'dense key agreement')
                        add_dense(key, vector)
                    else:
                        b, bp = bound(defect), bound(dn)
                        shift = intrinsic(ell, k) - vfact(k) + k * (ell * b - selected) + rest * (k + b - bp)
                        need(shift >= 0 and ell * k + rest == s, 'ordinary channel closes integrally')
                        value = (scale * data) << shift
                        prior = nxt.get(key, ('ordinary', 0))[1]
                        nxt[key] = ('ordinary', prior + value)
        states = {key: entry for key, entry in nxt.items()
                  if key[0] or (entry[1][0] if entry[0] == 'dense' else entry[1])}
    if reference is not None:
        audit_state(n, -1, states, reference[-1], core)
    need(set(states) <= {(0, ())}, 'terminal integral key')
    entry = states.get((0, ()), ('ordinary', 0))
    numerator = entry[1][0] if entry[0] == 'dense' else entry[1]
    quotient, remainder = divmod(numerator, oddfact(n))
    need(remainder == 0, ('exact terminal normalization', n, core))
    stats.update({'keys': sum(row['keys'] for row in rows),
                  'cells': sum(row['cells'] for row in rows), 'rows': rows})
    return quotient, stats


def partitions(n, low=1):
    if n == 0:
        yield ()
    for d in range(low, n + 1, 2):
        for tail in partitions(n - d, d):
            yield (d,) + tail


def ordinary_direct(n):
    total = Fraction(0)
    for part in partitions(n):
        multiplicity = Counter(part)
        exponent = sum((d - 1) // 2 for d in part) + sum(gcd(a, b) for a, b in combinations(part, 2))
        denominator = prod(d ** k * factorial(k) for d, k in multiplicity.items())
        total += Fraction(1 << exponent, denominator)
    return total


def lattice_controls():
    rng = Random(20261025)
    tests = 0
    negative = 0
    for ell in (3, 5, 7, 9, 11):
        for _ in range(30):
            s = ell + rng.randrange(7)
            defect = (0,) + tuple(rng.randrange(19) for _ in range((ell - 1) // 2))
            b = bound(defect)
            support = {b + j: rng.choice((-5, -3, -1, 1, 2, 4)) for j in range(5)}
            u = (ell + 2) * (b + 5)
            n = u + s
            source = sparse_vector(s, defect, support)
            arbitrary = [0] + [rng.randrange(-9, 10) for _ in range(s)]
            for k in range(s // ell + 1):
                for e in (range((s - ell * k) // ell + 1) if k else (0,)):
                    key, actual = dense_contribution(n, ell, s, defect, source, k, e)
                    rest, dn = key
                    dest = {}
                    for m, x in support.items():
                        shift = intrinsic(ell, k) - vfact(k) + k * (ell * m - defect[-1])
                        need(shift >= 0, 'relaxed-lattice sparse shift')
                        dest[m + k] = (transfer(u, ell, k, e) * x) << shift
                    need(actual == sparse_vector(rest, dn, dest), 'sparse/dense transfer square commutes')
                    _, arbitrary_result = dense_contribution(n, ell, s, defect, arbitrary, k, e)
                    need(all(type(x) is int for x in arbitrary_result), 'full product lattice closure')
                    negative += sum(intrinsic(ell, k) - vfact(k) - k * defect[-1] + r * k < 0
                                    for r in (range(1, rest + 1) if rest else (0,)))
                    tests += 1
    need(negative > 0, 'test universe includes negative ungauged exponents')

    # Genuine one-prefix hostile: exact half-moment, versus a malformed X_3=1.
    x = Fraction(4, 5) * oddfact(5)
    need(x == 12 and transfer(5, 3, 1, 0) == 7, '5-cycle normalization')
    source = sparse_vector(3, (0, 2), {1: int(x)})
    _, out = dense_contribution(8, 3, 3, (0, 2), source, 1, 0)
    need(out == [336], 'genuine negative-shift prefix repaired')
    need(Fraction(7, 2).denominator == 2, 'ungauged product lattice is not invariant')

    # Four actual prefix histories, same committed degree41 and D(3)=2.
    # The 9-cycle excess is 4-k; the final 5-cycle has no excess.
    actual_support = {}
    for k in range(1, 5):
        exponent = 9 * k * (k - 1) // 2 + 5 * k + 2
        w = Fraction(beta(k, 4 - k) << exponent, 9 ** k * factorial(k) * 5)
        scaled = w * oddfact(41)
        need(scaled.denominator == 1, 'genuine dense-prefix integrality')
        actual_support[k + 1] = int(scaled)
        need((scaled > 0) == ((41 - (k + 1)) % 2 == 0), 'genuine dense-prefix sign')
    need(len(actual_support) > 3, 'genuine adaptive dense fibre')
    source = sparse_vector(3, (0, 2), actual_support)
    _, out = dense_contribution(44, 3, 3, (0, 2), source, 1, 0)
    unshifted = sum(x << (3 * m) for m, x in actual_support.items())
    old = Fraction(transfer(41, 3, 1, 0) * unshifted, 2)
    need(old.denominator == 1 and out[0] == old, 'genuine compressed negative-shift repair')

    # Old horizon hostile: integrality does not recover missing moments.
    hidden = sparse_vector(3, (0,), {1: 8, 2: -6, 3: 1})
    need(hidden == [0, 0, 0, 192], 'horizon kernel survives integer gauge')
    return {'relaxed_lattice_transfers': tests, 'negative_ungauged_channels': negative,
            'genuine_dense_hostile': {'committed_degree': 41, 'remaining': 3,
                                      'cycle_counts': list(actual_support), 'b': 1,
                                      'ungauged_exponent': -1, 'gauged_exponent': 2},
            'genuine_dense_old_new_equal': out[0] == old, 'one_prefix_value': 336}


def main():
    old_path = Path(__file__).resolve().parents[2] / 'unprotected_ports_cancellation_20261025' / 'counting' / 'check.py'
    spec = spec_from_file_location('prior_rational_moment_controls', old_path)
    prior = module_from_spec(spec)
    spec.loader.exec_module(prior)
    samples = []
    stage_states = 0
    core_prefix = [1, 1, 0, 1, 2, 5, 26, 255, 4602, 147095, 8228360]
    ordinary_prefix = [1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056]
    for n in range(27):
        for core in (False, True):
            reference, rational = rational_stages(n, core)
            value, stats = integral_target(n, core, reference)
            need(value == rational, ('rational/integral totals', n, core))
            if core:
                baseline, old_stats = prior.adaptive_target(n)
                need(value == baseline, ('prior adaptive total', n))
                need(stats['cells'] == old_stats['scalar_cells'], ('adaptive cell graph unchanged', n))
                if n <= 22:
                    need(value == prior.direct(n), ('unmerged decorated partition oracle', n))
            else:
                need(value == ordinary_direct(n), ('unmerged ordinary partition oracle', n))
            if n < len(core_prefix):
                need(value == (core_prefix if core else ordinary_prefix)[n], ('known count prefix', n, core))
            stage_states += sum(row['keys'] for row in stats['rows'])
            if n in (10, 18, 26):
                samples.append({'n': n, 'mode': 'core' if core else 'ordinary', 'value': value,
                                'keys': stats['keys'], 'scalar_cells': stats['cells'],
                                'conversions': stats['conversions'],
                                'negative_ungauged_channels': stats['negative_ungauged_channels']})
    reference44, rational44 = rational_stages(44, True)
    value44, stats44 = integral_target(44, True, reference44)
    need(value44 == rational44 and stats44['negative_ungauged_channels'] > 0,
         'full order44 run exercises negative ungauged dense channels')
    for bad in (-1, True, 1.5, '3'):
        try:
            integral_target(bad)
        except ValueError:
            pass
        else:
            raise RuntimeError('invalid target accepted')
    lattice = lattice_controls()
    print(json.dumps({'status': 'PASS', 'ordinary_and_core_range': [0, 26],
                      'core_unmerged_oracle_range': [0, 22],
                      'intermediate_state_audits': stage_states,
                      'order44_control': {'same_intermediate_states': stats44['keys'],
                                          'scalar_cells': stats44['cells'],
                                          'negative_ungauged_channels': stats44['negative_ungauged_channels']},
                      'exact_transfer_cache_entries': transfer.cache_info().currsize,
                      'lattice_controls': lattice, 'samples': samples,
                      'scope': 'integer operations and unchanged adaptive scalar-cell counts; no runtime or positivity claim'}, indent=2))


if __name__ == '__main__':
    main()
