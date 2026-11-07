#!/usr/bin/env python3
"""Independent finite controls for the leaf-spectrum deductions. Stdlib only.

Checks raise explicitly and remain active under python -O. The all-order
classification and analytic limits are proved/inherited, not inferred here.
"""
from collections import Counter
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import permutations
from math import comb, factorial, log, sqrt
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def pal(s, k):
    return 0 if s % 2 and k % 2 else comb(s // 2 + k // 2, k // 2)


def spectrum(n):
    out = Counter()
    for k in range(n // 3 + 1):
        s = n - 3 * k
        number = comb(s + k, k) - pal(s, k)
        require(number % 2 == 0, 'nonintegral reversal orbit count')
        if number:
            out[3 ** k] += number // 2
    if n >= 4 and n % 2 == 0:
        for k in range((n - 4) // 6 + 1):
            s = (n - 4) // 2 - 3 * k
            out[5 * 9 ** k] += comb(s + k, k)
    return out


def words(n):
    if n == 0:
        yield ()
    for a in (1, 3):
        if n >= a:
            for w in words(n - a):
                yield (a,) + w


def word_spectrum(n):
    out = Counter()
    for w in words(n):
        if w < w[::-1]:
            out[3 ** w.count(3)] += 1
    if n >= 4 and n % 2 == 0:
        for w in words((n - 4) // 2):
            out[5 * 9 ** w.count(3)] += 1
    return out


def v3(h):
    k = 0
    while h % 3 == 0:
        h //= 3
        k += 1
    return k


def graph(n, code):
    out = [0] * n
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            u, v = (i, j) if code >> bit & 1 else (j, i)
            out[u] |= 1 << v
            bit += 1
    return out


def edge(out, u, v):
    return bool(out[u] >> v & 1)


def transform(out, p):
    n = len(out)
    inv = [0] * n
    for u, v in enumerate(p):
        inv[v] = u
    code = bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            code |= int(edge(out, inv[i], inv[j])) << bit
            bit += 1
    return code


def paths(out):
    return [p for p in permutations(range(len(out)))
            if all(edge(out, u, v) for u, v in zip(p, p[1:]))]


def h_dp(out):
    n = len(out)
    if not n:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v]:
                rest = out[v] & ~mask
                while rest:
                    low = rest & -rest
                    w = low.bit_length() - 1
                    dp[mask | low][w] += dp[mask][v]
                    rest -= low
    return sum(dp[-1])


def actual_census(n):
    """All labelled codes; identify ordinary/converse classes from permutations.

    Leaf decisions use literal path orbits and an anti-automorphism, not SCCs.
    """
    ps = list(permutations(range(n)))
    full = (1 << (n * (n - 1) // 2)) - 1
    canonical = {}
    for code in range(full + 1):
        if code in canonical:
            continue
        out = graph(n, code)
        orbit = {transform(out, p) for p in ps}
        representative = min(orbit)
        for c in orbit:
            canonical[c] = representative
    reps = sorted(set(canonical.values()))
    q_spectrum, labelled, details = Counter(), Counter(), []
    seen_q = set()
    for code in reps:
        out = graph(n, code)
        aut = [p for p in ps if transform(out, p) == code]
        anti = [p for p in ps if transform(out, p) == full ^ code]
        hp = paths(out)
        h = len(hp)
        require(h == h_dp(out), 'path DP disagrees with permutations')
        f = h // len(aut)
        require(h % len(aut) == 0, 'nonfree automorphism action')

        def path_rep(p):
            return min(tuple(g[v] for v in p) for g in aut)

        path_orbits = {path_rep(p) for p in hp}
        require(len(path_orbits) == f, 'path orbit count')
        blue = 0
        if anti:
            sigma = anti[0]
            blue = sum(path_rep(tuple(sigma[v] for v in p[::-1])) == p
                       for p in path_orbits)
            require((f - blue) % 2 == 0, 'black parity')
            degree = (f - blue) // 2
        else:
            degree = f
        qid = min(code, canonical[full ^ code])
        if degree == 1:
            labelled[h] += factorial(n) // len(aut)
            if qid not in seen_q:
                q_spectrum[h] += 1
                seen_q.add(qid)
        details.append(dict(code=code, H=h, aut=len(aut), self_converse=bool(anti),
                            B=blue, degree=degree))
    require(len(canonical) == full + 1, 'incomplete labelled universe')
    require(q_spectrum == spectrum(n), f'actual Q census differs n={n}')
    expected_labels = Counter()
    for h, count in spectrum(n).items():
        k = v3(h)
        ordinary_factor = 1 if h % 5 == 0 else 2
        weight = Fraction(ordinary_factor * factorial(n), 3 ** k)
        require(weight.denominator == 1, 'nonintegral labelled count')
        expected_labels[h] = count * weight.numerator
    require(labelled == expected_labels, 'labelled weighting fails')
    return dict(n=n, labelled_universe=full + 1, ordinary_classes=len(reps),
                Q_leaf_spectrum=dict(sorted(q_spectrum.items())),
                labelled_leaf_spectrum=dict(sorted(labelled.items())),
                all_class_details=details)


def join(blocks):
    out, offset = [], 0
    total = sum(len(b) for b in blocks)
    for block in blocks:
        later = ((1 << total) - 1) ^ ((1 << (offset + len(block))) - 1)
        out.extend((row << offset) | later for row in block)
        offset += len(block)
    return out


def parameters(a):
    with localcontext() as ctx:
        ctx.prec = 65
        a = Decimal(a.numerator) / Decimal(a.denominator)
        r = Decimal('0.7')
        for _ in range(30):
            r -= (r + a * r ** 3 - 1) / (1 + 3 * a * r ** 2)
        mu = (1 - r) / (3 - 2 * r)
        var = r * (1 - r) / (3 - 2 * r) ** 3
        return dict(rho=str(r), mu=str(mu), variance=str(var),
                    logH_rate=str(mu * Decimal(3).ln()))


def finite_moments(n, labelled):
    total = first = second = mixed = 0
    scale = n // 3
    for h, count in spectrum(n).items():
        k = v3(h)
        weight = count
        if labelled:
            weight *= (1 if h % 5 == 0 else 2) * 3 ** (scale - k)
        total += weight
        first += k * weight
        second += k * k * weight
        if h % 5 == 0:
            mixed += weight
    mean = Fraction(first, total)
    var = Fraction(second, total) - mean * mean
    return dict(n=n, mean_per_n=float(mean / n), variance_per_n=float(var / n),
                mixed_probability=float(Fraction(mixed, total)))


def main():
    checks = 0
    for n in range(23):
        require(spectrum(n) == word_spectrum(n), f'word census n={n}')
        checks += 1
    # Independent generating recurrences for A and the doubled half-word.
    a = [Counter() for _ in range(61)]
    a[0][0] = 1
    for n in range(1, 61):
        a[n].update(a[n - 1])
        if n >= 3:
            for k, count in a[n - 3].items():
                a[n][k + 1] += count
    for n in range(61):
        predicted = Counter()
        for k, count in a[n].items():
            p = 0
            for center_n, center_k in [(0, 0), (1, 0), (3, 1)]:
                if n >= center_n and (n - center_n) % 2 == 0 and k >= center_k and (k - center_k) % 2 == 0:
                    p += a[(n - center_n) // 2][(k - center_k) // 2]
            require((count - p) % 2 == 0, 'GF parity')
            if count > p:
                predicted[3 ** k] = (count - p) // 2
        if n >= 4 and n % 2 == 0:
            for k, count in a[(n - 4) // 2].items():
                predicted[5 * 9 ** k] = count
        require(predicted == spectrum(n), f'GF recurrence n={n}')
        if n >= 4:
            expected_max = (5 * 9 ** ((n - 4) // 6)
                            if n % 6 in (0, 4) else 3 ** ((n - 1) // 3))
            require(max(predicted) == expected_max, f'extremal n={n}')
        require(all(15 ** k not in predicted for k in range(1, 10)), '15 power leaf')
        checks += 1
    for k in range(1, 15):
        require(all(3 ** k not in spectrum(n) for n in range(3 * k + 1)), 'early pure leaf')
        require(spectrum(3 * k + 1)[3 ** k] > 0, 'missing minimal pure leaf')
    for k in range(10):
        require(all(5 * 9 ** k not in spectrum(n) for n in range(6 * k + 4)), 'early mixed leaf')
        require(spectrum(6 * k + 4)[5 * 9 ** k] > 0, 'missing minimal mixed leaf')
    censuses = [actual_census(n) for n in range(1, 6)]
    c3, s4, one = [2, 4, 1], [2, 4, 9, 3], [0]
    examples = []
    for name, blocks, expected in [
            ('minimal_H27_leaf', [one, c3, c3, c3], 27),
            ('mixed_H45_leaf', [c3, s4, c3], 45),
            ('order14_pure_max', [one, one, c3, c3, c3, c3], 81),
            ('order16_mixed_max', [c3, c3, s4, c3, c3], 405),
            ('nonleaf_balanced_H15', [c3, s4], 15)]:
        out = join(blocks)
        require(h_dp(out) == expected, f'positive example {name}')
        examples.append(dict(name=name, n=len(out), H=expected, rows=out))
    stats = {}
    for labelled, key, weight in [(False, 'uniform_Q_leaves', Fraction(1)),
                                   (True, 'uniform_labelled_preimages', Fraction(1, 3))]:
        exact = parameters(weight)
        rows = [finite_moments(n, labelled) for n in (60, 120, 300, 600)]
        require(abs(rows[-1]['mean_per_n'] - float(exact['mu'])) < .002, 'mean finite control')
        require(abs(rows[-1]['variance_per_n'] - float(exact['variance'])) < .002, 'variance finite control')
        stats[key] = dict(parameters=exact, finite_checks=rows)
    result = dict(status='PASS', formula_comparisons=checks, word_enumeration_max_order=22,
                  recurrence_max_order=60, actual_censuses=censuses,
                  spectra={n: dict(sorted(spectrum(n).items())) for n in (4, 6, 10, 14, 15, 16, 27)},
                  independent_DP_examples=examples, asymptotic_checks=stats)
    payload = json.dumps(result, indent=2, sort_keys=True) + '\n'
    (HERE / 'RESULT.json').write_text(payload)
    print('PASS: word, generating-function, min/max, balanced-power exclusion,')
    print('literal path-orbit census through n=5, labelled weights, larger DP, and moment controls.')
    print('RESULT.json sha256=' + hashlib.sha256(payload.encode()).hexdigest())


if __name__ == '__main__':
    main()
