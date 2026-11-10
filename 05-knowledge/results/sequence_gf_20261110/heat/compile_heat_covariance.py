#!/usr/bin/env python3
"""Exact new covariance recurrence and bounded direct controls; stdlib only."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from math import comb, factorial
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def covariance_polynomial(n):
    """(j,v) coefficient of y^j z^v in n! E[kappa N_sqrt(y) K_z kappa]."""
    result = {(0, 0): factorial(n)}
    for j in range(1, n // 2 + 1):
        k = n - 2 * j
        for p in range(1, min(j, k) + 1):
            result[j, 2 * j + p] = (factorial(k) * comb(k, p)
                                    * comb(j - 1, p - 1) * 2 ** p)
    return result


def evaluate(poly, y, z):
    return sum(Q(c) * y ** j * z ** v for (j, v), c in poly.items())


def recurrence(N, y, z):
    a = y * z ** 2
    b = (2 * z - 1) * a
    forcing = {0: 1, 2: 3 * b, 4: -a * a - 4 * a * b, 6: a * a * b}
    out = []
    for n in range(N + 1):
        coeffs = [n, -3 * b, (2 * b - a) * n + 4 * a - 2 * b,
                  a * a + 4 * a * b, (b * b - 2 * a * b) * n - 2 * b * b + 10 * a * b,
                  -a * a * b, -a * b * b * (n - 6)]
        out.append(Q(forcing.get(n, 0)) + sum(c * out[n - d]
                   for d, c in enumerate(coeffs, 1) if n >= d))
    return out


def path_count(n, mask):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pred = [0] * n
    for bit, (i, j) in enumerate(edges):
        u, v = (i, j) if mask >> bit & 1 else (j, i)
        pred[v] |= 1 << u
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for subset in range(1, 1 << n):
        for v in range(n):
            if not subset >> v & 1 or subset == 1 << v:
                continue
            rest = subset ^ (1 << v)
            p = rest & pred[v]
            while p:
                bit = p & -p
                dp[subset][v] += dp[rest][bit.bit_length() - 1]
                p -= bit
    return sum(dp[-1])


def fwht(a):
    a = a[:]
    width = 1
    while width < len(a):
        for i in range(0, len(a), 2 * width):
            for j in range(i, i + width):
                u, v = a[j], a[j + width]
                a[j], a[j + width] = u + v, u - v
        width *= 2
    return a


def direct_control(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    hs = [path_count(n, mask) for mask in range(1 << len(edges))]
    total = len(hs)
    coeffs = fwht(hs)
    energy = defaultdict(Q)
    # Squaring removes the FWHT's coordinate convention.
    for mask, c in enumerate(coeffs):
        if not c:
            continue
        e = mask.bit_count()
        require(e % 2 == 0, "Unexpected odd Fourier degree")
        vertices = set()
        for bit, pair in enumerate(edges):
            if mask >> bit & 1:
                vertices.update(pair)
        energy[e // 2, len(vertices)] += Q(c * c * 4 ** (n - 1), factorial(n) * total * total)
    require(dict(energy) == covariance_polynomial(n), "Covariance formula mismatch at n=%s" % n)
    raw = [sum(Q(h ** r, total) for h in hs) for r in range(5)]
    mean = raw[1]
    centered3 = sum((Q(h, 1) / mean - 1) ** 3 for h in hs) / total
    # Exact finite coefficient inversion on Y=(H-1)/2, using safe n! bound.
    M = (factorial(n) - 1) // 2
    beta = [sum(Q(comb((h - 1) // 2, r), total) for h in hs) for r in range(M + 1)]
    recovered = {2 * j + 1: sum((-1) ** (r - j) * comb(r, j) * beta[r]
                               for r in range(j, M + 1)) for j in range(M + 1)}
    histogram = Counter(hs)
    require({h: p for h, p in recovered.items() if p} == {h: Q(c, total) for h, c in histogram.items()},
            'Binomial support inversion mismatch')
    bonferroni_checks = 0
    for j in range(M + 1):
        actual = Q(histogram.get(2 * j + 1, 0), total)
        for length in set([0, 1, 2, 3, M - j, M - j + 1]):
            bound = sum((-1) ** (r - j) * comb(r, j) * beta[r]
                        for r in range(j, min(j + length, M) + 1))
            require(bound <= actual if length % 2 else bound >= actual,
                    'Alternating-moment bound failed at n=%s j=%s L=%s' % (n, j, length))
            if length >= M - j:
                require(bound == actual, 'Full truncation must be exact')
            bonferroni_checks += 1
    return {'n': n, 'labelled_universe': total, 'histogram': dict(sorted(histogram.items())),
            'moments_H_0_through_4': list(map(str, raw)), 'normalized_centered_third_moment': str(centered3),
            'walsh_covariance_formula_matches': True, 'binomial_support_inversion_matches': True,
            'odd_lower_even_upper_Bonferroni_checks': bonferroni_checks}


def birth_energy(N, a):
    """EGF exp(x+2 a x^3/(1-a x^2)); covariance of coherent birth martingales."""
    vals = [Q(1)]
    for n in range(N):
        fall2 = n * (n - 1)
        fall4 = n * (n - 1) * (n - 2) * (n - 3)
        get = lambda j: vals[j] if j >= 0 else Q(0)
        vals.append(vals[n] + 2 * a * fall2 * get(n - 1) + 4 * a * fall2 * get(n - 2)
                    - a * a * fall4 * (get(n - 3) + get(n - 4)))
    for n, value in enumerate(vals):
        exact = Q(1)
        for j in range(1, n // 2 + 1):
            for p in range(1, min(j, n - 2 * j) + 1):
                exact += Q(2 ** p * factorial(n) * comb(j - 1, p - 1),
                           factorial(n - 2 * j - p) * factorial(p)) * a ** j
        require(value == exact, "Coherent birth energy recurrence mismatch")
    return vals


def connected_triple_overlap(n):
    """C_3's z^n coefficients by two forest masks and the XOR third mask."""
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    hs = [path_count(n, mask) for mask in range(1 << len(edges))]
    ft = fwht(hs)
    weights = {}
    for mask, coeff in enumerate(ft):
        if coeff:
            e = mask.bit_count()
            w = Q(coeff * 2 ** (n - 1), len(hs) * factorial(n - e))
            require(w.denominator == 1, "Nonintegral forest weight")
            weights[mask] = w.numerator
    def connected(mask):
        if not mask:
            return False
        seen = {0}
        while True:
            old = len(seen)
            for bit, (i, j) in enumerate(edges):
                if mask >> bit & 1 and (i in seen or j in seen):
                    seen.update((i, j))
            if old == len(seen):
                return len(seen) == n
    coeffs = defaultdict(Q)
    for a, wa in weights.items():
        for b, wb in weights.items():
            c = a ^ b
            if c not in weights or not connected(a | b):
                continue
            exponents = (a.bit_count(), b.bit_count(), c.bit_count())
            coeffs[exponents] += Q(wa * wb * weights[c], factorial(n))
    return {' '.join('t%d^%d' % (i + 1, e) for i, e in enumerate(es)): str(c)
            for es, c in sorted(coeffs.items()) if c}


def main():
    weights = [(Q(1), Q(1)), (Q(1), Q(1, 2)), (Q(1, 4), Q(2, 3)),
               (Q(0), Q(1)), (Q(1), Q(0)), (Q(2), Q(3))]
    checks = []
    for y, z in weights:
        vals = recurrence(40, y, z)
        require(vals == [evaluate(covariance_polynomial(n), y, z) for n in range(41)],
                "Across-order covariance recurrence mismatch")
        checks.append({'y': str(y), 'z': str(z), 'checked_n': [0, 40],
                       'prefix_0_through_10': list(map(str, vals[:11]))})
    result = {'status': 'FINITE-EXACT new formulas only; inherited theorem claims not re-audited',
              'recurrence_controls': checks,
              'coherent_birth_energy_a_1_prefix_0_to_15': list(map(str, birth_energy(40, Q(1))[:16])),
              'coherent_birth_energy_a_1_over_4_prefix_0_to_15': list(map(str, birth_energy(40, Q(1, 4))[:16])),
              'covariance_polynomials': {str(n): {f'y^{j} z^{v}': c for (j, v), c in covariance_polynomial(n).items()}
                                        for n in range(11)},
              'direct_controls': [direct_control(n) for n in range(2, 6)],
              'C3_connected_kernel_coefficients_by_vertex_size': {str(n): connected_triple_overlap(n) for n in range(2, 6)}}
    path = Path(__file__).with_name('output.json')
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'output': str(path), 'recurrence_cases': len(checks), 'largest_recurrence_n': 40,
                      'direct_controls': result['direct_controls']}, indent=2))


if __name__ == '__main__':
    main()
