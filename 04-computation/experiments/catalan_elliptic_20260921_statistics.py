"""Exact finite controls for proper-divisor statistics; standard library only.

Run normally or with -O. Output contains exact counts, not evidence for an
unstated rate or an orbit law. All checks are explicit and survive -O.
"""
from array import array
from fractions import Fraction
from hashlib import sha256
from math import isqrt, log, pi
from pathlib import Path
import argparse
import json


def check(condition, message):
    if not condition:
        raise ArithmeticError(message)


def sieve(limit):
    spf = array('I', [0]) * (limit + 1)
    for p in range(2, limit + 1):
        if spf[p] == 0:
            for n in range(p, limit + 1, p):
                if spf[n] == 0:
                    spf[n] = p
    tau = array('I', [0]) * (limit + 1)
    omega = array('I', [0]) * (limit + 1)
    bigomega = array('I', [0]) * (limit + 1)
    exponent = array('I', [0]) * (limit + 1)
    mu = array('b', [0]) * (limit + 1)
    tau[1], mu[1] = 1, 1
    for n in range(2, limit + 1):
        p, m = spf[n], n // spf[n]
        bigomega[n] = bigomega[m] + 1
        if m % p == 0:
            exponent[n] = exponent[m] + 1
            omega[n] = omega[m]
            tau[n] = tau[m] // (exponent[m] + 1) * (exponent[n] + 1)
        else:
            exponent[n] = 1
            omega[n] = omega[m] + 1
            tau[n] = 2 * tau[m]
            mu[n] = -mu[m]
    return tau, omega, bigomega, mu


def trial_profile(n):
    result = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            a = 0
            while n % p == 0:
                n //= p
                a += 1
            result.append(a)
        p += 1
    if n > 1:
        result.append(1)
    return sorted(result)


def divisor_sum(x):
    r = isqrt(x)
    return 2 * sum(x // d for d in range(1, r + 1)) - r * r


def rational(n, d):
    f = Fraction(n, d)
    return {'numerator': f.numerator, 'denominator': f.denominator}


def run():
    limit, direct_limit = 1_000_000, 10_000
    tau, omega, bigomega, mu = sieve(limit)
    primes = [n for n in range(2, limit + 1) if tau[n] == 2]
    squarefree = [int(m != 0) for m in mu]
    prime_indicator = [int(t == 2) for t in tau]
    F = [tau[n] - 2 + int(n == 1) for n in range(limit + 1)]
    S = [(1 << omega[n]) - 1 - squarefree[n] + int(n == 1)
         for n in range(limit + 1)]
    U = [omega[n] - prime_indicator[n] for n in range(limit + 1)]
    D = [F[n] - S[n] - U[n] for n in range(limit + 1)]

    # Independent trial division and direct proper-divisor incidence.
    direct_F = [0] * (direct_limit + 1)
    direct_S = [0] * (direct_limit + 1)
    direct_U = [0] * (direct_limit + 1)
    for d in range(2, direct_limit + 1):
        prof = trial_profile(d)
        d_squarefree = all(a == 1 for a in prof)
        d_prime = prof == [1]
        for n in range(2 * d, direct_limit + 1, d):
            direct_F[n] += 1
            direct_S[n] += d_squarefree
            direct_U[n] += d_prime
    for n in range(1, direct_limit + 1):
        prof = trial_profile(n)
        check((F[n], S[n], U[n]) ==
              (direct_F[n], direct_S[n], direct_U[n]), ('direct divisors', n))
        check((omega[n], bigomega[n], squarefree[n]) ==
              (len(prof), sum(prof), int(all(a == 1 for a in prof))),
              ('independent factor profile', n))

    # The exact sign classification, including 1 and prime endpoints.
    for n in range(1, limit + 1):
        expected_zero = (n == 1 or prime_indicator[n] or
                         (omega[n] == 1 and bigomega[n] == 3) or
                         (omega[n] == 3 and bigomega[n] == 4))
        expected_negative = ((squarefree[n] and omega[n] >= 2) or
                             (omega[n] == 1 and bigomega[n] == 2) or
                             (omega[n] == 2 and bigomega[n] == 3))
        check((D[n] == 0) == bool(expected_zero), ('zero shapes', n))
        check((D[n] < 0) == bool(expected_negative), ('negative shapes', n))
        if omega[n] >= 4:
            check((D[n] < 0) == bool(squarefree[n]), ('sign detector', n))

    cutoffs = [1, 2, 3, 4, 8, 10, 100, 1000, 10_000, 100_000, limit]
    records = []
    for x in cutoffs:
        ps = [p for p in primes if p <= x]
        q = sum(squarefree[1:x + 1])
        all_divisors = sum(tau[1:x + 1])
        all_sf_divisors = sum(1 << omega[n] for n in range(1, x + 1))
        floor_sf = sum(squarefree[d] * (x // d) for d in range(1, x + 1))
        mobius_q = sum(mu[r] * (x // (r * r)) for r in range(1, isqrt(x) + 1))
        convolution = sum(mu[r] * divisor_sum(x // (r * r))
                          for r in range(1, isqrt(x) + 1))
        check(all_divisors == divisor_sum(x), ('hyperbola', x))
        check(all_sf_divisors == floor_sf == convolution, ('sf convolution', x))
        check(q == mobius_q, ('squarefree indicator', x))
        sums = {'F': sum(F[1:x + 1]), 'S': sum(S[1:x + 1]),
                'U': sum(U[1:x + 1]), 'D': sum(D[1:x + 1]),
                'tau': all_divisors, 'two_to_omega': all_sf_divisors,
                'omega': sum(omega[1:x + 1]), 'Omega': sum(bigomega[1:x + 1])}
        check(sums['F'] == all_divisors - 2*x + 1, ('F endpoint', x))
        check(sums['S'] == floor_sf - x - q + 1, ('S endpoint', x))
        check(sums['U'] == sum(x // p for p in ps) - len(ps), ('U endpoint', x))
        pp_sum = 0
        for p in ps:
            power = p
            while power <= x:
                pp_sum += x // power
                power *= p
        check(sums['Omega'] == pp_sum, ('prime power incidence', x))
        signs = {'negative': sum(d < 0 for d in D[1:x + 1]),
                 'zero': sum(d == 0 for d in D[1:x + 1]),
                 'positive': sum(d > 0 for d in D[1:x + 1])}
        check(sum(signs.values()) == x, ('sign partition', x))
        records.append({'X': x, 'prime_count': len(ps), 'squarefree_count': q,
                        'sums': sums, 'means_exact': {k: rational(v, x) for k, v in sums.items()},
                        'defect_sign_counts': signs,
                        'display_only': {'squarefree_density': q/x,
                                         'prime_density': len(ps)/x,
                                         'F_mean': sums['F']/x, 'S_mean': sums['S']/x,
                                         'U_mean': sums['U']/x, 'D_mean': sums['D']/x,
                                         'six_over_pi_squared': 6/pi**2,
                                         'log_X': log(x)}})

    rows = []
    for a in (1, 3, 5):
        ns = range(a, limit + 1, 6)
        count = len(ns)
        rows.append({'residue_mod_6': a, 'count': count,
                     'squarefree_count': sum(squarefree[n] for n in ns),
                     'negative_defect_count': sum(D[n] < 0 for n in ns),
                     'zero_defect_count': sum(D[n] == 0 for n in ns),
                     'positive_defect_count': sum(D[n] > 0 for n in ns),
                     'limiting_negative_density': '6/pi^2' if a == 3 else '9/pi^2'})

    crt_controls = []
    for qs in ([2, 3, 5], [2, 3, 5, 7], [2, 3, 5, 7, 11, 13]):
        polynomial, modulus = [1], 1
        for p in qs:
            nxt = [0] * (len(polynomial) + 1)
            for j, coeff in enumerate(polynomial):
                nxt[j] += (p - 1) * coeff
                nxt[j + 1] += coeff
            polynomial, modulus = nxt, modulus * p
        actual = [0] * len(polynomial)
        for n in range(modulus):
            actual[sum(n % p == 0 for p in qs)] += 1
        check(actual == polynomial, ('finite CRT law', qs))
        crt_controls.append({'primes': qs, 'modulus': modulus,
                             'hit_count_coefficients': polynomial,
                             'at_most_three_hits_fraction': rational(sum(polynomial[:4]), modulus)})

    controls = {str(n): {'F': F[n], 'S': S[n], 'U': U[n], 'D': D[n],
                         'tau': tau[n], 'omega': omega[n], 'Omega': bigomega[n]}
                for n in (1, 2, 4, 6, 8, 12, 30, 60, 210, 420)}
    return {'status': 'FINITE-EXACT controls; asymptotic claims require the note proofs',
            'universe': {'all_integers': [1, limit], 'direct_divisor_checks': [1, direct_limit],
                         'cutoffs': cutoffs, 'signed_extension': 'Functions evaluated at |n|; zero excluded'},
            'controls': controls, 'prefixes': records, 'odd_residue_rows': rows,
            'finite_prime_CRT': crt_controls,
            'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    output = run()
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + '\n',
                           encoding='utf-8', newline='\n')
    print(json.dumps({'output': str(args.output), 'source_sha256': output['source_sha256'],
                      'checks': 'PASS', 'max_X': 1_000_000}))
