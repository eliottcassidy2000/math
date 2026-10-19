#!/usr/bin/env python3
"""Exact finite controls and clearly labeled floating primorial diagnostics."""
from fractions import Fraction
from math import prod, log, log1p, exp, isqrt
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def primes(n):
    a = bytearray(b'\1') * (n+1)
    a[:2] = b'\0\0'
    for p in range(2, isqrt(n)+1):
        if a[p]:
            a[p*p:n+1:p] = b'\0' * ((n-p*p)//p+1)
    return [p for p in range(2, n+1) if a[p]]


def rho(p, m):
    t = m % p
    return 2 if t == 0 else 3 if (9*t*t-1) % p == 0 else 4


def is_clear(a, m, ps):
    return all(all(n % p for n in (6*a-1,6*a+1,6*(m-a)-1,6*(m-a)+1))
               for p in ps)


def main():
    pp = primes(97)
    tail_cases = 0
    for m in range(2, 202):
        for z in (11, 19, 43):
            ps = [p for p in pp if p >= z]
            actual = prod(Fraction(1,1)-Fraction(rho(p,m),p) for p in ps)
            den = prod((1-Fraction(1,p))**4 for p in ps)
            actual /= den
            generic = prod((1-Fraction(4,p))/(1-Fraction(1,p))**4 for p in ps)
            boost = prod(Fraction(p-rho(p,m),p-4) for p in ps)
            saturated = prod(Fraction(p-2,p-4) for p in ps)
            assert actual == generic*boost
            assert generic <= actual <= saturated
            tail_cases += 1

    intervals = 0
    targets = [103, 701, 5003] + [prod(p for p in pp if p <= y) for y in (13,29,43)]
    max_discrepancy = Fraction(0)
    for m in targets:
        for length in (5,17,101):
            if length >= m:
                continue
            for start in sorted({0, (m-length)//2, m-length-1}):
                for z in (5,7,11,13):
                    ps = [p for p in pp if 5 <= p < z]
                    w = prod(ps)
                    count_residues = prod(p-rho(p,m) for p in ps)
                    expected = Fraction(length*count_residues,w)
                    count = sum(is_clear(a,m,ps) for a in range(start+1,start+length+1))
                    error = abs(count-expected)
                    assert error <= count_residues
                    max_discrepancy = max(max_discrepancy,error/count_residues)
                    intervals += 1

    # These floating values illustrate convergence, not certified enclosures.
    # The analytic residual bound shown is rigorously proved separately:
    # -13/B <= log(product_{p>B} b_p(m)) <= 4 log N/(B log B).
    bigp = primes(1_000_000)
    diagnostics = []
    for y in (100,500,2000,10000):
        m = prod(p for p in bigp if p <= y)
        z = y**0.5
        B = 1_000_000
        log_t_finite = sum(log1p(-rho(p,m)/p)-4*log1p(-1/p)
                           for p in bigp if p >= z)
        N = m*(9*m*m-1)
        error_minus = 13/B
        error_plus = 4*log(N)/(B*log(B))
        diagnostics.append({
            'primorial_cutoff_y': y,
            'sieve_cutoff_sqrt_y': z,
            'finite_tail_product_through_1000000': exp(log_t_finite),
            'illustrative_lower_from_analytic_log_tail_bound': exp(log_t_finite-error_minus),
            'illustrative_upper_from_analytic_log_tail_bound': exp(log_t_finite+error_plus),
            'predicted_limit': 4,
            'rounding_status': 'ordinary floating arithmetic; not a certified numeric interval',
        })
    out = {
        'exact_status': 'all finite exact controls passed',
        'finite_tail_factorization_and_envelope_cases': tail_cases,
        'translated_interval_cases': intervals,
        'maximum_interval_error_divided_by_residue_count': str(max_discrepancy),
        'primorial_diagnostics': diagnostics,
        'scope': 'Finite checks support formulas; asymptotic uniformity and limits are proved in the note.',
    }
    (HERE/'short_interval_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    lines = [f'exact_status: {out["exact_status"]}', f'finite tail cases: {tail_cases}',
             f'translated interval cases: {intervals}',
             f'maximum discrepancy ratio: {max_discrepancy}',
             'Primorial diagnostics (floating, not certified intervals):']
    for row in diagnostics:
        lines.append(f'  y={row["primorial_cutoff_y"]}; finite product={row["finite_tail_product_through_1000000"]:.9f}; analytical-tail interval illustrated numerically=[{row["illustrative_lower_from_analytic_log_tail_bound"]:.9f}, {row["illustrative_upper_from_analytic_log_tail_bound"]:.9f}]; limit=4')
    text = '\n'.join(lines)+'\n'
    (HERE/'short_interval_checks.txt').write_text(text)
    print(text,end='')


if __name__ == '__main__':
    main()
