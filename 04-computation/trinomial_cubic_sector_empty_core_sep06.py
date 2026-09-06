#!/usr/bin/env python3
"""Exact boundary checks for the THM4440 cubic-core sector.

The all-parameter sector count is proved in the companion note. These are
named endpoint controls and exact polynomial identities, not a census.
"""
from fractions import Fraction as Q
from math import comb, gcd
from hashlib import sha256
import sympy as S

T, U, Z = S.symbols('tau u s')
BOUNDARY = -S.Rational(4, 27)
CHECKS = 0


def need(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def main():
    core = T+Z**2+Z**3
    need(S.factor(S.discriminant(core, Z)+T*(4+27*T)) == 0,
         'exact compressed-core discriminant')
    p2 = 6*T**2+20*U*T+U*(U-1)
    p3 = 72*T**3+504*U*T**2+84*U*(U-1)*T+U*(U-1)*(U-2)
    f2 = U**2-S.Rational(107, 27)*U+S.Rational(32, 243)
    f3 = U**3-S.Rational(139, 9)*U**2+S.Rational(2066, 81)*U-S.Rational(512, 2187)
    need(S.expand(p2.subs(T, BOUNDARY)-f2) == 0, 'quadratic sector endpoint')
    need(S.expand(p3.subs(T, BOUNDARY)-f3) == 0, 'cubic sector endpoint')
    need(f2.subs(U, 3) == -S.Rational(670, 243), 'quadratic g8 endpoint')
    need(f2.subs(U, 4) == S.Rational(68, 243), 'quadratic g9 endpoint')
    need(S.expand(27*S.diff(p3, T).subs(T, BOUNDARY).subs(U, Z+4)) ==
         2268*Z**2+11844*Z+11216, 'cubic monotonicity on sector')
    need(f3.subs(U, 4) == -S.Rational(177848, 2187), 'cubic u4 endpoint')
    need(f3.subs(U, 13) == -S.Rational(178820, 2187), 'cubic u13 endpoint')
    need(S.diff(f3, U).subs(U, 0) > 0 and S.diff(f3, U).subs(U, 4) < 0 and
         S.diff(f3, U).subs(U, 13) > 0, 'exact two-critical-point placement')
    shifted = S.Poly(f3.subs(U, Z+14).expand(), Z).all_coeffs()
    need(shifted == [1, S.Rational(239, 9), S.Rational(14666, 81),
                     S.Rational(161272, 2187)], 'positive cubic tail shift')
    print('core_real_rooted_sector -4/27 <= tau < 0')
    print('quadratic_threshold_polynomial', f2)
    print('cubic_threshold_polynomial', f3)
    print('cubic_positive_shift_u14', shifted)
    manifest = []
    for a, gs, polynomial, offset in ((15, (8, 11, 13), p2, 5),
                                      (21, (11, 13, 16, 17, 19, 20, 22, 23), p3, 7)):
        for g in gs:
            need(gcd(a, g) == 1, 'primitive named family control')
            p = S.Poly(polynomial.subs(U, g-offset), T)
            need(p.count_roots(-S.oo, 0) == p.degree(), 'all first roots negative')
            need(p.eval(BOUNDARY) != 0, 'no boundary equality in named controls')
            count = p.count_roots(BOUNDARY, 0)
            expected = int(g == 8) if a == 15 else int(g <= 20)
            need(count == expected, 'exact named sector-root count')
            print('named_control', (-a, 2*g-a, 3*g-a), 'g', g,
                  'first_degree', p.degree(), 'sector_roots', count)
            manifest.append((a, g, str(p), count))
    need([g for g in range(11, 21) if gcd(g, 21) == 1] == [11, 13, 16, 17, 19, 20],
         'complete finite primitive sector exception list')
    print('semantic_sha256', sha256(repr(manifest).encode()).hexdigest())
    print('explicit_checks', CHECKS)
    print('PASS: exact gauge-sector boundary and named Sturm controls')
    print('SCOPE: THM4440 sector overlap only; full doubled-row signs need separate proofs')


if __name__ == '__main__':
    main()
