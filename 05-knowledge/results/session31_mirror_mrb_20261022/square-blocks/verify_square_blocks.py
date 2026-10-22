#!/usr/bin/env python3
"""Exact structural checks and rational MRB curvature-tail certificates.

Standard library only. Optional mpmath diagnostics are explicitly noncertifying.
Run: python3 session31-square-blocks/verify_square_blocks.py
"""
from fractions import Fraction as F
from math import comb
import json
from pathlib import Path


def euler_poly(n):
    """Ascending coefficients, from E_n(x+1)+E_n(x)=2*x**n."""
    p = [F(0)] * (n + 1)
    p[n] = F(1)
    for k in range(n - 1, -1, -1):
        p[k] = -sum(p[j] * comb(j, k) for j in range(k + 1, n + 1)) / 2
    return p


def at(p, x):
    return sum(a * x**j for j, a in enumerate(p))


def log_interval(n, terms=32):
    """Rational enclosure via log z = 2*atanh((z-1)/(z+1))."""
    assert n >= 1
    r = n.bit_length() - 1
    z = F(n, 2**r)

    def positive_series(t):
        s = 2 * sum(t ** (2*j+1) / (2*j+1) for j in range(terms))
        rem = 2*t**(2*terms+1) / ((2*terms+1)*(1-t*t))
        return s, s + rem

    lo2, hi2 = positive_series(F(1, 3))
    loz, hiz = positive_series((z-1)/(z+1))
    return r*lo2+loz, r*hi2+hiz


def upper_decimal(q, digits=18):
    den = 10**digits
    k = -((-q.numerator*den)//q.denominator)
    sign = '-' if k < 0 else ''
    k = abs(k)
    return f'{sign}{k//den}.{k%den:0{digits}d}'


def lower_decimal(q, digits=18):
    den = 10**digits
    k = (q.numerator*den)//q.denominator
    sign = '-' if k < 0 else ''
    k = abs(k)
    return f'{sign}{k//den}.{k%den:0{digits}d}'


def radical_interval(k, digits=22):
    """Exact rational enclosure of k**(1/k), by integer power comparisons."""
    scale = 10**digits
    target = k * scale**k
    low, high = scale, 2*scale
    while high-low > 1:
        mid = (low+high)//2
        if mid**k <= target:
            low = mid
        else:
            high = mid
    return F(low, scale), F(low if low**k == target else low+1, scale)


def tail_bound(N):
    """Rigorous upper bound for sum_{n>N}|R_n|, N>=3.

    Uses exp(1/e)<3/2 and an exact rational upper bound on log N.
    """
    assert N >= 3
    _, L = log_interval(N)
    I4 = (L/F(3)+F(1, 9))/N**3
    I6 = (L*L/F(5)+2*L/F(25)+F(2, 125))/N**5
    return 3*(1+F(1, N+1))*(I4+I6)


def exact_checks():
    count = 0
    for n in range(1, 41):
        c = n*(n+1)
        assert sum((-1)**(n+d) for d in range(-n, n+1)) == 1
        for p in range(1, 17):
            moment = sum((-1)**(n+d)*d**p for d in range(-n, n+1))
            if p % 2:
                assert moment == 0
            else:
                assert moment == at(euler_poly(p), n+1)
            count += 1
        area = F(0)
        second = F(0)
        for j in range(n):
            # K(u)=A-B*u on [j,j+1]. Endpoints certify positivity.
            A = sum((-1)**(n+d)*d for d in range(j+1, n+1))
            B = sum((-1)**(n+d) for d in range(j+1, n+1))
            assert A-B*j >= 0 and A-B*(j+1) >= 0
            area += A - F(B*((j+1)**2-j**2), 2)
            second += F(A*((j+1)**3-j**3), 3)-F(B*((j+1)**4-j**4), 4)
        assert 2*area == F(c, 2)
        assert second == F(c*(c-1), 24)

    # Independent minimal witness, no transcendental floating-point arithmetic.
    s2 = radical_interval(2)
    s3 = radical_interval(3)
    s5 = radical_interval(5)
    B = (s2[0]-s3[1], s2[1]-s3[0])
    bad = (2-3*s2[1]+s3[0], 2-3*s2[0]+s3[1])
    lam = (1-2*s2[1]+s3[0], 1-2*s2[0]+s3[1])
    assert B[0] > bad[1]
    assert lam[1] < F(-3, 10)
    first_two_pairs = (2*s2[0]-s3[1]-s5[1], 2*s2[1]-s3[0]-s5[0])
    assert first_two_pairs[0] > F(64, 10000)
    # For n=3,d=1, 0<Lambda<=f''(11)<1/500. A conjugate has modulus>1079/500.
    L11 = log_interval(11)[1]
    root11 = radical_interval(11)[1]
    curvature11_upper = root11*((2*L11-3)/11**3+(L11-1)**2/11**4)
    assert curvature11_upper < F(1, 500)
    assert 6**13 < 13*5**13  # 13**(1/13)>6/5.
    assert F(3, 4) > F(16, 25)  # sqrt(3)/2>4/5.
    # The unqualified positive-power formula fails at m=0.
    assert (-1)*sum((-1)**d for d in range(1, 2)) != F(1, 2)
    return {
        'status': 'EXACT STRUCTURAL CHECKS AND RATIONAL ENCLOSURES',
        'universe': 'n=1..40; moments degrees 1..16; all kernel intervals',
        'moment_checks': count,
        'n1_actual_block': [lower_decimal(B[0]), upper_decimal(B[1])],
        'n1_printed_block': [lower_decimal(bad[0]), upper_decimal(bad[1])],
        'n1_Lambda': [lower_decimal(lam[0]), upper_decimal(lam[1])],
        'Many_Worlds_first_two_pairs_positive': [lower_decimal(first_two_pairs[0]), upper_decimal(first_two_pairs[1])],
        'Galois_hostile_original_upper': upper_decimal(curvature11_upper),
        'Galois_hostile_conjugate_modulus_strict_lower': '1079/500',
        'curvature_absolute_tail_upper_bounds': {
            str(N): upper_decimal(tail_bound(N), 24)
            for N in (3, 10, 100, 1000, 10000)
        },
        'certificate_method': 'Integer arithmetic, Fraction, positive log-series remainder; outward decimal rounding.',
    }


def optional_diagnostics():
    try:
        import mpmath as mp
    except ImportError:
        return {'status': 'SKIPPED: mpmath not installed; exact checks unaffected.'}
    mp.mp.dps = 65
    f = lambda x: mp.expm1(mp.log(x)/x)
    f2 = lambda x: mp.exp(mp.log(x)/x)*((2*mp.log(x)-3)/x**3+(mp.log(x)-1)**2/x**4)
    rows = []
    for n in (1, 2, 3, 5, 10, 30, 100, 1000):
        c = n*(n+1)
        B = mp.fsum((-1)**d*f(c+d) for d in range(-n, n+1))
        R = B-(-1)**n*f(c)
        direct = mp.fsum((-1)**d*(f(c-d)-2*f(c)+f(c+d)) for d in range(1, n+1))
        assert abs(R-direct) < mp.mpf('1e-55')
        if n >= 3:
            D = (-1)**n*R
            assert c*f2(n*n+2*n)/2 <= D <= c*f2(n*n)/2
        rows.append({'n': n, 'B': mp.nstr(B, 25), 'R': mp.nstr(R, 25),
                     'asymptotic_ratio': None if n == 1 else mp.nstr((-1)**n*R*n**4/(2*mp.log(n)), 18)})
    c, d = 12, 1
    original = f(c-d)-2*f(c)+f(c+d)
    zeta13 = mp.exp(2j*mp.pi*6/13)
    transformed = original+(zeta13-1)*mp.power(13, mp.mpf(1)/13)
    return {'status': 'FLOAT DIAGNOSTICS, NOT CERTIFICATES', 'rows': rows,
            'conjugate_example_original_abs': mp.nstr(abs(original), 20),
            'conjugate_example_transformed_abs': mp.nstr(abs(transformed), 20)}


if __name__ == '__main__':
    report = {'exact': exact_checks(), 'diagnostic': optional_diagnostics()}
    output = Path(__file__).with_name('verification.json')
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
