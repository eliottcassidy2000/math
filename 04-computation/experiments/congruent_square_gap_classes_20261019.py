#!/usr/bin/env python3
"""Exact controls for congruent classes of positive AP square-sum triples.

No third-party packages; all arithmetic is rational/integer. Finite tests do
not replace the proofs in congruent_square_gap_classes_20261019.md. Reproduce with:
  python3 04-computation/experiments/congruent_square_gap_classes_20261019.py
"""
from fractions import Fraction as F
from math import gcd, isqrt, lcm
import json


def point_add(P, Q, N):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    X, Y = Q
    if x == X and y == -Y:
        return None
    slope = (3*x*x-N*N)/(2*y) if P == Q else (Y-y)/(X-x)
    z = slope*slope-x-X
    return z, slope*(x-z)-y


def on_curve(P, N):
    return P is None or P[1]*P[1] == P[0]**3-N*N*P[0]


def val2(q):
    q = F(q)
    assert q
    def v(n):
        n = abs(n)
        return (n & -n).bit_length()-1
    return v(q.numerator)-v(q.denominator)


def squarefree_part(n):
    assert n > 0
    N = 1
    p = 2
    while p*p <= n:
        exponent = 0
        while n % p == 0:
            exponent += 1
            n //= p
        if exponent % 2:
            N *= p
        p += 1
    return N*n


def integer_square(n):
    return n >= 0 and isqrt(n)**2 == n


def triangle(P, N):
    assert P is not None and on_curve(P, N) and P[1]
    x, y = P
    a = abs((x*x-N*N)/y)
    b = abs(2*N*x/y)
    c = abs((x*x+N*N)/y)
    assert a*a+b*b == c*c and a*b/2 == N
    return a, b, c


def primitive_triangle(P, N):
    a, b, c = triangle(P, N)
    scale = lcm(a.denominator, b.denominator, c.denominator)
    A, B, C = int(scale*a), int(scale*b), int(scale*c)
    content = gcd(gcd(A, B), C)
    A, B, C = A//content, B//content, C//content
    A, B = min(A, B), max(A, B)
    assert gcd(A, B) == 1 and A*A+B*B == C*C
    assert A*B % 2 == 0
    assert (A*B//2) % N == 0 and integer_square(A*B//(2*N))
    return A, B, C


def eligible(P, N):
    a, b, c = triangle(P, N)
    return c*c > 4*a*b


def ap_triple(P, N):
    A, B, C = primitive_triangle(P, N)
    values = tuple(2*(C*C+i*4*A*B) for i in (-1, 0, 1))
    u, v, w = values
    assert v-u == w-v == 8*A*B
    assert gcd(gcd(u, v), w) == 2
    assert u+v == (2*(B-A))**2
    assert u+w == (2*C)**2
    assert v+w == (2*(A+B))**2
    assert (v-u) % N == 0 and integer_square((v-u)//N)
    return values


def ppt_point(A, B, C, N):
    area = A*B//2
    scale = isqrt(area//N)
    assert area == N*scale*scale
    # a=A/scale, b=B/scale, c=C/scale.
    x = F(N*(A+C), B)
    y = F(2*N*scale, B)*x
    assert x > N and on_curve((x, y), N)
    return x, y


def double_check(P, N):
    x, y = P
    Q = point_add(P, P, N)
    assert Q is not None and Q[1] and on_curve(Q, N)
    assert Q[0] == (x*x+N*N)**2/(4*y*y)
    assert val2(Q[0]) < 0
    if val2(x) < 0:
        assert val2(Q[0]) == val2(x)-2
    return Q


def main():
    count = 0
    first_pass = [0, 0, 0]
    strong_margin_pass = [0, 0, 0]
    curves = set()
    # Complete Euclidean-parameter box 2 <= m <= 100, 1 <= n < m,
    # gcd(m,n)=1, opposite parity. Both ordered-leg orientations are tested.
    for m in range(2, 101):
        for n in range(1, m):
            if gcd(m, n) != 1 or (m-n) % 2 == 0:
                continue
            A, B, C = m*m-n*n, 2*m*n, m*m+n*n
            N = squarefree_part(A*B//2)
            curves.add(N)
            for a, b in [(A, B), (B, A)]:
                P = ppt_point(a, b, C, N)
                P2 = double_check(P, N)
                P4 = double_check(P2, N)
                assert val2(P4[0]) == val2(P2[0])-2
                P3 = point_add(P2, P, N)
                candidates = [P, P2, P3]
                assert all(Q is not None and Q[1] and on_curve(Q, N) for Q in candidates)
                winners = [i for i, Q in enumerate(candidates) if eligible(Q, N)]
                assert winners, (N, P)
                margin_winners = []
                for i, Q in enumerate(candidates):
                    qa, qb, qc = triangle(Q, N)
                    if (qc*qc+2*qa*qb)**2 < 2*qc**4:
                        margin_winners.append(i)
                assert margin_winners, (N, P)
                strong_margin_pass[margin_winners[0]] += 1
                first_pass[winners[0]] += 1
                for Q in candidates:
                    values = ap_triple(Q, N)
                    assert (values[0] > 0) == eligible(Q, N)
                count += 1
    # Positive seed and non-positive double: not every multiple is admissible.
    P5 = (F(25, 4), F(75, 8))
    assert ap_triple(P5, 5) == (482, 3362, 6242)
    assert not eligible(point_add(P5, P5, 5), 5)
    Q = None
    successful_multiples = []
    distinct_ppts = set()
    for k in range(1, 25):
        Q = point_add(Q, P5, 5)
        assert on_curve(Q, 5)
        double_check(Q, 5)
        if eligible(Q, 5):
            successful_multiples.append(k)
            distinct_ppts.add(primitive_triangle(Q, 5))
            ap_triple(Q, 5)
    assert len(distinct_ppts) == len(successful_multiples)
    P15 = point_add(point_add(P5, P5, 5), P5, 5)
    assert P15 == (F(127351225,2439844), F(1430549626725,3811036328))
    assert primitive_triangle(P15, 5) == (124286848883560,642780541615689,654686219104361)
    assert ap_triple(P15, 5) == (
        218114746641670847960177053922,
        857228090970326756026018436642,
        1496341435298982664091859819362,
    )
    # Sharpness hostile for a first-two-multiples assertion.
    P30 = (F(45), F(225))
    Q30 = point_add(P30, P30, 30)
    R30 = point_add(P30, Q30, 30)
    assert Q30 == (F(169, 4), F(-1547, 8))
    assert R30 == (F(2790045, 121), F(-4660327575, 1331))
    assert [eligible(T,30) for T in [P30,Q30,R30]] == [False,False,True]
    assert primitive_triangle(R30,30) == (90025452,34597057445,34597174573)
    print(json.dumps({
        "status": "PASS",
        "euclidean_parameter_bound": 100,
        "ordered_seed_points_tested": count,
        "distinct_squarefree_congruent_classes": len(curves),
        "first_eligible_multiplier_counts_1_2_3": first_pass,
        "first_strong_margin_multiplier_counts_1_2_3": strong_margin_pass,
        "N5_eligible_multiples_through_24": successful_multiples,
        "N30_first_two_hostile": True,
        "N5_third_multiple_PPT": list(primitive_triangle(P15,5)),
        "N5_third_multiple_AP_triple": list(ap_triple(P15,5)),
    }, indent=2))


if __name__ == '__main__':
    main()
