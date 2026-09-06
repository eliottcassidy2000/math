#!/usr/bin/env python3
"""Exact controls for the NC2 root-rotation/correspondence obstruction.

Universe: coprime positive binomial widths M,N<=12; all primitive width-three
binomials with 4<=N<=20; and the sharp (M,N)=(3,4) formal roots through s^30.
The algebraic-degree theorem is proved in the companion note, not inferred
from this finite computation. No floating point, random tests, or asserts.
"""

from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd

import sympy as sp


checks = 0


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


def choose_fraction(a, n):
    value = F(1)
    for j in range(n):
        value *= (a - j) / (j + 1)
    return value


# Exact Q(omega), omega^2+omega+1=0, as coefficient pairs.
ZERO, ONE, OMEGA = (F(0), F(0)), (F(1), F(0)), (F(0), F(1))


def qa(a, b):
    return (a[0] + b[0], a[1] + b[1])


def qm(a, b):
    return (a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0] - a[1] * b[1])


def qp(a, n):
    value = ONE
    for _ in range(n):
        value = qm(value, a)
    return value


CAP = 30


def pa(a, b):
    result = dict(a)
    for degree, coefficient in b.items():
        result[degree] = qa(result.get(degree, ZERO), coefficient)
    return {d: c for d, c in result.items() if c != ZERO}


def pn(a):
    return {d: (-c[0], -c[1]) for d, c in a.items()}


def pm(a, b):
    result = {}
    for da, ca in a.items():
        for db, cb in b.items():
            if da + db <= CAP:
                result[da + db] = qa(result.get(da + db, ZERO), qm(ca, cb))
    return {d: c for d, c in result.items() if c != ZERO}


def pp(a, n):
    value = {0: ONE}
    for _ in range(n):
        value = pm(value, a)
    return value


def binomial_moments(M, N):
    polynomial = {-M: 1, N: 1}
    power = {0: 1}
    moments = []
    for _ in range(M + N):
        nxt = {}
        for a, ca in power.items():
            for b, cb in polynomial.items():
                nxt[a + b] = nxt.get(a + b, 0) + ca * cb
        power = nxt
        moments.append(power.get(0, 0))
    return moments


def main():
    primitive_rows = 0
    for M in range(1, 13):
        for N in range(1, 13):
            if gcd(M, N) != 1:
                continue
            primitive_rows += 1
            moments = binomial_moments(M, N)
            require(all(c == 0 for c in moments[:-1]), (M, N, "early return"))
            require(moments[-1] == comb(M + N, M), (M, N, "sharp return"))
            # alpha_j^M has distinct values for the d-th roots alpha_j.
            d = M + N
            require(len({(M * j) % d for j in range(d)}) == d,
                    (M, N, "critical value collision"))
    print(f"primitive_binomial_rows={primitive_rows}; widths=1..12")

    degrees = []
    for N in range(4, 21):
        if gcd(3, N) == 1:
            degrees.append((N, comb(N + 2, 2)))
    require(degrees[0] == (4, 15), "smallest interior width-three hostile")
    print(f"width_three_pair_degrees={degrees}")

    # Formal roots of u^3=s^3(1+u^7), from Lagrange inversion. Check the
    # defining equation independently rather than trusting that formula.
    roots = []
    for branch in range(3):
        root = {}
        for j in range((CAP - 1) // 7 + 1):
            degree = 7 * j + 1
            coefficient = choose_fraction(F(degree, 3), j) / degree
            root[degree] = qm((coefficient, F(0)), qp(OMEGA, branch * degree))
        roots.append(root)
        require(pa(pp(root, 3), pn(pm({3: ONE}, pa({0: ONE}, pp(root, 7))))) == {},
                (branch, "formal root equation"))

    trace_all = pa(pa(roots[0], roots[1]), roots[2])
    product_all = pm(pm(roots[0], roots[1]), roots[2])
    require(trace_all == {15: (F(2), F(0))}, "trace begins at s^15")
    require(product_all == {3: ONE, 24: (F(5), F(0))}, "product moment identity")
    norm = pm(roots[1], roots[2])
    norm_defect = pa(pm(norm, pa({0: ONE}, pp(roots[0], 7))), pn(pp(roots[0], 2)))
    require(min(norm_defect) == 23 and norm_defect[23] == (F(5), F(0)),
            "norm rational approximation begins at degree 23")
    print("sharp_hostile: first_return=7; CT=35; trace_defect_order=15; norm_defect_order=23")

    y, z, S, P = sp.symbols("y z S P")
    H = z**3 * y**7 - (1 + z**7) * y**3 + z**3
    Q = z**3 * y**6 + z**4 * y**5 + z**5 * y**4 + z**6 * y**3 - y**2 - z*y - z**2
    require(sp.expand(H - (y-z)*Q) == 0, "marked-root quotient")
    remainder = sp.rem(Q, y**2-S*y+P, y).expand()
    A = z**3*(S**5-4*S**3*P+3*S*P**2) + z**4*(S**4-3*S**2*P+P**2) + z**5*(S**3-2*S*P) + z**6*(S**2-P)-S-z
    B = z**3*(-S**4*P+3*S**2*P**2-P**3) + z**4*(-S**3*P+2*S*P**2) + z**5*(-S**2*P+P**2)-z**6*S*P+P-z**2
    require(sp.expand(remainder-A*y-B) == 0, "faithful pair correspondence")
    # Positive survivor: dilation by three reduces (-3,6) to width one.
    require(binomial_moments(3, 6)[2] == 3, "dilation survivor first return")
    # Hostile to dropping coprimality: critical values genuinely collide.
    require(len({3*j % 9 for j in range(9)}) == 3, "nonprimitive critical collision")
    print("quadratic_divisibility_equations=verified; dilation_survivor=(-3,6); collision_control=gcd(3,9)=3")

    # Entire degree-105 marked-pair monodromy universe, not a sample.
    configurations = {(i, pair) for i in range(7)
                      for pair in combinations([j for j in range(7) if j != i], 2)}
    require(len(configurations) == 105, "marked-pair degree over t")

    def cycle_lengths(permutation):
        remaining = set(configurations)
        lengths = []
        while remaining:
            current = min(remaining)
            length = 0
            while current in remaining:
                remaining.remove(current)
                length += 1
                current = (permutation[current[0]],
                           tuple(sorted(permutation[j] for j in current[1])))
            lengths.append(length)
        return sorted(lengths)

    simple = cycle_lengths((1, 0, 2, 3, 4, 5, 6))
    zero = cycle_lengths((1, 2, 0, 4, 5, 6, 3))
    require(simple == [1]*35+[2]*35, "simple branch marked-pair cycle type")
    require(zero == [3]+[4]*3+[6]+[12]*7, "zero branch marked-pair cycle type")
    ramification = 7*(105-len(simple)) + 105-len(zero)
    require(ramification == 338, "complete ramification sum")
    genus = (2-2*105+ramification)//2
    require(genus == 65, "Riemann-Hurwitz genus")
    print("marked_pair_curve: degree_over_t=105; simple_cycle_type=1^35,2^35; zero_cycle_type=3,4^3,6,12^7; genus=65")
    print(f"checks={checks}; status=PASS; algebraic_degree_proof=companion_note_not_finite_inference")


if __name__ == "__main__":
    main()
