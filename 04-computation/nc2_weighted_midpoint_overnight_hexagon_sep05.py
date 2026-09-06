#!/usr/bin/env python3
"""Exact controls for full-support midpoint deletion and its joint-sign boundary.

The unbounded proofs are in the companion note. No repository mathematical
producer is imported. Literal path DP, binomial convolution, Toeplitz minors,
exact root isolation and separate polynomial hostiles retain different maps.
"""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb

import sympy as sp


t, w = sp.symbols("t w")
GATES = 0
TRACE = sha256()


def require(condition, message):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(message)


def record(*data):
    TRACE.update((repr(data)+"\n").encode())


def path_counts(U, V, p, q, middle_U, middle_V):
    level = p*middle_U+q*middle_V
    table = {}
    for X in range(U+1):
        for Y in range(V+1):
            if X == Y == 0:
                pair = (1, 0)
            else:
                a, b = table.get((X-1, Y), (0, 0)), table.get((X, Y-1), (0, 0))
                pair = a[0]+b[0], a[1]+b[1]
            if p*X+q*Y == level and (middle_U-X) % q == 0:
                ell = (middle_U-X)//q
                if Y == middle_V+p*ell:
                    pair = 0, pair[0]+pair[1]
            table[X, Y] = pair
    missed, hit = table[U, V]
    return missed+hit, hit


def real_nonpositive(poly):
    if poly.is_zero:
        return True
    intervals = poly.intervals()
    return (sum(mult for _, mult in intervals) == poly.degree()
            and all(b <= 0 for (_, b), _ in intervals))


def determinant(matrix):
    if len(matrix) == 2:
        return matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]
    a, b, c = matrix
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])


def interval_value(poly, a, b):
    lower = upper = F(0)
    for coefficient in poly.all_coeffs():
        products = (lower*a, lower*b, upper*a, upper*b)
        lower, upper = min(products)+F(coefficient), max(products)+F(coefficient)
    return lower, upper


def nonpositive_at_roots(A, target):
    squarefree = A.sqf_part()
    shared = squarefree.gcd(target)
    remaining = squarefree.exquo(shared)
    remainder = target.rem(remaining) if remaining.degree() else sp.Poly(0, t)
    certificates = []
    for (a, b), mult in remaining.intervals():
        require(mult == 1 and b <= 0 and remaining.eval(0) != 0, "all distinct noncommon roots are negative")
        a, b = F(a), F(b)
        for _ in range(160):
            lo, hi = interval_value(remainder, a, b)
            if hi < 0:
                break
            if lo > 0:
                raise RuntimeError("positive endpoint-doubling value at first root")
            a, b = map(F, remaining.refine_root(a, b, eps=sp.Rational(b-a)/4))
        else:
            raise RuntimeError("root-sign refinement exhausted")
        certificates.append(tuple(map(str, (a, b, lo, hi))))
    require(len(certificates)+shared.degree() == squarefree.degree(), "all distinct first roots accounted for")
    return shared.degree(), len(certificates), certificates


def network_controls():
    cases = endpoints = weighted = zero_rows = minors = common = negative = 0
    for p in range(1, 4):
        for q in range(1, 4):
            for U in range(4):
                for V in range(4):
                    if U+V == 0:
                        continue
                    lo, hi = -(V//p), U//q
                    first = {j: comb(U+V+(p-q)*j, V+p*j) for j in range(lo, hi+1)}
                    lower, upper = -(2*V//p), 2*U//q
                    full, hit = {}, {}
                    for j in range(lower, upper+1):
                        X, Y = 2*U-q*j, 2*V+p*j
                        full[j] = comb(X+Y, Y)
                        hit[j] = sum(value*first.get(j-i, 0) for i, value in first.items())
                        dp_full, dp_hit = path_counts(X, Y, p, q, U, V)
                        require((dp_full, dp_hit) == (full[j], hit[j]), "independent full/selected path counts")
                        require(0 <= hit[j] <= full[j], "midpoint injection includes every signed support index")
                        endpoints += 1
                    for weight in (F(0), F(1, 2), F(1), F(2), F(7)):
                        coefficients = {j: F(full[j])+(weight-1)*hit[j] for j in full}
                        require(all(c >= 0 for c in coefficients.values()), "nonnegative weighted path coefficients")
                        poly = sp.Poly(sum(sp.Rational(c.numerator, c.denominator)*t**(j-lower) for j, c in coefficients.items()), t)
                        require(real_nonpositive(poly), "full weighted support is real-rooted with multiplicities retained")
                        zero_rows += int(poly.is_zero)
                        for size in (2, 3):
                            for I in combinations(range(4), size):
                                for J in combinations(range(4), size):
                                    minor = determinant([[coefficients.get(j-i, F(0)) for j in J] for i in I])
                                    require(minor >= 0, "independent finite Toeplitz minor control")
                                    minors += 1
                        record("weighted", p, q, U, V, str(weight), tuple(map(str, poly.all_coeffs())))
                        weighted += 1
                    A = sp.Poly(sum(c*t**(j-lo) for j, c in first.items()), t)
                    even_shift = 2*(V//p+1)
                    S = sp.Poly(sum(c*t**(j+even_shift//2) for j, c in first.items()), t)
                    B = sp.Poly(sum(c*t**(j+even_shift) for j, c in full.items()), t)
                    D = B-S*S
                    require(even_shift % 2 == 0 and min(j+even_shift for j in full) >= 0, "even Laurent clearing preserves the raw negative-root sign")
                    if A.degree():
                        shared, strict, cert = nonpositive_at_roots(A, B)
                        common += shared
                        negative += strict
                        record("one-factor sign", p, q, U, V, shared, cert)
                        for factor, multiplicity in A.factor_list()[1]:
                            if multiplicity >= 2:
                                require(D.rem(factor**(2*multiplicity-2)).is_zero, "multiple first roots retain the proved defect order")
                    cases += 1
    print(f"network universe: {cases} cases; {endpoints} literal endpoints; {weighted} weighted rows including {zero_rows} zero rows")
    print(f"independent Toeplitz controls: {minors} minors of sizes2,3 from index sets0..3")
    print(f"one-factor roots: {negative} strictly negative doubled values; {common} shared-root zeros; all multiplicities retained")


def degeneration_controls():
    rows = 0
    for r in range(1, 7):
        S = (t+1)**r
        for weight in (sp.Rational(1, 4), sp.Integer(1), sp.Integer(9)):
            for D in (sp.Integer(0), -(t+1)**(2*r-2), (t+1)**(2*r-1)):
                P = sp.Poly(D+weight*S*S, t)
                require(sum(mult for _, mult in P.intervals()) == P.degree(), "sharp allowed multiplicity degeneration is real-rooted")
                rows += 1
            wrong_sign = sp.Poly((t+1)**(2*r-2)+weight*S*S, t)
            require(sum(mult for _, mult in wrong_sign.intervals()) < wrong_sign.degree(), "wrong quadratic leading sign has nonreal roots")
            rows += 1
            if r >= 2:
                too_low = sp.Poly(-(t+1)**(2*r-3)+weight*S*S, t)
                require(sum(mult for _, mult in too_low.intervals()) < too_low.degree(), "losing three multiplicities creates a nonreal cubic pair")
                rows += 1
    print(f"square-pencil degeneration: {rows} exact allowed/boundary/hostile rows, r=1..6")


def joint_hostile():
    Fpoly = 1+t
    defect = 2*t*(2*t+3)**2
    actual = sp.Poly(Fpoly**2+defect, t)
    expected = sp.Poly(1+20*t+25*t*t+8*t**3, t)
    require(actual == expected, "exact cubic actual analogue")
    discriminant = sp.factor(sp.discriminant(defect+w*Fpoly**2, t))
    require(discriminant == 4*w*(2*w*w+9*w+432), "uniform nonnegative-pencil discriminant identity")
    virtual = sp.Poly(1+4*t+t*t, t)
    joint_actual = sp.Poly(sum(actual.nth(j)**2*t**j for j in range(4)), t)
    joint_defect = joint_actual-virtual
    require(virtual.eval(-1) == -2 and joint_actual.eval(-1) == 162 and joint_defect.eval(-1) == 164, "joint sign transfer is exactly refuted")
    require(joint_defect == sp.Poly(396*t+624*t*t+64*t**3, t), "joint defect retains nonnegative coefficients")
    require(real_nonpositive(joint_defect) and real_nonpositive(joint_actual), "joint real-rootedness does not rescue the sign")
    require(actual.eval(-1) == -2, "strict individual first-root sign also survives the hostile")
    quadratic_cases = 0
    for a0 in range(1, 4):
        for a2 in range(1, 4):
            for a1 in range(a0+a2, a0+a2+3):
                for b0 in range(1, 4):
                    for b2 in range(1, 4):
                        for b1 in range(b0+b2, b0+b2+3):
                            value = a0*b0-a1*b1+a2*b2
                            require(value <= -a0*b2-a2*b0 <= -2, "quadratic abstraction still forces joint negativity")
                            quadratic_cases += 1
    print(f"joint hostile: first root-1, virtual-2, actual162, defect164; both individual and joint PF predicates survive")
    print(f"narrow degree boundary: {quadratic_cases} exact quadratic controls; failure first exhibited in cubic degree")
    record("joint hostile", actual.all_coeffs(), str(discriminant), joint_defect.all_coeffs())


def main():
    print("WEIGHTED MIDPOINT DELETION AND THE JOINT-TRANSPORT BOUNDARY")
    print("scope=PROVED full-support weighted PF and all-root one-factor sign; abstract joint implication REFUTED")
    network_controls()
    degeneration_controls()
    joint_hostile()
    print("OPEN actual binomial joint signed transport; the cubic hostile is not a binomial-row counterexample")
    print("trace_sha256="+TRACE.hexdigest())
    print(f"PASS explicit_gates={GATES}")


if __name__ == "__main__":
    main()
