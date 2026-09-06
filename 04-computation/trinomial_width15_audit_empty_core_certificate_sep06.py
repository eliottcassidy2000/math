#!/usr/bin/env python3
"""Independent standard-library referee for the lower-carry width-15 family.

Nine rational evaluations certify polynomial identities of degree at most
eight. The bounds follow after removing (2g)_10=32K and explicitly cancelling
the two factors in the inverse-tau term; see the companion audit note.
No symbolic algebra package or concurrent producer is imported.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import comb, factorial, gcd, prod

CHECKS = 0
CPOLY = [-158894736, 91449662, -17466623, 1106459]
DPOLY = [-3266244982, 1990079478, -404127357, 27352253]
JPOLY = [-12588295720, 7685968926, -1564109559, 106089635]
HPOLY = [-853264618035013184, 1146316900083491136, -659512619763975524,
         210638390555236696, -40333487193583421, 4630220682540559,
         -295063248863031, 8051838961249]


def need(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def evaluate(polynomial, g):
    return sum(coefficient*g**j for j, coefficient in enumerate(polynomial))


def falling(g, length):
    return prod(g-j for j in range(length))


def normalized_remainder(g):
    """R/K from the canonical row, with tau^-1 retained explicitly."""
    u = Q(g-5)
    # q0/[K*u*(u-1)] has the following polynomial cancellation.
    q0_cancelled = Q(128*(2*g-11)*(2*g-13)*(2*g-14), factorial(15))
    constant, linear = -20*u*q0_cancelled, -6*q0_cancelled
    power = (Q(1), Q(0))  # tau^(j-1), starting with j=1.
    for j in range(1, 6):
        coefficient = Q(32*prod(2*g-k for k in range(10, 15-j)),
                        factorial(15-3*j)*factorial(2*j))
        constant += coefficient*power[0]
        linear += coefficient*power[1]
        a, b = power
        power = -b*u*(u-1)/6, a-Q(10, 3)*u*b
    return constant, linear


def literal_rows(g):
    a, b, c = 15, 2*g-15, 3*g-15
    state = {(0, 0): 1}
    answer = []
    for _ in range(2*g):
        following = {}
        for (charge, z), coefficient in state.items():
            for step, dz in ((-a, 0), (b, 0), (c, 1)):
                key = charge+step, z+dz
                following[key] = following.get(key, 0)+coefficient
        state = following
        answer.append({z: coefficient for (charge, z), coefficient in state.items()
                       if charge == 0})
    return answer


def shift(polynomial, offset):
    return [sum(polynomial[j]*comb(j, k)*offset**(j-k)
                for j in range(k, len(polynomial))) for k in range(len(polynomial))]


def main():
    manifest = []
    for g in range(8, 17):
        u = Q(g-5)
        constant, linear = normalized_remainder(g)
        need(constant == u*evaluate(CPOLY, g)/12259447200, ('constant identity', g))
        need(linear == Q(evaluate(DPOLY, g), 15324309000), ('linear identity', g))
        trace = 2*constant-Q(10, 3)*u*linear
        norm = constant*constant-Q(10, 3)*u*constant*linear+u*(u-1)*linear*linear/6
        need(trace == -u*evaluate(JPOLY, g)/18389170800, ('trace identity', g))
        need(norm == u*evaluate(HPOLY, g)/3757351141239696000000, ('norm identity', g))
        need(trace < 0 and norm > 0, ('positive identity controls', g))
        manifest.append((g, constant, linear, trace, norm))
    print('IDENTITY_AUDIT nine distinct parameters g=8,...,16')
    print('degree_bounds normalized_constant<=4, linear<=3, trace<=4, norm<=8')
    expected_shifts = {
        'J': [106089635, 982041681, 3029425902, 3114337032],
        'H': [8051838961249, 155839732966913, 1288856301033727,
              5903575385111259, 16172002313744184, 26489396415060828,
              24017452577936640, 9296533901880000]}
    for name, polynomial in (('J', JPOLY), ('H', HPOLY)):
        values = list(reversed(shift(polynomial, 8)))
        need(values == expected_shifts[name], 'independent binomial shift '+name)
        need(all(value > 0 for value in values), 'strict shifted positivity '+name)
        print('shifted_positive', name, values)
    for g in (8, 11, 13, 14):
        need(gcd(g, 15) == 1, 'primitive first-return control')
        native = literal_rows(g)
        need(all(not row for row in native[:g-1]), 'all earlier literal rows empty')
        first = [Q(falling(g, 5), 720)*(g-5)*(g-6),
                 Q(falling(g, 5), 720)*20*(g-5), Q(falling(g, 5), 720)*6]
        second = [Q(falling(2*g, 15-j), factorial(15-3*j)*factorial(2*j))
                  for j in range(6)]
        need(native[g-1] == {1+2*j: coefficient for j, coefficient in enumerate(first)},
             'full literal first factorial row')
        need(native[2*g-1] == {2*j: coefficient for j, coefficient in enumerate(second)},
             'full literal doubled row including lower carry')
        cbar, dbar = normalized_remainder(g)
        k = Q(falling(2*g, 10), 32)
        need(k > 0, 'positive common factor')
        print('native_control', (-15, 2*g-15, 3*g-15), 'first_mass', g,
              'first_channels', 3, 'second_channels', 6,
              'anchored_remainder', (k*cbar, k*dbar))
    print('semantic_sha256', sha256(repr(manifest).encode()).hexdigest())
    print('explicit_checks', CHECKS)
    print('PASS: independent all-parameter polynomial identities and positive shifts')
    print('SCOPE: (-15,2g-15,3g-15), integral g>=8,gcd(g,15)=1; not all width15 supports')


if __name__ == '__main__':
    main()
