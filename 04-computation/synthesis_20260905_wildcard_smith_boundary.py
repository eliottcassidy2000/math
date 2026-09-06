#!/usr/bin/env python3
"""Exact audit of the complete-residue two-jet Smith boundary.

Universe: primes 2,3,5,7,11,13; e=1,2,3; four complete residue
systems, including irregular positive/negative lifts and translated/reversed
nodes. Primary: exact rational DVR elimination. Independent path: modular
DVR elimination with inversions only of units, modulus exponent exceeding
the independently known total determinant valuation. All minors are also
enumerated in six rank<=6 cases by fraction-free determinants.

Controls: inherited s<p profiles; false uncorrected s=p profiles; the e=2
hostile against charging the geometric scale e for the missing line;
direct consecutive matrices in the newly closed final quadratic band.
No project modules, floating point, or Python assert statements are used.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import comb, factorial
import json
import sys


sys.stdout.reconfigure(newline="\n")
GATES = 0


def require(condition, label):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def vp(value, prime):
    require(value != 0, "valuation of nonzero input")
    value = Fraction(value)
    numerator, denominator = abs(value.numerator), value.denominator
    result = 0
    while numerator % prime == 0:
        numerator //= prime
        result += 1
    while denominator % prime == 0:
        denominator //= prime
        result -= 1
    return result


def matrix(nodes):
    n = 2 * len(nodes)
    return [[comb(q, r) * x ** (q-r) if q >= r else 0
             for q in range(n)] for x in nodes for r in (0, 1)]


def dvr_smith(entries, prime):
    a = [[Fraction(x) for x in row] for row in entries]
    size = len(a)
    result = []
    for k in range(size):
        i, j = min(((i, j) for i in range(k, size)
                    for j in range(k, size) if a[i][j]),
                   key=lambda pair: vp(a[pair[0]][pair[1]], prime))
        a[k], a[i] = a[i], a[k]
        for row in a:
            row[k], row[j] = row[j], row[k]
        pivot = a[k][k]
        result.append(vp(pivot, prime))
        for i in range(k+1, size):
            factor = a[i][k] / pivot
            require(not factor or vp(factor, prime) >= 0, "integral row operation")
            for j in range(k, size):
                a[i][j] -= factor*a[k][j]
        # Clearing the pivot row by column operations does not change the
        # remaining block, because the pivot column below it is now zero.
        for i in range(k+1, size):
            require(a[i][k] == 0, "pivot column cleared")
    require(result == sorted(result), "ordered DVR exponents")
    return tuple(result)


def modular_smith(nodes, prime):
    # The determinant is product_(i<j)(x_j-x_i)^4, so this precision
    # exceeds every possible invariant exponent without assuming the claim.
    total = 4 * sum(vp(y-x, prime) for i, x in enumerate(nodes)
                    for y in nodes[i+1:])
    exponent = total+1
    modulus = prime ** exponent
    size = 2*len(nodes)
    a = [[(pow(x, q, modulus) if r == 0 else
           (q*pow(x, q-1, modulus) % modulus if q else 0))
          for q in range(size)] for x in nodes for r in (0, 1)]

    def valuation(x):
        if x == 0:
            return exponent
        v = 0
        while x % prime == 0:
            v += 1
            x //= prime
        return v

    result = []
    for k in range(size):
        i, j = min(((i, j) for i in range(k, size)
                    for j in range(k, size)),
                   key=lambda pair: valuation(a[pair[0]][pair[1]]))
        a[k], a[i] = a[i], a[k]
        for row in a:
            row[k], row[j] = row[j], row[k]
        v = valuation(a[k][k])
        require(v < exponent, "independent modular precision sufficient")
        result.append(v)
        power = prime ** v
        residual_modulus = modulus // power
        inverse = pow(a[k][k] // power, -1, residual_modulus)
        for i in range(k+1, size):
            require(a[i][k] % power == 0, "modular pivot divides column")
            factor = ((a[i][k] // power) * inverse) % residual_modulus
            a[i] = [(left-factor*right) % modulus
                    for left, right in zip(a[i], a[k])]
            require(a[i][k] == 0, "modular column cleared")
    require(sum(result) == total, "independent modular determinant mass")
    return tuple(result)


def bareiss(entries):
    a = [list(row) for row in entries]
    n = len(a)
    if n == 0:
        return 1
    sign, previous = 1, 1
    for k in range(n-1):
        pivot_row = next((i for i in range(k, n) if a[i][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = pivot*a[i][j]-a[i][k]*a[k][j]
                require(numerator % previous == 0, "Bareiss exact division")
                a[i][j] = numerator // previous
        for i in range(k+1, n):
            a[i][k] = 0
        previous = pivot
    return sign*a[-1][-1]


def all_minor_divisors(nodes, prime):
    a = matrix(nodes)
    n = len(a)
    deltas = [0]
    for h in range(1, n+1):
        valuations = []
        for rows in combinations(range(n), h):
            for columns in combinations(range(n), h):
                value = bareiss([[a[i][j] for j in columns] for i in rows])
                if value:
                    valuations.append(vp(value, prime))
        require(bool(valuations), "nonzero minor at every rank")
        deltas.append(min(valuations))
    return tuple(deltas)


def old_profile(s, e):
    return (0, 0) + tuple(e*j for j in range(1, s)) + tuple(
        e*j for j in range(s+1, 2*s))


def profile(p, e):
    return ((0, 0) + tuple(e*j for j in range(1, p-1))
            + ((p-1)*e+1,) + tuple(e*j for j in range(p+1, 2*p-1))
            + ((2*p-1)*e-1,))


def divisor_formula(p, e):
    return ((0,) + tuple(e*(h-1)*(h-2)//2 for h in range(1, p+1))
            + tuple(e*(h*(h-1)//2-p)+1 for h in range(p+1, 2*p))
            + (2*e*p*(p-1),))


def rank_mod(entries, prime):
    a = [[x % prime for x in row] for row in entries]
    rows, columns, rank = len(a), len(a[0]), 0
    for column in range(columns):
        pivot = next((i for i in range(rank, rows) if a[i][column]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inverse = pow(a[rank][column], -1, prime)
        a[rank] = [x*inverse % prime for x in a[rank]]
        for i in range(rows):
            if i != rank:
                factor = a[i][column]
                a[i] = [(x-factor*y) % prime for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def consecutive_profile(n, p):
    result = []
    for residue in range(p):
        s = len(range(residue, n, p))
        require(s <= p, "quadratic consecutive universe")
        if s:
            result.extend(profile(p, 1) if s == p else old_profile(s, 1))
    return tuple(sorted(result))


def main():
    require(sys.argv[1:] == [], "no arguments")
    single, controls, residual, minors, consecutive = [], [], [], [], []
    for p in (2, 3, 5, 7, 11, 13):
        coordinates = (
            tuple(range(p)),
            tuple(u+p*(u*u-2*u-3) for u in range(p)),
            tuple(reversed(range(p))),
            tuple((3*u+1) % p + p*(u*u+5) for u in range(p))
            if p != 3 else tuple(u+3*(2*u*u-7) for u in range(p)),
        )
        for e in (1, 2, 3):
            expected = profile(p, e)
            require(expected == tuple(sorted(expected)), "claimed exponents ordered")
            require(sum(expected) == 2*e*p*(p-1), "claimed determinant mass")
            require(expected != old_profile(p, e), "uncorrected boundary hostile")
            running = [0]
            for value in expected:
                running.append(running[-1]+value)
            require(tuple(running) == divisor_formula(p, e), "consecutive differences")
            for variant, units in enumerate(coordinates):
                require(len(set(u % p for u in units)) == p, "complete residue sidecar")
                nodes = tuple((7-variant*p) + p**e*u for u in units)
                actual = dvr_smith(matrix(nodes), p)
                independent = modular_smith(nodes, p)
                require(actual == expected, ("boundary formula", p, e, variant))
                require(independent == actual, ("independent modular replay", p, e, variant))
                single.append((p, e, variant, actual))
        for s in sorted(set((1, max(1, p-1)))):
            nodes = tuple(p*u for u in range(s))
            actual = dvr_smith(matrix(nodes), p)
            require(actual == old_profile(s, 1), "inherited s<p positive control")
            controls.append((p, s, actual))
        for units in coordinates[:2]:
            witness = [[q*u**(q-1) for q in range(1, p+1)] for u in units]
            determinant = bareiss(witness)
            require(vp(determinant, p) == 1, "exact derivative saturation index p")
            vandermonde = 1
            for i, x in enumerate(units):
                for y in units[i+1:]:
                    vandermonde *= y-x
            require(determinant == factorial(p)*vandermonde, "saturation witness formula")
            for h in range(p+1, 2*p):
                derivatives = [[q*u**(q-1) if q else 0 for q in range(h)]
                               for u in units]
                full = [[comb(q, r)*u**(q-r) if q >= r else 0
                         for q in range(h)] for u in units for r in (0, 1)]
                require(rank_mod(derivatives, p) == p-1, "one missing derivative line")
                require(rank_mod(full, p) == h, "Hermite extension rank")
                residual.append((p, h, vp(determinant, p)))
    for p, e, units in ((2, 1, (0, 1)), (2, 2, (0, 1)),
                        (2, 3, (-4, 7)), (3, 1, (0, 1, 2)),
                        (3, 2, (0, 1, 2)), (3, 1, (-3, 7, 14))):
        nodes = tuple(5+p**e*u for u in units)
        actual = all_minor_divisors(nodes, p)
        require(actual == divisor_formula(p, e), "independent all-minor divisors")
        minors.append((p, e, units, actual))
    for p in (2, 3, 5):
        for n in range(p*(p-1), p*p+1):
            expected = consecutive_profile(n, p)
            actual = dvr_smith(matrix(tuple(range(n))), p)
            independent = modular_smith(tuple(range(n)), p)
            require(actual == expected, ("new consecutive band", p, n))
            require(actual == independent, ("independent consecutive replay", p, n))
            consecutive.append((p, n, actual))
    require(profile(2, 2) == (0, 0, 3, 5), "saturation cost is 1, not e")
    report = dict(single_scale=single, inherited_controls=controls,
                  derivative_ranks=residual, all_minor_divisors=minors,
                  consecutive_band=consecutive)
    digest = sha256(json.dumps(report, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    print("PASS: complete-residue two-jet Smith boundary")
    print("single-scale cases:", len(single), "; exact rational and independent modular paths")
    print("inherited s<p controls:", len(controls))
    print("derivative rank/saturation cases:", len(residual))
    print("all-minor determinantal-divisor cases:", len(minors))
    print("direct consecutive final-band cases:", len(consecutive))
    for p in (2, 3, 5, 7):
        print("boundary", p, "e=1", profile(p, 1), "e=2", profile(p, 2))
    for p, e, units, deltas in minors:
        print("all-minor", p, e, units, deltas)
    print("semantic_sha256:", digest)
    print("exact_gates:", GATES)
    print("SCOPE: p nodes at exact common scale; consecutive n<=p^2.")
    print("OPEN: multiscale clusters, n>p^2, higher jets, and all prize problems.")


if __name__ == "__main__":
    main()
