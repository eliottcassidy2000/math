#!/usr/bin/env python3
"""Exact audit of the full-residue two-jet Smith wall.

Universe: integer nodes a+p**e*u_i, with u_i running once through F_p,
e>=1. The unrestricted assertions are proved in the accompanying note;
this program checks finite controls by independent DVR and all-minor paths.
No assertions are used, so -O must produce identical output.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import gcd
import json
import sys


GATES = 0


def check(condition, label):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def valuation(x, p):
    x = Fraction(x)
    check(x != 0, "nonzero valuation")
    n, d = abs(x.numerator), x.denominator
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    while d % p == 0:
        d //= p
        v -= 1
    return v


def matrix(nodes):
    n = 2 * len(nodes)
    return [row for x in nodes for row in (
        [x**q for q in range(n)],
        [q*x**(q-1) if q else 0 for q in range(n)],
    )]


def dvr_matrix(a, p):
    a = [[Fraction(x) for x in row] for row in a]
    n = len(a)
    result = []
    for j in range(n):
        r, c = min(((r, c) for r in range(j, n) for c in range(j, n)
                    if a[r][c]), key=lambda rc: valuation(a[rc[0]][rc[1]], p))
        a[j], a[r] = a[r], a[j]
        for row in a:
            row[j], row[c] = row[c], row[j]
        pivot = a[j][j]
        result.append(valuation(pivot, p))
        for r in range(j+1, n):
            if a[r][j]:
                m = a[r][j] / pivot
                check(valuation(m, p) >= 0, "DVR row integrality")
                a[r] = [x-m*y for x, y in zip(a[r], a[j])]
        for c in range(j+1, n):
            if a[j][c]:
                m = a[j][c] / pivot
                check(valuation(m, p) >= 0, "DVR column integrality")
                for r in range(j, n):
                    a[r][c] -= m*a[r][j]
    check(result == sorted(result), "Smith divisibility order")
    return tuple(result)


def dvr_profile(nodes, p):
    result = dvr_matrix(matrix(nodes), p)
    check(sum(result) == 4*sum(valuation(y-x, p) for i, x in enumerate(nodes)
                              for y in nodes[i+1:]), "determinant conservation")
    return result


def predicted(s, p, e=1):
    check(1 <= s <= p and e >= 1, "single-scale domain")
    if s < p:
        return (0, 0) + tuple(e*j for j in range(1, s)) + tuple(
            e*j for j in range(s+1, 2*s))
    return ((0, 0) + tuple(e*j for j in range(1, p-1))
            + ((p-1)*e+1,) + tuple(e*j for j in range(p+1, 2*p-1))
            + ((2*p-1)*e-1,))


def determinant(a):
    a = [row[:] for row in a]
    n = len(a)
    if n == 0:
        return 1
    sign, previous = 1, 1
    for j in range(n-1):
        r = next((r for r in range(j, n) if a[r][j]), None)
        if r is None:
            return 0
        if r != j:
            a[r], a[j] = a[j], a[r]
            sign *= -1
        pivot = a[j][j]
        for r in range(j+1, n):
            for c in range(j+1, n):
                top = pivot*a[r][c]-a[r][j]*a[j][c]
                check(top % previous == 0, "Bareiss exact division")
                a[r][c] = top//previous
        for r in range(j+1, n):
            a[r][j] = 0
        previous = pivot
    return sign*a[-1][-1]


def all_minor_profile(nodes, p):
    a, divisors = matrix(nodes), [1]
    n = len(a)
    for h in range(1, n+1):
        g = 0
        for rr in combinations(range(n), h):
            for cc in combinations(range(n), h):
                minor = determinant([[a[r][c] for c in cc] for r in rr])
                g = gcd(g, abs(minor))
        check(g != 0, "nonzero determinantal divisor")
        check(g % divisors[-1] == 0, "determinantal divisor chain")
        divisors.append(g)
    return tuple(valuation(b//a, p) for a, b in zip(divisors, divisors[1:]))


def rank(a, p):
    a = [[x % p for x in row] for row in a]
    if not a:
        return 0
    j = 0
    for c in range(len(a[0])):
        r = next((r for r in range(j, len(a)) if a[r][c]), None)
        if r is None:
            continue
        a[j], a[r] = a[r], a[j]
        inv = pow(a[j][c], -1, p)
        a[j] = [(x*inv) % p for x in a[j]]
        for r in range(j+1, len(a)):
            m = a[r][c]
            a[r] = [(x-m*y) % p for x, y in zip(a[r], a[j])]
        j += 1
        if j == len(a):
            break
    return j


def wall_witnesses(p, units):
    certificates = []
    for h in range(p+1, 2*p):
        values = [[u**q for q in range(h)] for u in units]
        derivatives = [[q*u**(q-1) if q else 0 for q in range(h)]
                       for u in units]
        trace = [sum(row[q] for row in derivatives) for q in range(h)]
        check(all(t % p == 0 for t in trace), "derivative trace divisibility")
        saturated = derivatives[:-1] + [[t//p for t in trace]]
        check(rank(derivatives, p) == p-1, "missing derivative direction")
        check(rank(saturated, p) == p, "one divided trace repairs derivative rank")
        check(rank(derivatives+values, p) == h, "full confluent rank")
        selected = []
        work = saturated[:]
        for i, row in enumerate(values):
            if rank(work+[row], p) > len(work):
                work.append(row)
                selected.append(i)
            if len(work) == h:
                break
        check(len(work) == h, "basis extension with value rows")
        unscaled = derivatives+[values[i] for i in selected]
        v = valuation(determinant(unscaled), p)
        check(v == 1, "literal rank-h minor has exact residual valuation one")
        certificates.append((h, tuple(selected), v))
    return certificates


def consecutive_predicted(n, p):
    check(1 <= n <= p*p, "consecutive quadratic domain")
    return tuple(sorted(x for c in range(p) if c < n
                        for x in predicted((n-1-c)//p+1, p)))


def main():
    check(len(sys.argv) == 1, "no arguments")
    records, direct, minors, witnesses = [], [], [], []
    for p in (2, 3, 5, 7, 11):
        for e in (1, 2, 3):
            for variant in (0, 1):
                units = tuple(i if variant == 0 else i+p*(i*i-2*i+3)
                              for i in range(p))
                a = 0 if variant == 0 else 2-3*p
                nodes = tuple(a+p**e*u for u in units)
                actual = dvr_profile(nodes, p)
                check(actual == predicted(p, p, e), ("wall profile", p, e, variant))
                records.append((p, e, variant, actual))
                if e == 1:
                    witnesses.append((p, variant, wall_witnesses(p, units)))
        for s in (1, p-1):
            actual = dvr_profile(tuple(p*i for i in range(s)), p)
            check(actual == predicted(s, p), "inherited s<p control")
    for p in (2, 3, 5):
        for e in (1, 2):
            if p == 5:
                continue
            units = tuple(i+p*(i*i+1) for i in range(p))
            nodes = tuple(1+p**e*u for u in units)
            actual = all_minor_profile(nodes, p)
            check(actual == predicted(p, p, e), "independent all-minor wall")
            minors.append((p, e, actual))
    for p in (2, 3, 5):
        for n in range(p*(p-1)+1, p*p+1):
            actual = dvr_profile(tuple(range(n)), p)
            check(actual == consecutive_predicted(n, p), "new consecutive band")
            direct.append((p, n, actual))
    formula_count = 0
    for p in (2, 3, 5, 7, 11, 13, 17, 19):
        for n in range(1, p*p+1):
            profile = consecutive_predicted(n, p)
            expected = 4*sum((n-d)*valuation(d, p) for d in range(1, n))
            check(sum(profile) == expected, "consecutive index mass")
            check(len(profile) == 2*n, "consecutive rank")
            formula_count += 1
    hostile = dvr_profile((0, 2, 4), 2)
    check(hostile != (0, 0, 1, 2, 4, 5), "repeated-residue sidecar required")
    tensor_cases = []
    for left_e, right_e in ((1, 1), (1, 2), (2, 3)):
        left, right = matrix((0, 2**left_e)), matrix((0, 2**right_e))
        tensor = [[left[i][j]*right[k][ell] for j in range(4) for ell in range(4)]
                  for i in range(4) for k in range(4)]
        actual = dvr_matrix(tensor, 2)
        expected = tuple(sorted(a+b for a in predicted(2, 2, left_e)
                                for b in predicted(2, 2, right_e)))
        check(actual == expected, "tensor-grid full Smith partition")
        check(max(actual) == 3*(left_e+right_e)-2, "tensor-grid sharp precision")
        tensor_cases.append((left_e, right_e, actual))
    kernel_counts = []
    small = matrix((0, 2))
    for n in (1, 2, 3):
        modulus = 2**n
        count = sum(all(sum(a*b for a, b in zip(row, coeffs)) % modulus == 0
                        for row in small)
                    for coeffs in product(range(modulus), repeat=4))
        expected = 2**sum(min(n, a) for a in predicted(2, 2))
        check(count == expected, "complete finite-precision kernel count")
        kernel_counts.append((n, count))
    triple_cases = []
    for e in range(9):
        actual = dvr_profile((0, 2**e, 2**(e+1)), 2)
        delta = (0, min(e+1, 2*e), min(3*e+2, 4*e), 7*e+2, 12*e+4)
        expected = (0, 0)+tuple(b-a for a, b in zip(delta, delta[1:]))
        check(actual == expected, "dyadic equally spaced triple")
        triple_cases.append((e, actual))
    check(triple_cases[3][1][:4] != predicted(2, 2, 3),
          "adjoining same-residue node can change old Smith prefix")
    rows = ((1, 0), (1, 1), (2, 0), (2, 1))
    tariff_table = {}
    for h in range(1, 5):
        costs = {}
        for rr in combinations(range(4), h):
            for cc in combinations(range(2, 6), h):
                residual = [[q**r*u**(q-r) for q in cc]
                            for u, r in (rows[i] for i in rr)]
                det = determinant(residual)
                if det:
                    weight = sum(cc)-sum(rows[i][1] for i in rr)
                    costs[weight] = min(costs.get(weight, 10**9), valuation(det, 2))
        tariff_table[h] = costs
    check(tariff_table == {
        1: {1: 1, 2: 0, 3: 0, 4: 0, 5: 0},
        2: {3: 2, 4: 0, 5: 1, 6: 0, 7: 1, 8: 0, 9: 4},
        3: {7: 2, 8: 2, 9: 2, 10: 2, 11: 3},
        4: {12: 4},
    }, "all 69 symbolic triple residual minors")
    semantic = json.dumps((records, direct, minors, witnesses, hostile,
                           tensor_cases, kernel_counts, triple_cases, tariff_table),
                          separators=(",", ":"))
    print("CONFLUENT_TWOJET_PRIME_WALL_SYNTHESIS_SEP05")
    print("scope=p prime; e>=1; p integer nodes a+p^e*u_i; complete residue system u_i mod p")
    print("profile=(0,0,e,...,(p-2)e,(p-1)e+1,(p+1)e,...,(2p-2)e,(2p-1)e-1)")
    print("middle_determinantal_valuation=old_weighted_tariff+1")
    print("mechanism=divide_sum_of_derivative_rows_by_p;one_saturation_step")
    print(f"full_wall_DVR_cases={len(records)}")
    print(f"independent_all_minor_cases={len(minors)}")
    print(f"literal_divided_trace_witnesses={sum(len(w) for _, _, w in witnesses)}")
    print(f"new_band_direct_consecutive_cases={len(direct)}")
    print(f"consecutive_formula_rows_through_prime_19={formula_count}")
    print(f"direct_tensor_grid_cases={len(tensor_cases)}")
    print(f"complete_finite_precision_kernel_counts={kernel_counts}")
    print(f"dyadic_triple_DVR_cases={len(triple_cases)}")
    print(f"dyadic_triple_at_e3={triple_cases[3][1]}")
    print("symbolic_triple_residual_minors=69")
    for p, e, variant, profile in records:
        if variant == 0 and p <= 5:
            print(f"wall(p={p},e={e})={profile}")
    print(f"hostile_repeated_residue(nodes=(0,2,4),p=2)={hostile}")
    print(f"semantic_sha256={sha256(semantic.encode()).hexdigest()}")
    print(f"gates={GATES}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
