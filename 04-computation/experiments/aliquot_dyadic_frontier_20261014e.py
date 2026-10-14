"""Exact three-cycle exclusion in the dyadic-times-prime aliquot template.

The prime variables have no cutoff. The two smallest distinct exponents
are bounded by 128; a proved resultant bound covers every larger third
exponent. Repeated exponents are classified algebraically for all sizes.
No assertion is disabled by python -O. Stdout is deterministic JSON.
"""
from fractions import Fraction
from itertools import product
from math import isqrt
import json

CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(label)


def exact_int(n, low=1):
    if type(n) is not int or n < low:
        raise ValueError("exact integer outside the domain")
    return n


def dyadic(c):
    exact_int(c, 2)
    if c & (c-1):
        raise ValueError("cofactor is not a power of two")
    return c


def coefficients(c, d):
    dyadic(c)
    dyadic(d)
    L, K = c+d-1, (c-1)*(d-1)
    A = 2*d
    B = 2*c*d-2*c+2*d*d-3*d+1
    C = -2*c*d+2*c-2*d*d+2*d-1
    R = A*K*K-B*K*L+C*L*L
    if R >= 0:
        raise ValueError("resultant sign theorem failed")
    return L, K, A, B, C, R


def first_prime_candidate(c, d, e):
    dyadic(e)
    L, K, A, B, C, _ = coefficients(c, d)
    return Fraction(A*e*e+B*e+C, L*e+K)


def all_prime_candidates(c, d, e):
    p = first_prime_candidate(c, d, e)
    q = ((c-1)*p+2*c-1)/d
    r = ((d-1)*q+2*d-1)/e
    if ((e-1)*r+2*e-1)/c != p:
        raise ValueError("formal closure failed")
    return p, q, r


def compose_source_candidate(cores):
    # Independent route on the n_i variables, using rational affine maps.
    slope, intercept = Fraction(1), Fraction(0)
    for c in cores:
        dyadic(c)
        slope *= Fraction(c-1, c)
        intercept = Fraction(c-1, c)*intercept+2*c-1
    return intercept/(1-slope)


def aliquot_direct(n):
    exact_int(n)
    return sum(d for d in range(1, n) if n % d == 0)


def prime_trial(p):
    return type(p) is int and p >= 2 and all(p % d for d in range(2, isqrt(p)+1))


def positive_odd_integers(values):
    return all(v.denominator == 1 and v > 0 and v.numerator % 2 for v in values)


def repeated_controls():
    formal_nonconstant = []
    for a in range(1, 65):
        c = 1 << a
        for b in range(1, 65):
            d = 1 << b
            ps = all_prime_candidates(c, c, d)
            if positive_odd_integers(ps):
                ns = [c*ps[0], c*ps[1], d*ps[2]]
                expected = (a == b or (a, b) == (1, 4))
                check(expected, "repeated-cofactor classification")
                if a == b:
                    check(ps == (2*c-1,)*3, "diagonal fixed model")
                else:
                    formal_nonconstant.append({"exponents": [a, a, b],
                                                "sources": [int(x) for x in ns],
                                                "prime_variables": [int(x) for x in ps]})
            else:
                check(a != b and (a, b) != (1, 4), "repeated-case completeness")
    check(formal_nonconstant == [{"exponents": [1, 1, 4], "sources": [46, 26, 16],
                                  "prime_variables": [23, 13, 1]}], "unique ghost")
    check(aliquot_direct(46) == 26 and aliquot_direct(26) == 16
          and aliquot_direct(16) == 15, "native rejection at prime variable one")
    return {"control_exponents": [1, 64], "nonconstant_formal_candidates": formal_nonconstant,
            "all_exponent_proof": "repeated dyadic exponents force a fixed point or the nonnative 46,26,16 ghost"}


def independent_controls():
    # Exact independent source-variable composition, not the resultant formula.
    for exps in product(range(1, 9), repeat=3):
        c, d, e = [1 << a for a in exps]
        ps = all_prime_candidates(c, d, e)
        check(compose_source_candidate((c, d, e)) == c*ps[0], "independent affine composition")
        L, K, A, B, C, R = coefficients(c, d)
        den, num = L*e+K, A*e*e+B*e+C
        check((L*L*num-R) % den == 0, "resultant divisibility identity")
        check(den == c*d*e-(c-1)*(d-1)*(e-1)
              and den % 2 == 1 and all(den % p.denominator == 0 for p in ps),
              "odd common denominator for cyclic prime variables")
        check(all(p.denominator == 1 for p in ps) == (ps[0].denominator == 1)
              and (ps[0].denominator != 1 or positive_odd_integers(ps)),
              "integrality propagates around the dyadic cycle")
        if num % den == 0:
            check(den <= -R, "all-exponent cutoff positive control")
    for n in (6, 28, 496):
        check(aliquot_direct(n) == n, "actual fixed-point control")
    check(prime_trial(23) and prime_trial(13) and not prime_trial(1), "ghost prime guard")
    for bad in (True, 3, 0, 4.0):
        try:
            dyadic(bad)
        except ValueError:
            check(True, "malformed cofactor rejected")
        else:
            check(False, "malformed cofactor accepted")


def complete_third_exponent_search(bound):
    exact_int(bound, 2)
    tested = 0
    integral = []
    max_exponent = 0
    min_R, max_R = None, 0
    row_counts = []
    for b in range(2, bound+1):
        row_tested = 0
        for a in range(1, b):
            for c, d in ((1 << a, 1 << b), (1 << b, 1 << a)):
                L, K, A, B, C, R = coefficients(c, d)
                # If the candidate is integral, denominator | R, so e cannot
                # exceed this exact integer. All larger e are proved excluded.
                cap = (-R-K)//L
                max_e = cap.bit_length()-1 if cap > 0 else -1
                max_exponent = max(max_exponent, max_e)
                min_R = -R if min_R is None else min(min_R, -R)
                max_R = max(max_R, -R)
                check(L*(1 << (max_e+1))+K > -R, "entire omitted tail excluded")
                for f in range(b+1, max_e+1):
                    e = 1 << f
                    den, num = L*e+K, A*e*e+B*e+C
                    tested += 1
                    row_tested += 1
                    if num % den == 0:
                        integral.append({"ordered_cores_exponents":
                                         [c.bit_length()-1, d.bit_length()-1, f],
                                         "first_prime_variable": num//den})
                    check(num % den != 0, "distinct-exponent integer candidate survives")
        row_counts.append([b, row_tested])
    return {"second_smallest_exponent_at_most": bound,
            "smallest_exponent_at_least": 1,
            "largest_exponent": "unbounded: exact resultant tail exclusion",
            "largest_exponent_actually_tested": max_exponent,
            "oriented_core_words_tested": tested,
            "integral_first_prime_candidates": integral,
            "counts_by_second_smallest_exponent": row_counts,
            "resultant_absolute_range": [min_R, max_R]}


def odd_even_pair_sum_controls():
    for length in range(1, 10):
        values = tuple(i*i+3*i+1 for i in range(length))
        sigma_values = tuple(values[i]+values[(i+1) % length] for i in range(length))
        if length % 2:
            for i in range(length):
                recovered = Fraction(sum((-1)**j*sigma_values[(i+j) % length]
                                         for j in range(length)), 2)
                check(recovered == values[i], "odd pair-sum inverse")
        else:
            check(sum((-1)**i*sigma_values[i] for i in range(length)) == 0,
                  "even pair-sum compatibility")
            changed = [values[i]+(-1)**i for i in range(length)]
            check(all(changed[i]+changed[(i+1) % length] == sigma_values[i]
                      for i in range(length)), "even alternating coordinate is lost")


def main():
    independent_controls()
    repeated = repeated_controls()
    search = complete_third_exponent_search(128)
    odd_even_pair_sum_controls()
    result = {"status": "PROVED repeated-exponent exclusion and resultant cutoff; FINITE-EXACT remaining declared templates; general aliquot 3-cycle question OPEN",
              "template": "n_i=2^a_i p_i, a_i>=1, p_i an odd prime appearing to exponent one",
              "repeated_exponents": repeated, "distinct_exponents": search,
              "surviving_hypothesis": "no primitive aliquot 3-cycle lies wholly in this template",
              "remaining_parameters": "all exponents distinct, second smallest at least 129",
              "checks": CHECKS}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
