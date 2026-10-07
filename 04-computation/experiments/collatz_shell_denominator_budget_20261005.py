"""Exact coefficient denominators and damped integer shell selectors.

Integer coefficients are not an integrality certificate for a moment
readout.  No actual moment oracle or selected-source ROOT proof is read.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import gcd, lcm, prod
from itertools import product
import json


def need(condition, message):
    if not condition:
        raise ValueError(message)


def natural(value, minimum=0):
    need(type(value) is int and value >= minimum, "exact integer in range required")


def rational(value):
    need(type(value) in (int, F), "exact rational required")
    return F(value)


def shell(d):
    natural(d)
    return F(4*(1 << d), ((1 << d)+1)**2)


def multiply(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return tuple(out)


def evaluate(a, x):
    result = 0
    for coefficient in reversed(a):
        result = result*x+coefficient
    return result


def numerator(n):
    natural(n, 1)
    a = (1,)
    for d in range(1, n+1):
        a = multiply(a, (-(1 << (d+2)), ((1 << d)+1)**2))
    return a


def denominator(n):
    natural(n, 1)
    return prod(((1 << d)-1)**2 for d in range(1, n+1))


def denominator_of(vector):
    result = 1
    for value in vector:
        result = lcm(result, F(value).denominator)
    return result


def constant_exponent(n):
    natural(n, 1)
    return n*(n+5)//2


def damped_budget(n, k):
    natural(n, 1); natural(k)
    need(n % 2 == 1, "odd number of shell zeros required for a minorant")
    a = (0,)*k+numerator(n)
    d = denominator(n)
    bound = (1 << constant_exponent(n))*shell(n+1)**k
    return {"shell_zeros": n, "damping": k, "degree": n+k,
            "integer_coefficients": a, "normalization": d,
            "normalized_coefficients": tuple(F(c, d) for c in a),
            "integer_tail_bound": bound,
            "normalized_tail_bound": bound/d}


def universal_damping(n):
    natural(n, 3)
    need(n % 2 == 1, "odd n at least three required")
    return damped_budget(n, (n+9)//2)


@lru_cache(None)
def cyclotomic_at_two(n):
    natural(n, 1)
    previous = prod(cyclotomic_at_two(d) for d in range(1, n) if n % d == 0)
    value, remainder = divmod((1 << n)-1, previous)
    need(remainder == 0, "cyclotomic division failed")
    return value


def is_prime(p):
    if type(p) is not int or p < 2:
        return False
    d = 2
    while d*d <= p:
        if p % d == 0:
            return False
        d += 1
    return True


def valuation(value, prime):
    natural(value, 1)
    need(is_prime(prime), "prime required")
    result = 0
    while value % prime == 0:
        value //= prime
        result += 1
    return result


def order_two(prime):
    need(is_prime(prime) and prime != 2, "odd prime required")
    value, exponent = 2 % prime, 1
    while value != 1:
        value = 2*value % prime
        exponent += 1
    return exponent


def prime_budget(n, prime):
    natural(n, 1)
    e = order_two(prime)
    depth = valuation((1 << e)-1, prime)
    count = n//e
    factorial_depth, q = 0, count//prime
    while q:
        factorial_depth += q
        q //= prime
    return 2*(count*depth+factorial_depth)


def integrality_gate(readout, clearer=1):
    """Finite gate only; coefficients being integral are not an input substitute."""
    value = rational(readout)
    natural(clearer, 1)
    value *= clearer
    if value.denominator != 1:
        return "nonintegral_readout"
    if value == 0:
        return "zero_readout"
    return "nonzero_integer_not_smaller_than_one"


def main():
    checks = 0

    def check(condition, label):
        nonlocal checks
        need(condition, label)
        checks += 1

    import collatz_refinement_energy_dual_20261005 as old
    for n in range(1, 26):
        a, d = numerator(n), denominator(n)
        fractions = tuple(F(c, d) for c in a)
        check(denominator_of(fractions) == d, "least common coefficient denominator")
        check(gcd(*a) == 1, "primitive numerator")
        check(a[0] == (-1)**n*(1 << constant_exponent(n)), "constant coefficient")
        check(evaluate(a, 1) == d, "target normalization")
        for j in range(1, n+1):
            check(evaluate(a, shell(j)) == 0, "rational shell zeros")
        check(prod(cyclotomic_at_two(r)**(2*(n//r)) for r in range(1, n+1)) == d,
              "cyclotomic multiplicities")
        check(F(81, 1024)*(1 << (n*(n+1))) <= d < (1 << (n*(n+1))),
              "quadratic bit-height bounds")
        if n % 2:
            check(fractions == old.minorant(n), "independent rational product")
            check(d*old.error_bound(n) == (1 << constant_exponent(n)),
                  "denominator-cleared undamped error")
            check(sum(map(abs, a)) < 512*d, "norm versus height distinction")

    primes = tuple(p for p in range(3, 100, 2) if is_prime(p))
    for n in range(1, 51):
        for p in primes:
            check(prime_budget(n, p) == valuation(denominator(n), p), "prime clock budget")
    check(cyclotomic_at_two(2) == cyclotomic_at_two(6) == 3,
          "cyclotomic values need not be coprime")

    # Arbitrary extra degree does not reduce the normalization denominator.
    # All quadratics R with coefficients u/2,v/3,1-u/2-v/3 have R(1)=1.
    for n in (1, 3, 5, 7, 9):
        q = tuple(F(c, denominator(n)) for c in numerator(n))
        for u, v in product(range(-2, 3), repeat=2):
            extension = multiply(q, (F(u, 2), F(v, 3), 1-F(u, 2)-F(v, 3)))
            check(evaluate(extension, 1) == 1, "extended normalization")
            check(denominator_of(extension) % denominator(n) == 0,
                  "unavoidable denominator for added degree")
        # Elementary integer shears, with an integer inverse, preserve the
        # coefficient lattice.  A rational basis can relocate its denominator.
        for shear in (-3, -1, 1, 3):
            b = list(q)
            for i in range(len(b)-1):
                b[i] += shear*b[i+1]
            check(denominator_of(b) == denominator(n), "unimodular shear")

    damped_rows = []
    for n in range(3, 26, 2):
        budget = universal_damping(n)
        a, k, d = budget["integer_coefficients"], budget["damping"], budget["normalization"]
        upper = F(1, 1 << ((3*n-9)//2))
        check(budget["integer_tail_bound"] < upper <= 1, "small integer-coefficient envelope")
        check(denominator_of(budget["normalized_coefficients"]) == d,
              "damping retains the optimal denominator")
        check(sum(map(abs, budget["normalized_coefficients"])) == old.coefficient_norm(n),
              "damping preserves direct coefficient norm")
        for j in range(n+1, n+21):
            value = evaluate(a, shell(j))
            check(-budget["integer_tail_bound"] <= value < 0,
                  "nonzero small rational tail, not an integer")
        if n in (3, 5, 9, 15, 25):
            damped_rows.append({"N": n, "k": k, "degree": n+k,
                                "denominator_bits": d.bit_length(),
                                "undamped_cleared_error_bits": constant_exponent(n)+1,
                                "damped_bound_below_power_two": -(3*n-9)//2})

    # Three separate gates: a nonzero small rational form, a zero integral
    # form, and a nonzero integer that cannot have magnitude below one.
    hostile = evaluate((0, 0, -8, 9), shell(2))
    check(hostile == -F(14336, 15625) and 0 < abs(hostile) < 1,
          "integer coefficients do not imply integer readout")
    check(integrality_gate(hostile) == "nonintegral_readout", "readout lattice gate")
    check(integrality_gate(hostile, 15625) == "nonzero_integer_not_smaller_than_one",
          "honest readout clearing loses smallness")
    check(integrality_gate(evaluate(numerator(3), shell(2))) == "zero_readout",
          "integrality without nonvanishing")
    # An exactly supported head law supplies a legitimate integer lattice,
    # but its zero-target case makes the annihilated form exactly zero.
    mass = F(1, 7)
    readout = mass*evaluate(numerator(3), 1)+(1-mass)*evaluate(numerator(3), shell(2))
    check(7*readout == denominator(3), "genuine finite-head readout lattice")

    # The actual law has an unconditional infinite supply of positive shell
    # atoms: n_t=(4^(3t)-1)/3, one odd edge of valuation6t to ROOT.
    pole_controls = 0
    for m in (0, 2, 10, 100):
        previous_index = -1
        seen = set()
        for t in range(1, 13):
            source = ((1 << (6*t))-1)//3
            j = (source-3)//6
            check(source % 3 == 0 and source == 6*j+3, "explicit rooted leaf index")
            check(3*source+1 == (1 << (6*t)), "one-edge rooted leaf")
            weight = F(2, (3*t)*(3*t+1))
            check(weight > 0 and j > previous_index, "positive unbounded leaf family")
            previous_index = j
            if j > m:
                # Keep the shell index compressed: constructing2^j is not
                # needed to prove distinct poles or source-specific support.
                distance = j-m
                check(distance not in seen, "distinct positive pole addresses")
                seen.add(distance)
            pole_controls += 1

    hostiles = (lambda: numerator(True), lambda: denominator(0),
                lambda: damped_budget(2, 0), lambda: damped_budget(3, -1),
                lambda: universal_damping(1), lambda: prime_budget(3, 2),
                lambda: prime_budget(3, 9), lambda: integrality_gate(0.5),
                lambda: integrality_gate(1, True))
    for hostile_call in hostiles:
        try:
            hostile_call()
        except ValueError:
            check(True, "invalid input rejected")
        else:
            check(False, "invalid input accepted")

    print(json.dumps({"checks": checks,
        "status": "PASS; coefficient lattice is distinct from readout integrality",
        "denominator_universe": "N1..25",
        "prime_valuation_universe": "N1..50, odd primes<100",
        "extension_controls": "N1,3,5,7,9;25 normalized rational quadratics each",
        "damped_rows": damped_rows,
        "small_noninteger_hostile": str(hostile),
        "positive_pole_controls": pole_controls,
        "moment_oracle_inputs": 0, "type_hostiles": len(hostiles)}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
