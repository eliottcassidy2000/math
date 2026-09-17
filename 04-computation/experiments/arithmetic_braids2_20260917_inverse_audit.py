"""Independent brute-force audit of finite-word inverse completions.

Imports no primary completion code and uses no discrete logarithm routine.
Universe: every word of length 0..3 over {1,2,3}, ten b values and ten
targets, all 27 source-index classes. Exact Fraction arithmetic checks the
2-adic singular limit separately. All checks survive python -O.
"""
from fractions import Fraction
from itertools import product
from math import gcd
from pathlib import Path
import hashlib
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def valuation(n, p):
    require(n != 0, "finite valuation required")
    n = abs(n)
    result = 0
    while n % p == 0:
        n //= p
        result += 1
    return result


def rational_valuation(x, p):
    return valuation(x.numerator, p)-valuation(x.denominator, p)


def main():
    b_values = (-13, -11, -7, -5, -1, 1, 5, 7, 11, 13)
    targets = (-25, -17, -7, -5, -1, 1, 5, 7, 17, 25)
    cases = trajectory_checks = isometry_checks = singular_checks = 0
    maximum_bits = 0
    for length in range(4):
        for word in product((1, 2, 3), repeat=length):
            # Independent closed-form carry, rather than the primary recurrence.
            prefix_sums = [sum(word[:i]) for i in range(length)]
            total = sum(word)
            carry = sum(3**(length-1-i)*2**prefix_sums[i] for i in range(length))
            for b in b_values:
                for target in targets:
                    k0 = next(k for k in (1, 2) if (2**k*target-b) % 3 == 0)
                    leading = 2**(total+k0)*target
                    constant = b*(2**total+3*carry)
                    denominator = 3**(length+1)
                    period = 3**length
                    integral_t = [t for t in range(period)
                                  if (leading*4**t-constant) % denominator == 0]
                    require(len(integral_t) == 1, "brute integral t is not unique")
                    t0 = integral_t[0]
                    sources = [(leading*4**(t0+period*j)-constant)//denominator
                               for j in range(27)]
                    require(len({source % 27 for source in sources}) == 27,
                            "incomplete ternary cylinder image")
                    for j in range(1, 27):
                        require(valuation(sources[j]-sources[0], 3) == valuation(j, 3),
                                "3-adic isometry failed")
                        isometry_checks += 1

                    sign_witnesses = []
                    for j, source in enumerate(sources):
                        current = source
                        exponents = []
                        states = [source]
                        for step in range(length+1):
                            image = 3*current+b
                            k = valuation(image, 2)
                            current = image//2**k
                            exponents.append(k)
                            states.append(current)
                        require(tuple(exponents[:-1]) == word, "prefix mismatch")
                        require(exponents[-1] == k0+2*(t0+period*j), "final exponent mismatch")
                        require(current == target, "wrong target")
                        require(all(n % 2 for n in states), "even state")
                        require(all(gcd(n, b) == gcd(target, b) for n in states),
                                "gcd stratum not preserved")
                        if source*target > 0:
                            require(all(n*target > 0 for n in states), "intermediate sign change")
                            sign_witnesses.append(j)
                        maximum_bits = max(maximum_bits, abs(source).bit_length())
                        trajectory_checks += 1
                    require(sign_witnesses, "no same-sign source in the audit window")

                    limit = -Fraction(constant, denominator)
                    require(rational_valuation(limit, 2) == 0, "limit is not 2-adically odd")
                    require(rational_valuation(limit, 3) == -length-1, "limit 3-adic pole")
                    require((3**length*limit+b*carry)/2**total == -Fraction(b, 3),
                            "affine singular preimage identity")
                    current = limit
                    for k in word:
                        image = 3*current+b
                        require(rational_valuation(image, 2) == k, "limit prefix exponent mismatch")
                        current = image/2**k
                    require(current == -Fraction(b, 3) and 3*current+b == 0,
                            "limit does not reach singular point")
                    for j in (0, 1, 2, 26):
                        t = t0+period*j
                        difference = Fraction(sources[j])-limit
                        require(difference == Fraction(leading*4**t, denominator),
                                "singular difference identity")
                        require(rational_valuation(difference, 2) == total+k0+2*t,
                                "exact 2-adic convergence depth")
                        singular_checks += 1
                    cases += 1

    output = {
        "status": "FINITE-EXACT independent audit PASS; general proof supplied separately",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "universe": {
            "word_lengths": [0, 1, 2, 3], "word_alphabet": [1, 2, 3],
            "b_values": b_values, "targets": targets, "j_values": "0..26",
            "singular_limit_j_values": [0, 1, 2, 26],
        },
        "independence": "No primary imports; brute t search; closed-form carry; direct trajectories",
        "parameter_cases": cases,
        "direct_trajectory_checks": trajectory_checks,
        "isometry_valuation_checks": isometry_checks,
        "singular_limit_valuation_checks": singular_checks,
        "maximum_source_bits": maximum_bits,
        "checks": ["unique t modulo 3^L", "complete source image modulo27",
                   "exact prefix and last exponent", "same-sign source and intermediates",
                   "gcd(n,b) preservation", "limit maps to -b/3 after the prescribed word",
                   "v2(source-limit)=K+k0+2t", "v3(limit)=-(L+1)"],
        "limitations": ["No global convergence or natural density assertion",
                        "The singular limit is not 3-adically integral",
                        "Different finite words use different starting integers"],
    }
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(output, indent=2, sort_keys=True)+"\n", encoding="utf-8", newline="\n")
    print(f"independent inverse audit: {cases} parameter cases; PASS")
    print(f"direct trajectories: {trajectory_checks}; ternary isometry checks: {isometry_checks}; PASS")
    print(f"singular-limit checks: {singular_checks}; maximum source bits: {maximum_bits}; PASS")


if __name__ == "__main__":
    main()
