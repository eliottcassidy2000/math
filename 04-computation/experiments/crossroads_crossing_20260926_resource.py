"""Exact positive shadows of expanding rational Collatz cycles.

Tests an obstruction to ALL finite polynomial-valuation corrections to log height.
No floating arithmetic, randomness, or convergence assumption is used.
"""
from fractions import Fraction as F
from math import lcm
import json


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def vp(x, p):
    x = F(x)
    require(x != 0, "finite valuation required")
    a, b, v = abs(x.numerator), x.denominator, 0
    while a % p == 0:
        a //= p
        v += 1
    while b % p == 0:
        b //= p
        v -= 1
    return v


def poly(coefficients, x):
    out = 0
    for c in reversed(coefficients):
        out = out * x + c
    return out


def order(a, p):
    t, k = a % p, 1
    while t != 1:
        t = t * a % p
        k += 1
    return k


def rational_residue(r, modulus):
    return r.numerator * pow(r.denominator, -1, modulus) % modulus


def crt(congruences):
    residue, modulus = 0, 1
    for target, base in congruences:
        residue += modulus * ((target - residue) * pow(modulus, -1, base) % base)
        modulus *= base
    return residue, modulus


def choose_shadow(primes, polynomials):
    B, period = 1, 1
    for p in primes:
        if p >= 5:
            B = lcm(B, order(2, p))
            period = lcm(period, order(2, p), order(3, p))
    a = 1
    rejected_roots = 0
    while True:
        D = 3 ** a - 2 ** (a + B)
        if D > 0:
            r = F(-(3 ** a - 2 ** a), D)
            if all(poly(P, r) != 0 for P in polynomials):
                break
            rejected_roots += 1
        a += period
    require(r < -1, "negative expanding shadow")
    word = (1,) * a + (0,) * B
    state = r
    for bit in word:
        require(state.denominator % 2 == 1 and state.numerator % 2 == bit, "rational word is actually legal")
        state = (3 * state + 1) / 2 if bit else state / 2
    require(state == r, "rational periodic point closure")
    require(vp(r + 1, 2) == a and vp(r, 2) == 0, "direct word certificate")
    return a, B, r, word, rejected_roots


def family(primes, polynomials):
    a, B, r, word, rejected = choose_shadow(primes, polynomials)
    all_primes = sorted(set(primes) | {2})
    precision = {}
    for p in all_primes:
        require(r.denominator % p != 0, "shadow denominator is a unit at every selected prime")
        precision[p] = 1 + max(vp(poly(P, r), p) for P in polynomials)
    feature = [vp(poly(P, r), p) for p in primes for P in polynomials]
    multiplier = F(3 ** a, 2 ** len(word))
    controls = []
    cutoff = 1000
    for repeats in (1, 2, 5, 10, 20):
        length = repeats * len(word)
        congruences = []
        for p in all_primes:
            exponent = precision[p] + (length if p == 2 else 0)
            modulus = p ** exponent
            congruences.append((rational_residue(r, modulus), modulus))
        residue, modulus = crt(congruences)
        source = residue + modulus * (cutoff * 2 ** length + 1)
        require([vp(poly(P, source), p) for p in primes for P in polynomials] == feature, "source features match rational shadow")
        state = source
        minimum = source
        for bit in word * repeats:
            require(state > cutoff, "whole positive block stays outside finite core")
            require(state % 2 == bit, "actual positive parity prefix")
            require(all(poly(P, state) != 0 for P in polynomials), "no undefined feature on intermediate states")
            state = (3 * state + 1) // 2 if bit else state // 2
            minimum = min(minimum, state)
        target = state
        require(target == r + multiplier ** repeats * (source - r), "unwrapped affine endpoint")
        require(F(target, source) > multiplier ** repeats, "unbounded height gain")
        require([vp(poly(P, target), p) for p in primes for P in polynomials] == feature, "target features exactly equal source")
        require(minimum > cutoff, "last state also outside finite core")
        controls.append({"repeats": repeats, "steps": length, "source_bits": source.bit_length(),
                         "target_bits": target.bit_length(), "endpoint_features_equal": True})
    return {"primes": primes, "a": a, "B": B, "rational_shadow": str(r),
            "rejected_polynomial_roots": rejected, "polynomials_low_coefficient_first": polynomials,
            "feature_vector": feature, "precision": precision, "controls": controls}


def main():
    forms = [[1, 1], [-1, 1], [5, 3], [1, 2], [5, 1], [17, 1], [1, 0, 1], [1, 1, 1]]
    cases = [family(P, forms) for P in ([2], [2, 3], [2, 3, 11], [2, 3, 5, 11], [2, 3, 5, 7, 11, 13])]
    require(cases[0]["rejected_polynomial_roots"] >= 1, "hostile: n+5 excludes naive shadow -5")
    print(json.dumps({"status": "FINITE-EXACT controls of proved rational-shadow obstruction", "cases": cases}, indent=2))


if __name__ == "__main__":
    main()
