#!/usr/bin/env python3
"""Exact finite controls for the pointed parity-conjugacy integration.

No orbit search is used to certify an unbounded family. The infinite claims
have direct proofs in the companion note. Run normally and with python -O.
"""

from fractions import Fraction


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def step(n, a=3):
    return (a * n + 1) // 2 if n % 2 else n // 2


def itinerary(n, width, a=3):
    code = 0
    for j in range(width):
        code |= (n % 2) << j
        n = step(n, a)
    return code


def inverse_prefix(code, width, a=3):
    """Affine numerator method: T_a^m(x)=(a^r*x+b)/2^m."""
    b, r = 0, 0
    for j in range(width):
        if (code >> j) & 1:
            b = a * b + (1 << j)
            r += 1
    modulus = 1 << width
    return (-b * pow(pow(a, r, modulus), -1, modulus)) % modulus


def inverse_branches(code, width, a=3):
    """Independent reverse-branch method, arbitrary terminal tail zero."""
    modulus = 1 << width
    inverse_a = pow(a, -1, modulus)
    value = 0
    for j in reversed(range(width)):
        value = ((2 * value - 1) * inverse_a if (code >> j) & 1
                 else 2 * value) % modulus
    return value


def rational_step(x, a):
    if x.denominator % 2 == 0:
        raise ValueError("not a 2-adic integer")
    return (a * x + 1) / 2 if x.numerator % 2 else x / 2


def main():
    print("FINITE-EXACT: parity prefixes, pointed conjugacy and finite flips")
    print("Universe: a in {1,3,5,7,11}; all residues modulo 2^m, 1<=m<=12")
    total = 0
    for a in (1, 3, 5, 7, 11):
        for width in range(1, 13):
            modulus = 1 << width
            codes = [itinerary(n, width, a) for n in range(modulus)]
            require(len(set(codes)) == modulus, "prefix bijection")
            for n, code in enumerate(codes):
                total += 1
                require(inverse_prefix(code, width, a) == n, "affine inverse")
                require(inverse_branches(code, width, a) == n, "branch inverse")
                require(itinerary(step(n, a), width - 1, a) == code >> 1,
                        "shift compatibility")
                require(itinerary(n + modulus, width, a) == code,
                        "cylinder is independent of next source bit")
    print(f"Prefix samples {total}; all inverses and shift identities passed")

    width, k = 10, 4
    representatives = {}
    for n in range(1 << width):
        code = itinerary(n, width)
        orbit = {inverse_prefix(code ^ mask, width) for mask in range(1 << k)}
        require(len(orbit) == 1 << k, "free finite-flip action")
        require(all(itinerary(y, width) >> k == code >> k for y in orbit),
                "tail retained")
        representatives.setdefault(code >> k, orbit)
        require(representatives[code >> k] == orbit, "one finite tail class")
    require(len(representatives) == 1 << (width - k), "finite orbit count")
    print("Finite-flip model: 64 classes of size16 at width10, prefix4")
    print("The finite model has selectors; infinite no-measurable-selector uses the proof")

    for a in (1, 3, 5, 7, 11):
        x = Fraction(1, 4 - a)
        require(rational_step(x, a) == 2 * x, "odd root image")
        require(rational_step(2 * x, a) == x, "even root image")
        value, code = x, 0
        for j in range(24):
            code |= (value.numerator % 2) << j
            value = rational_step(value, a)
        alternating = sum(1 << j for j in range(0, 24, 2))
        require(code == alternating, "marked alternating itinerary")
        require(inverse_prefix(code, 24, a) ==
                x.numerator * pow(x.denominator, -1, 1 << 24) % (1 << 24),
                "rational and modular representations agree")
        print(f"T_3 ROOT phase1 maps under Q_{a}^-1 Q_3 to {x}")
    require(Fraction(1, 4 - 5) < 0, "positive-lattice hostile")

    # A short, independently evaluated supplied cycle for 5x+1.
    cycle = (1, 3, 8, 4, 2)
    require(all(step(cycle[j], 5) == cycle[(j + 1) % len(cycle)]
                for j in range(len(cycle))), "positive T_5 cycle")
    require(len(cycle) != 2, "T_5 positive ROOT differs from conjugate ROOT")
    print("Hostile: T_3 cycle(1,2) maps to T_5 cycle(-1,-2), not positive cycle(1,3,8,4,2)")
    print(f"PASS: {CHECKS} exact checks; no universal integer coverage claimed")


if __name__ == "__main__":
    main()
