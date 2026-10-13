"""Exact Tunnell coefficient receipts; cited arithmetic consequences are separate.

No numerical L-function, Sel4 computation, installation, or curve replacement.
Run with python -B or python -B -O. Every check uses explicit exceptions.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from math import isqrt
from collections import Counter
import hashlib
import json


CHECKS = 0


def check(value, label):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(label)


def positive(n):
    if type(n) is not int or n <= 0:
        raise ValueError("expected exact positive integer")
    return n


def odd_squarefree(n):
    positive(n)
    if n % 2 == 0:
        raise ValueError("Tunnell odd branch requires odd n")
    for p in range(3, isqrt(n) + 1, 2):
        if n % (p * p) == 0:
            raise ValueError("Tunnell receipt requires squarefree n")
    return n


def direct_count(n, c):
    """All signed x,z, then exact y-square test; finite bounds are exhaustive."""
    positive(n)
    if type(c) is not int or c not in (8, 32):
        raise ValueError("coefficient c must be exact 8 or 32")
    out = 0
    for x in range(-isqrt(n // 2), isqrt(n // 2) + 1):
        for z in range(-isqrt(n // c), isqrt(n // c) + 1):
            remainder = n - 2 * x * x - c * z * z
            if remainder < 0:
                continue
            y = isqrt(remainder)
            if y * y == remainder:
                out += 1 if y == 0 else 2
    return out


def theta_product(limit, c):
    """Independent polynomial convolution of three nonnegative-square series."""
    positive(limit)
    if type(c) is not int or c not in (8, 32):
        raise ValueError("coefficient c must be exact 8 or 32")
    poly = [1] + [0] * limit
    for weight in (1, 2, c):
        nxt = [0] * (limit + 1)
        for root in range(isqrt(limit // weight) + 1):
            degree, multiplicity = weight * root * root, (1 if root == 0 else 2)
            for j in range(limit - degree + 1):
                nxt[j + degree] += multiplicity * poly[j]
        poly = nxt
    return tuple(poly)


def eta_theta(limit):
    """q product_m(1-q^(8m))(1-q^(16m)) times sum_z q^(2z^2)."""
    positive(limit)
    poly = [0] * (limit + 1)
    poly[1] = 1
    for spacing in (8, 16):
        for degree in range(spacing, limit + 1, spacing):
            for j in range(limit, degree - 1, -1):
                poly[j] -= poly[j - degree]
    result = [0] * (limit + 1)
    for root in range(isqrt(limit // 2) + 1):
        degree, multiplicity = 2 * root * root, (1 if root == 0 else 2)
        for j in range(limit - degree + 1):
            result[j + degree] += multiplicity * poly[j]
    return tuple(result)


@dataclass(frozen=True)
class Receipt:
    n: int
    a32: int
    b8: int
    coefficient: int


def terminal(n):
    """Nonzero exact coefficient or None; None is no rank/corank conclusion."""
    odd_squarefree(n)
    a, b = direct_count(n, 32), direct_count(n, 8)
    if b % 2:
        raise ArithmeticError("unexpected parity")
    coefficient = a - b // 2
    return None if coefficient == 0 else Receipt(n, a, b, coefficient)


def audit(receipt):
    if type(receipt) is not Receipt:
        raise ValueError("expected Receipt")
    if any(type(x) is not int for x in
           (receipt.n, receipt.a32, receipt.b8, receipt.coefficient)):
        raise ValueError("receipt fields must be exact integers")
    if terminal(receipt.n) != receipt:
        raise ValueError("receipt does not match the source counts")
    return receipt


def census(limit=10000):
    """Same odd-squarefree sources; (n,pure-Sel2-dimension,Tunnell-a) rows.

    The arithmetic meaning of both gates is cited, not derived by this code.
    Zero a with dimension>=2 is unresolved by their union, not a negative BSD test.
    """
    import bsd_selmer_growth_20261013a as previous
    positive(limit)
    a32, b8, eta = theta_product(limit, 32), theta_product(limit, 8), eta_theta(limit)
    rows = []
    for n in range(1, limit + 1, 2):
        if 2 * eta[n] != 2 * a32[n] - b8[n]:
            raise ArithmeticError("independent coefficient engines disagree")
        try:
            primes, _, matrix = previous.monsky(n)
        except ValueError:
            continue
        width = 2 * len(primes)
        dimension = width - previous.rank(matrix, width)
        rows.append((n, dimension, eta[n]))
    return tuple(rows)


def main():
    import bsd_selmer_growth_20261013a as previous
    limit = 627
    a32, b8, eta = theta_product(limit, 32), theta_product(limit, 8), eta_theta(limit)
    for n in range(1, limit + 1):
        check(a32[n] == direct_count(n, 32), f"32-count {n}")
        check(b8[n] == direct_count(n, 8), f"8-count {n}")
        if n % 2:
            check(2 * eta[n] == 2 * a32[n] - b8[n], f"eta identity {n}")
    print("Exact count comparison: all n=1..627, both ternary forms.")
    print("Independent eta-theta coefficient comparison: all 314 odd n<=627.")
    for n, expected in ((17, (4, 16, -4)), (627, (40, 48, 16))):
        receipt = audit(terminal(n))
        check((receipt.a32, receipt.b8, receipt.coefficient) == expected, "target counts")
        primes, _, rows = previous.monsky(n)
        width = 2 * len(primes)
        pure_dimension = width - previous.rank(rows, width)
        check(pure_dimension == 2, "unchanged first-level obstruction")
        factor = Fraction(receipt.coefficient ** 2, 4)
        print(f"E{n}: A32={receipt.a32}; B8={receipt.b8}; a={receipt.coefficient}; "
              f"L(E,1)={factor}*beta/sqrt({n}); pure-Sel2-dimension={pure_dimension}.")
    check(terminal(5) is None, "zero coefficient is not a terminal receipt")
    x, y = Fraction(25, 4), Fraction(75, 8)
    check(y * y == x * x * x - 25 * x and y != 0, "E5 rational-point control")
    print("Zero control E5: a=0; exact point (25/4,75/8); no nonvanishing receipt.")
    bad = [True, 17.0, 0, -17, 2, 9, 45]
    for value in bad:
        try:
            terminal(value)
        except ValueError:
            check(True, "invalid source rejected")
        else:
            check(False, "invalid source accepted")
    good = terminal(17)
    for malformed in (replace(good, n=17.0), replace(good, a32=4.0),
                      replace(good, coefficient=4), replace(good, b8=8),
                      (17, 4, 16, -4)):
        try:
            audit(malformed)
        except ValueError:
            check(True, "forged receipt rejected")
        else:
            check(False, "forged receipt accepted")
    rows = census()
    baseline = tuple(n for n, d, a in rows if d <= 1)
    gained = tuple(n for n, d, a in rows if d >= 2 and a != 0)
    remaining = tuple(n for n, d, a in rows if d >= 2 and a == 0)
    check((len(rows), len(baseline), len(gained), len(remaining)) ==
          (4056, 2886, 724, 446), "fixed universe and union")
    gain_hist = dict(sorted(Counter(d for n, d, a in rows if d >= 2 and a).items()))
    remain_hist = dict(sorted(Counter(d for n, d, a in rows if d >= 2 and not a).items()))
    check(gain_hist == {2: 716, 4: 8}, "gain dimensions")
    check(remain_hist == {2: 272, 3: 164, 4: 10}, "remaining dimensions")
    check(set(baseline).isdisjoint(gained) and set(gained).isdisjoint(remaining), "disjoint sources")
    check(set(baseline) | set(gained) | set(remaining) == {n for n, _, _ in rows}, "complete source partition")
    check(17 in gained and 627 in gained and remaining[0] == 41 and 5 in baseline,
          "marked terminals and least residual")
    # Predeclared stratified controls, not an undisclosed full second census.
    selected = (set(gained[:12]) | set(gained[-12:]) | set(gained[::53]) |
                set(remaining[::37]) | {n for n, d, a in rows if d == 4 and a} | {17, 627, 5})
    coefficient_map = {n: a for n, d, a in rows}
    for n in sorted(selected):
        aa, bb = direct_count(n, 32), direct_count(n, 8)
        check(2 * coefficient_map[n] == 2 * aa - bb, f"bulk signed control {n}")
        receipt = terminal(n)
        check((receipt is None) == (coefficient_map[n] == 0), f"bulk gate control {n}")
        if receipt is not None:
            check(audit(receipt).coefficient == coefficient_map[n], f"bulk binding {n}")
    encoded = json.dumps(rows, separators=(",", ":")).encode("ascii")
    print(f"Declared census: {len(rows)} positive odd squarefree n<=10000; baseline={len(baseline)}.")
    print(f"Additional nonzero terminals={len(gained)}; union={len(baseline)+len(gained)}; remaining={len(remaining)}.")
    print(f"Gain histogram by pure-Sel2 dimension: {gain_hist}.")
    print(f"Remaining histogram by pure-Sel2 dimension: {remain_hist}.")
    print(f"Bulk coefficient cross-check: all 5000 odd indices; signed controls={len(selected)}.")
    print("All source/dimension/coefficient rows SHA256=" + hashlib.sha256(encoded).hexdigest())
    print("Remaining sources=" + json.dumps(remaining, separators=(",", ":")))
    print("CITED consequence: analytic rank0; CM rank0 and finite Sha; Sel(q^infinity) corank0.")
    print("Full BSD: conditional application of supplied accepted low-corank theorem.")
    print("No Sel4 group/order or numerical L-value is computed.")
    print(f"Checks: {CHECKS}")


if __name__ == "__main__":
    main()
