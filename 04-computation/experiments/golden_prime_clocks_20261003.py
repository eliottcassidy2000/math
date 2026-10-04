"""Exact golden-ratio residue clocks; see the companion proof note.

No floating point, third-party packages, or assertions disabled by python -O.
Public pair APIs use (a,b) for a+b*phi, phi**2=phi+1.  Matrix APIs use
row-major tuples of two rows.  An omitted modulus means exact integers.
"""

from math import gcd, lcm


IDENTITY = ((1, 0), (0, 1))
M = ((0, 1), (1, 1))


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def _modulus(modulus):
    if modulus is not None and (not isinstance(modulus, int) or modulus < 1):
        raise ValueError("modulus must be a positive integer or None")


def pair_mul(u, v, modulus=None):
    """Multiply two elements of Z[phi], optionally modulo a positive integer."""
    _modulus(modulus)
    a, b = u
    c, d = v
    result = (a * c + b * d, a * d + b * c + b * d)
    return result if modulus is None else tuple(x % modulus for x in result)


def pair_norm(u):
    a, b = u
    return a * a + a * b - b * b


def is_unit(u, modulus):
    _modulus(modulus)
    if modulus is None:
        raise ValueError("is_unit requires a positive integer modulus")
    return gcd(pair_norm(u), modulus) == 1


def pair_inverse(u, modulus):
    """Modular inverse, refusing nonunits rather than dividing by nonzero."""
    if not is_unit(u, modulus):
        raise ValueError("nonunit in Z[phi]/modulus")
    a, b = u
    c = pow(pair_norm(u), -1, modulus)
    return ((a + b) * c % modulus, -b * c % modulus)


def phi_power(k, modulus=None):
    """phi**k for any signed integer k; no modular division is required."""
    _modulus(modulus)
    if not isinstance(k, int):
        raise ValueError("exponent must be an integer")
    base = (0, 1) if k >= 0 else (-1, 1)
    result = (1, 0) if modulus is None else (1 % modulus, 0)
    k = abs(k)
    while k:
        if k & 1:
            result = pair_mul(result, base, modulus)
        base = pair_mul(base, base, modulus)
        k //= 2
    return result


def mat_mul(a, b, modulus=None):
    _modulus(modulus)
    result = tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                         for j in range(2)) for i in range(2))
    return result if modulus is None else tuple(
        tuple(x % modulus for x in row) for row in result)


def mat_power(a, k, modulus=None):
    """Nonnegative matrix power; signed shift powers use phi_power instead."""
    _modulus(modulus)
    if not isinstance(k, int) or k < 0:
        raise ValueError("matrix exponent must be a nonnegative integer")
    result = IDENTITY if modulus is None else ((1 % modulus, 0), (0, 1 % modulus))
    while k:
        if k & 1:
            result = mat_mul(result, a, modulus)
        a = mat_mul(a, a, modulus)
        k //= 2
    return result


def guard_clock(a, p):
    """Exact order of multiplication by phi modulo 2**a * 3**p."""
    if any(not isinstance(e, int) or e < 0 for e in (a, p)):
        raise ValueError("prime-power exponents must be nonnegative integers")
    period2 = 1 if a == 0 else 3 * 2 ** (a - 1)
    period3 = 1 if p == 0 else 8 * 3 ** (p - 1)
    return lcm(period2, period3)


def phi_order(modulus):
    """Independent finite enumeration, for small positive moduli only."""
    _modulus(modulus)
    if modulus is None:
        raise ValueError("finite order requires a modulus")
    state = (1 % modulus, 0)
    for k in range(1, modulus * modulus + 1):
        # Direct Fibonacci register step, not pair_mul or phi_power.
        state = (state[1], (state[0] + state[1]) % modulus)
        if state == (1 % modulus, 0):
            return k
    raise RuntimeError("invertible shift failed to return")


def prime_divisors(n):
    result = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            result.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        result.append(n)
    return result


def is_prime(n):
    return n >= 2 and prime_divisors(n) == [n]


def lucas(n):
    a, b = 2, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def poly_trim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def poly_mul(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    return poly_trim(result)


def poly_div_exact(a, b):
    """Monic exact division of coefficient lists in increasing degree."""
    need(b[-1] == 1, "divisor must be monic")
    a = a[:]
    q = [0] * (len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        k = len(a) - len(b)
        q[k] = a[-1]
        for j, x in enumerate(b):
            a[k + j] -= q[k] * x
        poly_trim(a)
    need(not any(a), "nonexact cyclotomic division")
    return poly_trim(q)


def cyclotomic_table(limit):
    result = {}
    for n in range(1, limit + 1):
        poly = [-1] + [0] * (n - 1) + [1]
        for d in range(1, n):
            if n % d == 0:
                poly = poly_div_exact(poly, result[d])
        result[n] = poly
    return result


def eval_phi(coefficients):
    result = (0, 0)
    for digit in reversed(coefficients):
        # Direct exact Horner register, independent of pair multiplication.
        result = (result[1] + digit, result[0] + result[1])
    return result


def check_exact_order(modulus, period):
    identity = (1 % modulus, 0)
    need(phi_power(period, modulus) == identity, "candidate period fails")
    for q in prime_divisors(period):
        need(phi_power(period // q, modulus) != identity,
             "candidate period is not minimal")


def main():
    print("golden_prime_clocks_20261003: exact integer checks")
    print("Universe: primes <=499; moduli 1..128; unit tables 1..64;")
    print("          2-exponents 0..12, 3-exponents 0..8; cyclotomics 1..105.")

    identity_checks = 0
    for n in range(1, 81):
        left = phi_power(2 * n)
        left = (left[0] + (-1) ** n, left[1])
        right = tuple(lucas(n) * x for x in phi_power(n))
        need(left == right, "Lucas identity")
        z = phi_power(n)
        need(pair_norm((z[0] - 1, z[1])) == 1 + (-1) ** n - lucas(n),
             "signed norm")
        need(pair_mul(z, phi_power(-n)) == (1, 0), "signed power inverse")
        a, b = z
        need(mat_power(M, n) == ((a, b), (b, a + b)), "matrix/pair power")
        identity_checks += 4
    need(phi_power(8) == (13, 21), "phi8")
    need(phi_power(10) == (34, 55), "phi10")
    need(pair_norm((1, 3)) == -5, "historical signed-norm hostile")
    print(f"Lucas, signed norm, signed powers, matrix identities: {identity_checks} checks")
    print("Signed norm hostile: N(phi^4-1)=N(1+3phi)=-5; quotient size=5.")

    prime_count = 0
    split_count = 0
    inert_count = 0
    for q in range(2, 500):
        if not is_prime(q):
            continue
        roots = [x for x in range(q) if (x * x - x - 1) % q == 0]
        order = phi_order(q)
        if q == 2:
            need(not roots and order == 3, "F4 case")
        elif q == 5:
            need(roots == [3] and order == 20, "ramified case")
        elif pow(5, (q - 1) // 2, q) == 1:
            need(len(roots) == 2 and (q - 1) % order == 0, "split clock")
            need(phi_power(q - 1, q) == (1, 0), "split Fermat law")
            split_count += 1
        else:
            need(not roots and 2 * (q + 1) % order == 0, "inert clock")
            need(phi_power(q + 1, q) == (q - 1, 0), "inert Frobenius law")
            need((q + 1) % order != 0, "minus identity is not identity")
            inert_count += 1
        prime_count += 1
    print(f"Prime laws: {prime_count} primes, {split_count} split, {inert_count} odd inert, 2 exceptions.")
    for q, expected in [(3, 8), (5, 20), (7, 16), (11, 10), (105, 80)]:
        need(phi_order(q) == expected, "named clock order")
        check_exact_order(q, expected)
        print(f"ord(phi mod {q})={expected}")
    need(phi_power(8, 7) == (6, 0), "7 half turn")
    need(phi_power(10, 11) == (1, 0), "11 full turn")
    need(phi_power(5, 11) != (1, 0), "11 not half turn")

    modulus_checks = 0
    for modulus in range(1, 129):
        period = phi_order(modulus)
        check_exact_order(modulus, period)
        need(mat_power(M, period, modulus) ==
             ((1 % modulus, 0), (0, 1 % modulus)), "matrix order agreement")
        modulus_checks += 1
    print(f"Independent Fibonacci-register/matrix clock comparisons: {modulus_checks}")

    unit_checks = 0
    direct_unit_checks = 0
    for modulus in range(1, 65):
        identity = (1 % modulus, 0)
        for a in range(modulus):
            for b in range(modulus):
                u = (a, b)
                unit = is_unit(u, modulus)
                if unit:
                    need(pair_mul(u, pair_inverse(u, modulus), modulus) == identity,
                         "unit inverse")
                else:
                    try:
                        pair_inverse(u, modulus)
                    except ValueError:
                        pass
                    else:
                        raise RuntimeError("accepted a nonunit inverse")
                if modulus <= 12:
                    direct = any(((a * c + b * d) % modulus,
                                  (a * d + b * c + b * d) % modulus) == identity
                                 for c in range(modulus) for d in range(modulus))
                    need(direct == unit, "direct inverse existence vs norm test")
                    direct_unit_checks += 1
                unit_checks += 1
    counts = {q: sum(is_unit((a, b), q) for a in range(q) for b in range(q))
              for q in (3, 5, 7, 11, 105)}
    need(counts == {3: 8, 5: 20, 7: 48, 11: 100, 105: 7680}, "unit counts")
    need(pair_mul((2, 1), (2, 1), 5) == (0, 0), "ramified nilpotent")
    need(pair_mul((102, 1), (63, 84), 105) == (0, 0), "105 zero divisor")
    print(f"Unit criterion: {unit_checks} pairs; {direct_unit_checks} independent brute inverse tests")
    print(f"Unit counts: {counts}; mod105 phi orbit covers 80/7680=1/96 of units.")
    print("Hostile: (phi-3)^2=0 mod5 although phi-3!=0; no nonzero-only division.")

    power_checks = 0
    for a in range(13):
        for p in range(9):
            modulus = 2 ** a * 3 ** p
            period = guard_clock(a, p)
            check_exact_order(modulus, period)
            need(mat_power(M, period, modulus) ==
                 ((1 % modulus, 0), (0, 1 % modulus)), "guard matrix clock")
            if a > 0 and p > 0:
                need(period == 2 ** max(a - 1, 3) * 3 ** max(1, p - 1),
                     "closed guard clock")
            power_checks += 1
    for a in range(1, 13):
        need(phi_order(2 ** a) == 3 * 2 ** (a - 1), "2-primary direct period")
    for p in range(1, 9):
        need(phi_order(3 ** p) == 8 * 3 ** (p - 1), "3-primary direct period")
    print(f"Guard clocks: {power_checks} mixed powers, plus 20 direct prime-power cycles")
    print("Example: modulus 2^8*3^3=6912, clock=" + str(guard_clock(8, 3)))

    # Test shift-clock usefulness and its loss of exact values.
    u = (105, 0)
    shifted = pair_mul(u, phi_power(80))
    need(shifted != u and pair_mul(u, phi_power(80, 105), 105) == (0, 0),
         "clock does not preserve exact integer")
    need(phi_power(80) != (1, 0), "modular identity is not real identity")
    # A tiny source guard requires the next binary digit to fix valuation.
    need((3 * 1 + 1) % 4 == (3 * 5 + 1) % 4 == 0,
         "coarse valuation class")
    need((3 * 1 + 1) % 8 == 4 and (3 * 5 + 1) % 8 == 0,
         "one-bit exact valuation hostile")
    print("Hostiles: shift by a full clock loses exact size; n=1,5 mod4 share v2>=2 but not v2=2.")

    table = cyclotomic_table(105)
    first = next(n for n, poly in table.items() if max(map(abs, poly)) > 1)
    need(first == 105, "first nonflat finite audit")
    c105 = table[105]
    need(len(c105) == 49 and c105 == c105[::-1], "Phi105 degree/reciprocity")
    need([(i, c) for i, c in enumerate(c105) if abs(c) > 1] == [(7, -2), (41, -2)],
         "Phi105 nonflat positions")
    # Independent rational product identity for squarefree 105.
    binomial = lambda n: [-1] + [0] * (n - 1) + [1]
    numerator = [1]
    denominator = [1]
    for n in (105, 3, 5, 7):
        numerator = poly_mul(numerator, binomial(n))
    for n in (35, 21, 15, 1):
        denominator = poly_mul(denominator, binomial(n))
    need(poly_mul(c105, denominator) == numerator, "independent Phi105 product")
    evaluation = eval_phi(c105)
    need(evaluation == (5224431949, 8453308464), "Phi105 evaluation")
    need(pair_norm(evaluation) == 16271615641 == 21211 * 767131,
         "Phi105 norm factorization")
    need(is_prime(21211) and is_prime(767131), "Phi105 norm prime factors")
    print("Cyclotomics: all 1..104 flat; Phi105 degree48; only nonflat coefficients (7,-2),(41,-2).")
    print(f"Phi105(phi)={evaluation}; norm=16271615641=21211*767131")

    repunit_count = 0
    for q in range(3, 106, 2):
        if not is_prime(q):
            continue
        repunit = eval_phi([1] * q)
        need(repunit == eval_phi(table[q]), "prime cyclotomic repunit")
        z = phi_power(q)
        need(repunit == pair_mul((0, 1), (z[0] - 1, z[1])), "repunit quotient")
        need(pair_norm(repunit) == lucas(q), "odd prime repunit norm")
        repunit_count += 1
    need(eval_phi([1, 1]) == phi_power(2), "11_phi=100_phi")
    p23 = eval_phi([1] * 23)
    need(p23 == pair_mul(phi_power(12), pair_mul((11, 2), (21, 1))),
         "prime-index repunit factorization")
    need(pair_norm((11, 2)) == 139 and pair_norm((21, 1)) == 461,
         "repunit factor norms")
    need(eval_phi([1] * 3) == (2, 2), "norm4 prime ideal exception")
    print(f"Prime-index repunits: {repunit_count} exact Lucas norms")
    print("Hostile: Phi23(phi)=phi^12*(11+2phi)*(21+phi), norm=64079=139*461.")
    print("Boundary: Phi3(phi)=2phi^2 is algebraically prime, quotient F4, norm4.")
    print("PASS: all exact checks completed; no global primality or Collatz-coverage claim.")


if __name__ == "__main__":
    main()
