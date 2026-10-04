"""Exact Laurent base-phi digits, carry certificates, and Fibonacci coordinates.

Run: python3 04-computation/experiments/golden_digit_carry_20261003.py
No floating point, external packages, or primality assumptions are used.
Finite audits are explicit in main(); the polynomial certificate is general.
Pairs (A, B) mean A+B*phi, where phi**2=phi+1.
"""

from collections import Counter
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from pathlib import Path


def pair_add(u, v, modulus=None):
    out = u[0] + v[0], u[1] + v[1]
    return out if modulus is None else tuple(x % modulus for x in out)


def pair_mul(u, v, modulus=None):
    a, b = u
    c, d = v
    out = a * c + b * d, a * d + b * c + b * d
    return out if modulus is None else tuple(x % modulus for x in out)


@lru_cache(maxsize=None)
def phi_power(k, modulus=None):
    if not isinstance(k, int):
        raise ValueError("exponent must be an integer")
    if modulus is not None and modulus < 1:
        raise ValueError("modulus must be positive")
    base = (0, 1) if k >= 0 else (-1, 1)
    result = (1, 0)
    exponent = abs(k)
    while exponent:
        if exponent & 1:
            result = pair_mul(result, base, modulus)
        base = pair_mul(base, base, modulus)
        exponent >>= 1
    return result if modulus is None else tuple(x % modulus for x in result)


def pair_sign(u):
    """Compare A+B*phi with zero using exact integer squares."""
    a, b = u
    c = 2 * a + b  # Twice the value is c+b*sqrt(5).
    if b == 0:
        return (c > 0) - (c < 0)
    if c == 0:
        return (b > 0) - (b < 0)
    if (c > 0) == (b > 0):
        return (c > 0) - (c < 0)
    if c > 0:
        return 1 if c * c > 5 * b * b else -1
    return 1 if 5 * b * b > c * c else -1


def clean_polynomial(coefficients):
    return {k: c for k, c in coefficients.items() if c}


def polynomial_add(a, b, scale=1):
    out = Counter(a)
    for k, c in b.items():
        out[k] += scale * c
    return clean_polynomial(out)


def convolve(a, b):
    out = Counter()
    for i, c in a.items():
        for j, d in b.items():
            out[i + j] += c * d
    return clean_polynomial(out)


def evaluate_laurent(coefficients, modulus=None):
    result = (0, 0)
    for exponent, coefficient in coefficients.items():
        result = pair_add(result, pair_mul((coefficient, 0),
                          phi_power(exponent, modulus), modulus), modulus)
    return result


def parse_phi_word(word):
    """Return the Laurent coefficient word and the recorded radix depth."""
    if word.count(".") > 1 or any(c not in "01." for c in word):
        raise ValueError("expected a binary positional word")
    whole, separator, fractional = word.partition(".")
    if not whole or (separator and not fractional):
        raise ValueError("digits required on each present side of the radix")
    digits = whole + fractional
    radix = len(fractional)
    coefficients = {len(digits) - 1 - i - radix: 1
                    for i, c in enumerate(digits) if c == "1"}
    return coefficients, radix


def read_phi_word(word, modulus=None):
    """Horner reader, followed by the invertible radix correction phi**(-R)."""
    _, radix = parse_phi_word(word)
    a = b = 0
    for symbol in word.replace(".", ""):
        a, b = b + int(symbol), a + b
        if modulus is not None:
            a, b = a % modulus, b % modulus
    return pair_mul((a, b), phi_power(-radix, modulus), modulus)


def support_word(support):
    support = frozenset(support)
    if not support:
        return "0"
    highest, lowest = max(0, max(support)), min(0, min(support))
    whole = "".join(str(int(k in support)) for k in range(highest, -1, -1))
    if lowest == 0:
        return whole
    fractional = "".join(str(int(k in support)) for k in range(-1, lowest - 1, -1))
    return whole + "." + fractional


@lru_cache(maxsize=None)
def canonical_phi_integer(n, max_fractional_places=512):
    """Exact greedy normal form, with a declared finite search bound.

    The finite audit checks termination in its stated universe. It does not
    silently use a decimal approximation or claim the bound works for all n.
    """
    if not isinstance(n, int) or n < 0:
        raise ValueError("source must be a nonnegative integer")
    if n == 0:
        return "0"
    highest = 0
    while pair_sign(pair_add((n, 0), pair_mul((-1, 0), phi_power(highest + 1)))) >= 0:
        highest += 1
    remainder = (n, 0)
    support = []
    for k in range(highest, -max_fractional_places - 1, -1):
        after = pair_add(remainder, pair_mul((-1, 0), phi_power(k)))
        if pair_sign(after) >= 0:
            support.append(k)
            remainder = after
        if remainder == (0, 0):
            word = support_word(support)
            assert "11" not in word.replace(".", "")
            return word
    raise ArithmeticError("declared finite fractional search bound exhausted")


def carry_certificate(raw, normalized, radix=None):
    """Unique Q with t**R*(raw-normalized)=(t**2-t-1)*Q.

    Both arguments are finite Laurent coefficient dictionaries. The chosen
    R clears all their negative exponents, including canceled terms. Q has
    integer coefficients and fixed R makes it unique. Reject unequal values.
    """
    raw, normalized = clean_polynomial(raw), clean_polynomial(normalized)
    required = max(0, -min(set(raw) | set(normalized) | {0}))
    if radix is None:
        radix = required
    if not isinstance(radix, int) or radix < required:
        raise ValueError("radix does not clear the Laurent inputs")
    difference = polynomial_add(raw, normalized, -1)
    remainder = {k + radix: c for k, c in difference.items()}
    quotient = {}
    while remainder and max(remainder) >= 2:
        degree = max(remainder)
        coefficient = remainder[degree]
        quotient[degree - 2] = coefficient
        remainder[degree] -= coefficient
        remainder[degree - 1] = remainder.get(degree - 1, 0) + coefficient
        remainder[degree - 2] = remainder.get(degree - 2, 0) + coefficient
        remainder = clean_polynomial(remainder)
    if remainder:
        raise ValueError(f"different exact phi values; nonzero remainder {remainder}")
    return radix, quotient


def recover_raw(normalized, radix, quotient):
    shifted = {k - radix: c for k, c in
               convolve({2: 1, 1: -1, 0: -1}, quotient).items()}
    return polynomial_add(normalized, shifted)


def carry_laurent(raw, normalized):
    """Return the radix-independent Laurent quotient (raw-normalized)/f."""
    radix, quotient = carry_certificate(raw, normalized)
    return {k - radix: c for k, c in quotient.items()}


def to_zeck_registers(pair):
    a, b = pair
    return a + 2 * b, 2 * a + 3 * b


def from_zeck_registers(registers):
    x, y = registers
    return 2 * y - 3 * x, 2 * x - y


def zeck_lift(n, zeck_module):
    """Read integer n's Fibonacci digits as phi powers, without a radix."""
    return read_phi_word("".join(map(str, zeck_module.zeckendorf_digits(n))))


def load_sibling_module(name):
    path = Path(__file__).with_name(name + ".py")
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    zeck = load_sibling_module("zeckendorf_guard_automaton_20261003")
    bound = 1000
    words = {n: canonical_phi_integer(n) for n in range(bound + 1)}
    moduli = (2, 3, 4, 5, 7, 8, 9, 11, 16, 27, 64, 81, 105)
    reader_checks = 0
    for n, word in words.items():
        coefficients, _ = parse_phi_word(word)
        assert evaluate_laurent(coefficients) == (n, 0) == read_phi_word(word)
        assert "11" not in word.replace(".", "")
        for modulus in moduli:
            assert read_phi_word(word, modulus) == (n % modulus, 0)
            assert evaluate_laurent(coefficients, modulus) == (n % modulus, 0)
            reader_checks += 1
    print(f"FINITE-EXACT greedy normal forms n=0..{bound}: {len(words)}")
    print(f"Independent Laurent / Horner modular checks: {reader_checks}; moduli={moduli}")
    for n in (2, 3, 5, 7, 11, 18, 105):
        print(f"{n} = {words[n]}_phi")

    # The same concatenated digits have different values if the radix is lost.
    assert read_phi_word("10.01") == (2, 0)
    assert read_phi_word("1001") == (2, 2)
    assert zeck.direct_digit_value((1, 0, 0, 1)) == 6
    print("Radix hostile: literal 10.01_phi=2; no-radix 1001_phi=2+2phi; Zeckendorf 1001=6")

    factors = tuple(parse_phi_word(words[n])[0] for n in (3, 5, 7))
    raw = {0: 1}
    for factor in factors:
        raw = convolve(raw, factor)
    normalized, _ = parse_phi_word(words[105])
    radix, quotient = carry_certificate(raw, normalized)
    assert evaluate_laurent(raw) == (105, 0)
    assert recover_raw(normalized, radix, quotient) == raw
    assert max(raw.values()) == 2
    factored_quotient = convolve({4: -1}, convolve({2: 1, 1: -1, 0: 1}, {8: 1, 4: 1, 0: 1}))
    assert quotient == factored_quotient
    qa, qb = evaluate_laurent(quotient)
    assert (qa, qb) == pair_mul((-16, 0), phi_power(8)) == (-208, -336)
    assert qa * qa + qa * qb - qb * qb == 256
    print("105 factor supports:", [tuple(sorted(factor, reverse=True)) for factor in factors])
    print("105 raw convolution (exponent, multiplicity):", sorted(raw.items(), reverse=True))
    print("105 normal support:", tuple(sorted(normalized, reverse=True)))
    print(f"105 carry certificate: R={radix}; Q(exponent, coefficient)={sorted(quotient.items(), reverse=True)}")
    print("105 exact factorization Q=-t^4*(t^2-t+1)*(t^8+t^4+1)=-t^4*Phi3*Phi6^2*Phi12")
    print("Carry specialization hostile: Q(phi)=-16phi^8=(-208,-336), norm256; this does not recover factors3,5,7")
    print("Carry recovery is exact; original factor list is an additional parse witness")

    # Fixed local rewrites remain valid at every tested positive/negative scale.
    local_checks = 0
    for k in range(-20, 21):
        local = ({k: 1, k + 1: 1}, {k + 2: 1})
        repeated = ({k: 2}, {k + 1: 1, k - 2: 1})
        for source, target in (local, repeated):
            r, q = carry_certificate(source, target)
            assert evaluate_laurent(source) == evaluate_laurent(target)
            assert recover_raw(target, r, q) == source
            local_checks += 1
    try:
        carry_certificate({0: 1}, {0: 2})
    except ValueError:
        pass
    else:
        raise AssertionError("unequal values accepted as a carry")

    # Independent product normalization and recovery, including composites.
    product_checks = 0
    for a, b in product(range(1, 31), repeat=2):
        source = convolve(parse_phi_word(words[a])[0], parse_phi_word(words[b])[0])
        target = parse_phi_word(words[a * b])[0]
        r, q = carry_certificate(source, target)
        assert recover_raw(target, r, q) == source
        assert evaluate_laurent(source) == (a * b, 0)
        product_checks += 1
    print(f"Carry audits: {local_checks} local rewrites; {product_checks} products a,b=1..30; unequal-value rejection")

    # Associativity transports carry witnesses through different parse trees.
    @lru_cache(maxsize=None)
    def sigma(n):
        return parse_phi_word(canonical_phi_integer(n))[0]

    @lru_cache(maxsize=None)
    def additive_carry(a, b):
        return carry_laurent(polynomial_add(sigma(a), sigma(b)), sigma(a + b))

    @lru_cache(maxsize=None)
    def multiplicative_carry(a, b):
        return carry_laurent(convolve(sigma(a), sigma(b)), sigma(a * b))

    composition_checks = 0
    for a, b, c in product(range(1, 11), repeat=3):
        additive_left = polynomial_add(additive_carry(a, b), additive_carry(a + b, c))
        additive_right = polynomial_add(additive_carry(b, c), additive_carry(a, b + c))
        assert additive_left == additive_right
        multiplicative_left = polynomial_add(multiplicative_carry(a * b, c),
                              convolve(multiplicative_carry(a, b), sigma(c)))
        multiplicative_right = polynomial_add(multiplicative_carry(a, b * c),
                               convolve(sigma(a), multiplicative_carry(b, c)))
        assert multiplicative_left == multiplicative_right
        composition_checks += 1
    print(f"Carry composition: {composition_checks} additive and {composition_checks} multiplicative associativity controls a,b,c=1..10")

    # Phi and Fibonacci readers are related by one unimodular coordinate map.
    bridge_checks = 0
    for n in range(bound + 1):
        digits = zeck.zeckendorf_digits(n)
        lifted = zeck_lift(n, zeck)
        shift = zeck.direct_digit_value(digits, shifted=True)
        assert to_zeck_registers(lifted) == (n, shift)
        assert from_zeck_registers((n, shift)) == lifted
        b = zeck.beatty_register(n)
        assert lifted == (n - 2 * b, b)
        for digit in (0, 1):
            updated_pair = pair_add(pair_mul(lifted, (0, 1)), (digit, 0))
            assert to_zeck_registers(updated_pair) == (shift + digit, n + shift + 2 * digit)
            bridge_checks += 1
    multiplication_checks = addition_checks = 0
    for a, b in product(range(101), repeat=2):
        u, v = zeck_lift(a, zeck), zeck_lift(b, zeck)
        projected_product = to_zeck_registers(pair_mul(u, v))[0]
        assert projected_product == a * b - u[1] * v[1]
        multiplication_checks += 1
        delta = zeck.beatty_register(a) + zeck.beatty_register(b) - zeck.beatty_register(a + b)
        defect = pair_add(pair_add(u, v), pair_mul((-1, 0), zeck_lift(a + b, zeck)))
        assert defect == (-2 * delta, delta)
        addition_checks += 1
    assert pair_mul(zeck_lift(2, zeck), zeck_lift(2, zeck)) == zeck_lift(3, zeck)
    print(f"Unimodular reader bridge: {bridge_checks} append checks n=0..{bound}")
    print(f"Zeckendorf-lift defects: {multiplication_checks} multiplication and {addition_checks} addition checks")
    print("Multiplication hostile: L(2)=phi; L(2)^2=phi^2=L(3), while 2*2=4")

    # A sparse two-unit formula describes Lucas numbers, including composites.
    lucas = [2, 1]
    for _ in range(2, 21):
        lucas.append(lucas[-1] + lucas[-2])
    for k in range(1, 21):
        assert pair_add(phi_power(2 * k), ((-1) ** k, 0)) == pair_mul((lucas[k], 0), phi_power(k))
    assert lucas[4:7] == [7, 11, 18]
    assert evaluate_laurent({6: 1, -6: 1}) == (18, 0)
    print("Lucas controls k=1..20: L4=7, L5=11, L6=18; 18=phi^6+phi^-6 is sparse and composite")

    # The reader computes the odd Collatz numerator by two shifts and one unit.
    collatz_checks = 0
    for n in range(1, 334):
        source = parse_phi_word(words[n])[0]
        raw_numerator = polynomial_add(convolve(source, {2: 1, -2: 1}), {0: 1})
        target = parse_phi_word(canonical_phi_integer(3 * n + 1))[0]
        r, q = carry_certificate(raw_numerator, target)
        assert recover_raw(target, r, q) == raw_numerator
        assert evaluate_laurent(raw_numerator) == (3 * n + 1, 0)
        collatz_checks += 1
    print(f"Collatz numerator: {collatz_checks} exact two-shift-plus-unit carry checks; division guards are separate")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
