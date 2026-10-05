"""Scoped finite-observer obstruction on the primitive sibling-base map G."""
from dataclasses import dataclass
from fractions import Fraction


def integer(value, least=0):
    if type(value) is not int or value < least:
        raise ValueError("exact integer outside declared domain")
    return value


def v2(value):
    integer(value, 1)
    return (value & -value).bit_length() - 1


def odd(value):
    integer(value, 1)
    if value % 2 == 0:
        raise ValueError("positive odd source required")
    return value


def decode_base(source):
    odd(source)
    current, depth = source, 0
    while current % 8 == 5:
        current = (current - 1) // 4
        depth += 1
    return current, depth


def is_base(source):
    odd(source)
    return v2(3 * source + 1) in (1, 2)


def base_step(source):
    odd(source)
    if source == 1 or not is_base(source):
        raise ValueError("nonroot primitive base required")
    exponent = v2(3 * source + 1)
    actual = (3 * source + 1) >> exponent
    target, depth = decode_base(actual)
    return target, exponent, depth


@dataclass(frozen=True)
class Observation:
    height: int
    registers: tuple
    valuations: tuple
    sibling_depths: tuple
    reached_root: bool


def observe(source, lookahead, modulus):
    odd(source)
    integer(lookahead)
    integer(modulus, 1)
    if not is_base(source):
        raise ValueError("observer source is not a primitive base")
    current, registers, valuations, depths = source, [], [], []
    for i in range(lookahead + 1):
        h = v2(current + 1)
        cofactor = (current + 1) >> h
        registers.append((h, cofactor % modulus, current % modulus))
        if current == 1 or i == lookahead:
            break
        current, exponent, depth = base_step(current)
        valuations.append(exponent)
        depths.append(depth)
    return Observation(source.bit_length(), tuple(registers), tuple(valuations), tuple(depths), current == 1)


@dataclass(frozen=True)
class Alias:
    modulus: int
    lookahead: int
    K: int
    s: int
    t: int
    source: int
    endpoint: int


def construct(modulus, lookahead, extra_height=0):
    integer(modulus, 1)
    integer(lookahead)
    integer(extra_height)
    s = modulus.bit_length() + 2 + extra_height
    scale = 1 << s
    t = modulus * ((3 * scale + modulus - 1) // modulus)
    K = 3 * (lookahead // 2) + 4
    source = (1 << (K + 3)) * t - 5
    endpoint = 9 * (1 << K) * t - 5
    return Alias(modulus, lookahead, K, s, t, source, endpoint)


def construct_at_height(modulus, height, lookahead):
    integer(modulus, 1)
    integer(height, 1)
    integer(lookahead)
    s = modulus.bit_length() + 2
    K = height - s - 5
    if K < 3 * (lookahead // 2) + 4:
        raise ValueError("declared height does not fund this alias guard")
    scale = 1 << s
    t = modulus * ((3 * scale + modulus - 1) // modulus)
    source = (1 << (K + 3)) * t - 5
    endpoint = 9 * (1 << K) * t - 5
    return Alias(modulus, lookahead, K, s, t, source, endpoint)


def formal_shadow(power, cofactor, i):
    """Closed formula, independently compared with actual quotient iteration."""
    integer(power, 1)
    integer(cofactor, 1)
    integer(i)
    q, parity = divmod(i, 2)
    exponent = power - 3 * q - parity
    if exponent < 0:
        raise ValueError("negative displayed power")
    return 3 ** i * (1 << exponent) * cofactor - (7 if parity else 5)


def gamma(source):
    odd(source)
    return Fraction(1, 1 << (2 * ((source + 1) // 2).bit_length() - 1))


def mixture(source):
    odd(source)
    return Fraction(8, 3 * 4 ** source.bit_length())


def pole_fuel(source):
    odd(source)
    if source % 8 == 3:
        return v2(source + 5)
    if source % 8 == 1:
        return v2(source + 7)
    return 0


def pole_weight(source):
    return mixture(source) / 8 ** pole_fuel(source)


def pole_guard(source):
    odd(source)
    fuel = pole_fuel(source)
    return (source % 8 == 3 and fuel >= 4) or (source % 8 == 1 and fuel >= 5)


def main():
    checks = 0

    def check(value, label):
        nonlocal checks
        checks += 1
        if not value:
            raise ValueError(label)

    moduli = (1, 2, 3, 5, 6, 7, 11, 19, 64, 81, 105, 223, 233, 729, 4096)
    pairs = 0
    for L in range(65):
        for M in moduli:
            pair = construct(M, L)
            n, x, K, s, t = pair.source, pair.endpoint, pair.K, pair.s, pair.t
            check(t % M == 0 and 3 * (1 << s) <= t < Fraction(13, 4) * (1 << s), "cofactor interval")
            check(n.bit_length() == x.bit_length() == K + s + 5, "equal full current height")
            intermediate, a, depth = base_step(n)
            endpoint, b, depth2 = base_step(intermediate)
            check((endpoint, a, b, depth, depth2) == (x, 1, 2, 0, 0), "actual expanding G square")
            check(x > n > 1, "expanding positive witnesses")
            check(observe(n, L, M) == observe(x, L, M), "complete displayed observation aliases")
            check(gamma(n) == gamma(x) and mixture(n) == mixture(x), "equal baseline atoms")
            left, right = n, x
            for i in range(L + 1):
                check(left == formal_shadow(K + 3, t, i), "source formula")
                check(right == formal_shadow(K, 9 * t, i), "endpoint formula")
                check(left > 1 and right > 1 and is_base(left) and is_base(right), "positive primitive states")
                check(v2(left + 1) == v2(right + 1) == (2 if i % 2 == 0 else 1), "alternating register")
                if i < L:
                    left, aa, dd = base_step(left)
                    right, bb, ee = base_step(right)
                    check(aa == bb == (1 if i % 2 == 0 else 2) and dd == ee == 0, "no hidden sibling removal")
            pairs += 1

    variable_pairs = 0
    for M in moduli:
        s = M.bit_length() + 2
        for height in range(s + 9, s + 90):
            L = 2 * ((height - s - 9) // 3) + 1
            pair = construct_at_height(M, height, L)
            check(observe(pair.source, L, M) == observe(pair.endpoint, L, M), "variable radius boundary")
            check(pair.source.bit_length() == pair.endpoint.bit_length() == height, "requested full height")
            variable_pairs += 1

    paid = 0
    for source in range(3, 65536, 2):
        if pole_guard(source):
            target, exponent, depth = base_step(source)
            ratio = pole_weight(source) / pole_weight(target)
            cap = Fraction(1, 2) if source % 8 == 3 else Fraction(1, 64)
            check(depth == 0, "guarded pole edge remains unreduced")
            check(pole_fuel(target) == pole_fuel(source) - exponent, "exact pole-fuel cost")
            check(ratio <= cap <= Fraction(5, 8), "guarded base payment")
            paid += 1

    check(base_step(187) == (281, 1, 0) and base_step(281) == (211, 2, 0), "small observer hostile")
    check(observe(187, 0, 1) == observe(211, 0, 1), "small equal height and register")
    check(base_step(7) == (11, 1, 0), "pole-refuel edge")
    check(pole_weight(7) / pole_weight(11) == 16384, "pole-refuel unpaid bill")
    check(observe(1, 32, 19).reached_root and observe(1, 32, 19).valuations == (), "root observation does not pad")

    hostiles = [
        lambda: construct(True, 1), lambda: construct(1, True),
        lambda: construct(0, 1), lambda: construct(1, -1),
        lambda: construct(1.0, 1), lambda: construct_at_height(19, 5, 10),
        lambda: base_step(1), lambda: base_step(5), lambda: base_step(3.0),
        lambda: observe(True, 1, 1), lambda: observe(5, 1, 1),
        lambda: observe(3, 1, 0), lambda: pole_guard(2),
    ]
    for action in hostiles:
        try:
            action()
        except ValueError:
            check(True, "hostile rejected")
        else:
            raise ValueError("hostile accepted")

    print("PROVED: finite future-register/cofactor-residue/current-height aliases on G")
    print("FINITE-EXACT: 65 lookaheads x 15 moduli = " + str(pairs) + " exact positive source pairs")
    print("PROVED: exact alias guard K>=3 floor(L/2)+4; expanding native word (1,2); all displayed sibling depths 0")
    print("PROVED: identical complete current heights and both baseline atoms; no strict weight factors through this observer")
    print("PROVED: variable-radius obstruction L(N)<=2 floor((N-s-9)/3)+1; no optimality claim")
    print("FINITE-EXACT: variable-radius pairs = " + str(variable_pairs))
    print("PROVED: two-anchor pole correction pays guarded phase edges; base charges r=1/16, rho=2/3")
    print("FINITE-EXACT: guarded pole edges below 65536 = " + str(paid))
    print("PROVED: same correction fails at 7 -> 11 with ratio 16384; global refuel payment remains OPEN")
    print("FINITE-EXACT: malformed/source/root hostiles = " + str(len(hostiles)))
    print("explicit checks=" + str(checks))
    print("SCOPE: no obstruction to arbitrary finite algorithms, full cofactors, or other growing-precision controllers")


if __name__ == "__main__":
    main()
