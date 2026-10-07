"""Explicit finite-bank pair-chain absorption cylinders; no orbit conjecture.

All algorithms use exact integers/Fractions. The finite universe is signed
translations +/-1..512, with three direct positive-integer lifts per cylinder.
"""
from fractions import Fraction as F


def need(ok, message):
    if not ok:
        raise ValueError(message)


def exact_integer(n):
    need(type(n) is int, "exact integer required")
    return n


def valuation(n):
    exact_integer(n)
    need(n > 0, "positive integer required")
    return (n & -n).bit_length() - 1


def bits_checked(bits):
    need(type(bits) is tuple and all(type(b) is int and b in (0, 1) for b in bits),
         "tuple of exact parity bits required")
    return bits


def three(k):
    exact_integer(k)
    return F(3 ** k) if k >= 0 else F(1, 3 ** (-k))


def transition(k, c, bit):
    exact_integer(k)
    need(type(c) in (int, F), "exact rational translation required")
    c = F(c)
    need(c.denominator % 2 == 1, "translation must be 2-adically integral")
    need(type(bit) is int and bit in (0, 1), "exact parity bit required")
    sigma = c.numerator % 2
    if not bit:
        return (k, c / 2) if not sigma else (k + 1, (3 * c + 1) / 2)
    return ((k, (3 * c + 1 - three(k)) / 2) if not sigma
            else (k - 1, (c - three(k - 1)) / 2))


def positive_certificate(e):
    exact_integer(e)
    need(e > 0, "positive translation required")
    bits = []
    while e:
        if e % 2 == 0:
            bits.append(0)
            e //= 2
        else:
            a = valuation(3 * e + 1)
            bits.extend([0] * a + [1])
            e = (3 * e + 1 - (1 << a)) // (1 << (a + 1))
    return tuple(bits)


def certificate(e):
    """Driving-y bits that force T^L(y+e)=T^L(y). Empty for e=0."""
    exact_integer(e)
    if e == 0:
        return ()
    positive = positive_certificate(abs(e))
    if e > 0:
        return positive
    k, c = 0, F(-e)
    mirrored = []
    for bit in positive:
        mirrored.append(bit ^ (c.numerator % 2))
        k, c = transition(k, c, bit)
    need((k, c) == (0, 0), "positive construction must absorb")
    return tuple(mirrored)


def source_residue(bits):
    """Unique driving-source class modulo 2^len(bits)."""
    bits_checked(bits)
    r, p, carry = 0, 1, 0
    for i, bit in enumerate(bits):
        current = (p * r + carry) // (1 << i)
        r += (bit ^ (current % 2)) * (1 << i)
        if bit:
            p, carry = 3 * p, 3 * carry + (1 << i)
    return r


def pair_receipt(e):
    exact_integer(e)
    bits = certificate(e)
    k, c = 0, F(e)
    for bit in bits:
        k, c = transition(k, c, bit)
    need((k, c) == (0, 0), "receipt failed absorption")
    return bits, source_residue(bits), F(1, 1 << len(bits))


def terras(n):
    exact_integer(n)
    return (3 * n + 1) // 2 if n % 2 else n // 2


def direct_pair(e, y, bits):
    exact_integer(e)
    exact_integer(y)
    bits_checked(bits)
    need(y > 0 and y + e > 0, "both integer orbits must be positive")
    need(y % (1 << len(bits)) == source_residue(bits), "incompatible source phase")
    x = y + e
    for bit in bits:
        need(x != y, "receipt must have no earlier equal-time merge")
        need(y % 2 == bit, "literal driver bit differs")
        x, y = terras(x), terras(y)
    need(x == y, "literal pair did not merge")
    return x


def integer_bank_bound(m):
    exact_integer(m)
    need(m >= 1, "positive bank radius required")
    return 5 * (m - 1).bit_length() + 3


def main():
    checks = 0

    def check(ok, message):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(message)

    largest = (0, 0)
    direct = 0
    for e in list(range(-512, 0)) + list(range(1, 513)):
        word, residue, mass = pair_receipt(e)
        length, radius = len(word), abs(e)
        if length > largest[1]:
            largest = (e, length)
        check(length == len(positive_certificate(radius)), "mirror preserves cost")
        check(length <= integer_bank_bound(radius), "rational looser bank bound")
        # Equivalent to L <= kappa*log2(e)+3, without computing any logarithm.
        q = max(0, length - 3)
        check(4 ** q <= radius ** 2 * 3 ** q, "exact sharp logarithmic inequality")
        check(mass >= F(1, 1 << integer_bank_bound(radius)), "explicit bank probability")
        k, c = 0, F(e)
        for i, bit in enumerate(word):
            check((k, c) != (0, 0), "constructed chain has no early absorption")
            k, c = transition(k, c, bit)
        check((k, c) == (0, 0), "chain absorption")
        modulus = 1 << length
        lift = max(1, (max(0, -e) - residue) // modulus + 1)
        for j in range(3):
            y = residue + (lift + j) * modulus
            endpoint = direct_pair(e, y, word)
            check(endpoint > 0, "positive direct common endpoint")
            direct += 1
        # Swapping the positive driver shifts the source class by radius.
        if e < 0:
            rp = source_residue(positive_certificate(radius))
            check(residue == (rp + radius) % modulus, "signed source translation")

    # Independently enumerate every parity class through length 9.
    phase_controls = 0
    for length in range(10):
        residues = set()
        for mask in range(1 << length):
            word = tuple((mask >> i) & 1 for i in range(length))
            r = source_residue(word)
            residues.add(r)
            x = r
            for bit in word:
                check(x % 2 == bit, "direct parity inversion")
                x = terras(x)
            phase_controls += 1
        check(len(residues) == 1 << length, "complete phase bijection")

    check(pair_receipt(0) == ((), 0, F(1)), "identity empty boundary")
    w, r, mass = pair_receipt(1)
    check((w, r, mass) == ((0, 0, 1), 4, F(1, 8)), "sharp e=1 control")
    check(direct_pair(1, 4, w) == 2, "literal 4/5 merge")
    check((terras(2), terras(1)) == (1, 2)
          and (terras(1), terras(2)) == (2, 1), "exact closed hostile two-cycle")
    x, y = 2, 1
    for _ in range(20):
        check(x != y and {x, y} == {1, 2}, "same-time nonmerge hostile")
        x, y = terras(x), terras(y)

    # Equal Terras time can precede the first common odd state by many steps.
    x, y, k, c = 2417, 805, 1, F(2)
    equal_times, odd_equal_times, identity_times = [], [], []
    for time in range(17):
        if x == y:
            equal_times.append((time, x))
            if x % 2:
                odd_equal_times.append((time, x))
        if (k, c) == (0, 0):
            identity_times.append(time)
        k, c = transition(k, c, y % 2)
        x, y = terras(x), terras(y)
    check(equal_times[0] == (9, 128) and identity_times[0] == 9,
          "Mersenne equal-T merge clock")
    check(odd_equal_times[0] == (16, 1), "Mersenne first common odd clock")
    recovered = []
    for source in (2417, 805):
        word, n = [], source
        while n != 1:
            a = valuation(3 * n + 1)
            word.append(a)
            n = (3 * n + 1) >> a
        recovered.append(tuple(word))
    check(recovered == [(2, 6, 8), (4, 1, 1, 10)], "literal odd words")
    check(sum(recovered[0]) == sum(recovered[1]) == 9 + valuation(128) == 16,
          "odd-clock delay is common-state valuation")

    hostiles = [lambda: certificate(True), lambda: certificate(1.0),
                lambda: integer_bank_bound(0), lambda: integer_bank_bound(True),
                lambda: source_residue((True,)), lambda: source_residue([0]),
                lambda: direct_pair(1, 1, w), lambda: direct_pair(-1, 1, ()),
                lambda: transition(0, F(1, 2), 0), lambda: positive_certificate(0)]
    for action in hostiles:
        try:
            action()
        except ValueError:
            check(True, "typed/phase hostile")
        else:
            check(False, "hostile accepted")

    check(3 ** 5 < 4 ** 4, "kappa strictly below five")
    print("PROVED: explicit signed finite-bank absorption cylinder; no individual-integer theorem")
    print("universe: e=-512..-1,1..512; translations=1024")
    print("direct_positive_integer_lifts=" + str(direct))
    print("all_parity_words_lengths0..9=" + str(phase_controls))
    print("maximum_constructed_length=" + str(largest[1]) + "; first_translation=" + str(largest[0]))
    print("L(e)<=kappa*log2(abs(e))+3; kappa=2/log2(4/3)<5")
    print("bank_radius512: exact_integer_deadline=" + str(integer_bank_bound(512)))
    print("sharp_boundary: e=1, bits001, source4mod8, mass1/8")
    print("hostile: source1 with e1 stays the out-of-phase pair(1,2)")
    print("clock_hostile: pair2417/805 firstTerrasMerge9at128; firstOddMerge16at1")
    print("typed_and_phase_hostiles=" + str(len(hostiles)))
    print("checks=" + str(checks))


if __name__ == "__main__":
    main()
