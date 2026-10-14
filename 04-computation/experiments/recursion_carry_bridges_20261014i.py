"""Exact triangular quotient of the Mahler ceiling map; no Z-number oracle.

All public finite readers retain their precision and exact input types.
Run with python -B (and -O); stdout is deterministic JSON.
"""
from fractions import Fraction
from itertools import product
from math import isqrt
import json


def integer(x, minimum=None):
    if type(x) is not int or (minimum is not None and x < minimum):
        raise ValueError("exact integer outside domain")
    return x


def rational(x):
    if type(x) not in (int, Fraction):
        raise ValueError("exact rational required")
    x = Fraction(x)
    if x.denominator % 2 == 0:
        raise ValueError("odd denominator required")
    return x


def bit(x):
    integer(x, 0)
    if x > 1:
        raise ValueError("bit required")
    return x


def binary_word(w):
    if type(w) is not tuple:
        raise ValueError("tuple required")
    for x in w:
        bit(x)
    return w


def parity(x):
    return rational(x).numerator % 2


def triangular(x):
    x = rational(x)
    return x * (x + 1) / 2


def mahler(x):
    x = rational(x)
    return (3 * x + parity(x)) / 2


def fold(x):
    """Return the exact even-representative coordinate and retained sheet."""
    x = rational(x)
    s = parity(x)
    return ((-x - 1) if s else x) / 2, s


def unfold(z, sheet):
    z = rational(z)
    bit(sheet)
    return -2 * z - 1 if sheet else 2 * z


def qchart(z):
    z = rational(z)
    return z * (2 * z + 1)


def folded_step(z):
    z = rational(z)
    return -(3 * z + 1) / 2 if parity(z) else 3 * z / 2


def marked_step(z, sheet):
    rational(z)
    bit(sheet)
    return folded_step(z), sheet ^ parity(z)


def inverse_q_residue(a, bits):
    """Unique z mod 2^bits with z(2z+1)=a mod 2^bits."""
    integer(bits, 0)
    integer(a, 0)
    if a >= 1 << bits:
        raise ValueError("noncanonical residue")
    z = 0
    for j in range(bits):
        if (z * (2 * z + 1) - a) % (1 << (j + 1)):
            z += 1 << j
    return z


def phi_residue(a, output_bits):
    """Read Phi(a) mod 2^k from canonical a mod 2^(k+1)."""
    integer(output_bits, 0)
    z = inverse_q_residue(a, output_bits + 1)
    z1 = -(3 * z + 1) // 2 if z % 2 else 3 * z // 2
    return z1 * (2 * z1 + 1) % (1 << output_bits)


def parity_cylinder(a, depth):
    """Read exactly depth triangular parity bits from a mod 2^depth."""
    integer(depth, 0)
    integer(a, 0)
    if a >= 1 << depth:
        raise ValueError("noncanonical source residue")
    out = []
    for remaining in range(depth, 0, -1):
        out.append(a % 2)
        a = phi_residue(a, remaining - 1)
    return tuple(out)


def lift_derivative(w, initial_sheet):
    binary_word(w)
    bit(initial_sheet)
    c = [initial_sheet]
    for d in w:
        c.append(c[-1] ^ d)
    return tuple(c)


def primitive_period(w):
    binary_word(w)
    if not w:
        raise ValueError("nonempty periodic word required")
    return next(p for p in range(1, len(w) + 1)
                if len(w) % p == 0 and w == w[:p] * (len(w) // p))


def periodic_lift(w, initial_sheet):
    binary_word(w)
    bit(initial_sheet)
    if not w:
        raise ValueError("nonempty periodic word required")
    full = w * (2 if sum(w) % 2 else 1)
    c = lift_derivative(full, initial_sheet)
    if c[-1] != c[0]:
        raise ArithmeticError("closure failure")
    return c[:-1]


def tail(w, shift=0):
    """Exact weighted tail of a supplied nonempty periodic carry word."""
    binary_word(w)
    integer(shift, 0)
    if not w:
        raise ValueError("nonempty word required")
    m = len(w)
    base = sum(Fraction(w[(j + shift) % m]) * Fraction(2, 3) ** (j + 1)
               for j in range(m))
    return base / (1 - Fraction(2, 3) ** m)


def affine_matrix(e):
    """Pullback on the column (1,n,T(n)); branch must remain authenticated."""
    bit(e)
    return ((Fraction(1), Fraction(0), Fraction(0)),
            (Fraction(e, 2), Fraction(3, 2), Fraction(0)),
            (Fraction(3 * e, 8), Fraction(6 * e - 3, 8), Fraction(9, 4)))


def matvec(a, v):
    return tuple(sum(x * y for x, y in zip(row, v)) for row in a)


def divisor_counts(n):
    integer(n, 2)
    ds = [d for d in range(2, n) if n % d == 0]
    prime = lambda d: all(d % k for k in range(2, isqrt(d) + 1))
    sf = lambda d: all(d % (k * k) for k in range(2, isqrt(d) + 1))
    return len(ds), sum(prime(d) for d in ds), sum(sf(d) for d in ds)


def encode(x):
    if type(x) is Fraction:
        return {"numerator": x.numerator, "denominator": x.denominator}
    if type(x) is dict:
        return {k: encode(v) for k, v in x.items()}
    if type(x) in (tuple, list):
        return [encode(v) for v in x]
    return x


def main():
    checks = 0

    def check(condition):
        nonlocal checks
        checks += 1
        if not condition:
            raise ArithmeticError("control failed")

    def rejects(f):
        try:
            f()
        except ValueError:
            check(True)
        else:
            check(False)

    for den in (1, 3, 5, 7, 19):
        for num in range(-128, 129):
            x = Fraction(num, den)
            z, s = fold(x)
            check(unfold(z, s) == x)
            check(triangular(x) == qchart(z))
            check(mahler(-1 - x) == -1 - mahler(x))
            check(unfold(*marked_step(z, s)) == mahler(x))
            check(qchart(folded_step(z)) == triangular(mahler(x)))
            check(matvec(affine_matrix(parity(x)), (1, x, triangular(x))) ==
                  (1, mahler(x), triangular(mahler(x))))
            check(mahler(2 * x) - 2 * mahler(x) == -parity(x))
            check(triangular(2 * mahler(x)) - triangular(mahler(2 * x)) ==
                  parity(x) * (3 * x + 1))

    for bits in range(0, 12):
        m = 1 << bits
        qimages = [z * (2 * z + 1) % m for z in range(m)]
        check(len(set(qimages)) == m)
        for z, a in enumerate(qimages):
            check(inverse_q_residue(a, bits) == z)
        counts = [0] * m
        for a in range(2 * m):
            counts[phi_residue(a, bits)] += 1
        check(counts == [2] * m)
        check(len({parity_cylinder(a, bits) for a in range(m)}) == m)

    for k in range(1, 41):
        n = 1 << (k + 1)
        check(triangular(n) % (1 << k) == 0)
        check(triangular(mahler(n)) % (1 << k) == 1 << (k - 1))
    for n in range(0, 257):
        a = triangular(n)
        c, d = [], []
        x = Fraction(n)
        for _ in range(13):
            c.append(parity(x))
            d.append(parity(triangular(x)))
            x = mahler(x)
        check(tuple(d[:-1]) == tuple(c[j] ^ c[j + 1] for j in range(12)))
        check(lift_derivative(tuple(d[:-1]), c[0]) == tuple(c))
        check(parity_cylinder(int(a) % 4096, 12) == tuple(d[:-1]))

    primitive_controls = 0
    for m in range(1, 9):
        for w in product((0, 1), repeat=m):
            if primitive_period(w) != m:
                continue
            primitive_controls += 1
            for s in (0, 1):
                c = periodic_lift(w, s)
                check(primitive_period(c) == m * (2 if sum(w) % 2 else 1))
                check(tuple(c[j] ^ c[(j + 1) % len(c)] for j in range(len(c))) ==
                      w * (len(c) // m))

    cycles = []
    for x0 in (Fraction(-9, 19), Fraction(-10, 19)):
        x, c, ts = x0, [], []
        for _ in range(3):
            c.append(parity(x))
            ts.append(triangular(x))
            x = mahler(x)
        check(x == x0)
        cycles.append({"source": x0, "carries": c, "triangular_orbit": ts,
                       "tails": [tail(tuple(c), j) for j in range(3)]})
    check(cycles[0]["triangular_orbit"] == cycles[1]["triangular_orbit"])
    check(all(t < 1 for t in cycles[0]["tails"]))
    check(all(t > 1 for t in cycles[1]["tails"]))
    check(all(a + b == 2 for a, b in zip(cycles[0]["tails"], cycles[1]["tails"])))
    check(mahler(Fraction(-2, 5)) == Fraction(-3, 5))
    check(mahler(Fraction(-3, 5)) == Fraction(-2, 5))
    check(triangular(Fraction(-2, 5)) == triangular(Fraction(-3, 5)) == Fraction(-3, 25))

    examples = {n: divisor_counts(n) for n in (2, 4, 8, 60, 240)}
    check(examples == {2: (0, 0, 0), 4: (1, 1, 1), 8: (2, 1, 1),
                       60: (10, 3, 7), 240: (18, 3, 7)})
    fs, f, g = [], 0, 1
    for _ in range(61):
        fs.append(f)
        f, g = g, f + g
    states = [(fs[r] % 10, fs[r + 1] % 10) for r in range(60)]
    tri = [r * (r + 1) // 2 % 30 for r in range(60)]
    check(len(set(states)) == 60)
    check(len({(fs[r] % 10, tri[r]) for r in range(60)}) == 48)
    check(len({(fs[r] % 10, tri[r], r % 2) for r in range(60)}) == 60)
    common = sorted({v for v in fs[:41] if v > 0 and isqrt(8 * v + 1) ** 2 == 8 * v + 1})
    check(common == [1, 3, 21, 55])
    check(7 * 15 == 105 and 105 not in fs)
    check(next(p for p in range(1, 25) if all((r+p)*(r+p+1)//2 % 8 == r*(r+1)//2 % 8
                                             for r in range(32))) == 16)

    for f in (lambda: mahler(True), lambda: triangular(1.0),
              lambda: fold(Fraction(1, 2)), lambda: unfold(1, True),
              lambda: inverse_q_residue(1, False), lambda: inverse_q_residue(4, 2),
              lambda: phi_residue(4, 1), lambda: parity_cylinder(2, 1),
              lambda: lift_derivative([0, 1], 0), lambda: lift_derivative((True,), 0),
              lambda: periodic_lift((), 0), lambda: tail((0, 2)),
              lambda: affine_matrix(1.0), lambda: divisor_counts(1)):
        rejects(f)

    return {"status": "PROVED scoped + FINITE-EXACT controls",
            "checks": checks, "primitive_derivative_words": primitive_controls,
            "finite_fibonacci_indices": [0, 40], "positive_common_values": common,
            "divisor_examples_F_U_S": examples,
            "sheet_loss_hostile": cycles,
            "triangular_fixed_points": [Fraction(0), Fraction(-3, 25)],
            "precision": "Phi mod 2^k needs a mod 2^(k+1), sharply; no wrapped fixed-level graph",
            "universes": {"rational_numerators": [-128, 128], "odd_denominators": [1, 3, 5, 7, 19],
                          "finite_dyadic_bits": [0, 11], "periodic_word_lengths": [1, 8],
                          "integer_itinerary_sources": [0, 256]},
            "scope": "The marked quotient is exact; strict real safety and positive integer stabilization remain separate. No Mahler or Collatz closure."}


if __name__ == "__main__":
    print(json.dumps(encode(main()), indent=2, sort_keys=True))
