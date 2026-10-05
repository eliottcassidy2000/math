"""Exact atomic measure on Collatz renewal/parity words. Standard library only.

Run: python -X utf8 -B 04-computation/experiments/algorithmic_collatz_measure_20261005.py
No Kolmogorov complexities are estimated; no universal termination is assumed.
"""
from fractions import Fraction


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def shortcut(n):
    return (3 * n + 1) // 2 if n & 1 else n // 2


def parity_after_odd(m, length):
    """Works for signed m too; source is 2m-1. Continue through root1."""
    n = shortcut(2 * m - 1)
    bits = []
    for _ in range(length):
        bits.append(n & 1)
        n = shortcut(n)
    return tuple(bits)


def weight(m):
    if type(m) is not int or m < 1:
        raise ValueError("positive integer source index required")
    return Fraction(1, 1 << (2 * m.bit_length() - 1))


def prefix_residue(bits):
    """Affine parity lifting, independent of forward source enumeration.

    After j shortcut steps the current value is (P*n+C)/2**j.
    """
    if any(type(e) is not int or e not in (0, 1) for e in bits):
        raise ValueError("binary integer tuple required")
    p, c, n_residue = 1, 0, 0
    for j, e in enumerate((1,) + tuple(bits)):
        modulus = 1 << (j + 1)
        n_residue = ((e * (1 << j) - c) * pow(p, -1, modulus)) % modulus
        if e:
            p, c = 3 * p, 3 * c + (1 << j)
    return ((n_residue + 1) // 2) % (1 << len(bits))


def cylinder(bits):
    r = prefix_residue(bits)
    return Fraction(1, 4 ** len(bits)) + (weight(r) if r else 0)


def interval_count(lo, hi, residue, modulus):
    return (hi - residue) // modulus - (lo - 1 - residue) // modulus


def independent_residue_mass(residue, length, blocks):
    """Count a congruence in gamma blocks; return sum and whole omitted tail."""
    mass = Fraction(0)
    for k in range(blocks):
        count = interval_count(1 << k, (1 << (k + 1)) - 1, residue, 1 << length)
        mass += count * Fraction(1, 1 << (2 * k + 1))
    return mass, Fraction(1, 1 << blocks)


def run():
    checks = 0
    def check(condition, message):
        nonlocal checks
        require(condition, message)
        checks += 1

    print("Atomic Collatz word measure; exact finite audit; universal coverage OPEN")
    print("mu[prefix] = 4^-L + weight(r) for r>0; otherwise 4^-L")
    for length in range(11):
        row = {}
        total = Fraction(0)
        for m in range(1, (1 << length) + 1):
            bits = parity_after_odd(m, length)
            check(bits not in row, "finite parity map not injective")
            row[bits] = m % (1 << length)
            check(prefix_residue(bits) == row[bits], "affine inverse disagrees")
            exact = cylinder(bits)
            total += exact
            # Equal split of every omitted block after k>=L.
            partial, tail = independent_residue_mass(row[bits], length, length + 5)
            check(exact == partial + tail / (1 << length), "block tail formula")
            check(exact == cylinder(bits + (0,)) + cylinder(bits + (1,)), "measure split")
            check(exact >= Fraction(1, 4 ** length), "full support lower bound")
        check(total == 1, "partition mass")
        print("length", length, "cylinders", len(row), "total", total)

    for m in range(1, 258):
        length = m.bit_length() + 9
        bits = parity_after_odd(m, length)
        check(cylinder(bits) == weight(m) + Fraction(1, 4 ** length), "atom plateau")
        check(cylinder(bits) >= weight(m), "atom disappears")
    for length in (1, 2, 5, 20, 100, 512):
        minus_one = parity_after_odd(0, length)
        check(minus_one == (1,) * length, "negative fixed point")
        check(cylinder(minus_one) == Fraction(1, 4 ** length), "zero-address hostile")
        minus_three = parity_after_odd(-1, length)
        check(cylinder(minus_three) == 3 * Fraction(1, 4 ** length), "negative-index control")

    # Same renewal coordinate, checked directly against accelerated valuations.
    for source in range(1, 256, 2):
        n, bits = source, []
        for _ in range(9):
            z = 3 * n + 1
            a = (z & -z).bit_length() - 1
            bits.extend([0] * (a - 1) + [1])
            n = z >> a
        check(tuple(bits) == parity_after_odd((source + 1) // 2, len(bits)), "renewal map")

    # Infinite sparse address with ones in positions 1,2,4,8,...: liminf b/L=1/2.
    # Finite probes only; the limiting local-dimension proof is in the note.
    sparse = sum(1 << (1 << j) for j in range(10))
    print("sparse-address near-gap (L, bitlength(r), 2*bitlength/L):")
    for j in range(3, 10):
        length = 1 << j
        r = sparse % (1 << length)
        b = r.bit_length()
        check(b == (1 << (j - 1)) + 1, "sparse coordinate")
        print(length, b, Fraction(2 * b, length))

    for m in (1, 2, 14, 128, 257):
        length = m.bit_length() + 3
        bits = parity_after_odd(m, length)
        print("positive source", 2 * m - 1, "atom", weight(m), "L", length,
              "cylinder", cylinder(bits))
    print("negative source -1: cylinder mass4^-L, local dimension2, computable sequence")
    print("positive integer paths: atoms, measure-relative random, effective/local dimension0")
    print("Checks:", checks)
    print("PASS; computation tests finite formulas, not infinite coverage or complexity values")


if __name__ == "__main__":
    run()
