"""Standalone exact controls for Bernoulli endpoint and Faulhaber claims.

Run: python bernoulli_endpoint_em_audit_20260925.py
Uses only Python standard library; no floating-point arithmetic.
Finite universe is printed. These controls do not prove Collatz or Erdos-Moser.
"""
from fractions import Fraction
from math import comb, factorial


def psi(x):
    x = Fraction(x)
    return x - (x.numerator // x.denominator) - Fraction(1, 2)


def jump(z, m):
    return Fraction(1, m) + psi(Fraction(z - 1, m)) - psi(Fraction(z, m))


def v2(z):
    assert z != 0
    z = abs(z)
    return (z & -z).bit_length() - 1


def bernoulli_up_to(nmax):
    result = [Fraction(1)]
    for n in range(1, nmax + 1):
        result.append(-sum(Fraction(comb(n + 1, j)) * result[j]
                           for j in range(n)) / (n + 1))
    return result


jump_checks = 0
for z in range(-100, 101):
    for m in range(2, 101):
        actual = jump(z, m)
        expected = int(z % m == 0)
        assert actual == expected, (z, m, actual, expected)
        jump_checks += 1

valuation_checks = 0
for z in range(-100, 101):
    partial = Fraction(0)
    for cutoff in range(1, 11):
        partial += jump(z, 2 ** cutoff)
        expected = cutoff if z == 0 else min(v2(z), cutoff)
        assert partial == expected, (z, cutoff, partial, expected)
        valuation_checks += 1

# Hostile endpoint control: setting the sawtooth to zero at integers
# produces half a mass on both sides of the boundary, not a divisor indicator.
def psi_mid(x):
    x = Fraction(x)
    return Fraction(0) if x.denominator == 1 else psi(x)


def jump_mid(z, m):
    return Fraction(1, m) + psi_mid(Fraction(z - 1, m)) - psi_mid(Fraction(z, m))


assert jump_mid(0, 2) == Fraction(1, 2)
assert jump_mid(1, 2) == Fraction(1, 2)
assert jump(0, 2) == 1 and jump(1, 2) == 0

# Compute Bernoulli numbers by their defining binomial recurrence;
# compare Faulhaber expression with an independently computed direct power sum.
B = bernoulli_up_to(30)
assert B[1] == Fraction(-1, 2)
assert B[2] == Fraction(1, 6)
assert all(B[n] == 0 for n in range(3, 31, 2))
faulhaber_checks = 0
for m in range(2, 40):
    for k in range(2, 30):
        lhs = Fraction(sum(j ** k for j in range(1, m)), m ** k)
        rhs = Fraction(m, k + 1) - Fraction(1, 2)
        for r in range(1, k // 2 + 1):
            falling = factorial(k) // factorial(k - 2 * r + 1)
            rhs += B[2 * r] * Fraction(falling, factorial(2 * r) * m ** (2 * r - 1))
        assert lhs == rhs, (m, k, lhs, rhs)
        faulhaber_checks += 1

print(f'PASS: {jump_checks} exact indicator checks; z=-100..100, m=2..100.')
print(f'PASS: {valuation_checks} exact partial-valuation checks; z=-100..100, cutoff=1..10.')
print('PASS: endpoint hostile controls; midpoint convention fails as expected.')
print(f'PASS: {faulhaber_checks} exact Faulhaber checks; m=2..39, k=2..29.')
print('Bernoulli convention: B1=-1/2. Universe checks are finite; identities have separate algebraic proofs.')
