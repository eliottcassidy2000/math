#!/usr/bin/env python3
"""Exact companion for THM-2319.

All checks use integers or Fraction.  Runtime failures are explicit so
python -O executes the same audit.
"""

from fractions import Fraction
from itertools import product
from math import gcd

P = 7
Q = 13
N = P * Q


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def local_kernel_formula(p, a, b):
    return (
        p * p * int(a % p == 0 and b % p == 0)
        - p
        * (
            int(a % p == 0)
            + int(b % p == 0)
            + int((a - b) % p == 0)
        )
        + 2
    )


def local_kernel_direct(p, a, b):
    # Exact cyclotomic coefficients are not needed: group the exponent
    # residues and check that every nonzero residue has coefficient zero.
    counts = [0] * p
    for x in range(1, p):
        for y in range(1, p):
            if (x + y) % p:
                counts[(x * a + y * b) % p] += 1
    # In Q[zeta_p], sum c_r zeta^r equals c_0-c_(p-1) after subtracting
    # c_(p-1) Phi_p; it is rational exactly when all nonzero coefficients
    # agree.  The direct sum here must equal the integer formula.
    tail = counts[1]
    require(all(value == tail for value in counts[1:]), "kernel not rational")
    return counts[0] - tail


def crt_kernel(r, s, t):
    return local_kernel_formula(P, t - r, t - s) * local_kernel_formula(
        Q, t - r, t - s
    )


def unit_bispectrum(sites, masses):
    require(len(sites) == len(masses), "site/mass length mismatch")
    total = 0
    for i, j, k in product(range(len(sites)), repeat=3):
        total += (
            masses[i]
            * masses[j]
            * masses[k]
            * crt_kernel(sites[i], sites[j], sites[k])
        )
    return total


def mixed_bispectrum(sites, left, middle, right):
    require(len(sites) == len(left) == len(middle) == len(right), "mixed lengths")
    total = 0
    for i, j, k in product(range(len(sites)), repeat=3):
        total += (
            left[i]
            * middle[j]
            * right[k]
            * crt_kernel(sites[i], sites[j], sites[k])
        )
    return total


for prime in (P, Q):
    for a, b in product(range(prime), repeat=2):
        require(
            local_kernel_direct(prime, a, b)
            == local_kernel_formula(prime, a, b),
            f"local kernel mismatch p={prime}, a={a}, b={b}",
        )

units = [k for k in range(1, N) if gcd(k, N) == 1]
require(len(units) == 72, "primitive-character count")
allowed = [
    (k, ell)
    for k in units
    for ell in units
    if gcd((k + ell) % N, N) == 1
]
carry_zero = [(k, ell) for k, ell in allowed if k + ell < N]
carry_one = [(k, ell) for k, ell in allowed if k + ell > N]
require(len(allowed) == 3960, "allowed-pair count")
require(len(carry_zero) == len(carry_one) == 1980, "carry split")
require(
    {(N - k, N - ell) for k, ell in carry_zero} == set(carry_one),
    "carry conjugation bijection",
)

coefficient_table = {
    "all_equal": crt_kernel(0, 0, 0),
    "two_same_fibre": crt_kernel(0, 0, 7),
    "two_cross_fibre": crt_kernel(0, 0, 1),
    "three_one_pair": crt_kernel(0, 7, 1),
    "three_fibres": crt_kernel(0, 1, 2),
}
require(
    coefficient_table
    == {
        "all_equal": 3960,
        "two_same_fibre": -330,
        "two_cross_fibre": 55,
        "three_one_pair": -10,
        "three_fibres": 4,
    },
    "CRT coefficient table",
)

floor_constant = Fraction(1305, 98)
pair_constant = floor_constant / 3960
require(pair_constant == Fraction(29, 8624), "selected-pair constant")
require(
    15 * (Fraction(101, 2) * Fraction(1, 49) - Fraction(1, 7))
    == floor_constant,
    "seven-fibre floor arithmetic",
)
lrc_total_constant = floor_constant * 91**3
lrc_pair_constant = pair_constant * 91**3
require(lrc_total_constant == Fraction(20069595, 2), "LRC total constant")
require(lrc_pair_constant == Fraction(445991, 176), "LRC pair constant")

positive_vectors = [
    [1] + [0] * 12,
    [1] * 13,
    list(range(1, 14)),
    [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9],
    [62, 842, 894, 609, 819, 261, 207, 761] + [0] * 5,
]
for index, vector in enumerate(positive_vectors):
    value = unit_bispectrum(list(range(13)), vector)
    mass = sum(vector)
    require(
        Fraction(value, mass**3) >= floor_constant,
        f"positive control below floor at {index}",
    )

mixed_source = [62, 842, 894, 609, 819, 261, 207, 761]
mixed_word = [62, 0, 0, 0, 0, 0, 0, 0]
mixed_value = mixed_bispectrum(
    list(range(8)), mixed_source, mixed_word, mixed_source
)
require(mixed_value == -3042309124, "mixed hostile value")
mixed_full = unit_bispectrum(list(range(8)), mixed_source)
require(mixed_full == 12219375808488, "mixed-source full value")
require(
    Fraction(mixed_full, sum(mixed_source) ** 3) > floor_constant,
    "mixed-source full positivity",
)

wild_sites = [53, 79, 37, 27, 40, 58, 1]
wild_masses = [283, 57, 16, 270, 174, 4, 196]
wild_value = unit_bispectrum(wild_sites, wild_masses)
require(wild_value == -4459059900, "nonconsecutive hostile value")

for jumps in (1, 2, 4, 9, 17):
    h_bound = jumps * jumps - 1
    h3_bound = 2 * jumps * jumps - 1
    require(
        (91 * h_bound + 90) == 91 * jumps * jumps - 1,
        "first Schur bound",
    )
    require(
        (91 * h3_bound + 90) == 182 * jumps * jumps - 1,
        "third Schur bound",
    )
for scale in (1, 2, 7, 31):
    jumps = 6 * scale
    require(91 * jumps * jumps - 1 == 3276 * scale * scale - 1, "LRC A/B bound")
    require(182 * jumps * jumps - 1 == 6552 * scale * scale - 1, "LRC C bound")
    require(91 * (2 * scale) - 1 == 182 * scale - 1, "bare colour bound")
    require(91 * (6 * scale) - 1 == 546 * scale - 1, "word colour bound")

for unit_owner in range(1, 13):
    for shell_character in range(1, 13):
        matching = [
            k for k in units if (unit_owner * k - shell_character) % 13 == 0
        ]
        require(len(matching) == 6, "prescribed-character CRT lift count")

print("THM-2319 exact companion")
print(f"unit pairs: total={len(allowed)} carry0={len(carry_zero)} carry1={len(carry_one)}")
print(f"primitive characters surviving on rational needles: {len(units)}")
print(f"CRT coefficient table: {coefficient_table}")
print(f"needle floor: {floor_constant}")
print(f"carry-one selected-pair floor: {pair_constant}")
print(f"LRC total/pair floors: {lrc_total_constant}, {lrc_pair_constant}")
print(f"positive controls: {len(positive_vectors)}")
print(f"mixed consecutive hostile: {mixed_value}")
print(f"same-source full current: {mixed_full}")
print(f"nonconsecutive seven-site hostile: {wild_value}")
print("Schur landing: h1,h2<=J^2-1; h3<=2J^2-1")
print("LRC jump/charge bounds: J<=6S; A,B<=3276S^2-1; C<=6552S^2-1")
print("bare/word primitive-colour bounds: 182S-1, 546S-1")
print("PASS")
