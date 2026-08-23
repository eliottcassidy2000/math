#!/usr/bin/env python3
"""Independent hostile checker for fixed-order inert-prime amplification.

The checker does not import the candidate.  It separately tests the
ordered-collision union bound, empty/small banks, the THM-3793 primitive
shell typing, complete small singleton fibres, the row mass/cutoff, and the
liminf constants.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import comb, factorial, gcd
from pathlib import Path


GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    if condition is not True:
        raise RuntimeError(label)
    GATES += 1


def factor_integer(number: int) -> dict[int, int]:
    factors: dict[int, int] = {}
    trial = 2
    remaining = number
    while trial * trial <= remaining:
        while remaining % trial == 0:
            factors[trial] = factors.get(trial, 0) + 1
            remaining //= trial
        trial += 1
    if remaining > 1:
        factors[remaining] = factors.get(remaining, 0) + 1
    return factors


def exact_cube_root(number: int) -> int | None:
    low, high = 0, 1
    while high**3 < number:
        high *= 2
    while low + 1 < high:
        middle = (low + high) // 2
        if middle**3 < number:
            low = middle
        else:
            high = middle
    if high**3 == number:
        return high
    if low**3 == number:
        return low
    return None


def positive_distinct_representations(value: int) -> tuple[tuple[int, int], ...]:
    representations: list[tuple[int, int]] = []
    x = 1
    while 2 * x**3 < value:
        y = exact_cube_root(value - x**3)
        if y is not None and x < y:
            representations.append((x, y))
        x += 1
    return tuple(representations)


def elementary(weights: tuple[Fraction, ...], degree: int) -> Fraction:
    return sum(
        (prod(values) for values in combinations(weights, degree)),
        Fraction(0),
    )


def prod(values) -> Fraction | int:
    result: Fraction | int = 1
    for value in values:
        result *= value
    return result


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "04-computation" / "two_cube_all_fixed_prime_product_amplification_thm3793.py"
CANDIDATE_SHA256 = "3b6b8f4118ad4b383328cab7d5adf551fd43dda26e7d65adb4b6870acd47bd2a"
candidate_bytes = CANDIDATE.read_bytes()
gate(sha256(candidate_bytes).hexdigest() == CANDIDATE_SHA256, "candidate source hash")
gate(not any(isinstance(node, ast.Assert)
             for node in ast.walk(ast.parse(candidate_bytes.decode("utf-8")))),
     "candidate is assertion-free")


# Ordered versus unordered elementary-symmetric mass, including empty banks
# and every j larger than the bank.  The collision union tax is exact for
# each chosen pair of positions but may double-count higher collisions, so
# its direction is collision_mass <= C(j,2) B A^(j-2).
weight_banks = tuple(
    tuple(Fraction(1, p) for p in (5, 11, 17, 23)[:size])
    for size in range(5)
)
collision_cases = 0
for weights in weight_banks:
    count = len(weights)
    A = sum(weights, Fraction(0))
    B = sum((weight**2 for weight in weights), Fraction(0))
    for j in range(1, 7):
        e_j = elementary(weights, j)
        ordered = sum(
            (prod(weights[index] for index in indices)
             for indices in permutations(range(count), j)),
            Fraction(0),
        )
        gate(ordered == factorial(j) * e_j,
             f"ordered/unordered factor bank={count},j={j}")
        gate((e_j > 0) == (count >= j),
             f"empty/small-bank positivity bank={count},j={j}")
        if j == 1:
            gate(ordered == A, f"j=1 separate identity bank={count}")
        else:
            all_tuple_mass = A**j
            collision_mass = all_tuple_mass - ordered
            union_tax = comb(j, 2) * B * A ** (j - 2)
            gate(collision_mass >= 0,
                 f"collision mass nonnegative bank={count},j={j}")
            gate(collision_mass <= union_tax,
                 f"ordered collision union bound bank={count},j={j}")
            gate(ordered >= all_tuple_mass - union_tax,
                 f"elementary lower bound bank={count},j={j}")
            collision_cases += 1
        if j == 2:
            gate(ordered == A**2 - B,
                 f"j=2 exact collision identity bank={count}")


# Complete finite rows for products of distinct inert primes.  The complete
# coordinate-fibre enumeration is independent of THM-3793 and verifies both
# within-row and cross-row singleton/disjointness on this hostile bank.
prime_bank = (5, 11, 17)
seen_values: dict[int, tuple[int, int, int]] = {}
row_count = 0
value_count = 0
primitive_rows = 0
nonprimitive_rows = 0
for j in range(1, len(prime_bank) + 1):
    subsets = tuple(combinations(prime_bank, j))
    gate(len(subsets) == comb(len(prime_bank), j),
         f"unordered product count j={j}")
    for subset in subsets:
        d = int(prod(subset))
        factors_d = factor_integer(d)
        gate(factors_d == {p: 1 for p in subset}, f"squarefree product d={d}")
        gate(all(p % 3 == 2 for p in factors_d), f"inert factors d={d}")
        Z = max(subset)
        gate(d <= Z**j, f"product cutoff d={d},j={j}")
        row_values: set[int] = set()
        for x in range(1, (d + 1) // 2):
            y = d - x
            common = gcd(x, y)
            gate(common == gcd(x, d), f"gcd shell d={d},x={x}")
            gate(d % common == 0, f"common scale divides row sum d={d},x={x}")
            primitive_sum = d // common
            primitive_factors = factor_integer(primitive_sum)
            gate(all(p in factors_d and exponent <= 1
                     for p, exponent in primitive_factors.items()),
                 f"primitive quotient squarefree d={d},x={x}")
            gate(gcd(x // common, y // common) == 1
                 and x // common + y // common == primitive_sum,
                 f"primitive shell reconstruction d={d},x={x}")
            if common == 1:
                primitive_rows += 1
            else:
                nonprimitive_rows += 1

            value = x**3 + y**3
            gate(value < d**3 <= Z ** (3 * j),
                 f"strict exponent-three height d={d},x={x}")
            fibre = positive_distinct_representations(value)
            gate(fibre == ((x, y),), f"complete singleton fibre d={d},x={x}")
            gate(value not in row_values, f"within-row disjoint d={d},x={x}")
            gate(value not in seen_values, f"cross-row disjoint d={d},x={x}")
            row_values.add(value)
            seen_values[value] = (d, x, y)
        number = (d - 1) // 2
        gate(len(row_values) == number, f"odd row cardinality d={d}")
        mass_floor = Fraction(number, d**2)
        gate(mass_floor == Fraction(d - 1, 2 * d**2), f"row mass floor d={d}")
        gate(mass_floor >= Fraction(2, 5 * d), f"uniform row mass d={d}")
        row_count += 1
        value_count += number

gate(primitive_rows > 0 and nonprimitive_rows > 0,
     "primitive and nonprimitive shell controls both occur")


# The exponent-two primitive boundary is safe; exponent three is genuinely
# outside THM-3793 and admits the standard collision.
for x in range(1, 13):
    y = 25 - x
    common = gcd(x, y)
    primitive_sum = 25 // common
    gate(max(factor_integer(primitive_sum).values(), default=0) <= 2,
         f"inert exponent-two boundary x={x}")
    gate(positive_distinct_representations(x**3 + y**3) == ((x, y),),
         f"exponent-two singleton x={x}")

hostile = 515375
hostile_fibre = set(positive_distinct_representations(hostile))
gate(hostile_fibre == {(15, 80), (54, 71)}, "complete exponent-three hostile fibre")
gate(gcd(54, 71) == 1 and 54 + 71 == 5**3,
     "hostile primitive sum has inert exponent three")
gate(15 + 80 == 5 * 19 and (5 % 3, 19 % 3) == (2, 1),
     "competing row leaves the all-inert shell")


# Exact constants.  THM-3730 supplies A(Z)>=1/2 log log Z-C_P and B(Z)<1/4.
# The row factor is 2/5 and the unordered correction is 1/j!.
constants: dict[int, Fraction] = {}
for j in range(1, 13):
    derived = Fraction(2, 5) * Fraction(1, factorial(j)) * Fraction(1, 2) ** j
    expected = Fraction(2, 5 * 2**j * factorial(j))
    gate(derived == expected, f"liminf constant j={j}")
    constants[j] = derived
gate(constants[1] == Fraction(1, 5), "j=1 constant")
gate(constants[2] == Fraction(1, 20), "j=2 recovers THM-3793")
gate(Fraction(6, 25) < Fraction(1, 4), "elementary global B<1/4 envelope")

# A fixed real power R is beaten by first choosing one fixed integer j>R;
# the choice is made before the X-limit.  Rational controls include negative,
# integral, and nonintegral exponents.
real_power_controls = (
    Fraction(-7, 3), Fraction(0), Fraction(1, 2), Fraction(1),
    Fraction(7, 3), Fraction(10), Fraction(101, 7),
)
for power in real_power_controls:
    j = max(1, power.numerator // power.denominator + 1)
    gate(Fraction(j) > power, f"fixed-j quantifier for real power {power}")
    gate(constants.get(j, Fraction(2, 5 * 2**j * factorial(j))) > 0,
         f"positive fixed-j constant for real power {power}")

semantic = {
    "candidate_sha256": CANDIDATE_SHA256,
    "row": "d is an unordered product of j distinct inert primes; g=gcd(x,d); s=d/g is squarefree",
    "singleton": "THM-3793 applies to every row entry; singleton fibres force within/cross-row disjointness",
    "mass": "row > (d-1)/(2d^2) >= 2/(5d); m<d^3<=Z^(3j)",
    "collision": "j!e_j >= A^j-C(j,2)BA^(j-2), separately j=1: e_1=A",
    "passage": "Z=X^(1/(3j)); loglogZ=loglogX-log(3j)",
    "liminf": "2/(5*2^j*j!) for each fixed integer j>=1",
    "growth": "for every fixed real R, H(X)/(loglogX)^R tends to infinity",
    "scope": "explicit critical-mass lower bound only; no support asymptotic, residue, or fibre converse",
    "finite": {
        "collision_cases": collision_cases,
        "rows": row_count,
        "values": value_count,
        "primitive_entries": primitive_rows,
        "nonprimitive_entries": nonprimitive_rows,
    },
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert)
             for node in ast.walk(ast.parse(source))),
     "independent checker is assertion-free")

print("audit=two_cube_all_fixed_prime_product_amplification")
print(f"candidate_sha256={CANDIDATE_SHA256}")
print(f"finite_rows={row_count};finite_values={value_count};primitive={primitive_rows};nonprimitive={nonprimitive_rows}")
print(f"collision_cases={collision_cases};banks=empty_through_four;orders=1..6")
print("boundary=j1_separate;j2_exact;exponent2_safe;exponent3_hostile")
print("liminf_constant=2/(5*2^j*j!);fixed_j_before_X_limit")
print("growth=for_each_real_R:H(X)/(loglogX)^R_to_infinity")
print("scope=no_support_asymptotic_or_residue")
print(f"semantic_sha256={sha256(semantic_blob).hexdigest()}")
print(f"GATES={GATES}")
print("RESULT=PASS;J1_EMPTY_BANK_FORMULA_REPAIRED")
