#!/usr/bin/env python3
"""Exact Berggren descendant-angle language and harmonic-density audit.

THM-3334 gives the two irrational slope walls.  THM-2596 gives the three
Möbius branch maps.  This companion proves that both walls have period four
under the inverse branch map, derives the exact level-count recurrences, and
checks the resulting 16/41, 9/41, 16/41 density split.

All arithmetic is exact.  Runtime gates survive ``python -O``.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
D = 145
EXPECTED_SEMANTIC_DIGEST = "b0a83446f06d4836ebb9125a3e88fbba8fe41fe44ed79017f85d3c8cfd229129"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


class ExactDigest:
    def __init__(self):
        self._hash = sha256()

    def add(self, value):
        self._hash.update(repr(value).encode("ascii"))
        self._hash.update(bytes((10,)))

    def hexdigest(self):
        return self._hash.hexdigest()


@dataclass(frozen=True)
class Quad:
    """Element a+b*sqrt(145), with a,b rational."""

    a: Fraction
    b: Fraction = Fraction(0)

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Quad) else Quad(Fraction(value))

    def __add__(self, other):
        other = Quad.coerce(other)
        return Quad(self.a + other.a, self.b + other.b)

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return Quad(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-Quad.coerce(other))

    def __rsub__(self, other):
        return Quad.coerce(other) - self

    def __mul__(self, other):
        other = Quad.coerce(other)
        return Quad(
            self.a * other.a + D * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    def __rmul__(self, other):
        return self * other

    def inverse(self):
        denominator = self.a * self.a - D * self.b * self.b
        require(denominator != 0, ("zero quadratic denominator", self))
        return Quad(self.a / denominator, -self.b / denominator)

    def __truediv__(self, other):
        return self * Quad.coerce(other).inverse()

    def __rtruediv__(self, other):
        return Quad.coerce(other) / self

    def sign(self):
        if self.b == 0:
            return (self.a > 0) - (self.a < 0)
        if self.a == 0:
            return (self.b > 0) - (self.b < 0)
        sign_a = (self.a > 0) - (self.a < 0)
        sign_b = (self.b > 0) - (self.b < 0)
        if sign_a == sign_b:
            return sign_a
        comparison = self.a * self.a - D * self.b * self.b
        if self.a > 0:
            return (comparison > 0) - (comparison < 0)
        return -((comparison > 0) - (comparison < 0))

    def __lt__(self, other):
        return (self - other).sign() < 0

    def __repr__(self):
        return f"Quad({self.a!r},{self.b!r})"


ZERO = Quad(Fraction(0))
ONE = Quad(Fraction(1))
SQRT_D = Quad(Fraction(0), Fraction(1))
ALPHA = (SQRT_D - 9) / 8
BETA = (SQRT_D - 8) / 9


def inverse_branch(value):
    if value < Fraction(1, 3):
        return "C", value / (1 - 2 * value)
    if value < Fraction(1, 2):
        return "B", 1 / value - 2
    return "A", 2 - 1 / value


def wall_orbit(value):
    orbit = []
    current = value
    for _ in range(4):
        letter, following = inverse_branch(current)
        orbit.append((letter, current))
        current = following
    require(current == value, ("wall does not close", value, orbit, current))
    return tuple(orbit)


AFFINE_CDF = {
    "A": (Fraction(2, 3), Fraction(1, 3)),
    "B": (Fraction(2, 3), Fraction(-1, 3)),
    "C": (Fraction(0), Fraction(1, 3)),
}


def compose_affine(outer, inner):
    c_out, s_out = outer
    c_in, s_in = inner
    return c_out + s_out * c_in, s_out * s_in


def orbit_cdf_return(orbit):
    result = (Fraction(0), Fraction(1))
    for letter, _ in reversed(orbit):
        result = compose_affine(AFFINE_CDF[letter], result)
    return result


def child(letter, pair):
    m, n = pair
    if letter == "A":
        return n, 2 * n - m
    if letter == "B":
        return n, 2 * n + m
    if letter == "C":
        return m, n + 2 * m
    raise ValueError(letter)


def rational_less_wall(pair, wall):
    m, n = pair
    return Quad(Fraction(m, n)) < wall


def direct_levels(max_depth):
    level = ((1, 2),)
    records = []
    all_seen = set()
    for depth in range(max_depth + 1):
        require(len(level) == 3**depth, ("level size", depth, len(level)))
        require(len(set(level)) == len(level), ("duplicate in level", depth))
        require(not (set(level) & all_seen), ("cross-level duplicate", depth))
        all_seen.update(level)
        for m, n in level:
            require(0 < m < n, ("positive slope", depth, m, n))
            require(gcd(m, n) == 1, ("primitive slope", depth, m, n))
            require((m - n) % 2 != 0, ("opposite parity", depth, m, n))

        below_alpha = sum(rational_less_wall(pair, ALPHA) for pair in level)
        below_beta = sum(rational_less_wall(pair, BETA) for pair in level)
        records.append(
            (
                depth,
                below_alpha,
                below_beta - below_alpha,
                len(level) - below_beta,
            )
        )
        level = tuple(child(letter, pair) for pair in level for letter in "ABC")
    return tuple(records), len(all_seen)


def recurrence_counts(max_depth):
    u = [0, 1, 4, 10]
    v = [0, 2, 6, 16]
    for depth in range(max_depth - 3):
        u.append(32 * 3**depth - u[depth])
        v.append(50 * 3**depth - v[depth])
    records = []
    for depth in range(max_depth + 1):
        records.append((depth, u[depth], v[depth] - u[depth], 3**depth - v[depth]))
    return tuple(records)


def suffix_counts(start, depth):
    level = (start,)
    for _ in range(depth):
        level = tuple(child(letter, pair) for pair in level for letter in "ABC")
    u = sum(rational_less_wall(pair, ALPHA) for pair in level)
    v = sum(rational_less_wall(pair, BETA) for pair in level)
    return u, v - u, len(level) - v


def main():
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "floating literal",
    )

    require(4 * ALPHA * ALPHA + 9 * ALPHA - 4 == ZERO, "alpha wall")
    require(9 * BETA * BETA + 16 * BETA - 9 == ZERO, "beta wall")
    require(ZERO < ALPHA and ALPHA < BETA and BETA < ONE, "wall order")

    alpha_orbit = wall_orbit(ALPHA)
    beta_orbit = wall_orbit(BETA)
    require(tuple(letter for letter, _ in alpha_orbit) == tuple("BABB"), alpha_orbit)
    require(tuple(letter for letter, _ in beta_orbit) == tuple("BCBB"), beta_orbit)
    alpha_return = orbit_cdf_return(alpha_orbit)
    beta_return = orbit_cdf_return(beta_orbit)
    require(alpha_return == (Fraction(32, 81), Fraction(-1, 81)), alpha_return)
    require(beta_return == (Fraction(50, 81), Fraction(-1, 81)), beta_return)

    alpha_limit = alpha_return[0] / (1 - alpha_return[1])
    beta_limit = beta_return[0] / (1 - beta_return[1])
    densities = (alpha_limit, beta_limit - alpha_limit, 1 - beta_limit)
    require(densities == (Fraction(16, 41), Fraction(9, 41), Fraction(16, 41)), densities)

    direct, distinct_nodes = direct_levels(10)
    recurrence = recurrence_counts(40)
    require(direct == recurrence[:11], (direct, recurrence[:11]))

    expected_residuals = (
        (-16, -9, 25),
        (-7, 14, -7),
        (20, 1, -21),
        (-22, 3, 19),
        (16, 9, -25),
        (7, -14, 7),
        (-20, -1, 21),
        (22, -3, -19),
    )
    residuals = []
    for depth, u_count, acute_count, d_count in recurrence:
        residual = (
            41 * u_count - 16 * 3**depth,
            41 * acute_count - 9 * 3**depth,
            41 * d_count - 16 * 3**depth,
        )
        residuals.append(residual)
        require(residual == expected_residuals[depth % 8], (depth, residual))

    # The four-step contraction is independent of the starting rational
    # state.  Exhaust a nontrivial bank of prefix states and suffix depths.
    prefix_level = ((1, 2),)
    prefix_states = []
    for _ in range(5):
        prefix_states.extend(prefix_level)
        prefix_level = tuple(
            child(letter, pair) for pair in prefix_level for letter in "ABC"
        )
    suffix_checks = 0
    for start in prefix_states:
        values = [suffix_counts(start, depth) for depth in range(8)]
        for depth in range(4):
            u0, a0, d0 = values[depth]
            u4, a4, d4 = values[depth + 4]
            require(u4 == 32 * 3**depth - u0, ("prefix U", start, depth))
            require(a4 == 18 * 3**depth - a0, ("prefix acute", start, depth))
            require(d4 == 32 * 3**depth - d0, ("prefix D", start, depth))
            suffix_checks += 1

    semantic = ExactDigest()
    semantic.add(("walls", ALPHA, BETA))
    semantic.add(("orbits", alpha_orbit, beta_orbit))
    semantic.add(("returns", alpha_return, beta_return))
    semantic.add(("densities", densities))
    semantic.add(("direct", direct, distinct_nodes))
    semantic.add(("recurrence", recurrence))
    semantic.add(("residuals", tuple(residuals)))
    semantic.add(("suffix_checks", len(prefix_states), suffix_checks))
    digest = semantic.hexdigest()
    if EXPECTED_SEMANTIC_DIGEST != "TO_BE_FILLED":
        require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    print("BERGGREN DESCENDANT-ANGLE LANGUAGE / DENSITY EXACT AUDIT")
    print(f"source_sha256_lf={lf_hash(source)}")
    print("status=PROOF_COMPLETE_ELEMENTARY_CANDIDATE+VERIFIED_EXACT;not_yet_canon")
    print("walls=alpha=(sqrt(145)-9)/8;beta=(sqrt(145)-8)/9")
    print(
        "inverse_wall_itineraries="
        f"alpha:{''.join(letter for letter, _ in alpha_orbit)};"
        f"beta:{''.join(letter for letter, _ in beta_orbit)}"
    )
    print(f"four_step_cdf_returns=alpha:{alpha_return};beta:{beta_return}")
    print("count_recursions=U[n+4]=32*3^n-U[n];A[n+4]=18*3^n-A[n];D[n+4]=32*3^n-D[n]")
    print(f"level_densities_U_acute_D={densities}")
    print(f"residual_period8={expected_residuals}")
    print(f"direct_tree_depth=10;distinct_nodes={distinct_nodes};counts={direct}")
    print(f"uniform_prefix_bank={len(prefix_states)};four_step_suffix_checks={suffix_checks}")
    print("language=regular_by_periodic_wall_itinerary;shortlex_discrepancy=O(log_N)_proof_side")
    print("harmonic=each_shortlex_type_has_logarithmic_coefficient_equal_to_its_density")
    print("scope=no_RXTX_tensor_map_no_tournament_current_no_LRC_or_JC_consequence")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
