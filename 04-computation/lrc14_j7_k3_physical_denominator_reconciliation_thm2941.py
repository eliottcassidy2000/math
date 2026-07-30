#!/usr/bin/env python3
"""Physical ``(E,z1)`` to conditional denominator-GF bridge for LRC14 k=3.

This exact-computation instrument intersects two different
necessary relaxations without identifying either one with a realized cover.

1.  The projected suffix ledger of THM-2941 supplies a *physical* row
    ``(E,z1)``.  Its distinguished first denominator is

        d1 = L / gcd(L,z1),       L = 14 lcm(E).

2.  If the other three drift denominators form a multiset ``M``, the resolving
    denominator is ``D=lcm(d1,M)``.  For fixed ``(D,d1)`` we count the later
    three-multisets by divisor-lattice Mobius inversion.  The retained feature
    is

        (multiplicities of d=2,3,4; large-mask capacity; uniform count c),

    where ``c`` counts all four denominators not dividing ``q=D/7``.

The row first has to admit a resolving ``D`` passing THM-2928's support
transfer.  It is then tested by two upper relaxations:

* expected spike: with ``N_c=#{q-fibres of S_D having load>c}``, exact literal
  phase averaging and compact/open strictness require

      N_c=0  or  55*N_c < 13*(4-c)*q;

* support status: d=2,3,4 receive their exact largest ``S_D`` residue-class
  load when active, their activity marginals are d/7, and d>4 receive full
  ambient ceiling capacity.  The upward-event fractional-cover theorem gives
  the exact real one-threshold status allowance.

The conjunction is still only necessary.  It forgets later numerators, order,
distinctness, common-table compatibility across thresholds, literal phases,
and the projected high-wall constraint after the inherited physical row.

The physical ledger is reconstructed independently from the canonical carrier
primitive using an exact vectorized integer formula, then checked against the
canonical count 376,020 and selected canonical per-body profiles.  No finite
label horizon is treated as exhaustive: the inherited THM-2941 omitted-label
tail is retained exactly.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import multiprocessing as mp
import os
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import cmp_to_key, lru_cache
from itertools import combinations, combinations_with_replacement
from math import comb, gcd, lcm
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PHYSICAL_PATH = (
    ROOT / "04-computation" / "lrc14_j7_aligned_projected_arc_suffix_thm2941.py"
)
COMBINED_PATH = (
    ROOT / "04-computation" / "lrc14_three_drift_body_projection_fiber_thm2928.py"
)
Z380_CLOSURE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_frontier_fibre_closure_thm2941.py"
)
Z380_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" / "lrc14_j7_k3_frontier_fibre_closure_thm2941.out"
)
Z250_CLOSURE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_uniform_ray_status_closure_thm2941.py"
)
Z250_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" / "lrc14_j7_k3_uniform_ray_status_closure_thm2941.out"
)
Z378_CLOSURE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
)
Z378_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" / "lrc14_j7_k3_z378_ray_status_closure_thm2941.out"
)
Z364_CLOSURE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z364_ray_status_closure_thm2941.py"
)
Z364_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" / "lrc14_j7_k3_z364_ray_status_closure_thm2941.out"
)
Z350_CLOSURE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z350_ray_status_closure_thm2941.py"
)
Z350_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" / "lrc14_j7_k3_z350_ray_status_closure_thm2941.out"
)
Z336_CLOSURE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z336_ray_status_closure_thm2941.py"
)
Z336_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" / "lrc14_j7_k3_z336_ray_status_closure_thm2941.out"
)
Z330_SCALAR_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z330_scalar_slice_thm2941.py"
)
Z330_SCALAR_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" / "lrc14_j7_k3_z330_scalar_slice_thm2941.out"
)
Z330_CLOSURE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z330_ray_status_closure_thm2941.py"
)
Z330_OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k3_z330_ray_status_closure_thm2941.out"
)
Z328_SCALAR_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z328_scalar_slice_thm2941.py"
)
Z328_SCALAR_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" / "lrc14_j7_k3_z328_scalar_slice_thm2941.out"
)
Z328_CLOSURE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z328_ray_status_closure_thm2941.py"
)
Z328_OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k3_z328_ray_status_closure_thm2941.out"
)
EXPECTED_PHYSICAL_SHA256 = (
    "a003d287f618eb301edf6974d0b67dc128c4f380a169e7809ed5b5754e8b8303"
)
EXPECTED_COMBINED_SHA256 = (
    "42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4"
)
EXPECTED_Z380_CLOSURE_SHA256 = (
    "64f98439f677668c82045e7f9107cbfdff467afd8f16975c7e37d8ae5c5c9f26"
)
EXPECTED_Z380_OUTPUT_SHA256 = (
    "a1c77b24488240f1ee0295e427ee4583b7d8215caf6615f424bf325350fb56b6"
)
EXPECTED_Z250_CLOSURE_SHA256 = (
    "dfa4788297b8c31fc9b5dce1afadf29d20b267cb4159fa95dadb9346b1980b36"
)
EXPECTED_Z250_OUTPUT_SHA256 = (
    "5abccb7ef700cec83b9989e8abcd83bc24f51c0a35f7f9054522da0dd62109fe"
)
EXPECTED_Z378_CLOSURE_SHA256 = (
    "6ff8676255d51d818d7c24102a8fc755e673544f0ac6b99be4bfc262c892df1e"
)
EXPECTED_Z378_OUTPUT_SHA256 = (
    "a48a7ca3c3142476deac9b81a48c7f6937bc3312220579520d2440f5c76133a9"
)
EXPECTED_Z364_CLOSURE_SHA256 = (
    "b76f82921920970f2f9423d83a145d928d4ab7ca97af66484ec5ef0935832c16"
)
EXPECTED_Z364_OUTPUT_SHA256 = (
    "8e052b92f1896e3f2a0b69e8bb9bb519700c445001958572871fa05a82f28f65"
)
EXPECTED_Z350_CLOSURE_SHA256 = (
    "e3e99d6fa63f966518b7b07f7c2a984faee5d6c0b0f69874e9688c5b449d2c9c"
)
EXPECTED_Z350_OUTPUT_SHA256 = (
    "ff6ae15ed608e588d9200a91622e43e50601c0ab2e62f55d8d86ab6bae6589d6"
)
EXPECTED_Z336_CLOSURE_SHA256 = (
    "6e991432a1ce5ec9c1cbe97199cb5fd647c5edc89bdcc60ce560382a70835fcf"
)
EXPECTED_Z336_OUTPUT_SHA256 = (
    "e257da724128208a9c80fc5f3e8f0cd4151b2073d3a8afa4814ca1e274f168ac"
)
EXPECTED_Z330_SCALAR_SHA256 = (
    "5eb30248f2beac09b29ebccc0dc8b203361f97c4b8552021130290566f8de0d1"
)
EXPECTED_Z330_SCALAR_OUTPUT_SHA256 = (
    "c2bf51b8b661d98d1ce2a8114c6187a6eefb17856e9684a46704a62ccecc81ed"
)
EXPECTED_Z330_CLOSURE_SHA256 = (
    "255aa2522d1ed6bc493ed397eff490434da9cba736ebc0e14bb5514857a22657"
)
EXPECTED_Z330_OUTPUT_SHA256 = (
    "3c62328f51c043a953ea20b788d612313522108c3489614300b0a286c1cfd97b"
)
EXPECTED_Z328_SCALAR_SHA256 = (
    "75fd846d9070b267003c70e175a9af1225815c3567a5a2558891f5af7a61f8f3"
)
EXPECTED_Z328_SCALAR_OUTPUT_SHA256 = (
    "f5c11f364a626141af181d84f39d48030ef91a8ddf7d74b9602bb15cd7eb626e"
)
EXPECTED_Z328_CLOSURE_SHA256 = (
    "1869024359822ac811d79721ba5e3a53fe58e516ec2e5bec2ebd932bd56c3d7a"
)
EXPECTED_Z328_OUTPUT_SHA256 = (
    "40e70a0edf567c302955298b4d4b434544537084e8eaa19524103fe4d82bb1f2"
)

EXPECTED_PHYSICAL_ROWS = 376_020
EXPECTED_SUPPORT_ROWS = 26_970
EXPECTED_SUPPORT_DIVISORS = 217
SUPPORT_CUTOFF = Q(125, 143)
THREE_SAFE_FLOOR = Q(55, 91)
ETA3 = Q(3, 91)
PROJECTED_RATIO3 = Q(13, 132)
K = 3
DRIFTS = 4
LATER_ARITY = 3


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


INPUT_SPECS = (
    ("physical_source", PHYSICAL_PATH, EXPECTED_PHYSICAL_SHA256),
    ("projection_source", COMBINED_PATH, EXPECTED_COMBINED_SHA256),
    ("z380_closure_source", Z380_CLOSURE_PATH, EXPECTED_Z380_CLOSURE_SHA256),
    ("z380_closure_output", Z380_OUTPUT_PATH, EXPECTED_Z380_OUTPUT_SHA256),
    ("z250_closure_source", Z250_CLOSURE_PATH, EXPECTED_Z250_CLOSURE_SHA256),
    ("z250_closure_output", Z250_OUTPUT_PATH, EXPECTED_Z250_OUTPUT_SHA256),
    ("z378_closure_source", Z378_CLOSURE_PATH, EXPECTED_Z378_CLOSURE_SHA256),
    ("z378_closure_output", Z378_OUTPUT_PATH, EXPECTED_Z378_OUTPUT_SHA256),
    ("z364_closure_source", Z364_CLOSURE_PATH, EXPECTED_Z364_CLOSURE_SHA256),
    ("z364_closure_output", Z364_OUTPUT_PATH, EXPECTED_Z364_OUTPUT_SHA256),
    ("z350_closure_source", Z350_CLOSURE_PATH, EXPECTED_Z350_CLOSURE_SHA256),
    ("z350_closure_output", Z350_OUTPUT_PATH, EXPECTED_Z350_OUTPUT_SHA256),
    ("z336_closure_source", Z336_CLOSURE_PATH, EXPECTED_Z336_CLOSURE_SHA256),
    ("z336_closure_output", Z336_OUTPUT_PATH, EXPECTED_Z336_OUTPUT_SHA256),
    ("z330_scalar_source", Z330_SCALAR_PATH, EXPECTED_Z330_SCALAR_SHA256),
    (
        "z330_scalar_output",
        Z330_SCALAR_OUTPUT_PATH,
        EXPECTED_Z330_SCALAR_OUTPUT_SHA256,
    ),
    ("z330_closure_source", Z330_CLOSURE_PATH, EXPECTED_Z330_CLOSURE_SHA256),
    ("z330_closure_output", Z330_OUTPUT_PATH, EXPECTED_Z330_OUTPUT_SHA256),
    ("z328_scalar_source", Z328_SCALAR_PATH, EXPECTED_Z328_SCALAR_SHA256),
    (
        "z328_scalar_output",
        Z328_SCALAR_OUTPUT_PATH,
        EXPECTED_Z328_SCALAR_OUTPUT_SHA256,
    ),
    ("z328_closure_source", Z328_CLOSURE_PATH, EXPECTED_Z328_CLOSURE_SHA256),
    ("z328_closure_output", Z328_OUTPUT_PATH, EXPECTED_Z328_OUTPUT_SHA256),
)
INPUT_SHA256 = {}
for _label, _path, _expected in INPUT_SPECS:
    _observed = file_sha256(_path)
    require(_observed == _expected, ("frozen input changed", _label, _observed))
    INPUT_SHA256[_label] = _observed


def require_inputs_unchanged(stage: str) -> None:
    """Reject a concurrent checkout mutation instead of printing a mixed run."""
    for label, path, _expected in INPUT_SPECS:
        observed = file_sha256(path)
        require(observed == INPUT_SHA256[label], ("input drift", stage, label, observed))


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot load", path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


physical = load_module("lrc14_k3_physical_source", PHYSICAL_PATH)
combined = load_module("lrc14_k3_projection_source", COMBINED_PATH)
support = combined.support_module
require_inputs_unchanged("after_module_load")


def mobius(number: int) -> int:
    result = 1
    remaining = number
    prime = 2
    while prime * prime <= remaining:
        if remaining % prime:
            prime += 1
            continue
        remaining //= prime
        if remaining % prime == 0:
            return 0
        result = -result
        while remaining % prime == 0:
            remaining //= prime
        prime += 1
    if remaining > 1:
        result = -result
    return result


@lru_cache(maxsize=None)
def divisors_of(number: int) -> tuple[int, ...]:
    return tuple(support.divisors(number))


def multichoose(types: int, copies: int) -> int:
    if copies == 0:
        return 1
    if types == 0:
        return 0
    return comb(types + copies - 1, copies)


def denominator_feature(D: int, d: int) -> tuple[int, int, int, int, int]:
    """Return the additive proof feature of one denominator symbol."""
    require(D % d == 0 and D % 7 == 0, ("bad feature denominator", D, d))
    q = D // 7
    return (
        int(d == 2),
        int(d == 3),
        int(d == 4),
        (D // d) * ((d + 6) // 7) if d > 4 else 0,
        int(q % d != 0),
    )


@lru_cache(maxsize=None)
def conditional_unweighted_count(D: int, d1: int) -> int:
    """Closed Mobius formula for fixed d1 plus three later symbols."""
    require(D % d1 == 0, ("d1 does not divide D", D, d1))
    total = 0
    for E in divisors_of(D):
        if E % d1:
            continue
        total += mobius(D // E) * multichoose(len(divisors_of(E)) - 1, 3)
    require(total > 0, ("empty conditional exact-lcm alphabet", D, d1))
    return total


@lru_cache(maxsize=None)
def conditional_c_distribution(D: int, d1: int) -> Counter:
    """Closed Mobius formula by total uniform-one count c."""
    require(D % 7 == 0 and D % d1 == 0, ("bad conditional c request", D, d1))
    q = D // 7
    fixed_c = int(q % d1 != 0)
    result = Counter()
    for E in divisors_of(D):
        if E % d1:
            continue
        sign = mobius(D // E)
        if not sign:
            continue
        alphabet = tuple(d for d in divisors_of(E) if d > 1)
        uniform = sum(q % d != 0 for d in alphabet)
        spike = len(alphabet) - uniform
        for later_uniform in range(4):
            result[fixed_c + later_uniform] += sign * multichoose(
                uniform, later_uniform
            ) * multichoose(spike, 3 - later_uniform)
    require(
        all(value >= 0 for value in result.values()),
        ("negative conditional c coefficient", D, d1, result),
    )
    result += Counter()
    require(
        sum(result.values()) == conditional_unweighted_count(D, d1),
        ("conditional c formula lost shapes", D, d1),
    )
    return result


@lru_cache(maxsize=128)
def conditional_feature_distribution(D: int, d1: int) -> Counter:
    """Exact-lcm GF for fixed d1 and a later denominator 3-multiset."""
    require(D % 7 == 0 and D % d1 == 0, ("bad conditional GF request", D, d1))
    result = Counter()
    fixed = denominator_feature(D, d1)
    for E in divisors_of(D):
        if E % d1:
            continue
        sign = mobius(D // E)
        if not sign:
            continue
        groups = Counter(
            denominator_feature(D, d)
            for d in divisors_of(E)
            if d > 1
        )
        # State: later_used, m2,m3,m4,large_capacity,uniform_count.
        states = {(0, *fixed): 1}
        for unit, alphabet_size in groups.items():
            additions = Counter()
            for state, multiplicity in tuple(states.items()):
                used, m2, m3, m4, large, uniform = state
                for copies in range(1, 4 - used):
                    additions[
                        (
                            used + copies,
                            m2 + copies * unit[0],
                            m3 + copies * unit[1],
                            m4 + copies * unit[2],
                            large + copies * unit[3],
                            uniform + copies * unit[4],
                        )
                    ] += multiplicity * multichoose(alphabet_size, copies)
            for state, multiplicity in additions.items():
                states[state] = states.get(state, 0) + multiplicity
        for state, multiplicity in states.items():
            used, m2, m3, m4, large, uniform = state
            if used == 3:
                result[((m2, m3, m4), large, uniform)] += sign * multiplicity
    require(
        all(value >= 0 for value in result.values()),
        ("negative conditional feature coefficient", D, d1),
    )
    result += Counter()
    collapsed = Counter()
    for (_pattern, _large, c), multiplicity in result.items():
        collapsed[c] += multiplicity
    require(
        collapsed == conditional_c_distribution(D, d1),
        ("conditional feature GF disagrees with c formula", D, d1),
    )
    return result


def brute_conditional_controls(max_D: int = 140) -> tuple[int, int]:
    """Literal three-multiset controls, independent of Mobius inversion."""
    c_cases = 0
    feature_cases = 0
    for D in range(7, max_D + 1, 7):
        alphabet = tuple(d for d in divisors_of(D) if d > 1)
        for d1 in alphabet:
            brute = Counter()
            for later in combinations_with_replacement(alphabet, 3):
                if lcm(d1, *later) != D:
                    continue
                units = (d1, *later)
                pattern = (
                    units.count(2),
                    units.count(3),
                    units.count(4),
                )
                large = sum(
                    (D // d) * ((d + 6) // 7)
                    for d in units
                    if d > 4
                )
                c = sum((D // 7) % d != 0 for d in units)
                brute[(pattern, large, c)] += 1
            require(
                brute == conditional_feature_distribution(D, d1),
                ("literal conditional feature mismatch", D, d1),
            )
            brute_c = Counter()
            for (_pattern, _large, c), multiplicity in brute.items():
                brute_c[c] += multiplicity
            require(
                brute_c == conditional_c_distribution(D, d1),
                ("literal conditional c mismatch", D, d1),
            )
            c_cases += 1
            feature_cases += 1
    return c_cases, feature_cases


def solve_square(rows: tuple, rhs: tuple) -> tuple[Q, ...] | None:
    size = len(rows)
    matrix = [
        [Q(value) for value in row] + [Q(rhs[index])]
        for index, row in enumerate(rows)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if matrix[row][column]),
            None,
        )
        if pivot is None:
            return None
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        scale = matrix[column][column]
        matrix[column] = [value / scale for value in matrix[column]]
        for row in range(size):
            if row == column or not matrix[row][column]:
                continue
            scale = matrix[row][column]
            matrix[row] = [
                value - scale * pivot_value
                for value, pivot_value in zip(matrix[row], matrix[column])
            ]
    return tuple(matrix[index][-1] for index in range(size))


@lru_cache(maxsize=None)
def maximum_upward_mass(
    weights: tuple[int, ...], marginals: tuple[Q, ...], need: int
) -> Q:
    count = len(weights)
    require(len(marginals) == count, "status arity mismatch")
    if need <= 0:
        return Q(1)
    if need > sum(weights):
        return Q(0)
    minimal = []
    for mask in range(1 << count):
        weight = sum(
            weights[index] for index in range(count) if (mask >> index) & 1
        )
        if weight < need:
            continue
        if any(
            weight - weights[index] >= need
            for index in range(count)
            if (mask >> index) & 1
        ):
            continue
        minimal.append(tuple((mask >> index) & 1 for index in range(count)))
    require(minimal, ("missing minimal true status", weights, need))
    constraints = [(row, Q(1)) for row in minimal]
    constraints.extend(
        (
            tuple(int(index == coordinate) for index in range(count)),
            Q(0),
        )
        for coordinate in range(count)
    )
    optimum = None
    for chosen in combinations(range(len(constraints)), count):
        point = solve_square(
            tuple(constraints[index][0] for index in chosen),
            tuple(constraints[index][1] for index in chosen),
        )
        if point is None or any(value < 0 for value in point):
            continue
        if any(
            sum(value * coefficient for value, coefficient in zip(point, row)) < 1
            for row in minimal
        ):
            continue
        objective = sum(
            marginal * value for marginal, value in zip(marginals, point)
        )
        if optimum is None or objective < optimum:
            optimum = objective
    require(optimum is not None, "fractional-cover vertex set empty")
    return min(Q(1), optimum)


@lru_cache(maxsize=None)
def status_need_limit(weights: tuple[int, ...], marginals: tuple[Q, ...]) -> int:
    if not weights:
        return 0
    thresholds = sorted(
        {
            sum(
                weights[index]
                for index in range(len(weights))
                if (mask >> index) & 1
            )
            for mask in range(1, 1 << len(weights))
        }
    )
    low = 0
    for threshold in thresholds:
        if maximum_upward_mass(weights, marginals, threshold) > THREE_SAFE_FLOOR:
            low = threshold
        else:
            break
    require(
        maximum_upward_mass(weights, marginals, low) > THREE_SAFE_FLOOR,
        "status allowance is not admissible",
    )
    if low < sum(weights):
        require(
            maximum_upward_mass(weights, marginals, low + 1)
            <= THREE_SAFE_FLOOR,
            "status allowance is not maximal",
        )
    return low


def small_vectors(
    pattern: tuple[int, int, int], small_loads: dict[int, int]
) -> tuple[tuple[int, ...], tuple[Q, ...]]:
    weights = []
    marginals = []
    for d, copies in zip((2, 3, 4), pattern):
        if not copies:
            continue
        require(d in small_loads and small_loads[d] > 0, ("missing small load", d))
        weights.extend([small_loads[d]] * copies)
        marginals.extend([Q(d, 7)] * copies)
    return tuple(weights), tuple(marginals)


def _rational_compare(left: tuple[int, int, int], right: tuple[int, int, int]) -> int:
    """Sort descending by value, then ascending by label."""
    cross_left = left[0] * right[1]
    cross_right = right[0] * left[1]
    if cross_left != cross_right:
        return -1 if cross_left > cross_right else 1
    return -1 if left[2] < right[2] else (1 if left[2] > right[2] else 0)


RATIONAL_KEY = cmp_to_key(_rational_compare)


def top_insert_raw(
    rows: list[tuple[int, int, int]], item: tuple[int, int, int], limit: int = 3
) -> None:
    rows.append(item)
    index = len(rows) - 1
    while index and _rational_compare(rows[index], rows[index - 1]) < 0:
        rows[index], rows[index - 1] = rows[index - 1], rows[index]
        index -= 1
    del rows[limit:]


def rational_value(item: tuple[int, int, int]) -> Q:
    return Q(item[0], item[1])


def physical_rows_for_body(body: tuple[int, ...]) -> tuple:
    """Exact vectorized replay of THM-2941's k=3 projected suffix row."""
    carrier = physical.A.carrier_for(body)
    ruler = physical.A.RULER
    h_numerator = sum(right - left for left, right in carrier)
    components = len(carrier)
    canonical_L = 14 * lcm(*body)
    require(h_numerator > 0 and components > 0, (body, "empty carrier"))

    labels = np.arange(physical.A.BASE_LABEL, physical.HORIZON + 1, dtype=np.int64)
    coverage_numerators = np.zeros(labels.shape, dtype=np.int64)
    for left, right in carrier:
        for sign, endpoint in ((-1, left), (1, right)):
            argument = labels * endpoint
            cycles = argument // ruler
            remainder = argument % ruler
            primitive = (
                cycles * physical.A.ONE_SEVENTH
                + np.minimum(remainder, physical.A.ONE_FOURTEENTH)
                + np.maximum(0, remainder - physical.A.THIRTEEN_FOURTEENTHS)
            )
            coverage_numerators += sign * primitive
    delta_numerators = 7 * coverage_numerators - h_numerator * labels
    delta_denominators = 7 * ruler * labels

    h = Q(h_numerator, ruler)
    analytic_bound = Q(24 * components, 49) / (h * ETA3)
    analytic_cap = analytic_bound.numerator // analytic_bound.denominator
    wall = PROJECTED_RATIO3 * canonical_L
    high_floor = max(15, wall.numerator // wall.denominator + 1)
    lower = h * ETA3
    ordinary_tail = (6 * components, 49 * (physical.HORIZON + 1), physical.HORIZON + 1)
    high_tail = (
        6 * components,
        49 * max(physical.HORIZON + 1, high_floor),
        physical.HORIZON + 1,
    )

    arbitrary_top: list[tuple[int, int, int]] = []
    high_top: list[tuple[int, int, int]] = []
    survivors = []
    for offset in range(len(labels) - 1, -1, -1):
        first = int(labels[offset])
        item = (
            int(delta_numerators[offset]),
            int(delta_denominators[offset]),
            first,
        )
        if first % canonical_L == 0:
            continue
        if first <= analytic_cap:
            arbitrary = sorted(arbitrary_top + [ordinary_tail] * 3, key=RATIONAL_KEY)
            if first < high_floor:
                chosen_high = sorted(high_top + [high_tail], key=RATIONAL_KEY)[0]
                if chosen_high[2] == physical.HORIZON + 1:
                    rest = arbitrary[:2]
                else:
                    rest = [
                        candidate
                        for candidate in arbitrary
                        if candidate[2] != chosen_high[2]
                    ][:2]
                chosen = (chosen_high, *rest)
            else:
                chosen = tuple(arbitrary[:3])
            upper = rational_value(item) + sum(
                (rational_value(candidate) for candidate in chosen), Q(0)
            )
            if upper >= lower:
                survivors.append(first)
        top_insert_raw(arbitrary_top, item)
        if first >= high_floor:
            top_insert_raw(high_top, item)
    survivors.sort()
    return body, canonical_L, tuple(survivors), analytic_cap, high_floor


def reconstruct_physical_rows(workers: int) -> tuple[list[tuple], list[tuple]]:
    bodies = tuple(combinations(range(1, 15), 6))
    if workers == 1:
        profiles = [physical_rows_for_body(body) for body in bodies]
    else:
        with mp.get_context("spawn").Pool(workers) as pool:
            profiles = list(pool.imap(physical_rows_for_body, bodies, chunksize=1))
    profiles.sort(key=lambda row: row[0])
    rows = [
        (body, L, z1, L // gcd(L, z1))
        for body, L, survivors, _cap, _floor in profiles
        for z1 in survivors
    ]
    require(len(rows) == EXPECTED_PHYSICAL_ROWS, "physical row count changed")
    require(len(set((body, z1) for body, _L, z1, _d1 in rows)) == len(rows), "duplicate row")
    require(max(z1 for _body, _L, z1, _d1 in rows) == 380, "physical frontier changed")
    return profiles, rows


def selected_physical_controls(profiles: list[tuple]) -> int:
    by_body = {row[0]: row for row in profiles}
    controls = (
        (1, 2, 3, 4, 5, 6),
        (1, 4, 8, 10, 12, 14),
        (2, 4, 8, 10, 12, 14),
        (1, 5, 9, 11, 12, 14),
        (6, 9, 10, 12, 13, 14),
    )
    for body in controls:
        canonical = physical.profile(body)
        expected = canonical["constrained_counts"]["projected"][3]
        require(len(by_body[body][2]) == expected, ("canonical body count mismatch", body))
        canonical_frontier = canonical["constrained_frontier"]["projected"][3]
        observed_max = max(by_body[body][2], default=None)
        expected_max = None if canonical_frontier is None else canonical_frontier["first"]
        require(observed_max == expected_max, ("canonical body frontier mismatch", body))
    return len(controls)


def build_support_rows() -> tuple[dict, int, int]:
    by_body = defaultdict(list)
    body_divisor_rows = 0
    for body in combinations(range(1, 15), 6):
        L, ranges = support.safe_cell_ranges(body)
        for D in divisors_of(L):
            body_divisor_rows += 1
            support_count = support.support_size_bitset(D, ranges)
            if Q(support_count, D) > SUPPORT_CUTOFF:
                continue
            require(D % 7 == 0, ("nonseptimal support row", body, D))
            arcs = combined.projected_support_arcs(D, ranges)
            require(
                sum(right - left for left, right in arcs) == support_count,
                ("projected support size mismatch", body, D),
            )
            q = D // 7
            histogram_q = combined.residue_load_histogram(arcs, q)
            N_by_c = tuple(
                sum(count for load, count in histogram_q if load > c)
                for c in range(5)
            )
            small_loads = {}
            for d in (2, 3, 4):
                if D % d == 0:
                    histogram_d = combined.residue_load_histogram(arcs, d)
                    small_loads[d] = combined.top_class_load(histogram_d, 1)
            by_body[body].append((D, support_count, N_by_c, small_loads))
        by_body[body].sort()
    require(body_divisor_rows == 251_536, "body/divisor universe changed")
    require(sum(map(len, by_body.values())) == EXPECTED_SUPPORT_ROWS, "support row count changed")
    require(
        len({D for records in by_body.values() for D, *_rest in records})
        == EXPECTED_SUPPORT_DIVISORS,
        "support divisor universe changed",
    )
    return by_body, body_divisor_rows, len(by_body)


def expected_spike_passes(N_c: int, c: int, q: int) -> bool:
    return N_c == 0 or 55 * N_c < 13 * (4 - c) * q


def suffix_counter(counter: Counter) -> tuple[tuple[int, ...], tuple[int, ...]]:
    keys = tuple(sorted(counter))
    suffix = [0] * (len(keys) + 1)
    for index in range(len(keys) - 1, -1, -1):
        suffix[index] = suffix[index + 1] + counter[keys[index]]
    return keys, tuple(suffix)


def count_at_least(
    keys: tuple[int, ...], suffix: tuple[int, ...], threshold: int
) -> int:
    import bisect

    return suffix[bisect.bisect_left(keys, threshold)]


def distribution_summary(distribution: Counter) -> tuple:
    """Compact all row queries before the large feature table is discarded."""
    by_c = Counter()
    by_pattern = defaultdict(Counter)
    by_pattern_c = defaultdict(Counter)
    for (pattern, large_capacity, c), multiplicity in distribution.items():
        by_c[c] += multiplicity
        by_pattern[pattern][large_capacity] += multiplicity
        by_pattern_c[(pattern, c)][large_capacity] += multiplicity
    return (
        sum(distribution.values()),
        by_c,
        {pattern: suffix_counter(counter) for pattern, counter in by_pattern.items()},
        {key: suffix_counter(counter) for key, counter in by_pattern_c.items()},
    )


def classify_all_groups(grouped_rows: dict, support_rows: dict) -> tuple:
    """Process by D, reusing each conditional GF and then releasing it."""
    stages = ("raw", "expected_spike", "support_status", "joint")
    group_counts = {key: Counter() for key in grouped_rows}
    group_divisors = {
        key: {stage: set() for stage in stages} for key in grouped_rows
    }
    d1_by_body = defaultdict(set)
    for body, d1 in grouped_rows:
        d1_by_body[body].add(d1)
    support_by_D = defaultdict(list)
    for body, records in support_rows.items():
        for D, support_count, N_by_c, small_loads in records:
            support_by_D[D].append((body, support_count, N_by_c, small_loads))

    coefficient_semantic = hashlib.sha256()
    pair_count = 0
    for D in sorted(support_by_D):
        records_by_d1 = defaultdict(list)
        for body, support_count, N_by_c, small_loads in support_by_D[D]:
            for d1 in d1_by_body[body]:
                if D % d1 == 0:
                    records_by_d1[d1].append(
                        (body, support_count, N_by_c, small_loads)
                    )
        for d1 in sorted(records_by_d1):
            distribution = conditional_feature_distribution(D, d1)
            pair_count += 1
            for feature, multiplicity in sorted(distribution.items()):
                coefficient_semantic.update(
                    f"{D}|{d1}|{feature}|{multiplicity}\n".encode()
                )
            total, by_c, by_pattern, by_pattern_c = distribution_summary(distribution)
            for body, support_count, N_by_c, small_loads in records_by_d1[d1]:
                key = (body, d1)
                expected_cs = {
                    c
                    for c in by_c
                    if expected_spike_passes(N_by_c[c], c, D // 7)
                }
                expected_count = sum(by_c[c] for c in expected_cs)
                support_count_here = 0
                joint_count = 0
                for pattern, (keys, suffix) in by_pattern.items():
                    weights, marginals = small_vectors(pattern, small_loads)
                    threshold = support_count - status_need_limit(weights, marginals)
                    support_count_here += count_at_least(keys, suffix, threshold)
                    for c in expected_cs:
                        pattern_c = by_pattern_c.get((pattern, c))
                        if pattern_c is not None:
                            joint_count += count_at_least(
                                pattern_c[0], pattern_c[1], threshold
                            )
                for stage, count in (
                    ("raw", total),
                    ("expected_spike", expected_count),
                    ("support_status", support_count_here),
                    ("joint", joint_count),
                ):
                    group_counts[key][stage] += count
                    if count:
                        group_divisors[key][stage].add(D)
            # The next pair does not need the full feature table.  Clearing
            # prevents the all-physical bridge from retaining gigabytes of
            # exact-lcm coefficient dictionaries.
            conditional_feature_distribution.cache_clear()

    classifications = {
        key: (
            tuple(group_counts[key][stage] for stage in stages),
            tuple(tuple(sorted(group_divisors[key][stage])) for stage in stages),
        )
        for key in grouped_rows
    }
    return classifications, pair_count, coefficient_semantic.hexdigest()


def row_record(row: tuple, classification: tuple) -> tuple:
    body, L, z1, d1 = row
    counts, divisors = classification
    return (z1, body, L, d1, counts, tuple(len(values) for values in divisors))


def summarize(
    rows: list[tuple], group_classifications: dict, excluded: set[tuple]
) -> dict:
    stages = ("physical", "raw", "expected_spike", "support_status", "joint")
    row_counts = Counter()
    profile_occurrences = Counter()
    max_rows = {stage: None for stage in stages}
    min_rows = {stage: None for stage in stages}
    killed_min = {stage: None for stage in stages[1:]}
    semantic = hashlib.sha256()
    by_d1_joint = Counter()
    frontier_z = {stage: None for stage in stages}
    frontier_rows = {stage: [] for stage in stages}
    relation_names = (
        "raw_not_expected",
        "raw_not_support",
        "expected_not_joint",
        "support_not_joint",
    )
    relation_counts = Counter()
    relation_min = {name: None for name in relation_names}
    relation_max = {name: None for name in relation_names}
    relation_small = {name: [] for name in relation_names}
    relation_semantic = {name: hashlib.sha256() for name in relation_names}
    for row in rows:
        body, L, z1, d1 = row
        if (body, z1) in excluded:
            continue
        classification = group_classifications[(body, d1)]
        counts = classification[0]
        record = row_record(row, classification)
        semantic.update(f"{record}\n".encode())
        row_counts["physical"] += 1
        physical_key = (z1, body, L, d1)
        if min_rows["physical"] is None or physical_key < min_rows["physical"]:
            min_rows["physical"] = physical_key
        if max_rows["physical"] is None or physical_key > max_rows["physical"]:
            max_rows["physical"] = physical_key
        if frontier_z["physical"] is None or z1 > frontier_z["physical"]:
            frontier_z["physical"] = z1
            frontier_rows["physical"] = [physical_key]
        elif z1 == frontier_z["physical"]:
            frontier_rows["physical"].append(physical_key)
        for index, stage in enumerate(stages[1:]):
            profile_occurrences[stage] += counts[index]
            if counts[index]:
                row_counts[stage] += 1
                if min_rows[stage] is None or physical_key < min_rows[stage]:
                    min_rows[stage] = physical_key
                if max_rows[stage] is None or physical_key > max_rows[stage]:
                    max_rows[stage] = physical_key
                if frontier_z[stage] is None or z1 > frontier_z[stage]:
                    frontier_z[stage] = z1
                    frontier_rows[stage] = [physical_key]
                elif z1 == frontier_z[stage]:
                    frontier_rows[stage].append(physical_key)
                if stage == "joint":
                    by_d1_joint[d1] += 1
            elif killed_min[stage] is None or physical_key < killed_min[stage]:
                killed_min[stage] = physical_key
        predicates = {
            "raw_not_expected": bool(counts[0] and not counts[1]),
            "raw_not_support": bool(counts[0] and not counts[2]),
            "expected_not_joint": bool(counts[1] and not counts[3]),
            "support_not_joint": bool(counts[2] and not counts[3]),
        }
        for name, present in predicates.items():
            if not present:
                continue
            relation_counts[name] += 1
            relation_semantic[name].update(f"{record}\n".encode())
            if relation_min[name] is None or physical_key < relation_min[name]:
                relation_min[name] = physical_key
            if relation_max[name] is None or physical_key > relation_max[name]:
                relation_max[name] = physical_key
            if relation_small[name] is not None:
                relation_small[name].append(record)
                if len(relation_small[name]) > 20:
                    relation_small[name] = None
    return {
        "row_counts": row_counts,
        "profile_occurrences": profile_occurrences,
        "min_rows": min_rows,
        "max_rows": max_rows,
        "killed_min": killed_min,
        "semantic": semantic.hexdigest(),
        "joint_d1_top": by_d1_joint.most_common(30),
        "frontier_z": frontier_z,
        "frontier_rows": {
            stage: tuple(sorted(values)) for stage, values in frontier_rows.items()
        },
        "relation_counts": relation_counts,
        "relation_min": relation_min,
        "relation_max": relation_max,
        "relation_small": relation_small,
        "relation_semantic": {
            name: digest.hexdigest() for name, digest in relation_semantic.items()
        },
    }


def render_summary(name: str, summary: dict) -> list[str]:
    lines = [f"universe={name}"]
    lines.append(f"  row_counts={summary['row_counts']}")
    lines.append(f"  profile_occurrences={summary['profile_occurrences']}")
    lines.append(f"  min_rows={summary['min_rows']}")
    lines.append(f"  max_rows={summary['max_rows']}")
    lines.append(f"  first_killed_rows={summary['killed_min']}")
    lines.append(f"  frontier_z={summary['frontier_z']}")
    lines.append(f"  frontier_rows={summary['frontier_rows']}")
    lines.append(f"  relation_counts={summary['relation_counts']}")
    lines.append(f"  relation_min={summary['relation_min']}")
    lines.append(f"  relation_max={summary['relation_max']}")
    lines.append(f"  relation_small={summary['relation_small']}")
    lines.append(f"  relation_semantic_sha256={summary['relation_semantic']}")
    lines.append(f"  joint_d1_top30={summary['joint_d1_top']}")
    lines.append(f"  semantic_sha256={summary['semantic']}")
    return lines


def main(workers: int, output: Path | None) -> None:
    require(workers >= 1, "workers must be positive")
    c_controls, feature_controls = brute_conditional_controls()
    profiles, rows = reconstruct_physical_rows(workers)
    canonical_body_controls = selected_physical_controls(profiles)
    support_rows, body_divisor_rows, support_body_count = build_support_rows()

    grouped_rows = defaultdict(list)
    for row in rows:
        grouped_rows[(row[0], row[3])].append(row)
    group_classifications, distribution_pairs, coefficient_semantic = (
        classify_all_groups(grouped_rows, support_rows)
    )

    closed_380 = {((1, 4, 8, 10, 12, 14), 380)}
    closed_250 = {((1, 4, 8, 10, 12, 14), 250)}
    closed_378 = {(body, z1) for body, _L, z1, _d1 in rows if z1 == 378}
    closed_364 = {(body, z1) for body, _L, z1, _d1 in rows if z1 == 364}
    closed_350 = {(body, z1) for body, _L, z1, _d1 in rows if z1 == 350}
    closed_336 = {(body, z1) for body, _L, z1, _d1 in rows if z1 == 336}
    closed_330 = {(body, z1) for body, _L, z1, _d1 in rows if z1 == 330}
    closed_328 = {(body, z1) for body, _L, z1, _d1 in rows if z1 == 328}
    require(closed_380 <= {(body, z1) for body, _L, z1, _d1 in rows}, "z380 row missing")
    require(closed_250 <= {(body, z1) for body, _L, z1, _d1 in rows}, "z250 row missing")
    require(len(closed_378) == 9, "z378 physical row count changed")
    require(len(closed_364) == 25, "z364 physical row count changed")
    require(len(closed_350) == 53, "z350 physical row count changed")
    require(len(closed_336) == 8, "z336 physical row count changed")
    require(len(closed_330) == 1, "z330 physical row count changed")
    require(len(closed_328) == 9, "z328 physical row count changed")

    z328_rows = [row for row in rows if row[2] == 328]
    z328_summary = summarize(z328_rows, group_classifications, set())
    z328_expected_survivors = tuple(
        sorted(
            body
            for body, _L, _z1, d1 in z328_rows
            if group_classifications[(body, d1)][0][1]
        )
    )
    z328_expected_kills = tuple(
        sorted(body for body, _L, _z1, _d1 in z328_rows)
    )
    z328_expected_kills = tuple(
        body for body in z328_expected_kills if body not in z328_expected_survivors
    )
    require(
        z328_summary["row_counts"]
        == Counter(
            physical=9,
            raw=9,
            support_status=9,
            expected_spike=4,
            joint=4,
        ),
        "z328 necessary-state split changed",
    )
    require(
        z328_summary["profile_occurrences"]
        == Counter(
            raw=481_415,
            support_status=376_039,
            expected_spike=107_748,
            joint=107_740,
        ),
        "z328 profile split changed",
    )
    require(
        z328_expected_kills
        == (
            (1, 2, 6, 8, 12, 14),
            (1, 4, 6, 8, 12, 14),
            (1, 4, 8, 10, 12, 14),
            (2, 4, 6, 8, 12, 14),
            (2, 4, 8, 10, 12, 14),
        ),
        "z328 expected-spike kill set changed",
    )
    require(
        z328_expected_survivors
        == (
            (1, 2, 8, 10, 12, 14),
            (1, 6, 8, 10, 12, 14),
            (2, 6, 8, 10, 12, 14),
            (2, 8, 10, 11, 12, 14),
        ),
        "z328 expected-spike survivor set changed",
    )

    exclusions = (
        ("old_376020", set()),
        ("minus_closed_z380", closed_380),
        ("minus_closed_z380_z250", closed_380 | closed_250),
        (
            "pre_z364_control_minus_z380_z250_z378",
            closed_380 | closed_250 | closed_378,
        ),
        (
            "frozen_exact_minus_z380_z250_z378_z364",
            closed_380 | closed_250 | closed_378 | closed_364,
        ),
        (
            "pre_z336_control_minus_z380_z250_z378_z364_z350",
            closed_380 | closed_250 | closed_378 | closed_364 | closed_350,
        ),
        (
            "newest_exact_minus_z380_z250_z378_z364_z350_z336",
            closed_380
            | closed_250
            | closed_378
            | closed_364
            | closed_350
            | closed_336,
        ),
        (
            "minus_z380_z250_z378_z364_z350_z336_z330",
            closed_380
            | closed_250
            | closed_378
            | closed_364
            | closed_350
            | closed_336
            | closed_330,
        ),
        (
            "current_exact_minus_z380_z250_z378_z364_z350_z336_z330_z328",
            closed_380
            | closed_250
            | closed_378
            | closed_364
            | closed_350
            | closed_336
            | closed_330
            | closed_328,
        ),
    )
    summaries = [
        (name, summarize(rows, group_classifications, excluded))
        for name, excluded in exclusions
    ]
    require(summaries[0][1]["row_counts"]["physical"] == 376_020, "old ledger changed")
    require(summaries[1][1]["row_counts"]["physical"] == 376_019, "z380 exclusion changed")
    require(summaries[2][1]["row_counts"]["physical"] == 376_018, "z250 exclusion changed")
    require(summaries[3][1]["row_counts"]["physical"] == 376_009, "current ledger changed")
    require(
        summaries[4][1]["row_counts"]["physical"] == 375_984,
        "incoming z364 exclusion changed",
    )
    require(
        summaries[5][1]["row_counts"]
        == Counter(
            physical=375_931,
            raw=375_931,
            support_status=375_931,
            expected_spike=247_572,
            joint=247_571,
        ),
        "incoming z350 exclusion changed",
    )
    require(
        summaries[6][1]["row_counts"]
        == Counter(
            physical=375_923,
            raw=375_923,
            support_status=375_923,
            expected_spike=247_570,
            joint=247_569,
        ),
        "incoming z336 exclusion changed",
    )
    require(
        summaries[7][1]["row_counts"]
        == Counter(
            physical=375_922,
            raw=375_922,
            support_status=375_922,
            expected_spike=247_570,
            joint=247_569,
        ),
        "z330 exclusion changed",
    )
    require(
        summaries[7][1]["profile_occurrences"]
        == Counter(
            raw=75_423_449_780,
            support_status=62_057_393_714,
            expected_spike=18_784_011_176,
            joint=18_778_518_180,
        ),
        "z330 profile ledger changed",
    )
    require(
        summaries[8][1]["row_counts"]
        == Counter(
            physical=375_913,
            raw=375_913,
            support_status=375_913,
            expected_spike=247_566,
            joint=247_565,
        ),
        "z328 exclusion changed",
    )
    require(
        summaries[8][1]["profile_occurrences"]
        == Counter(
            raw=75_422_968_365,
            support_status=62_057_017_675,
            expected_spike=18_783_903_428,
            joint=18_778_410_440,
        ),
        "z328 profile ledger changed",
    )
    require(
        summaries[8][1]["frontier_z"]["physical"] == 324,
        "post-z328 physical frontier changed",
    )

    for _name, summary in summaries:
        require(
            summary["row_counts"]["raw"]
            <= summary["row_counts"]["physical"],
            "raw resolving rows exceed physical rows",
        )
        require(
            summary["row_counts"]["expected_spike"]
            <= summary["row_counts"]["raw"],
            "expected-spike rows exceed raw resolving rows",
        )
        require(
            summary["row_counts"]["support_status"]
            <= summary["row_counts"]["raw"],
            "support-status rows exceed raw resolving rows",
        )
        require(
            summary["row_counts"]["joint"]
            <= min(
                summary["row_counts"]["expected_spike"],
                summary["row_counts"]["support_status"],
            ),
            "joint rows exceed a parent screen",
        )

    feature_semantic = hashlib.sha256()
    for (body, d1), classification in sorted(group_classifications.items()):
        counts, divisors = classification
        feature_semantic.update(
            f"{body}|{d1}|{counts}|{tuple(map(len, divisors))}\n".encode()
        )

    lines = [
        "LRC14 k=3 physical-row / conditional-denominator reconciliation",
        f"physical_source_sha256={INPUT_SHA256['physical_source']}",
        f"projection_source_sha256={INPUT_SHA256['projection_source']}",
        f"closed_z380_source_sha256={INPUT_SHA256['z380_closure_source']}",
        f"closed_z380_output_sha256={INPUT_SHA256['z380_closure_output']}",
        f"closed_z250_source_sha256={INPUT_SHA256['z250_closure_source']}",
        f"closed_z250_output_sha256={INPUT_SHA256['z250_closure_output']}",
        f"closed_z378_source_sha256={INPUT_SHA256['z378_closure_source']}",
        f"closed_z378_output_sha256={INPUT_SHA256['z378_closure_output']}",
        f"incoming_z364_closure_source_sha256={INPUT_SHA256['z364_closure_source']}",
        f"incoming_z364_closure_output_sha256={INPUT_SHA256['z364_closure_output']}",
        f"incoming_z350_closure_source_sha256={INPUT_SHA256['z350_closure_source']}",
        f"incoming_z350_closure_output_sha256={INPUT_SHA256['z350_closure_output']}",
        f"incoming_z336_closure_source_sha256={INPUT_SHA256['z336_closure_source']}",
        f"incoming_z336_closure_output_sha256={INPUT_SHA256['z336_closure_output']}",
        f"closed_z330_scalar_source_sha256={INPUT_SHA256['z330_scalar_source']}",
        f"closed_z330_scalar_output_sha256={INPUT_SHA256['z330_scalar_output']}",
        f"closed_z330_source_sha256={INPUT_SHA256['z330_closure_source']}",
        f"closed_z330_output_sha256={INPUT_SHA256['z330_closure_output']}",
        f"closed_z328_scalar_source_sha256={INPUT_SHA256['z328_scalar_source']}",
        f"closed_z328_scalar_output_sha256={INPUT_SHA256['z328_scalar_output']}",
        f"closed_z328_source_sha256={INPUT_SHA256['z328_closure_source']}",
        f"closed_z328_output_sha256={INPUT_SHA256['z328_closure_output']}",
        "physical_universe=(six-body E subset {1,...,14}, external z1>=15,"
        "THM-2941 projected suffix with exact labels through 7000 and rigorous omitted tail)",
        "conditional_universe=(fixed d1=L/gcd(L,z1), later denominator 3-multiset,"
        "D=exact lcm, support-transfer resolving rows only)",
        "scope=necessary upper relaxations; later numerators/order/distinctness/phase/common-table/high-wall not realized",
        f"mobius_literal_c_control_cases={c_controls}",
        f"mobius_literal_feature_control_cases={feature_controls}",
        f"canonical_physical_body_controls={canonical_body_controls}",
        f"physical_rows={len(rows)}",
        f"physical_bodies={len({body for body, _L, _z1, _d1 in rows})}",
        f"physical_body_d1_groups={len(grouped_rows)}",
        f"body_divisor_rows={body_divisor_rows}",
        f"support_body_count={support_body_count}",
        f"support_rows={sum(map(len, support_rows.values()))}",
        f"support_divisors={len({D for records in support_rows.values() for D, *_rest in records})}",
        f"conditional_GF_pairs={distribution_pairs}",
        f"conditional_coefficient_semantic_sha256={coefficient_semantic}",
        f"conditional_group_semantic_sha256={feature_semantic.hexdigest()}",
        f"closed_z380={sorted(closed_380)}",
        f"closed_z250={sorted(closed_250)}",
        f"closed_z378={sorted(closed_378)}",
        f"incoming_exact_closed_z364={sorted(closed_364)}",
        f"incoming_exact_closed_z350={sorted(closed_350)}",
        f"incoming_exact_closed_z336={sorted(closed_336)}",
        f"exact_status_closed_z330={sorted(closed_330)}",
        f"exact_status_closed_z328={sorted(closed_328)}",
        f"z328_expected_spike_kills_5={z328_expected_kills}",
        f"z328_expected_spike_survivors_4={z328_expected_survivors}",
    ]
    lines.extend(render_summary("z328_before_exact_status", z328_summary))
    for name, summary in summaries:
        lines.extend(render_summary(name, summary))
    lines.extend(
        [
            f"fractional_cover_cache={maximum_upward_mass.cache_info()}",
            f"status_limit_cache={status_need_limit.cache_info()}",
            "expected_spike_predicate=N_c==0 OR 55*N_c<13*(4-c)*q",
            "input_snapshot_guard=start/after_module_load/end:PASS",
            "all_exact_controls=PASS",
        ]
    )
    require_inputs_unchanged("end")
    text = "\n".join(lines) + "\n"
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=None)
    arguments = parser.parse_args()
    os.environ["PYTHONHASHSEED"] = str(arguments.hash_seed)
    main(arguments.workers, arguments.output)
