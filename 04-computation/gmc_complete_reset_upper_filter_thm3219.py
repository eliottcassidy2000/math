#!/usr/bin/env python3
"""Exact upper-filter controls for THM-3219."""

from __future__ import annotations

import ast
import hashlib
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = ROOT / "04-computation" / "gmc_pole_prefix_hasse_current_scout.py"
THM3209 = (
    ROOT / "01-canon" / "theorems"
    / "THM-3209-depth-eight-complete-quotient-reset-and-negative-singleton-tangent.md"
)
DEPENDENCIES = (
    (UPSTREAM,
     "151edb9b8cee4807d3f08dc17af32e45420021ba30dfd116c04d9fcaf8bbd5b7"),
    (THM3209,
     "52b97aa1bdf8a0e0984fd2668a7602836b8ac2c43d4950eab1856b547081d0a2"),
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(data).hexdigest()


def load_upstream_prefix(maximum_degree: int) -> dict[str, object]:
    tree = ast.parse(UPSTREAM.read_text(encoding="utf-8"))
    prefix = []
    for node in tree.body:
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name)
                        and target.id == "MAXIMUM_DEGREE"
                        for target in node.targets)):
            node.value = ast.Constant(maximum_degree)
        prefix.append(node)
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name)
                        and target.id == "UNIVERSE"
                        for target in node.targets)):
            break
    module = ast.fix_missing_locations(ast.Module(body=prefix, type_ignores=[]))
    namespace = {"__file__": str(UPSTREAM)}
    exec(compile(module, str(UPSTREAM), "exec"), namespace)
    return namespace


dependency_hashes = tuple((path.name, lf_sha256(path)) for path, _ in DEPENDENCIES)
for (path, expected), (_, actual) in zip(DEPENDENCIES, dependency_hashes):
    require(actual == expected, ("dependency hash drift", str(path), actual))

UP = load_upstream_prefix(5)
UP["residual_roots"] = lru_cache(maxsize=None)(UP["residual_roots"])
UP["all_monomial_values"] = lru_cache(maxsize=None)(UP["all_monomial_values"])
all_monomial_values = UP["all_monomial_values"]
dominant_row = UP["dominant_row"]
partitions = UP["partitions"]
reduced_poles = UP["reduced_poles"]
residual_roots = UP["residual_roots"]
BANK = UP["BANKS"][1]


def phi_complete(degree: int, removed: tuple[int, ...]) -> Fraction:
    return sum(
        coefficient * sum(
            all_monomial_values(
                residual_roots(1, row, 1, 3), removed
            )[shape]
            for shape in partitions(degree)
        )
        for coefficient, row in BANK
    )


def phi_power(degree: int, removed: tuple[int, ...]) -> Fraction:
    return sum(
        coefficient * all_monomial_values(
            residual_roots(1, row, 1, 3), removed
        )[(degree,)]
        for coefficient, row in BANK
    )


POLES, _ = reduced_poles(1, BANK, 1, 3)
RESET = tuple(residual_roots(1, dominant_row(1), 1, 3))
EXPECTED_RESET = (1, 3, 3, 4, 5, 6, 7, 8)
EXPECTED_POLES = (8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1)
require(RESET == EXPECTED_RESET, ("reset drift", RESET))
require(POLES == EXPECTED_POLES, ("pole bank drift", POLES))

remainder_counter = Counter(POLES) - Counter(RESET)
REMAINDER = tuple(sorted(remainder_counter.elements()))
require(REMAINDER == (1, 1, 1, 2, 2, 2, 4, 5),
        ("remainder drift", REMAINDER))


# The first live complete row and the zero power marginals, before any prefix.
BASE_COMPLETE = tuple(phi_complete(degree, ()) for degree in range(6))
BASE_POWER = tuple(phi_power(degree, ()) for degree in range(1, 6))
require(BASE_COMPLETE == (0, 0, 0, 0, 0, 1440),
        ("base complete rows", BASE_COMPLETE))
require(BASE_POWER == (0, 0, 0, 0, 0),
        ("base power rows", BASE_POWER))


specification = tuple(sorted(remainder_counter.items()))
rows = []
prefix_invariance_checks = 0
nonzero_negative_complete = 0
for multiplicities in product(*(range(bound + 1) for _, bound in specification)):
    tau = tuple(
        value
        for (value, _), multiplicity in zip(specification, multiplicities)
        for _ in range(multiplicity)
    )
    removed = tuple(sorted(RESET + tau))

    # Directly verify the algebraic prefix-invariance inputs on every physical
    # prefix, including the reset itself.
    complete = phi_complete(5, removed)
    power = phi_power(5, removed)
    require(complete == 1440 and power == 0,
            ("prefix invariance", tau, complete, power))
    prefix_invariance_checks += 2

    quotient = all_monomial_values(RESET, removed)
    quotient_power = quotient[(5,)]
    quotient_complete = sum(quotient[shape] for shape in partitions(5))
    expected_power = -sum(value**5 for value in tau)
    require(quotient_power == expected_power,
            ("negative alphabet power", tau, quotient_power, expected_power))
    response = complete * quotient_power - power * quotient_complete
    expected_response = -1440 * sum(value**5 for value in tau)
    require(response == expected_response,
            ("principal response", tau, response, expected_response))
    if tau:
        require(response < 0, ("nonempty completion not negative", tau))
        if quotient_complete:
            nonzero_negative_complete += 1
        rows.append((len(tau), tau, response))
    else:
        require(response == 0, "reset lost zero response")


DEPTH_CENSUS = tuple(Counter(depth for depth, _, _ in rows)[depth]
                     for depth in range(1, 9))
require(len(rows) == 63, ("nonempty completion count", len(rows)))
require(DEPTH_CENSUS == (4, 8, 12, 14, 12, 8, 4, 1),
        ("depth census", DEPTH_CENSUS))
require(nonzero_negative_complete == 25,
        ("multi-negative complete hostile count", nonzero_negative_complete))

sharp = min((abs(response), tau, response) for _, tau, response in rows)
require(sharp == (1440, (1,), -1440), ("sharp minimum", sharp))
full = next(response for depth, tau, response in rows
            if depth == 8 and tau == REMAINDER)
require(full == -6117120, ("full completion", full))


source = Path(__file__).read_text(encoding="utf-8")
require(not any(isinstance(node, ast.Assert) for node in ast.parse(source).body),
        "assert-sensitive top-level test")

print("THM-3219 complete-reset upper-filter exact control")
print("dependency_hashes=" + repr(dependency_hashes))
print("poles=" + repr(POLES))
print("reset=" + repr(RESET))
print("remainder=" + repr(REMAINDER))
print("base_complete_degrees0_through5=" + repr(BASE_COMPLETE))
print("base_power_degrees1_through5=" + repr(BASE_POWER))
print("physical_filter_states_including_reset=64")
print("nonempty_completion_states=" + repr(len(rows)))
print("completion_depth_census=" + repr(DEPTH_CENSUS))
print("prefix_invariance_checks=" + repr(prefix_invariance_checks))
print("nonzero_h5_negative_alphabet_hostiles="
      + repr(nonzero_negative_complete))
print("sharp_nonzero_response=" + repr(sharp))
print("full_completion_response=" + repr(full))
print("formula=G5_Q_plus_tau_(5)=-1440*sum_tau_r^5")
print("principal_filter_intersection=delta_reset")
print("scope=complete_upper_filter_only_not_unrelated_deeper_states")
print("all_exact_checks=PASS")
