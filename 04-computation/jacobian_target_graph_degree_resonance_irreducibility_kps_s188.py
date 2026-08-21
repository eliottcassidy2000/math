#!/usr/bin/env python3
"""Exact controls for the THM-3564 target-graph resonance gate.

The finite tables below are hostile controls.  The all-degree conclusion is
the displayed valuation comparison, implemented by ``resonance_candidate``.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"failed truth gate: {label}")


a, b, x = sp.symbols("a b x")


def term_valuations(d: int, k: int) -> tuple[int, int, int]:
    """Valuations at a=infinity of the three core-cubic terms."""
    require(d >= 0, "nonnegative a-degree")
    if d == 0:
        return -2 + 3 * k, k, 0
    return -(2 + 2 * d) + 3 * k, -d + k, -d


def resonance_candidate(d: int) -> int | None:
    """The only k not excluded by a unique-minimum valuation."""
    if d >= 1 and d % 3 == 1:
        return (d + 2) // 3
    return None


# Direct symbolic degree controls.  The test polynomial has exact a-degree d;
# its lower terms prevent accidental monomial simplifications.
degree_rows = []
for d in range(0, 13):
    if d == 0:
        phi = b + 1
    else:
        phi = (b + 1) * a**d + (b**2 + 1) * a ** (d - 1) + b
    L_phi = sp.expand(27 * a**2 * phi**2 + 18 * a * b * phi + 16 * a - b**3 * phi - b**2)
    middle = sp.expand(4 + 3 * b * phi)
    constant = sp.expand(2 * phi)
    observed = (sp.degree(L_phi, a), sp.degree(middle, a), sp.degree(constant, a))
    expected = (2 + 2 * d, d, d) if d >= 1 else (2, 0, 0)
    require(observed == expected, f"coefficient degree row d={d}")
    degree_rows.append((d, *observed))


# Exhaust a broad two-dimensional window as an independent hostile control.
# Exactly one non-unique row occurs for each resonant d, at k=(d+2)/3.
valuation_rows = []
for d in range(0, 31):
    exceptional = []
    for k in range(-40, 41):
        values = term_valuations(d, k)
        if values.count(min(values)) != 1:
            exceptional.append(k)
    candidate = resonance_candidate(d)
    expected = [] if candidate is None else [candidate]
    require(exceptional == expected, f"valuation exception row d={d}")
    valuation_rows.append((d, candidate))


# Verify the all-integer inequalities underlying the case split at symbolic
# sample endpoints.  The proof in THM-3564 handles the unbounded tails.
for d in range(1, 31):
    for k in range(-40, 1):
        A, B, C = term_valuations(d, k)
        require(A < B and A < C, f"k<=0 unique A, d={d}, k={k}")
    for k in range(1, 41):
        A, B, C = term_valuations(d, k)
        require(B > C, f"k>0 excludes B, d={d}, k={k}")
        if 3 * k != d + 2:
            require(A != C, f"nonresonant A/C separation, d={d}, k={k}")

for k in range(-40, 41):
    values = term_valuations(0, k)
    require(values.count(min(values)) == 1, f"d=0 unique minimum k={k}")


print("THM-3564 target-graph degree-resonance audit")
print("coefficient a-degrees: d=0 -> (2,0,0); d>=1 -> (2+2d,d,d)")
print("d k_candidate")
for d, candidate in valuation_rows[:16]:
    print(f"{d} {'none' if candidate is None else candidate}")
print("for d=0, valuations are (-2+3k,k,0): unique minimum for every k in Z")
print("for d>=1, valuations are (-(2+2d)+3k,-d+k,-d)")
print("k<=0: first term uniquely minimal")
print("k>0: middle term exceeds constant; a tie can only satisfy 3k=d+2")
print("resonant degrees: d=1 mod 3; all other nonzero target graphs have irreducible core cubic")
print("all active truth gates passed")
