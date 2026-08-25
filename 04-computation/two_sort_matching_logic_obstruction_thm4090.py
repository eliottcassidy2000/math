#!/usr/bin/env python3
"""Finite semantic boundary audit for THM-4090.

The symbolic theorem is proved in the companion canon record.  This script
exhausts small set-valued interpretations of ``f : b -> a`` and checks the
semantic formulas used in that proof.  It deliberately uses no imports whose
semantics can change under ``python -O``.
"""

from __future__ import annotations

import sys
from itertools import product

sys.dont_write_bytecode = True


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def gamma_holds(b_size: int, a_size: int, images: tuple[int, ...]) -> bool:
    """Totality of forall x,y:b. f(x and y), using subset bit masks."""

    full_a = (1 << a_size) - 1
    for x in range(b_size):
        for y in range(b_size):
            value = images[x] if x == y else 0
            if value != full_a:
                return False
    return True


def phi_holds(b_size: int) -> bool:
    """Totality of forall x:b. x."""

    intersection = (1 << b_size) - 1
    for x in range(b_size):
        intersection &= 1 << x
    return intersection == (1 << b_size) - 1


def feed_closure() -> set[tuple[str, str]]:
    """Reflexive-transitive closure for the sole feed edge b -> a."""

    closure = {("a", "a"), ("b", "b"), ("b", "a")}
    changed = True
    while changed:
        changed = False
        for left, middle in tuple(closure):
            for middle2, right in tuple(closure):
                if middle == middle2 and (left, right) not in closure:
                    closure.add((left, right))
                    changed = True
    return closure


def main() -> None:
    model_count = 0
    gamma_model_count = 0
    size_rows: list[tuple[int, int, int]] = []

    for b_size in range(1, 5):
        for a_size in range(1, 4):
            count = 0
            subsets_a = range(1 << a_size)
            for images in product(subsets_a, repeat=b_size):
                model_count += 1
                if gamma_holds(b_size, a_size, images):
                    gamma_model_count += 1
                    count += 1
                    require(phi_holds(b_size), "semantic entailment failed")
            require(count == (1 if b_size == 1 else 0), "wrong Gamma-model count")
            size_rows.append((b_size, a_size, count))

    require(model_count == 5050, "finite universe drift")
    require(gamma_model_count == 3, "Gamma-model total drift")
    require(phi_holds(1) and not phi_holds(2), "singleton boundary drift")
    feeds = feed_closure()
    require(("a", "b") not in feeds, "hypothesis sort unexpectedly feeds target")

    print("THM-4090 finite semantic audit")
    print(f"set-valued models checked: {model_count}")
    print(f"Gamma-models checked: {gamma_model_count}")
    print("Gamma-model counts by (|b|,|a|):")
    for b_size, a_size, count in size_rows:
        print(f"  ({b_size},{a_size}): {count}")
    print("feed closure: a->a, b->a, b->b; hostile a->b absent")
    print("all finite semantic gates PASS")


if __name__ == "__main__":
    main()
