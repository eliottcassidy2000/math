#!/usr/bin/env python3
"""Finite semantic controls for THM-4090.

The proof of non-derivability is symbolic (the sort-flow induction in the
theorem file).  This companion exhausts every finite set-valued unary model
with 1 <= |M_b| <= 4 and 1 <= |M_a| <= 3 and checks the semantic half of the witness.
It also checks the directed sort-flow closure and a hostile formula obtained
by deleting the load-bearing intersection of two independently valued
variables.

Reproduction:
  python -B 04-computation/matching_logic_two_sort_obstruction_thm4090.py
  python -B -O 04-computation/matching_logic_two_sort_obstruction_thm4090.py
"""

from itertools import product


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def image_of_subset(operation: tuple[int, ...], subset_mask: int) -> int:
    """Pointwise extension of f : M_b -> P(M_a)."""
    image = 0
    for point, value_mask in enumerate(operation):
        if subset_mask & (1 << point):
            image |= value_mask
    return image


def gamma_holds(operation: tuple[int, ...], nb: int, na: int) -> bool:
    """Totality of forall x:b forall y:b. f(x and y)."""
    full_a = (1 << na) - 1
    for x in range(nb):
        for y in range(nb):
            intersection = (1 << x) & (1 << y)
            if image_of_subset(operation, intersection) != full_a:
                return False
    return True


def hostile_gamma_holds(operation: tuple[int, ...], nb: int, na: int) -> bool:
    """Totality of forall x:b. f(x), with the intersection removed."""
    full_a = (1 << na) - 1
    return all(image_of_subset(operation, 1 << x) == full_a for x in range(nb))


def phi_holds(nb: int) -> bool:
    """Totality of forall x:b. x."""
    intersection = (1 << nb) - 1
    for x in range(nb):
        intersection &= 1 << x
    return intersection == (1 << nb) - 1


def transitive_closure(vertices: tuple[str, ...], edges: set[tuple[str, str]]) -> set[tuple[str, str]]:
    closure = set(edges)
    closure.update((v, v) for v in vertices)
    changed = True
    while changed:
        changed = False
        for u in vertices:
            for v in vertices:
                for w in vertices:
                    if (u, v) in closure and (v, w) in closure and (u, w) not in closure:
                        closure.add((u, w))
                        changed = True
    return closure


def main() -> None:
    print("THM-4090 TWO-SORT MATCHING-LOGIC SEMANTIC AUDIT")
    print("signature: sorts b,a; one symbol f:b->a")
    print("Gamma = {forall x:b forall y:b. f(x and y)}")
    print("phi   = forall x:b. x")
    print()

    total_models = 0
    gamma_models = 0
    entailment_failures = 0
    hostile_models = 0
    hostile_entailment_failures = 0

    print("|M_b| |M_a| operations Gamma-models Gamma-and-not-phi hostile-models hostile-and-not-phi")
    for nb in range(1, 5):
        for na in range(1, 4):
            operations = (1 << na) ** nb
            row_gamma = 0
            row_bad = 0
            row_hostile = 0
            row_hostile_bad = 0
            for operation in product(range(1 << na), repeat=nb):
                if gamma_holds(operation, nb, na):
                    row_gamma += 1
                    if not phi_holds(nb):
                        row_bad += 1
                if hostile_gamma_holds(operation, nb, na):
                    row_hostile += 1
                    if not phi_holds(nb):
                        row_hostile_bad += 1
            total_models += operations
            gamma_models += row_gamma
            entailment_failures += row_bad
            hostile_models += row_hostile
            hostile_entailment_failures += row_hostile_bad
            print(
                f"{nb:5d} {na:5d} {operations:10d} {row_gamma:12d} "
                f"{row_bad:17d} {row_hostile:14d} {row_hostile_bad:19d}"
            )
            require(row_gamma == (1 if nb == 1 else 0), "unexpected Gamma model count")
            require(row_bad == 0, "finite semantic countermodel to Gamma |= phi")
            require(row_hostile == 1, "unexpected hostile model count")
            require(row_hostile_bad == (0 if nb == 1 else 1), "hostile did not expose the intersection")

    feed = transitive_closure(("b", "a"), {("b", "a")})
    require(("a", "b") not in feed, "hypothesis sort unexpectedly feeds conclusion sort")
    require(not phi_holds(2), "phi must be invalid on a two-element b-carrier")

    print()
    print(f"finite models exhausted = {total_models}")
    print(f"Gamma models            = {gamma_models}")
    print(f"Gamma-and-not-phi        = {entailment_failures}")
    print(f"hostile models           = {hostile_models}")
    print(f"hostile-and-not-phi      = {hostile_entailment_failures}")
    print("feed closure             =", sorted(feed))
    print("a feeds b                = False")
    print("phi valid at |M_b|=2     = False")
    print("FINITE SEMANTIC CONTROLS: PASS")


if __name__ == "__main__":
    main()
