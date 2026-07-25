#!/usr/bin/env python3
"""Exact referee for THM-2169's deletion-adapted cap profiles.

All Selberg/Kraft checks use Fraction arithmetic.  The last two searches
also verify that 1247 is the best M=29 bound inside the symmetric profile
templates used by the proof; the theorem itself does not claim global
optimality.
"""

from fractions import Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def kraft(terms: tuple[tuple[int, int], ...]) -> Fraction:
    """Return sum multiplicity*14/(6*cap+13)."""
    return sum(
        (Fraction(14 * multiplicity, 6 * cap + 13) for cap, multiplicity in terms),
        Fraction(),
    )


PROFILES = {
    "max-low": ((1, 1), (105, 12)),
    "max-high": ((7, 1), (36, 12)),
    "other-low": ((1, 1), (35, 1), (126, 11)),
    "other-mid": ((7, 1), (9, 1), (46, 11)),
    "other-high": ((22, 1), (8, 1), (36, 11)),
    "other-top": ((27, 1), (9, 1), (34, 11)),
}


def deletion_bound(m: int, is_max: bool) -> tuple[str, int]:
    if is_max:
        if m <= 7:
            return "max-low", 106 * m
        return "max-high", 43 * m
    if m <= 7:
        return "other-low", 161 * m
    if m <= 22:
        return "other-mid", 55 * m
    if m <= 27:
        return "other-high", 44 * m
    return "other-top", 43 * m


def best_symmetric_max_deletion(m: int) -> tuple[int, int, int]:
    best = (10**9, -1, -1)
    for special in range(1, m):
        for common in range(1, 200):
            if kraft(((special, 1), (common, 12))) < 1:
                best = min(best, (m * (special + common), special, common))
                break
    return best


def best_symmetric_other_deletion(m: int) -> tuple[int, int, int, int]:
    best = (10**9, -1, -1, -1)
    for max_cap in range(1, m):
        for delete_cap in range(1, 200):
            for common in range(1, 200):
                if kraft(((max_cap, 1), (delete_cap, 1), (common, 11))) < 1:
                    bound = m * (delete_cap + max(max_cap, common))
                    best = min(best, (bound, max_cap, delete_cap, common))
                    break
    return best


def main() -> None:
    print("THM-2169 DELETION-ADAPTED ANISOTROPIC PROFILE REFEREE")
    print("all Kraft sums are exact rational numbers")
    print()
    for name, terms in PROFILES.items():
        value = kraft(terms)
        require(value < 1, f"profile is not strictly admissible: {name}")
        print(f"{name:>10}: {value.numerator}/{value.denominator} < 1")

    worst = (68, 1, "signed-unit")
    for m in range(2, 30):
        for is_max in (True, False):
            name, bound = deletion_bound(m, is_max)
            require(kraft(PROFILES[name]) < 1, f"bad selected profile at M={m}")
            max_cap = PROFILES[name][0][0]
            require(max_cap < m, f"independence cap failed at M={m}, profile={name}")
            require(bound <= 1247, f"bound exceeds 1247 at M={m}, profile={name}")
            worst = max(worst, (bound, m, name))

    print()
    print(f"exhaustive height audit M=1..29: worst={worst[0]} at M={worst[1]} ({worst[2]})")
    require(worst[0] == 1247, "unexpected universal deletion bound")

    best_max = best_symmetric_max_deletion(29)
    best_other = best_symmetric_other_deletion(29)
    print(
        "M=29 symmetric-template optimum, delete max:"
        f" bound={best_max[0]}, caps=({best_max[1]}; {best_max[2]}x12)"
    )
    print(
        "M=29 symmetric-template optimum, delete other:"
        f" bound={best_other[0]}, caps=({best_other[1]},"
        f" {best_other[2]}; {best_other[3]}x11)"
    )
    require(best_max[0] == best_other[0] == 1247, "template optimum changed")

    primitive_l1 = 11 * 1247 + 1246
    positive_carries = primitive_l1 - 1
    require(primitive_l1 == 14963 and positive_carries == 14962, "carry arithmetic failed")
    print(f"primitive twelve-coordinate l1 bound: {primitive_l1}")
    print(f"sign-sharp positive-row carry alphabet: {positive_carries}")
    print("all assertions passed")


if __name__ == "__main__":
    main()
