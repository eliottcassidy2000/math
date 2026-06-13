#!/usr/bin/env python3
"""
lonely_runner_k13_microstaircase_s363.py

codex-2026-05-31 S363

Probe for the next Lonely Runner frontier after the public k<=12 result.
Here k=13, so n=k+1=14 and the conjecture is the fourteen-runner case in
the common "one runner stationary" reduction.

The recent k=10 and k=12 proof avoids an expensive lift by proving an
analytic tight-class lemma when n=k+1 is an odd prime.  This file checks what
survives when n=14.  The exact r/14 analogue fails, but actual p-grid times
have finer "micro-staircase" cells that can resolve the obstruction.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


N = 14
K = 13

# Found by deterministic local search, then verified exhaustively below.
GRID14_OBSTRUCTION = (7, 4, 9, 6, 7, 8, 5, 0, 1, 12, 13, 12, 7)

# A small right-hand perturbation cell: floor(14*i*alpha), i=1..13.
MICRO_STAIR = (0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4)

VERIFIER_LEDGER = [
    # p, |I(13,p,1)| up to equivalence, find-cover seconds, final set size
    (17, 99, 0.001502, 0),
    (19, 55, 0.000907, 0),
    (23, 6, 0.000298, 0),
    (29, 15511, 0.115932, 0),
    (31, 11831, 0.093718, 0),
    (37, 4048, 0.043737, 0),
    (41, 1001, 0.013994, 0),
    (43, 163507, 1.335027, 0),
    (47, 84023, 0.559538, 0),
    (53, 14414, 0.094802, 0),
    (59, 1135845, 8.499790, 0),
    (61, 845384, 6.618051, 0),
    (67, 240365, 2.082710, 0),
    (71, 3302243, 25.440498, 0),
    (73, 2611797, 19.905382, 0),
    (79, 515605, 3.868331, 0),
    (83, 115903, 0.887726, 0),
    (89, 11804159, 91.812971, 0),
    (97, 1619090, 12.470085, 0),
    (101, 12697411, 113.117640, 0),
]


@dataclass(frozen=True)
class Witness:
    s: int
    r: int
    bins: tuple[int, ...]
    residues: tuple[int, ...]


def rbf(n: int, k: int, r: int, denominator: int) -> tuple[int, ...]:
    """Return floor(n*{i*r/denominator}) for i=1..k using integers."""

    return tuple((n * ((i * r) % denominator)) // denominator for i in range(1, k + 1))


def residues_after_shift(
    vector: tuple[int, ...], n: int, s: int, bins: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple((s * value + bin_value) % n for value, bin_value in zip(vector, bins))


def is_good_residue(residue: int, n: int) -> bool:
    return residue not in (0, n - 1)


def resolve_on_grid(vector: tuple[int, ...], n: int, denominator: int) -> Witness | None:
    """Search s/n + r/denominator cells for a witness to the tight lift."""

    k = len(vector)
    for s in range(n):
        for r in range(denominator):
            bins = rbf(n, k, r, denominator)
            residues = residues_after_shift(vector, n, s, bins)
            if all(is_good_residue(residue, n) for residue in residues):
                return Witness(s=s, r=r, bins=bins, residues=residues)
    return None


def count_blocked_candidates(vector: tuple[int, ...], n: int, denominator: int) -> int:
    k = len(vector)
    blocked = 0
    for s in range(n):
        for r in range(denominator):
            bins = rbf(n, k, r, denominator)
            residues = residues_after_shift(vector, n, s, bins)
            if any(not is_good_residue(residue, n) for residue in residues):
                blocked += 1
    return blocked


def pattern_interval(n: int, pattern: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    lo = Fraction(0)
    hi = Fraction(1)
    for i, value in enumerate(pattern, start=1):
        lo = max(lo, Fraction(value, n * i))
        hi = min(hi, Fraction(value + 1, n * i))
    return lo, hi


def first_grid_point_in_interval(
    denominator: int, lo: Fraction, hi: Fraction
) -> int | None:
    for r in range(denominator):
        x = Fraction(r, denominator)
        if lo <= x < hi:
            return r
    return None


def fmt_frac(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def main() -> None:
    print("Lonely Runner k=13 micro-staircase probe (codex-2026-05-31 S363)")
    print("Reduced form: k=13 moving runners, threshold=1/14.\n")

    obstruction = GRID14_OBSTRUCTION
    print("Exact r/14 tight-lift audit")
    print(f"  obstruction_vector_mod_14={obstruction}")
    blocked = count_blocked_candidates(obstruction, N, N)
    print(f"  blocked_candidates={blocked}/{N * N}")
    print(f"  resolved_on_r_over_14={resolve_on_grid(obstruction, N, N) is not None}")
    print("  interpretation=the prime-field tight-class lemma does not copy directly to n=14")
    print()

    lo, hi = pattern_interval(N, MICRO_STAIR)
    residues = residues_after_shift(obstruction, N, 1, MICRO_STAIR)
    print("Micro-staircase cell")
    print(f"  bins={MICRO_STAIR}")
    print(f"  alpha_interval=[{fmt_frac(lo)}, {fmt_frac(hi)})")
    print(f"  interval_width={fmt_frac(hi - lo)}")
    print(f"  residues_with_s_1={residues}")
    print(f"  resolves_obstruction={all(is_good_residue(r, N) for r in residues)}")
    print()

    primes = [17, 19, 29, 43, 97, 181, 191, 199, 211, 251, 307, 401, 701]
    print("Prime-grid witnesses for the same obstruction")
    for p in primes:
        witness = resolve_on_grid(obstruction, N, p)
        if witness is None:
            print(f"  p={p}: unresolved")
            continue
        in_cell = first_grid_point_in_interval(p, lo, hi)
        cell_note = f" first_micro_r={in_cell}" if in_cell is not None else ""
        print(
            "  "
            f"p={p}: s={witness.s} r={witness.r} "
            f"bins={witness.bins} residues={witness.residues}{cell_note}"
        )
    print()

    print("Published-verifier k=13 small-prime ledger")
    print("  source=external vzsky/13-lonely-runners code, run locally in temp")
    print("  config=LrcVerifier<13>::Config = Squeeze<2>, Squeeze<3>, Print")
    total_find = sum(row[2] for row in VERIFIER_LEDGER)
    max_row = max(VERIFIER_LEDGER, key=lambda row: row[1])
    print(f"  checked_primes={len(VERIFIER_LEDGER)}")
    print(f"  all_final_sets_empty={all(row[3] == 0 for row in VERIFIER_LEDGER)}")
    print(f"  total_find_cover_seconds={total_find:.6f}")
    print(
        "  "
        f"largest_initial_bad_set=p{max_row[0]} "
        f"size={max_row[1]} find_cover_seconds={max_row[2]:.6f}"
    )
    for p, size, seconds, final_size in VERIFIER_LEDGER:
        print(
            "  "
            f"p={p:3d} I13p1_size={size:8d} "
            f"find_cover_s={seconds:10.6f} final_size={final_size}"
        )

    print()
    print("Takeaway")
    print("  The next proof does not yet follow from these data.")
    print("  However, the hard part is sharply isolated: computing I(13,p,1).")
    print("  Once that set is found for tested primes, a c=2 squeeze already empties it.")
    print("  The analytic replacement to chase is a micro-staircase lemma for n=14.")


if __name__ == "__main__":
    main()
