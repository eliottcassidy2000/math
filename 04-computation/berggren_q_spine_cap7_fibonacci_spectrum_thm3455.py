#!/usr/bin/env python3
"""Exact deterministic companion for THM-3455.

The script intersects THM-3453's fifteen cap-seven divisor atoms with the
parabolic Berggren labels q_t=4t(t+1)+3.  It verifies the surviving roots,
the exact rank-4/6/7 priority spectrum over its complete period, the complete
Pisano-period pullback at t=F_n, and the resulting natural/harmonic densities.

The script does not search for a physical LRC time or current.  All
truth-bearing arithmetic is integral or rational, and explicit exceptions
remain active under ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from math import lcm
from pathlib import Path


EXPECTED_SEMANTIC_SHA256 = "330499ea2bcbf3d2a0da6d3512870ebfae83c9e8268c78e7291ed00f0d95e652"
NO_CAP7 = 99
T_PERIOD = 1683
FIBONACCI_PERIOD = 360

ATOM_RANKS = {
    8: 4,
    9: 4,
    10: 5,
    11: 6,
    12: 5,
    13: 7,
    14: 7,
    15: 6,
    23: 6,
    25: 6,
    29: 7,
    38: 7,
    51: 7,
    68: 7,
    148: 7,
}


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def q_label(t: int) -> int:
    require(t >= 0, ("negative q index", t))
    return 4 * t * (t + 1) + 3


def cap7_rank_from_atoms(modulus: int) -> int:
    ranks = [rank for atom, rank in ATOM_RANKS.items() if modulus % atom == 0]
    return min(ranks, default=NO_CAP7)


def cap7_rank_from_residues(t: int) -> int:
    value = q_label(t)
    if value % 9 == 0:
        return 4
    if value % 11 == 0:
        return 6
    if value % 51 == 0:
        return 7
    return NO_CAP7


def divisors(number: int) -> tuple[int, ...]:
    return tuple(value for value in range(1, number + 1) if number % value == 0)


def minimal_period(sequence: tuple[int, ...]) -> int:
    size = len(sequence)
    return next(
        period
        for period in divisors(size)
        if all(sequence[index] == sequence[index % period] for index in range(size))
    )


atoms = tuple(sorted(ATOM_RANKS))
atom_roots = {
    atom: tuple(t for t in range(atom) if q_label(t) % atom == 0)
    for atom in atoms
}
expected_atom_roots = {
    8: (),
    9: (2, 6),
    10: (),
    11: (1, 9),
    12: (),
    13: (),
    14: (),
    15: (),
    23: (),
    25: (),
    29: (),
    38: (),
    51: (3, 20, 30, 47),
    68: (),
    148: (),
}
require(atom_roots == expected_atom_roots, ("atom roots", atom_roots))

for t in range(T_PERIOD):
    require(
        cap7_rank_from_atoms(q_label(t)) == cap7_rank_from_residues(t),
        ("atom sieve", t, q_label(t)),
    )

t_rank_sequence = tuple(cap7_rank_from_residues(t) for t in range(T_PERIOD))
t_rank_counts = {rank: t_rank_sequence.count(rank) for rank in (4, 6, 7, NO_CAP7)}
require(t_rank_counts == {4: 374, 6: 238, 7: 72, NO_CAP7: 999}, t_rank_counts)
require(minimal_period(t_rank_sequence) == T_PERIOD, "nonminimal t period")

t_densities = {
    4: Fraction(2, 9),
    6: Fraction(14, 99),
    7: Fraction(8, 187),
    NO_CAP7: Fraction(111, 187),
}
for rank, count in t_rank_counts.items():
    require(Fraction(count, T_PERIOD) == t_densities[rank], ("t density", rank, count))
require(sum((t_densities[rank] for rank in (4, 6, 7)), Fraction(0)) == Fraction(76, 187), "t support density")

# The first four geometric spine labels realize all four outcomes.
require(
    tuple((t, q_label(t), cap7_rank_from_residues(t)) for t in range(1, 5))
    == ((1, 11, 6), (2, 27, 4), (3, 51, 7), (4, 83, NO_CAP7)),
    "initial controls",
)
require(q_label(20) == 1683 and cap7_rank_from_residues(20) == 4, "overlap priority")
require(all(cap7_rank_from_residues(t) != 5 for t in range(T_PERIOD)), "rank five appeared")


# Complete Fibonacci pullback modulo lcm(9,11,51).
require(lcm(9, 11, 51) == T_PERIOD, "wrong joint modulus")
fibonacci = [0, 1]
for _ in range(FIBONACCI_PERIOD + 1):
    fibonacci.append((fibonacci[-1] + fibonacci[-2]) % T_PERIOD)
require((fibonacci[FIBONACCI_PERIOD], fibonacci[FIBONACCI_PERIOD + 1]) == (0, 1), "Pisano return")
for proper in divisors(FIBONACCI_PERIOD)[:-1]:
    require((fibonacci[proper], fibonacci[proper + 1]) != (0, 1), ("early Pisano return", proper))

fibonacci_rank_sequence = tuple(
    cap7_rank_from_residues(fibonacci[n])
    for n in range(FIBONACCI_PERIOD)
)
fibonacci_rank_counts = {
    rank: fibonacci_rank_sequence.count(rank)
    for rank in (4, 6, 7, NO_CAP7)
}
require(
    fibonacci_rank_counts == {4: 60, 6: 90, 7: 22, NO_CAP7: 188},
    fibonacci_rank_counts,
)
require(minimal_period(fibonacci_rank_sequence) == FIBONACCI_PERIOD, "nonminimal Fibonacci rank period")

rank4_residues_mod24 = tuple(index for index in range(24) if fibonacci_rank_sequence[index] == 4)
rank6_residues_mod120 = tuple(index for index in range(120) if fibonacci_rank_sequence[index] == 6)
rank7_residues_mod360 = tuple(index for index in range(360) if fibonacci_rank_sequence[index] == 7)
require(rank4_residues_mod24 == (3, 16, 20, 21), rank4_residues_mod24)
require(
    rank6_residues_mod120
    == (1, 2, 9, 11, 12, 19, 22, 29, 31, 32, 39, 41, 42, 49, 52,
        59, 61, 62, 71, 72, 79, 81, 82, 89, 91, 101, 102, 109, 111, 119),
    rank6_residues_mod120,
)
require(
    rank7_residues_mod360
    == (4, 14, 28, 43, 76, 86, 100, 115, 134, 148, 158,
        173, 187, 206, 220, 230, 244, 245, 278, 316, 317, 350),
    rank7_residues_mod360,
)
require(minimal_period(tuple(rank == 4 for rank in fibonacci_rank_sequence)) == 24, "rank4 period")
require(minimal_period(tuple(rank == 6 for rank in fibonacci_rank_sequence)) == 120, "rank6 period")
require(minimal_period(tuple(rank == 7 for rank in fibonacci_rank_sequence)) == 360, "rank7 period")

fibonacci_densities = {
    4: Fraction(1, 6),
    6: Fraction(1, 4),
    7: Fraction(11, 180),
    NO_CAP7: Fraction(47, 90),
}
for rank, count in fibonacci_rank_counts.items():
    require(Fraction(count, FIBONACCI_PERIOD) == fibonacci_densities[rank], ("Fibonacci density", rank))
require(
    sum((fibonacci_densities[rank] for rank in (4, 6, 7)), Fraction(0)) == Fraction(43, 90),
    "Fibonacci support density",
)


repo_root = Path(__file__).resolve().parents[1]
dependency_paths = (
    Path("01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md"),
    Path("01-canon/theorems/THM-3415-zero-mode-cochain-global-rank-five-support.md"),
    Path("01-canon/theorems/THM-3416-zero-mode-cochain-global-rank-six-support.md"),
    Path("01-canon/theorems/THM-3453-global-literal-half-twist-cap-seven-support-classification.md"),
    Path("01-canon/theorems/THM-3454-fibonacci-selected-u-spine-farey-lorentz-isometry-and-one-tie-edge-order.md"),
)
dependency_hashes = tuple((path.as_posix(), lf_sha256(repo_root / path)) for path in dependency_paths)

semantic_payload = {
    "atom_roots": atom_roots,
    "t_period": T_PERIOD,
    "t_rank_counts": t_rank_counts,
    "t_densities": {rank: str(value) for rank, value in t_densities.items()},
    "fibonacci_period": FIBONACCI_PERIOD,
    "fibonacci_rank_counts": fibonacci_rank_counts,
    "fibonacci_densities": {rank: str(value) for rank, value in fibonacci_densities.items()},
    "rank4_residues_mod24": rank4_residues_mod24,
    "rank6_residues_mod120": rank6_residues_mod120,
    "rank7_residues_mod360": rank7_residues_mod360,
    "dependency_hashes": dependency_hashes,
}
semantic_sha256 = sha256(
    json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
).hexdigest()
if EXPECTED_SEMANTIC_SHA256:
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, (semantic_sha256, EXPECTED_SEMANTIC_SHA256))


print("THM-3455 exact companion")
print("status=FINITE-EXACT SUPPORT FOR AN ELEMENTARY COROLLARY OF THM-3453")
print("spine=q_t=4t(t+1)+3=(2t+1)^2+2=Q_(t-1);t>=1")
print(f"atom_roots={atom_roots}")
print("surviving_atoms=9,11,51;rank_priority=4,6,7")
print("rank4_iff=t mod 9 in {2,6}")
print("rank6_iff=not_rank4 and t mod 11 in {1,9}")
print("rank7_iff=not_rank4_or_6 and t mod 51 in {3,20,30,47}")
print("rank_gt7=otherwise;rank5_is_empty")
print(f"t_rank_period={minimal_period(t_rank_sequence)};counts={t_rank_counts}")
print("t_rank_densities=rank4:2/9;rank6:14/99;rank7:8/187;gt7:111/187")
print("t_cap7_density=76/187;harmonic_coefficient=76/187")
print("q_label_counting=cap7_q_t_up_to_X=(38/187)*sqrt(X)+O(1);reciprocal_q_subseries_converges")
print(f"Pisano_mod_1683={FIBONACCI_PERIOD};rank_period={minimal_period(fibonacci_rank_sequence)}")
print(f"Fibonacci_rank_counts={fibonacci_rank_counts}")
print(f"Fibonacci_rank4_mod24={rank4_residues_mod24}")
print(f"Fibonacci_rank6_mod120={rank6_residues_mod120}")
print(f"Fibonacci_rank7_mod360={rank7_residues_mod360}")
print("Fibonacci_index_densities=rank4:1/6;rank6:1/4;rank7:11/180;gt7:47/90")
print("Fibonacci_index_cap7_density=43/90;harmonic_coefficient=43/90")
print("Fibonacci_value_boundary=labelled_index_harmonic_series_diverges_but_reciprocal_F_n,positive_depth_F_n_minus_1,and_q_(F_n)_series_converge")
print("controls=(t,q_t,rank)=(1,11,6),(2,27,4),(3,51,7),(4,83,>7);t=20 gives q=1683 rank4 priority")
for path, value in dependency_hashes:
    print(f"dependency_sha256[{path}]={value}")
print(f"semantic_sha256={semantic_sha256}")
print("scope=literal_half_twist_has_fixed_common_zero-cochain_centre;no_other-centre,nonzero-current,decrement,LRC14,or_JC_consequence")
