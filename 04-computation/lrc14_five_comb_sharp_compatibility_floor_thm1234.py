#!/usr/bin/env python3
"""Independent exact referee and bank launcher for THM-1234.

The Python half checks the folded pair formula against interval geometry,
the heavy-ratio cutoff/list, the disconnected-partition ledger, the witness,
the five-to-six averaging constants, and the tournament fingerprint.  It then
compiles and runs the companion C++ bank.  The bank uses exact integer lower
rounding only; no floating-point value participates in its certificate.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
from math import gcd
from pathlib import Path
from shutil import which
from subprocess import check_output
from tempfile import TemporaryDirectory


TARGET = F(44, 273)
THETA = F(29, 30576)


def require(condition: bool, message: str) -> None:
    """Optimization-stable certificate check (unlike Python ``assert``)."""
    if not condition:
        raise RuntimeError(message)


def folded(z: int) -> int:
    r = z % 14
    return r * (14 - r)


def rho(a: int, b: int) -> F:
    common = gcd(a, b)
    a //= common
    b //= common
    if a > b:
        a, b = b, a
    return F(4 * a * b + folded(a + b) - folded(b - a), 196 * a * b)


def danger_pieces(speed: int) -> tuple[tuple[F, F], ...]:
    radius = F(1, 14 * speed)
    pieces: list[tuple[F, F]] = []
    for tooth in range(speed):
        center = F(tooth, speed)
        left, right = center - radius, center + radius
        if left < 0:
            pieces.extend(((F(0), right), (F(1) + left, F(1))))
        elif right > 1:
            pieces.extend(((left, F(1)), (F(0), right - 1)))
        else:
            pieces.append((left, right))
    return tuple(pieces)


def rho_geometry(a: int, b: int) -> F:
    total = F(0)
    pieces_a = sorted(danger_pieces(a))
    pieces_b = sorted(danger_pieces(b))
    index_a = index_b = 0
    while index_a < len(pieces_a) and index_b < len(pieces_b):
        left_a, right_a = pieces_a[index_a]
        left_b, right_b = pieces_b[index_b]
        left, right = max(left_a, left_b), min(right_a, right_b)
        if left < right:
            total += right - left
        if right_a <= right_b:
            index_a += 1
        if right_b <= right_a:
            index_b += 1
    return total


def tournament_fingerprint(speeds: tuple[int, ...]) -> tuple:
    size = len(speeds)
    adjacency = [[False] * size for _ in range(size)]
    flips = 0
    ties: list[tuple[int, int]] = []
    tie_path = (0, 1, 3, 2, 4)
    tie_rank = {vertex: rank for rank, vertex in enumerate(tie_path)}
    for i, j in combinations(range(size), 2):
        defect = F(1, 49) - rho(speeds[i], speeds[j])
        if defect > 0:
            source, target = i, j
        elif defect < 0:
            source, target = j, i
            flips += 1
        else:
            # The declared tie Hamiltonian path supplies the zero gauge.
            source, target = sorted((i, j), key=tie_rank.__getitem__)
            ties.append((i, j))
        adjacency[source][target] = True

    scores = tuple(sum(row) for row in adjacency)
    triangles = tuple(
        vertices
        for vertices in combinations(range(size), 3)
        if all(
            sum(adjacency[i][j] for j in vertices if j != i) == 1
            for i in vertices
        )
    )
    reach = [
        [i == j or adjacency[i][j] for j in range(size)]
        for i in range(size)
    ]
    for k in range(size):
        for i in range(size):
            for j in range(size):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    seen: set[int] = set()
    scc_sizes: list[int] = []
    for i in range(size):
        if i not in seen:
            component = {
                j for j in range(size) if reach[i][j] and reach[j][i]
            }
            seen.update(component)
            scc_sizes.append(len(component))
    hamiltonian_paths = tuple(
        path
        for path in permutations(range(size))
        if all(adjacency[path[i]][path[i + 1]] for i in range(size - 1))
    )
    return scores, triangles, tuple(scc_sizes), len(hamiltonian_paths), flips, tuple(ties), hamiltonian_paths[0]


print("THM-1234 exact referee: sharp compatibility-sensitive five-comb floor")

geometry_checks = 0
for a in range(1, 151):
    for b in range(a + 1, 151):
        if gcd(a, b) != 1:
            continue
        require(rho(a, b) == rho_geometry(a, b), f"pair formula failed at {a}:{b}")
        geometry_checks += 1
require(geometry_checks == 6857, "geometry check count changed")
print(f"pair_formula_geometry_checks={geometry_checks}")

# If a reduced a<b is heavy, THM-1166's defect bound gives
# theta < 1/(7b), so b < 1/(7theta)=4368/29 and hence b<=150.
require(F(1, 7) / THETA == F(4368, 29), "heavy denominator cutoff changed")
heavy_types = tuple(
    (a, b)
    for b in range(2, 151)
    for a in range(1, b)
    if gcd(a, b) == 1 and F(1, 49) - rho(a, b) > THETA
)
require(len(heavy_types) == 72, "heavy ratio count changed")
require(max(b for _, b in heavy_types) == 125, "heavy ratio maximum denominator changed")
heavy_digest = sha256(
    " ".join(f"{a}:{b}" for a, b in heavy_types).encode()
).hexdigest()
require(
    heavy_digest == "785e3b41577df477c6c04a284efaeb99a319382cea130ae128dc4d6ce87520c3",
    "heavy ratio digest changed",
)
print(
    "heavy_ratio_bank: cutoff=150 count=72 actual_max_denominator=125 "
    f"sha256={heavy_digest}"
)
require(196 * 125**8 < 2**64, "pair denominator no longer fits uint64")
print(f"uint64_pair_denominator_ceiling={196 * 125**8} < {2**64}")

partition_rows = []
for partition in (
    (4, 1),
    (3, 2),
    (3, 1, 1),
    (2, 2, 1),
    (2, 1, 1, 1),
    (1, 1, 1, 1, 1),
):
    cross_edges = (25 - sum(part * part for part in partition)) // 2
    internal = sum(
        (
            F(1, 12)
            if part == 4
            else F(1, 24)
            if part == 3
            else F(1, 91)
            if part == 2
            else F(0)
        )
        for part in partition
    )
    lower = internal + cross_edges * (F(1, 49) - THETA)
    require(lower >= TARGET, f"disconnected partition failed: {partition}")
    partition_rows.append((partition, cross_edges, lower, lower - TARGET))
require(partition_rows[0][2] == TARGET, "4+1 partition is no longer sharp")
print("disconnected_partition_rows=" + repr(partition_rows))

witness = (1, 12, 13, 27, 156)
pair_word = sorted(rho(a, b) for a, b in combinations(witness, 2))
require(
    pair_word
    == sorted(2 * [F(1, 63), F(17, 819), F(1, 84), F(1, 91), F(23, 1092)]),
    "telemetry witness pair word changed",
)
witness_mass = sum(pair_word, F(0))
require(witness_mass == TARGET, "sharp witness mass changed")
require(witness_mass - F(13, 81) == F(5, 7371), "13/81 witness gap changed")
print("sharp_witness=(1,12,13,27,156) R=44/273 gap_over_13/81=5/7371")

fallback_rows = (
    (1, 12, 13, 27, 156),
    (1, 12, 13, 156, 351),
    (4, 9, 108, 117, 1404),
    (9, 52, 108, 117, 1404),
)
for row in fallback_rows:
    require(
        sum((rho(a, b) for a, b in combinations(row, 2)), F(0)) == TARGET,
        f"exact fallback row changed: {row}",
    )
print(f"connected_heavy_exact_fallback_rows={fallback_rows}")

# Six five-subsets contain each of the fifteen pairs four times.
six_floor = F(6, 4) * TARGET
require(six_floor == F(22, 91), "six-speed average changed")
require(F(6, 7) * six_floor == F(132, 637), "slow-gap RHS changed")
require(F(196, 13) * F(132, 637) == F(528, 169), "gcd RHS changed")
union_ceiling = F(6, 7) - F(1, 3) * six_floor
require(union_ceiling == F(212, 273), "union ceiling changed")
require(F(7, 6) * union_ceiling == F(106, 117), "common-period ratio changed")
print(
    "six_speed_consequences: R6>=22/91 mean=132/637 "
    "gcd_RHS=528/169-(72/13)delta union<=212/273 G0/a<106/117"
)

seven_floor = F(21, 10) * TARGET
require(seven_floor == F(22, 65), "seven-speed average changed")
require(F(2, 7) * seven_floor == F(44, 455), "seven-speed uncovered mass changed")
require(F(7) - 2 * seven_floor == F(411, 65), "seven-speed period bound changed")
fano_sum = F(256, 539) * seven_floor
require(fano_sum == F(512, 3185), "Fano line-period sum changed")
require(F(7) / fano_sum == F(22295, 512), "Fano pigeonhole changed")
print(
    "seven_speed_consequences: R7>=22/65 uncovered>=44/455 "
    "common_G/m<=411/65 fano_sum_m/G>=512/3185 "
    "some_fano_line_G/m<=22295/512"
)

fingerprint = tournament_fingerprint(witness)
require(
    fingerprint
    == (
        (3, 2, 2, 1, 2),
        ((0, 1, 4), (0, 2, 4), (1, 2, 3), (2, 3, 4)),
        (5,),
        13,
        4,
        (),
        (0, 1, 3, 2, 4),
    ),
    "tournament fingerprint changed",
)
print(
    "tournament: scores=(3,2,2,1,2) cycles=4 SCC=(5,) "
    "hamiltonian_paths=13 flips=4 ties=0 tie_path=(0,1,3,2,4)"
)

cpp_source = Path(__file__).with_name(
    "lrc14_five_comb_sharp_compatibility_floor_thm1234_bank.cpp"
)
compiler = which("clang++") or which("g++")
if compiler is None:
    raise RuntimeError("a C++17 compiler is required for the exact bank")
with TemporaryDirectory(prefix="thm1234-") as temporary:
    binary = Path(temporary) / "thm1234_bank"
    check_output(
        [compiler, "-O3", "-std=c++17", str(cpp_source), "-o", str(binary)],
        text=True,
    )
    bank_output = check_output([str(binary)], text=True)
print(bank_output, end="")
expected_bank_output = """heavy_threshold=29/30576 denominator_cutoff=150 types=72
level=2 states=72 attempts=144 pruned=0
level=3 states=10048 attempts=20592 pruned=180
level=4 states=1189651 attempts=4300386 pruned=1524044
five_extensions=678030510 fixed_certified=678030462 exact_fallback_occurrences=48
minimum_fixed_lower_numerator=161172154 / 1000000000 witness=1,12,13,27,156
exact_fallback_unique=4
fallback_row=1,12,13,27,156 exact_R=44/273
fallback_row=1,12,13,156,351 exact_R=44/273
fallback_row=4,9,108,117,1404 exact_R=44/273
fallback_row=9,52,108,117,1404 exact_R=44/273
fallback_pair_word=2*(1/63+17/819+1/84+1/91+23/1092)=44/273
fallback_grid_deficit_numerator=1958 / 273000000000
"""
require(bank_output == expected_bank_output, "exact C++ bank output changed")
print("PASS: sharp R5>=44/273 and sharpened R6/R7/Fano consequences")
