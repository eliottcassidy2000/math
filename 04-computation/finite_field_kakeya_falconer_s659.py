#!/usr/bin/env python3
"""S659: finite-field Kakeya/Falconer carrier ledger.

This scout specializes the HYP-2234 frontier-carrier lane to finite-field
Kakeya/Falconer.  The point is not to claim progress on the external
conjectures.  The point is to import the correct small finite carrier:

  Kakeya   -> one affine line in every projective direction over F_p^2
  Falconer -> distance fibers Delta(E) from the quadratic form x^2 + y^2
  Repo use -> keep concurrency, pinned distances, unit directions, characters,
              and LRC shell/owner analogues before collapsing to cardinality.

Tournament Analysis:
  Vertices are proof carriers, not points or lines.
  Pairwise observable is side-channel retention for LRC/unit-distance/Falconer
  transfer, with external theorem usefulness as a secondary gauge.
  The tie Hamiltonian path is the first ranking in the majority profile.

Assumption challenge:
  Candidate vertices considered: points, lines, directions, distances,
  pinned centers, quadratic characters, polynomial certificates, shell labels,
  and proof obligations.  This script chooses proof carriers because the
  preserved predicate is "which side channel survives a scalar quotient?"
  It destroys exact embedding data and must reattach it for geometry proofs.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations, product
from random import Random
from typing import Iterable


Point = tuple[int, int]
Direction = int | str


def inv(a: int, p: int) -> int:
    return pow(a % p, -1, p)


def legendre(a: int, p: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def directions(p: int) -> list[Direction]:
    return list(range(p)) + ["v"]


def line_points(p: int, direction: Direction, intercept: int) -> frozenset[Point]:
    if direction == "v":
        return frozenset((intercept % p, y) for y in range(p))
    slope = int(direction)
    return frozenset((x, (slope * x + intercept) % p) for x in range(p))


def line_mask(p: int, direction: Direction, intercept: int) -> int:
    mask = 0
    for x, y in line_points(p, direction, intercept):
        mask |= 1 << (x * p + y)
    return mask


def points_from_mask(p: int, mask: int) -> frozenset[Point]:
    pts: set[Point] = set()
    for x in range(p):
        for y in range(p):
            if mask & (1 << (x * p + y)):
                pts.add((x, y))
    return frozenset(pts)


def quadratic_distance(a: Point, b: Point, p: int) -> int:
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) % p


def direction_class(dx: int, dy: int, p: int) -> Direction:
    dx %= p
    dy %= p
    if dx == 0:
        return "v"
    return (dy * inv(dx, p)) % p


def parabola_line_family(p: int) -> list[tuple[Direction, int]]:
    """Odd-prime Besicovitch/Kakeya parabola carrier in the affine plane."""
    family = [(m, (-m * m) % p) for m in range(p)]
    family.append(("v", 0))
    return family


def union_mask_for_family(p: int, family: Iterable[tuple[Direction, int]]) -> int:
    mask = 0
    for direction, intercept in family:
        mask |= line_mask(p, direction, intercept)
    return mask


def multiplicity_histogram(p: int, family: Iterable[tuple[Direction, int]]) -> Counter[int]:
    multiplicity: Counter[Point] = Counter()
    for direction, intercept in family:
        for point in line_points(p, direction, intercept):
            multiplicity[point] += 1
    return Counter(multiplicity.values())


@dataclass(frozen=True)
class CarrierProfile:
    p: int
    label: str
    size: int
    formula_size: int
    multiplicity_hist: Counter[int]
    distance_count: int
    pinned_min: int
    pinned_max: int
    pinned_full_points: int
    unit_edges: int
    unit_direction_count: int
    unit_degree_hist: Counter[int]
    paley_distance_hist: Counter[int]


def profile_family(p: int, label: str, family: list[tuple[Direction, int]]) -> CarrierProfile:
    mask = union_mask_for_family(p, family)
    pts = points_from_mask(p, mask)
    distances: set[int] = set()
    pinned: list[int] = []
    paley_distance_hist: Counter[int] = Counter()
    unit_edges = 0
    unit_dirs: set[Direction] = set()
    unit_degree: Counter[Point] = Counter()

    for a in pts:
        seen = {quadratic_distance(a, b, p) for b in pts}
        pinned.append(len(seen))
        distances |= seen

    pts_list = sorted(pts)
    for a, b in combinations(pts_list, 2):
        d = quadratic_distance(a, b, p)
        paley_distance_hist[legendre(d, p)] += 1
        if d == 1:
            unit_edges += 1
            unit_degree[a] += 1
            unit_degree[b] += 1
            unit_dirs.add(direction_class(a[0] - b[0], a[1] - b[1], p))

    formula_size = (p * p + 2 * p - 1) // 2
    return CarrierProfile(
        p=p,
        label=label,
        size=len(pts),
        formula_size=formula_size,
        multiplicity_hist=multiplicity_histogram(p, family),
        distance_count=len(distances),
        pinned_min=min(pinned),
        pinned_max=max(pinned),
        pinned_full_points=sum(1 for x in pinned if x == p),
        unit_edges=unit_edges,
        unit_direction_count=len(unit_dirs),
        unit_degree_hist=Counter(unit_degree.values()),
        paley_distance_hist=paley_distance_hist,
    )


def exact_kakeya_min(p: int) -> tuple[int, int, int]:
    """Exact branch-and-bound over one affine line in every direction."""
    dirs = directions(p)
    masks = [[line_mask(p, direction, b) for b in range(p)] for direction in dirs]
    best = (p * p + 2 * p - 1) // 2
    best_count = 0
    nodes = 0
    total = len(dirs)

    def rec(i: int, mask: int) -> None:
        nonlocal best, best_count, nodes
        nodes += 1
        current = mask.bit_count()
        if current > best:
            return
        if i == total:
            if current < best:
                best = current
                best_count = 1
            elif current == best:
                best_count += 1
            return

        # A future line can share at most one point with each already chosen
        # direction.  This is deliberately simple, but prunes p=7 comfortably.
        lower = current + sum(max(0, p - (i + j)) for j in range(total - i))
        if lower > best:
            return

        options = sorted(masks[i], key=lambda lm: ((mask | lm).bit_count(), lm))
        for line in options:
            rec(i + 1, mask | line)

    rec(0, 0)
    return best, best_count, nodes


def exhaustive_line_family_distribution(p: int) -> tuple[Counter[int], Counter[int], Counter[int], Counter[tuple[int, int, int, int]]]:
    line_sets = [
        [line_points(p, direction, b) for b in range(p)]
        for direction in directions(p)
    ]
    union_hist: Counter[int] = Counter()
    distance_hist: Counter[int] = Counter()
    pinned_min_hist: Counter[int] = Counter()
    joint_hist: Counter[tuple[int, int, int, int]] = Counter()
    for choices in product(range(p), repeat=p + 1):
        pts: set[Point] = set()
        for i, b in enumerate(choices):
            pts |= set(line_sets[i][b])
        distance_set = {quadratic_distance(a, b, p) for a in pts for b in pts}
        pinned = [
            len({quadratic_distance(a, b, p) for b in pts})
            for a in pts
        ]
        key = (len(pts), len(distance_set), min(pinned), max(pinned))
        union_hist[len(pts)] += 1
        distance_hist[len(distance_set)] += 1
        pinned_min_hist[min(pinned)] += 1
        joint_hist[key] += 1
    return union_hist, distance_hist, pinned_min_hist, joint_hist


def random_line_family_distribution(p: int, trials: int) -> Counter[tuple[int, int, int]]:
    rng = Random(659 + p)
    hist: Counter[tuple[int, int, int]] = Counter()
    for _ in range(trials):
        family = [(direction, rng.randrange(p)) for direction in directions(p)]
        pts = points_from_mask(p, union_mask_for_family(p, family))
        distance_set = {quadratic_distance(a, b, p) for a in pts for b in pts}
        pinned_min = min(
            len({quadratic_distance(a, b, p) for b in pts})
            for a in pts
        )
        hist[(len(pts), len(distance_set), pinned_min)] += 1
    return hist


def directed_3cycles(adj: list[list[int]]) -> int:
    total = 0
    for i, j, k in combinations(range(len(adj)), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (
            adj[i][k] and adj[k][j] and adj[j][i]
        ):
            total += 1
    return total


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach(start: int, reverse: bool = False) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for u in range(n):
                edge = adj[u][v] if reverse else adj[v][u]
                if edge and u not in seen:
                    seen.add(u)
                    stack.append(u)
        return seen

    remaining = set(range(n))
    sizes: list[int] = []
    while remaining:
        v = next(iter(remaining))
        comp = reach(v) & reach(v, reverse=True)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            if not dp[mask][v]:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[-1])


def majority_tournament(vertices: list[str], rankings: list[list[str]]) -> list[list[int]]:
    positions = [{v: i for i, v in enumerate(ranking)} for ranking in rankings]
    n = len(vertices)
    adj = [[0] * n for _ in range(n)]
    for i, u in enumerate(vertices):
        for j, v in enumerate(vertices):
            if i == j:
                continue
            wins = sum(1 for pos in positions if pos[u] < pos[v])
            if wins > len(rankings) / 2 or (
                wins == len(rankings) / 2 and positions[0][u] < positions[0][v]
            ):
                adj[i][j] = 1
    return adj


def print_counter(counter: Counter[object]) -> str:
    return "{" + ", ".join(f"{k}: {v}" for k, v in sorted(counter.items(), key=lambda kv: str(kv[0]))) + "}"


def main() -> None:
    print("S659 finite-field Kakeya/Falconer carrier ledger")
    print("=" * 78)
    print()

    print("A. External theorem anchors")
    print("-" * 78)
    print("Dvir 2009: finite-field Kakeya set = contains a line in every direction;")
    print("  polynomial method gives |K| >= C_n q^n in F_q^n.")
    print("Iosevich-Rudnev 2005: finite-field Falconer distance sets use Fourier")
    print("  machinery and Kloosterman sums to bound |Delta(E)| from |E|.")
    print("Hart-Iosevich-Koh-Senger-Uriarte-Tuero 2008: finite-field distance")
    print("  graphs expose the unit-distance and pseudorandom side channels.")

    print("\nB. Exact small Kakeya line-family minima")
    print("-" * 78)
    for p in [3, 5, 7]:
        best, count, nodes = exact_kakeya_min(p)
        formula = (p * p + 2 * p - 1) // 2
        print(
            f"p={p}: exact_min={best}, parabola_formula={formula}, "
            f"min_family_count={count}, branch_nodes={nodes}"
        )
    print("For p=11,13 below we use the same parabola carrier; exact search is not")
    print("needed for the repo transfer, which is side-channel retention.")

    print("\nC. Parabola Kakeya carriers as Falconer/unit-distance packets")
    print("-" * 78)
    profiles = [
        profile_family(p, "parabola", parabola_line_family(p))
        for p in [3, 5, 7, 11, 13]
    ]
    for prof in profiles:
        print(
            f"p={prof.p}: |K|={prof.size} formula={prof.formula_size} "
            f"density={prof.size}/{prof.p ** 2}, distances={prof.distance_count}/{prof.p}"
        )
        print(
            f"  multiplicity_hist={print_counter(prof.multiplicity_hist)}, "
            f"pinned_min/max/full={prof.pinned_min}/{prof.pinned_max}/{prof.pinned_full_points}"
        )
        print(
            f"  unit_edges={prof.unit_edges}, unit_direction_count={prof.unit_direction_count}, "
            f"unit_degree_hist={print_counter(prof.unit_degree_hist)}"
        )
        print(f"  paley_distance_hist={print_counter(prof.paley_distance_hist)}")

    print("\nD. Exhaustive scalar-leak audit over all p=5 line families")
    print("-" * 78)
    union_hist, distance_hist, pinned_min_hist, joint_hist = exhaustive_line_family_distribution(5)
    print(f"union_size_hist={print_counter(union_hist)}")
    print(f"distance_count_hist={print_counter(distance_hist)}")
    print(f"pinned_min_hist={print_counter(pinned_min_hist)}")
    print("top joint states (union_size, distance_count, pinned_min, pinned_max):")
    for key, value in joint_hist.most_common(8):
        print(f"  {key}: {value}")
    print("Interpretation: over F_5^2 every one-line-per-direction family already")
    print("has full distance support, but pinned fibers and union/concurrency still vary.")

    print("\nE. Random larger-prime line-family audit")
    print("-" * 78)
    for p in [7, 11]:
        hist = random_line_family_distribution(p, 2000)
        print(f"p={p}, trials=2000, top states (union_size, distance_count, pinned_min):")
        for key, value in hist.most_common(8):
            print(f"  {key}: {value}")

    print("\nF. Tournament Analysis over proof carriers")
    print("-" * 78)
    vertices = [
        "direction_cover_kakeya",
        "concurrency_multiplicity_packet",
        "falconer_distance_support",
        "pinned_distance_fibers",
        "unit_distance_direction_spine",
        "paley_character_orientation",
        "polynomial_vanishing_certificate",
        "raw_cardinality_density",
        "LRC_C27_shell_owner_lift",
    ]
    rankings = [
        [
            "LRC_C27_shell_owner_lift",
            "concurrency_multiplicity_packet",
            "pinned_distance_fibers",
            "unit_distance_direction_spine",
            "direction_cover_kakeya",
            "falconer_distance_support",
            "polynomial_vanishing_certificate",
            "paley_character_orientation",
            "raw_cardinality_density",
        ],
        [
            "polynomial_vanishing_certificate",
            "direction_cover_kakeya",
            "concurrency_multiplicity_packet",
            "raw_cardinality_density",
            "pinned_distance_fibers",
            "falconer_distance_support",
            "paley_character_orientation",
            "unit_distance_direction_spine",
            "LRC_C27_shell_owner_lift",
        ],
        [
            "falconer_distance_support",
            "pinned_distance_fibers",
            "paley_character_orientation",
            "unit_distance_direction_spine",
            "direction_cover_kakeya",
            "concurrency_multiplicity_packet",
            "polynomial_vanishing_certificate",
            "LRC_C27_shell_owner_lift",
            "raw_cardinality_density",
        ],
        [
            "unit_distance_direction_spine",
            "direction_cover_kakeya",
            "pinned_distance_fibers",
            "paley_character_orientation",
            "falconer_distance_support",
            "concurrency_multiplicity_packet",
            "LRC_C27_shell_owner_lift",
            "polynomial_vanishing_certificate",
            "raw_cardinality_density",
        ],
        [
            "LRC_C27_shell_owner_lift",
            "concurrency_multiplicity_packet",
            "pinned_distance_fibers",
            "unit_distance_direction_spine",
            "direction_cover_kakeya",
            "falconer_distance_support",
            "polynomial_vanishing_certificate",
            "paley_character_orientation",
            "raw_cardinality_density",
        ],
    ]
    adj = majority_tournament(vertices, rankings)
    score_hist = Counter(sum(row) for row in adj)
    print(f"vertices={vertices}")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={directed_3cycles(adj)}")
    print(f"scc_sizes={scc_sizes(adj)}")
    print(f"hamiltonian_paths={hamiltonian_paths(adj)}")
    print("ranking:")
    for score, vertex in sorted(((sum(adj[i]), v) for i, v in enumerate(vertices)), reverse=True):
        print(f"  {score}: {vertex}")

    print("\nG. Repo synthesis")
    print("-" * 78)
    print("1. Kakeya direction coverage is the finite-field analogue of an LRC clock")
    print("   cover: every projective direction is hit, but cardinality alone forgets")
    print("   which concurrency point owns the overlap.")
    print("2. Falconer distance support saturates too early in the small fields tested.")
    print("   The stronger transferable carrier is pinned distances plus unit-direction")
    print("   spine, parallel to S655's odd wall pairs and C=27 shell labels.")
    print("3. The LRC14 move is a Kakeya-style no-leak lemma: after the scalar wall is")
    print("   hit, retain owner/concurrency labels.  In finite fields this is line")
    print("   multiplicity; in LRC14 it is Res_27 gcd shell plus carry/owner.")


if __name__ == "__main__":
    main()
