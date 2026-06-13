#!/usr/bin/env python3
"""S649: Euler/Heegner prime windows as a tournament boundary carrier.

For the Euler-Rabinowitsch quadratics

    f_p(x) = x^2 + x + p,

the familiar lucky primes are

    p in {2, 3, 5, 11, 17, 41}.

They are the p-values with d=4p-1 in the Heegner class-number-one list
{1,2,3,7,11,19,43,67,163}, apart from the small d=1,2,3 exceptions.

This script records the exact finite boundary:

    f_p(x) is prime for x=0..p-2,
    f_p(p-1) = p^2.

The tournament-facing observation is the indexing split:

    zero-based prime run length = p-1 = fixed Hamiltonian spine length;
    interior prime run length   = p-2 = Moon strong-tournament c3 floor.

The script does not claim that tournaments prove the Heegner theorem.  It
names a carrier dictionary and tests the finite scout around it.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations


EULER_LUCKY_PRIMES = (2, 3, 5, 11, 17, 41)
HEEGNER_NUMBERS = (1, 2, 3, 7, 11, 19, 43, 67, 163)


def choose2(n: int) -> int:
    return n * (n - 1) // 2


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    d = 5
    step = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += step
        step = 6 - step
    return True


def primes_upto(limit: int) -> list[int]:
    return [n for n in range(2, limit + 1) if is_prime(n)]


def f(p: int, x: int) -> int:
    return x * x + x + p


def first_composite_index(p: int) -> tuple[int, int]:
    """Return (x, f_p(x)) for the first composite value."""
    x = 0
    while is_prime(f(p, x)):
        x += 1
    return x, f(p, x)


def lucky_prime_search(limit: int) -> list[int]:
    """Search p where f_p is prime exactly through x=p-2."""
    out: list[int] = []
    for p in primes_upto(limit):
        fail_x, fail_value = first_composite_index(p)
        if fail_x == p - 1 and fail_value == p * p:
            out.append(p)
    return out


def heegner_projection_rows() -> list[tuple[int, str, str]]:
    rows = []
    for d in HEEGNER_NUMBERS:
        if (d + 1) % 4:
            rows.append((d, "-", "outside d=4p-1 Euler-family shape"))
            continue
        p = (d + 1) // 4
        if is_prime(p):
            rows.append((d, str(p), "Euler lucky prime"))
        else:
            rows.append((d, str(p), "degenerate p, not prime"))
    return rows


def tournament_boundary_rows() -> list[dict[str, int | str | bool]]:
    rows: list[dict[str, int | str | bool]] = []
    for p in EULER_LUCKY_PRIMES:
        d = 4 * p - 1
        fail_x, fail_value = first_composite_index(p)
        rows.append(
            {
                "p": p,
                "d": d,
                "heegner": d in HEEGNER_NUMBERS,
                "prime_run_zero_based": fail_x,
                "prime_inputs": f"0..{fail_x - 1}",
                "interior_prime_inputs": max(p - 2, 0),
                "first_failure_x": fail_x,
                "first_failure_value": fail_value,
                "spine_edges": p - 1,
                "moon_c3_floor": p - 2 if p >= 3 else 0,
                "off_path_fiber": choose2(p - 1),
            }
        )
    return rows


@dataclass(frozen=True)
class Lens:
    name: str
    arithmetic_exactness: int
    tournament_use: int
    side_channel_retained: int
    proof_transfer_value: int
    scalar_risk: int


LENSES = (
    Lens("norm_line_d_equals_4p_minus_1", 5, 3, 5, 5, 1),
    Lens("endpoint_square_failure_p2", 5, 4, 4, 5, 1),
    Lens("spine_length_p_minus_1", 4, 5, 4, 4, 1),
    Lens("interior_length_p_minus_2_moon_floor", 4, 5, 4, 4, 1),
    Lens("heegner_class_number_one_side_channel", 5, 3, 5, 5, 1),
    Lens("off_path_fiber_choose_p_minus_1_2", 3, 5, 4, 3, 2),
    Lens("input_positions_as_vertices", 3, 4, 3, 3, 2),
    Lens("raw_lucky_prime_list", 2, 1, 1, 1, 5),
)


def lens_score(lens: Lens) -> int:
    return (
        3 * lens.arithmetic_exactness
        + 3 * lens.tournament_use
        + 3 * lens.side_channel_retained
        + 2 * lens.proof_transfer_value
        - 4 * lens.scalar_risk
    )


def majority_tournament(lenses: tuple[Lens, ...]) -> dict[tuple[str, str], str]:
    arcs: dict[tuple[str, str], str] = {}
    for a, b in combinations(lenses, 2):
        sa = lens_score(a)
        sb = lens_score(b)
        if sa > sb or (sa == sb and a.name < b.name):
            arcs[(a.name, b.name)] = a.name
        else:
            arcs[(a.name, b.name)] = b.name
    return arcs


def tournament_fingerprints(lenses: tuple[Lens, ...]) -> dict[str, object]:
    arcs = majority_tournament(lenses)
    out_neighbors: dict[str, set[str]] = {lens.name: set() for lens in lenses}
    for (a, b), winner in arcs.items():
        loser = b if winner == a else a
        out_neighbors[winner].add(loser)

    score_hist = Counter(len(out_neighbors[lens.name]) for lens in lenses)
    c3 = 0
    for a, b, c in combinations([lens.name for lens in lenses], 3):
        ab = b in out_neighbors[a]
        bc = c in out_neighbors[b]
        ca = a in out_neighbors[c]
        ba = a in out_neighbors[b]
        cb = b in out_neighbors[c]
        ac = c in out_neighbors[a]
        if (ab and bc and ca) or (ba and cb and ac):
            c3 += 1

    # Count Hamiltonian paths by DP over this lens tournament.
    names = [lens.name for lens in lenses]
    idx = {name: i for i, name in enumerate(names)}
    adj = [0] * len(names)
    for name, outs in out_neighbors.items():
        i = idx[name]
        for out in outs:
            adj[i] |= 1 << idx[out]
    dp = [[0] * len(names) for _ in range(1 << len(names))]
    for i in range(len(names)):
        dp[1 << i][i] = 1
    for mask in range(1 << len(names)):
        for last in range(len(names)):
            ways = dp[mask][last]
            if not ways:
                continue
            avail = adj[last] & ~mask
            while avail:
                bit = avail & -avail
                nxt = bit.bit_length() - 1
                avail -= bit
                dp[mask | bit][nxt] += ways
    h_paths = sum(dp[-1])
    ranking = sorted(lenses, key=lambda lens: (-lens_score(lens), lens.name))
    return {
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3_cycles": c3,
        "hamiltonian_paths": h_paths,
        "ranking": [lens.name for lens in ranking],
    }


def print_table(headers: list[str], rows: list[list[object]]) -> None:
    widths = [len(header) for header in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(str(cell)))
    print(" | ".join(header.rjust(widths[i]) for i, header in enumerate(headers)))
    print("-+-".join("-" * width for width in widths))
    for row in rows:
        print(" | ".join(str(cell).rjust(widths[i]) for i, cell in enumerate(row)))


def main() -> None:
    print("Euler/Heegner prime-window boundary scout (S649)")
    print()
    print("Exact identity:")
    print("  f_p(x) = x^2 + x + p = ((2x+1)^2 + d)/4 with d=4p-1")
    print("  f_p(p-1) = p^2, so x=p-1 is the forced square endpoint")
    print()

    search_limit = 500
    found = lucky_prime_search(search_limit)
    print(f"Lucky-prime search p<={search_limit}: {found}")
    print(f"Matches expected set: {tuple(found) == EULER_LUCKY_PRIMES}")
    print()

    print("Heegner projection through d=4p-1")
    print_table(
        ["d", "p=(d+1)/4", "role"],
        [[d, p, role] for d, p, role in heegner_projection_rows()],
    )
    print()

    print("Prime window and tournament boundary dictionary")
    rows = tournament_boundary_rows()
    print_table(
        [
            "p",
            "d",
            "prime x",
            "run",
            "interior",
            "fail x",
            "fail value",
            "spine",
            "Moon floor",
            "off-path fiber",
        ],
        [
            [
                row["p"],
                row["d"],
                row["prime_inputs"],
                row["prime_run_zero_based"],
                row["interior_prime_inputs"],
                row["first_failure_x"],
                row["first_failure_value"],
                row["spine_edges"],
                row["moon_c3_floor"],
                row["off_path_fiber"],
            ]
            for row in rows
        ],
    )
    print()

    print("Boundary reading")
    print("  zero-based run length p-1 == Hamiltonian spine length on p vertices")
    print("  interior inputs x=1..p-2 have count p-2 == Moon c3 floor for strong p-tournaments")
    print("  endpoint x=p-1 folds to p^2, a square sink rather than another prime")
    print("  after fixing the spine, the tournament deformation budget is C(p-1,2)")
    print()

    print("Tournament Analysis")
    print("  challenged vertices: primes, Heegner discriminants, input positions, residues,")
    print("    Hamiltonian-path arcs, directed 3-cycles, off-path arcs, and proof obligations")
    print("  chosen vertices: carrier lenses; quotient preserves the p-1/p-2 boundary")
    print("    and destroys the class-number-one proof channel if scalarized")
    fp = tournament_fingerprints(LENSES)
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3_cycles={fp['directed_3_cycles']}")
    print(f"  hamiltonian_paths={fp['hamiltonian_paths']}")
    print("  majority ranking:")
    for i, name in enumerate(fp["ranking"], start=1):
        print(f"    {i}. {name}")
    print()

    print("Hypothesis generated")
    print("  HYP-2225: The lucky-prime window is a spine/interior/end-square carrier.")
    print("  The exact arithmetic is Heegner/Rabinowitsch; the repo-useful transfer is")
    print("  to keep both side channels: class-number-one norms and tournament spine/fiber labels.")


if __name__ == "__main__":
    main()
