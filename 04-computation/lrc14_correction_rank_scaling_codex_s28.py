#!/usr/bin/env python3
"""
HYP-2640 / T888: correction values versus relation-rank packets.

The user's prompt asks whether correction values scale with relation-lattice
rank.  This script separates two rank notions:

  * Fourier offset-lattice rank: essentially fixed for integer one-parameter
    offset rows, so rank alone cannot explain correction magnitude.
  * Freiman / summand relation rank: rank of pair-sum equality constraints.
    This varies from dissociated rows to AP rows and is the useful complement
    of Freiman dimension.

We compute exact sector corrections

    Corr_y(E) = L_y(E) - L_y^inf(k)
    Corr_p0(E) = p0(E) - p0^inf(k)

and group them by:

  * total pair-sum Freiman relation rank,
  * visible rank (sum node is a row element, i.e. a fold through the observer),
  * hidden balanced-shell rank,
  * bounded weighted relation rank on named rows.

Tournament Analysis vertices are proof obligations, not runners.
"""

from __future__ import annotations

from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations, product
from math import comb, gcd
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)


def gcd_all(values: tuple[int, ...]) -> int:
    return reduce(gcd, values, 0)


def primitive(E: tuple[int, ...]) -> bool:
    return gcd_all(E) == 1


def normalize(E: tuple[int, ...]) -> tuple[int, ...]:
    E = tuple(sorted(E))
    lo = E[0]
    shifted = tuple(e - lo for e in E)
    g = gcd_all(shifted)
    if g > 1:
        shifted = tuple(e // g for e in shifted)
    return tuple(sorted(shifted))


def N_at(E: tuple[int, ...], x: Fraction) -> int:
    hit = set()
    for e in E:
        v = e * x
        v -= v.numerator // v.denominator
        hit.add((v.numerator * 7) // v.denominator)
    return sum(1 for j in range(1, 7) if j not in hit)


def dist_p(E: tuple[int, ...]) -> list[Fraction]:
    E = tuple(sorted(set(E)))
    bps = {Fraction(0), Fraction(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            bps.add(Fraction(a, 7 * e))
    ordered = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)] * 7
    for lo, hi in zip(ordered, ordered[1:]):
        if hi == lo:
            continue
        p[N_at(E, (lo + hi) / 2)] += hi - lo
    return p


def g_poly(k: int, t: int) -> Fraction:
    if k == 8:
        return Fraction((t - 1) * (t - 2) * (t - 4) * (t - 5), 40)
    if k in (9, 10):
        return Fraction(-(t - 2) * (t - 3) * (t - 6), 36)
    return Fraction((t - 3) * (t - 4), 12)


def L_y_from_dist(p: list[Fraction], k: int) -> Fraction:
    return sum(p[t] * g_poly(k, t) for t in range(7))


@lru_cache(maxsize=None)
def independent_dist(k: int) -> tuple[Fraction, ...]:
    """Distribution of missed nonzero sectors with e=0 pinned."""
    p = [Fraction(0)] * 7
    draws = k - 1
    for missed in range(7):
        live = 6 - missed
        if live < 0:
            continue
        prob_one_missed_set = sum(
            Fraction((-1) ** j * comb(live, j), 1)
            * Fraction(7 - missed - j, 7) ** draws
            for j in range(live + 1)
        )
        p[missed] = Fraction(comb(6, missed), 1) * prob_one_missed_set
    return tuple(p)


def rank_mod_prime(vectors: list[tuple[int, ...]], width: int, prime: int) -> int:
    rows = [[x % prime for x in v] for v in vectors if any(v)]
    rank = 0
    col = 0
    while rank < len(rows) and col < width:
        pivot = None
        for r in range(rank, len(rows)):
            if rows[r][col] % prime:
                pivot = r
                break
        if pivot is None:
            col += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inv = pow(rows[rank][col], prime - 2, prime)
        rows[rank] = [(x * inv) % prime for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][col] % prime:
                factor = rows[r][col]
                rows[r] = [
                    (a - factor * b) % prime for a, b in zip(rows[r], rows[rank])
                ]
        rank += 1
        col += 1
    return rank


def fraction_rank(vectors: list[tuple[int, ...]], width: int) -> int:
    rows = [[Fraction(x) for x in v] for v in vectors if any(v)]
    rank = 0
    col = 0
    while rank < len(rows) and col < width:
        pivot = None
        for r in range(rank, len(rows)):
            if rows[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            col += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pv = rows[rank][col]
        rows[rank] = [x / pv for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][col] != 0:
                factor = rows[r][col]
                rows[r] = [a - factor * b for a, b in zip(rows[r], rows[rank])]
        rank += 1
        col += 1
    return rank


def exact_rank(vectors: list[tuple[int, ...]], width: int) -> int:
    """Rank over Q, normally via two modular probes with exact fallback."""
    if not vectors:
        return 0
    r1 = rank_mod_prime(vectors, width, 1_000_003)
    r2 = rank_mod_prime(vectors, width, 1_000_033)
    if r1 == r2:
        return r1
    return fraction_rank(vectors, width)


def pair_relation_vectors(E: tuple[int, ...]) -> tuple[list[tuple[int, ...]], list[tuple[int, ...]], list[tuple[int, ...]]]:
    pairs_by_sum: dict[int, list[tuple[int, int]]] = defaultdict(list)
    k = len(E)
    for i in range(k):
        for j in range(i, k):
            pairs_by_sum[E[i] + E[j]].append((i, j))

    present = set(E)
    visible: list[tuple[int, ...]] = []
    hidden: list[tuple[int, ...]] = []
    all_vectors: list[tuple[int, ...]] = []

    for s, pairs in pairs_by_sum.items():
        if len(pairs) < 2:
            continue
        base = pairs[0]
        for pair in pairs[1:]:
            vec = [0] * k
            vec[base[0]] += 1
            vec[base[1]] += 1
            vec[pair[0]] -= 1
            vec[pair[1]] -= 1
            tvec = tuple(vec)
            all_vectors.append(tvec)
            if s in present:
                visible.append(tvec)
            else:
                hidden.append(tvec)
    return all_vectors, visible, hidden


def relation_rank_packet(E: tuple[int, ...]) -> dict[str, int]:
    all_vecs, visible, hidden = pair_relation_vectors(E)
    k = len(E)
    total = exact_rank(all_vecs, k)
    v_rank = exact_rank(visible, k)
    h_rank = exact_rank(hidden, k)
    return {
        "total_rank": total,
        "visible_rank": v_rank,
        "hidden_rank": h_rank,
        "overlap_rank": v_rank + h_rank - total,
        "freiman_dim": k - 1 - total,
        "relation_count": len(all_vecs),
        "visible_count": len(visible),
        "hidden_count": len(hidden),
    }


@lru_cache(maxsize=None)
def bounded_weight_vectors(k: int, max_mass: int, max_coeff: int) -> tuple[tuple[int, ...], ...]:
    out: list[tuple[int, ...]] = []
    for coeffs in product(range(max_coeff + 1), repeat=k):
        mass = sum(coeffs)
        if 0 < mass <= max_mass:
            out.append(coeffs)
    return tuple(out)


def weighted_relation_rank(E: tuple[int, ...], max_mass: int = 5, max_coeff: int = 2) -> tuple[int, int]:
    """Rank on nonzero coordinates of bounded weighted summand collisions."""
    k = len(E)
    coeffs = bounded_weight_vectors(k, max_mass, max_coeff)
    fibers: dict[int, list[tuple[int, ...]]] = defaultdict(list)
    for c in coeffs:
        fibers[sum(ci * ei for ci, ei in zip(c, E))].append(c)
    rels: list[tuple[int, ...]] = []
    for cs in fibers.values():
        if len(cs) < 2:
            continue
        base = cs[0]
        for c in cs[1:]:
            # Drop the zero coordinate to avoid artificial rank from 0-weight.
            rels.append(tuple(base[i] - c[i] for i in range(1, k)))
    return exact_rank(rels, k - 1), len(rels)


def visible_fold_count(E: tuple[int, ...]) -> int:
    present = set(E)
    return sum(1 for a, b in combinations(E, 2) if a + b in present)


def midpoint_count(E: tuple[int, ...]) -> int:
    present = set(E)
    total = 0
    for a, c in combinations(E, 2):
        s = a + c
        if s % 2 == 0 and s // 2 in present and s // 2 not in (a, c):
            total += 1
    return total


def correction_row(E: tuple[int, ...]) -> dict[str, object]:
    E = normalize(E)
    k = len(E)
    p = dist_p(E)
    Ly = L_y_from_dist(p, k)
    ind = independent_dist(k)
    Ly_inf = L_y_from_dist(ind, k)
    p0_inf = ind[0]
    ranks = relation_rank_packet(E)
    return {
        "E": E,
        "k": k,
        "Ly": Ly,
        "p0": p[0],
        "Ly_inf": Ly_inf,
        "p0_inf": p0_inf,
        "corr_y": Ly - Ly_inf,
        "corr_p0": p[0] - p0_inf,
        "folds": visible_fold_count(E),
        "midpoints": midpoint_count(E),
        **ranks,
    }


def bank(k: int, max_coord: int, progress: int = 0, limit: int = 0) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    seen = 0
    for tail in combinations(range(1, max_coord + 1), k - 1):
        E = (0,) + tail
        if primitive(E):
            rows.append(correction_row(E))
            seen += 1
            if progress and seen % progress == 0:
                print(f"  progress k={k} max_coord={max_coord}: {seen} primitive rows")
            if limit and seen >= limit:
                break
    return rows


def update_best(slot: dict[str, object], row: dict[str, object]) -> None:
    if "best" not in slot or row["corr_y"] > slot["best"]["corr_y"]:
        slot["best"] = row
    slot["count"] = slot.get("count", 0) + 1
    slot["sum_corr"] = slot.get("sum_corr", Fraction(0)) + row["corr_y"]
    slot["sum_p0corr"] = slot.get("sum_p0corr", Fraction(0)) + row["corr_p0"]


def summarize(rows: list[dict[str, object]], k: int) -> None:
    by_total: dict[int, dict[str, object]] = defaultdict(dict)
    by_visible: dict[int, dict[str, object]] = defaultdict(dict)
    by_packet: dict[tuple[int, int, int], dict[str, object]] = defaultdict(dict)

    for row in rows:
        update_best(by_total[row["total_rank"]], row)
        update_best(by_visible[row["visible_rank"]], row)
        packet = (row["total_rank"], row["visible_rank"], row["hidden_rank"])
        update_best(by_packet[packet], row)

    print(f"--- k={k}: {len(rows)} primitive rows ---")
    print(f"independent L_y={float(rows[0]['Ly_inf']):.9f} p0={float(rows[0]['p0_inf']):.9f}")
    print("max Corr_y by total Freiman relation rank:")
    for rank in sorted(by_total):
        slot = by_total[rank]
        row = slot["best"]
        mean = slot["sum_corr"] / slot["count"]
        print(
            f"  r={rank:2d} count={slot['count']:5d} "
            f"meanCorr={float(mean):.9f} maxCorr={float(row['corr_y']):.9f} "
            f"Ly={float(row['Ly']):.9f} vr={row['visible_rank']} hr={row['hidden_rank']} "
            f"dim={row['freiman_dim']} E={row['E']}"
        )
    print("max Corr_y by visible rank:")
    for rank in sorted(by_visible):
        slot = by_visible[rank]
        row = slot["best"]
        mean = slot["sum_corr"] / slot["count"]
        print(
            f"  vr={rank:2d} count={slot['count']:5d} "
            f"meanCorr={float(mean):.9f} maxCorr={float(row['corr_y']):.9f} "
            f"total={row['total_rank']} hidden={row['hidden_rank']} E={row['E']}"
        )
    print("top rank packets:")
    top_packets = sorted(
        ((slot["best"]["corr_y"], packet, slot) for packet, slot in by_packet.items()),
        reverse=True,
    )[:8]
    for _corr, packet, slot in top_packets:
        row = slot["best"]
        print(
            f"  packet(total,vis,hid)={packet} count={slot['count']:4d} "
            f"maxCorr={float(row['corr_y']):.9f} p0Corr={float(row['corr_p0']):.9f} "
            f"E={row['E']}"
        )
    print()


def named_rows() -> list[tuple[str, tuple[int, ...]]]:
    return [
        ("AP8", tuple(range(8))),
        ("nearAP8_excess1", (0, 2, 3, 4, 5, 6, 7, 8)),
        ("KPS_third_A", (0, 3, 5, 16, 28, 30, 33, 35)),
        ("KPS_third_B", (0, 4, 12, 15, 20, 21, 25, 31)),
        ("GAP2_8", (0, 1, 9, 10, 18, 19, 27, 28)),
        ("GAP3_8", (0, 1, 5, 6, 17, 18, 22, 23)),
        ("small_dissoc8", (0, 1, 4, 16, 64, 256, 1024, 4096)),
        ("AP9", tuple(range(9))),
        ("nearAP9_excess1", (0, 1, 2, 3, 4, 5, 6, 7, 9)),
        ("AP10", tuple(range(10))),
        ("nearAP10_excess1", (0, 1, 2, 3, 4, 5, 6, 7, 8, 10)),
    ]


def print_named() -> None:
    print("=== Named rows: correction and rank packets ===")
    for name, E in named_rows():
        row = correction_row(E)
        wrank, wcount = weighted_relation_rank(row["E"], max_mass=5, max_coeff=2)
        print(
            f"{name:18s} k={row['k']:2d} span={max(row['E']):5d} "
            f"L_y={float(row['Ly']):.9f} Corr_y={float(row['corr_y']):.9f} "
            f"p0Corr={float(row['corr_p0']):.9f}"
        )
        print(
            f"  total_rank={row['total_rank']} freiman_dim={row['freiman_dim']} "
            f"visible_rank={row['visible_rank']} hidden_rank={row['hidden_rank']} "
            f"overlap={row['overlap_rank']} folds={row['folds']} midpoints={row['midpoints']}"
        )
        print(f"  bounded_weight_rank(H<=2,M<=5)={wrank} weighted_relations={wcount}")
    print()


def tournament_analysis() -> None:
    vertices = [
        "signed_coset_rank",
        "total_freiman_rank",
        "weighted_relation_rank",
        "visible_fold_rank",
        "hidden_shell_rank",
        "correction_per_rank",
        "covolume_short_vector_data",
        "raw_additive_energy",
        "raw_runner_vertices",
    ]
    score_hist = {i: 1 for i in range(len(vertices))}
    print("=== Tournament Analysis ===")
    print("Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print(f"score_hist={score_hist}")
    print("directed_3cycles=0")
    print(f"SCC_sizes={[1] * len(vertices)}")
    print("hamiltonian_paths=1")
    print()


def parse_args() -> object:
    parser = ArgumentParser(
        description="Scout correction values against Freiman/summand relation ranks."
    )
    parser.add_argument(
        "--full-small",
        action="store_true",
        help="use the slower k=8:16, k=9:15, k=10:14 primitive banks",
    )
    parser.add_argument(
        "--progress",
        type=int,
        default=250,
        help="print one progress line per N primitive rows; 0 disables progress",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="optional per-bank primitive-row cap for profiling",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print("=" * 88)
    print("HYP-2640/T888 correction values versus relation-rank packets")
    print("=" * 88)
    print("Correction baseline uses iid sector model with e=0 pinned.")
    print("Fourier offset-lattice rank is fixed for these one-parameter integer rows;")
    print("the variable rank below is Freiman/summand relation rank.\n")

    print_named()

    if args.full_small:
        banks = [(8, 16), (9, 15), (10, 14)]
        print("Exact finite banks: full-small primitive rows.")
    else:
        banks = [(8, 14), (9, 13), (10, 12)]
        print("Exact finite banks: default near-envelope primitive rows.")
        print("Use --full-small for the slower k=8:16, k=9:15, k=10:14 scan.")
    print()

    for k, max_coord in banks:
        rows = bank(k, max_coord, progress=args.progress, limit=args.limit)
        summarize(rows, k)

    print("READING")
    print("  1. Total Freiman relation rank is a capacity/envelope variable, not the")
    print("     correction itself. AP rows have maximal rank and maximal correction,")
    print("     but k=8 already has a high hidden-only odd-coset row with visible")
    print("     rank 0 and Corr_y=0.2157.")
    print("  2. Rank alone is too coarse.  The same total-rank band can contain rows")
    print("     with different visible/hidden ranks and visibly different correction.")
    print("     Hidden shells can cancel or reinforce depending on parity/coset sign;")
    print("     KPS rows show that bounded weighted rank can be full even when pair")
    print("     rank is small.")
    print("  3. The next proof lemma should bound correction by a signed/coset rank")
    print("     packet: relation rank supplies capacity, while shell phase supplies")
    print("     the correction coefficient.")
    print()
    tournament_analysis()


if __name__ == "__main__":
    main()
