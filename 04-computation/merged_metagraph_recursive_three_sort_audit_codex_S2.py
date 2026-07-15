#!/usr/bin/env python3
"""Exact cross-n audit of nodes, tilings, and blue/black complement-lines.

The merged metagraph is treated as a three-sorted incidence object:

    X_n  = fixed-Hamiltonian-path staircase tilings,
    M_n  = converse-orbits of tournament isomorphism classes,
    L_n  = orbits of X_n under all-tile complementation.

There is a quotient q_n:X_n->M_n and a two-to-one quotient X_n->L_n.
Each line [t] in L_n has the unordered node endpoints
{q_n(t),q_n(complement(t))}; these endpoints may coincide.  The line is
blue exactly when t is fixed by anti-diagonal grid reflection, and black
otherwise.

Deleting either endpoint of the distinguished Hamiltonian path gives maps
d_top,d_bottom:X_n->X_{n-1}.  They commute with complementation and therefore
descend to L_n, but they do not descend to functions M_n->M_{n-1}.  The exact
node recursion is a weighted correspondence through X_n.  This script checks
that distinction, line multiplicities, self-loops, and the full endpoint-
coupled line transport tensor in the committed exact atlases for 3 <= n <= 7.
The tensor quotients only by simultaneous endpoint swap; an independently
sorted endpoint-pair tensor is retained separately as a strict marginal.

The symmetric two-end (Mode-B) restriction is audited as a second recursion.
It tracks the quotient-code reflection defect, inherited versus fresh
blackness, all four loop transitions, and the same-class/converse-sheet
holonomy of merged loops.

The general blue/black transition law is also checked.  With

    m_n = binom(n-1,2),  h_n = floor((n-1)^2/4),
    B_n = 2^(h_n-1),     K_n = 2^(m_n-1)-B_n,

the number of high-n lines of each colour that top-delete to each low colour
is

    BB = 2^(n-3),
    BK = B_n-BB,
    KB = 2^(n-2) B_(n-1)-BB,
    KK = 2^(n-2) K_(n-1)-BK.

The BB masks are exactly the masks constant on every difference diagonal
x-y=d, 2 <= d <= n-1.  This supplies the elementary proof: high and low grid
reflection together generate unit translation along each difference diagonal.

Tournament Analysis is deliberately not imposed on this audit.  The native
pair relation is a symmetric involutive multigraph relation with multiplicity
and self-loops; turning the three carrier sorts into a binary tournament would
erase the data being tested.  The challenged vertex-set assumption is instead
made explicit: vertices here may be tiling presentations, quotient nodes, or
line instances, and the maps between those sorts are part of the object.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SMALL = ROOT / "05-knowledge/results/tournament_tiling_metagraph_address_codex_S4.json"
DEFAULT_N7 = ROOT / "05-knowledge/results/tournament_tiling_metagraph_address_n7_codex_S4.json"
DEFAULT_OUTPUT = ROOT / "05-knowledge/results/merged_metagraph_recursive_three_sort_audit_codex_S2.out"


def tile_schema(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1))


def reflection_permutation(n: int) -> tuple[int, ...]:
    tiles = tile_schema(n)
    index = {tile: bit for bit, tile in enumerate(tiles)}
    return tuple(index[(n - y + 1, n - x + 1)] for x, y in tiles)


def permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    return sum(((mask >> source) & 1) << target for source, target in enumerate(permutation))


def is_grid_symmetric(mask: int, n: int) -> bool:
    return mask == permute_mask(mask, reflection_permutation(n))


def reflection_defect(mask: int, n: int) -> int:
    """The complement-invariant line defect mask (1+reflection)mask."""
    return mask ^ permute_mask(mask, reflection_permutation(n))


def is_difference_striped(mask: int, n: int) -> bool:
    values: dict[int, int] = {}
    for bit, (x, y) in enumerate(tile_schema(n)):
        value = (mask >> bit) & 1
        diagonal = x - y
        if diagonal in values and values[diagonal] != value:
            return False
        values[diagonal] = value
    return True


def endpoint_delete(mask: int, n: int, side: str) -> int:
    """Delete first/top n or last/bottom 1 from the fixed path."""
    target_index = {tile: bit for bit, tile in enumerate(tile_schema(n - 1))}
    answer = 0
    for source_bit, (x, y) in enumerate(tile_schema(n)):
        if side == "top":
            if x == n:
                continue
            target_tile = (x, y)
        elif side == "bottom":
            if y == 1:
                continue
            target_tile = (x - 1, y - 1)
        else:
            raise ValueError(f"unknown side {side!r}")
        answer |= ((mask >> source_bit) & 1) << target_index[target_tile]
    return answer


def load_atlases(small_path: Path, n7_path: Path) -> dict[int, dict]:
    small = json.loads(small_path.read_text())
    n7 = json.loads(n7_path.read_text())
    assert small["schema_version"] == 1
    assert n7["schema_version"] == 1

    atlases: dict[int, dict] = {}
    for record in small["sizes"]:
        q = [None] * record["tilings"]
        class_codes = [None] * record["tilings"]
        for tiling in record["tiling_map"]:
            q[tiling["mask"]] = tiling["node_rank"]
            class_codes[tiling["mask"]] = int(tiling["class_code"], 16)
        assert all(value is not None for value in q)
        assert all(value is not None for value in class_codes)
        atlases[record["n"]] = {"record": record, "q": q, "class_codes": class_codes}

    atlases[7] = {
        "record": n7,
        "q": n7["node_rank_by_mask"],
        "class_codes": n7["class_code_by_mask"],
    }
    assert sorted(atlases) == [3, 4, 5, 6, 7]
    return atlases


def line_census(q: list[int], class_codes: list[int], n: int) -> dict:
    m = len(tile_schema(n))
    full = (1 << m) - 1
    instances: Counter[tuple[str, bool]] = Counter()
    support: set[tuple[int, int]] = set()
    colored_support: set[tuple[str, int, int]] = set()
    loop_holonomy: Counter[tuple[str, str]] = Counter()
    records = []
    for mask in range(1 << m):
        other = mask ^ full
        if mask > other:
            continue
        u, v = sorted((q[mask], q[other]))
        color = "B" if is_grid_symmetric(mask, n) else "K"
        instances[(color, u == v)] += 1
        if u == v:
            sheet = "same" if class_codes[mask] == class_codes[other] else "swap"
            loop_holonomy[(color, sheet)] += 1
        support.add((u, v))
        colored_support.add((color, u, v))
        records.append((mask, other, u, v, color))
    return {
        "instances": instances,
        "support": support,
        "colored_support": colored_support,
        "loop_holonomy": loop_holonomy,
        "records": records,
    }


def color_formula(n: int) -> Counter[tuple[str, str]]:
    m = (n - 1) * (n - 2) // 2
    low_m = (n - 2) * (n - 3) // 2
    h = (n - 1) ** 2 // 4
    low_h = (n - 2) ** 2 // 4
    blue = 1 << (h - 1)
    black = (1 << (m - 1)) - blue
    low_blue = 1 << (low_h - 1)
    low_black = (1 << (low_m - 1)) - low_blue if low_m else 0
    bb = 1 << (n - 3)
    bk = blue - bb
    kb = (1 << (n - 2)) * low_blue - bb
    kk = (1 << (n - 2)) * low_black - bk
    assert black == kb + kk
    return Counter({("B", "B"): bb, ("B", "K"): bk, ("K", "B"): kb, ("K", "K"): kk})


def sorted_counter(counter: Counter) -> str:
    return "{" + ", ".join(f"{key}:{counter[key]}" for key in sorted(counter, key=repr)) + "}"


def format_fraction(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def mode_b_delete(mask: int, n: int) -> int:
    """Delete both distinguished path endpoints, retaining the interior."""
    return endpoint_delete(endpoint_delete(mask, n, "top"), n - 1, "bottom")


def line_sheet_status(q: list[int], class_codes: list[int], mask: int, other: int) -> str:
    """E for a cross-edge, S/C for same/converse-sheet merged loops."""
    if q[mask] != q[other]:
        return "E"
    return "S" if class_codes[mask] == class_codes[other] else "C"


def analyze(atlases: dict[int, dict]) -> str:
    lines = [
        "MERGED METAGRAPH RECURSIVE THREE-SORT AUDIT",
        "exact committed atlases; all counts include line multiplicity and self-loops",
        "",
        "SIZE CENSUS",
        "n m tilings classes merged fibre[min,max] B[cross,loop] K[cross,loop] "
        "support[colored,plain,loops] loop_hol[Bsame,Bswap,Ksame,Kswap]",
    ]

    censuses: dict[int, dict] = {}
    for n in sorted(atlases):
        record = atlases[n]["record"]
        q = list(map(int, atlases[n]["q"]))
        class_codes = list(map(int, atlases[n]["class_codes"]))
        m = len(tile_schema(n))
        assert len(q) == 1 << m == record["tilings"]
        assert len(set(q)) == record["merged_nodes"]

        reflection = reflection_permutation(n)
        full = (1 << m) - 1
        assert all(q[mask] == q[permute_mask(mask, reflection)] for mask in range(1 << m))
        assert all(
            is_grid_symmetric(mask, n) == is_grid_symmetric(mask ^ full, n)
            for mask in range(1 << m)
        )

        fibres = Counter(q)
        census = line_census(q, class_codes, n)
        censuses[n] = census
        instances = census["instances"]
        blue = instances[("B", False)] + instances[("B", True)]
        black = instances[("K", False)] + instances[("K", True)]
        assert blue == record["blue_lines"]
        assert black == record["black_lines"]

        node_sc = {int(node["rank"]): bool(node["self_converse"]) for node in record["nodes"]}
        assert set(node_sc) == set(q)
        black_pair_multiplicity: Counter[tuple[int, int]] = Counter()
        blue_half_degree: Counter[int] = Counter()
        black_half_degree: Counter[int] = Counter()
        for mask, other, u, v, color in census["records"]:
            reflected_mask = permute_mask(mask, reflection)
            reflected_line = min(reflected_mask, reflected_mask ^ full)
            assert (reflected_line == mask) == (color == "B")
            if color == "K":
                assert min(
                    permute_mask(reflected_line, reflection),
                    permute_mask(reflected_line, reflection) ^ full,
                ) == mask
                black_pair_multiplicity[(u, v)] += 1
            degree = blue_half_degree if color == "B" else black_half_degree
            if u == v:
                degree[u] += 2
            else:
                degree[u] += 1
                degree[v] += 1
        assert all(multiplicity % 2 == 0 for multiplicity in black_pair_multiplicity.values())
        assert all(black_half_degree[node] % 2 == 0 for node in node_sc)
        assert all(
            (blue_half_degree[node] % 2 == 1) if self_converse else (blue_half_degree[node] == 0)
            for node, self_converse in node_sc.items()
        )

        loop_support = sum(u == v for u, v in census["support"])
        holonomy = census["loop_holonomy"]
        assert holonomy[("B", "swap")] == 0
        lines.append(
            f"{n} {m} {1 << m} {record['classes']} {record['merged_nodes']} "
            f"[{min(fibres.values())},{max(fibres.values())}] "
            f"B[{instances[('B', False)]},{instances[('B', True)]}] "
            f"K[{instances[('K', False)]},{instances[('K', True)]}] "
            f"[{len(census['colored_support'])},{len(census['support'])},{loop_support}] "
            f"[{holonomy[('B', 'same')]},{holonomy[('B', 'swap')]},"
            f"{holonomy[('K', 'same')]},{holonomy[('K', 'swap')]}]"
        )

    lines.extend(
        [
            "",
            "CROSS-n RECURSION",
            "n P_support parent-target histogram line-degree colors[BB,BK,KB,KK] "
            "tensor_support[coupled,projected] edge-images[high->low,low<-high] loop-transport",
        ]
    )

    for n in range(4, 8):
        high_q = list(map(int, atlases[n]["q"]))
        low_q = list(map(int, atlases[n - 1]["q"]))
        m = len(tile_schema(n))
        low_m = len(tile_schema(n - 1))
        full = (1 << m) - 1
        low_full = (1 << low_m) - 1
        reflection = reflection_permutation(n)
        low_reflection = reflection_permutation(n - 1)
        lift_degree = 1 << (n - 2)

        top = [endpoint_delete(mask, n, "top") for mask in range(1 << m)]
        bottom = [endpoint_delete(mask, n, "bottom") for mask in range(1 << m)]

        # Cube identities: regular endpoint extensions, complement naturality,
        # and exchange of the two faces by grid reflection.
        assert set(Counter(top).values()) == {lift_degree}
        assert set(Counter(bottom).values()) == {lift_degree}
        for mask in range(1 << m):
            assert top[mask ^ full] == (top[mask] ^ low_full)
            assert bottom[mask ^ full] == (bottom[mask] ^ low_full)
            reflected = permute_mask(mask, reflection)
            assert top[reflected] == permute_mask(bottom[mask], low_reflection)
            if n >= 5:
                assert endpoint_delete(top[mask], n - 1, "bottom") == endpoint_delete(
                    bottom[mask], n - 1, "top"
                )

        transfers = {}
        for side, images in (("top", top), ("bottom", bottom)):
            transfer: Counter[tuple[int, int]] = Counter()
            for mask, high_node in enumerate(high_q):
                transfer[(high_node, low_q[images[mask]])] += 1
            transfers[side] = transfer
        assert transfers["top"] == transfers["bottom"]
        transfer = transfers["top"]

        high_fibres = Counter(high_q)
        low_fibres = Counter(low_q)
        assert all(sum(transfer[(u, v)] for v in low_fibres) == size for u, size in high_fibres.items())
        assert all(
            sum(transfer[(u, v)] for u in high_fibres) == lift_degree * size
            for v, size in low_fibres.items()
        )
        parent_targets: dict[int, set[int]] = defaultdict(set)
        for (u, v), count in transfer.items():
            assert count > 0
            parent_targets[u].add(v)
        parent_histogram = Counter(len(targets) for targets in parent_targets.values())

        color_transition: Counter[tuple[str, str]] = Counter()
        projected_low_lines: Counter[tuple[int, int, str]] = Counter()
        coupled_tensor_top: Counter[tuple[tuple[int, int, int, int], str, str]] = Counter()
        coupled_tensor_bottom: Counter[tuple[tuple[int, int, int, int], str, str]] = Counter()
        projected_tensor_top: Counter[tuple[int, int, int, int, str, str]] = Counter()
        projected_tensor_bottom: Counter[tuple[int, int, int, int, str, str]] = Counter()
        loop_transport: Counter[tuple[str, str, bool, bool]] = Counter()
        high_edge_images: dict[tuple[int, int], set[tuple[int, int]]] = defaultdict(set)
        low_edge_sources: dict[tuple[int, int], set[tuple[int, int]]] = defaultdict(set)
        blue_lifts_by_low_line: Counter[int] = Counter()

        for mask, other, high_u, high_v, high_color in censuses[n]["records"]:
            high_left, high_right = high_q[mask], high_q[other]
            low_mask = top[mask]
            low_other = top[other]
            assert low_other == (low_mask ^ low_full)
            low_left, low_right = low_q[low_mask], low_q[low_other]
            low_u, low_v = sorted((low_q[low_mask], low_q[low_other]))
            low_color = "B" if is_grid_symmetric(low_mask, n - 1) else "K"
            color_transition[(high_color, low_color)] += 1
            projected_low_lines[(low_u, low_v, low_color)] += 1
            coupled = min(
                (high_left, high_right, low_left, low_right),
                (high_right, high_left, low_right, low_left),
            )
            coupled_tensor_top[(coupled, high_color, low_color)] += 1
            projected_tensor_top[(high_u, high_v, low_u, low_v, high_color, low_color)] += 1
            loop_transport[(high_color, low_color, high_u == high_v, low_u == low_v)] += 1
            high_edge_images[(high_u, high_v)].add((low_u, low_v))
            low_edge_sources[(low_u, low_v)].add((high_u, high_v))
            if high_color == "B":
                blue_lifts_by_low_line[min(low_mask, low_other)] += 1

            bottom_mask = bottom[mask]
            bottom_other = bottom[other]
            bottom_left, bottom_right = low_q[bottom_mask], low_q[bottom_other]
            bottom_u, bottom_v = sorted((low_q[bottom_mask], low_q[bottom_other]))
            bottom_color = "B" if is_grid_symmetric(bottom_mask, n - 1) else "K"
            bottom_coupled = min(
                (high_left, high_right, bottom_left, bottom_right),
                (high_right, high_left, bottom_right, bottom_left),
            )
            coupled_tensor_bottom[(bottom_coupled, high_color, bottom_color)] += 1
            projected_tensor_bottom[
                (high_u, high_v, bottom_u, bottom_v, high_color, bottom_color)
            ] += 1

            # The BB intersection is exactly the difference-striped subcube.
            assert (
                high_color == "B" and low_color == "B"
            ) == is_difference_striped(mask, n)
            if high_color == "B":
                assert low_color == (
                    "B" if is_grid_symmetric(bottom_mask, n - 1) else "K"
                )

        assert coupled_tensor_top == coupled_tensor_bottom
        assert projected_tensor_top == projected_tensor_bottom
        projected_from_coupled: Counter[tuple[int, int, int, int, str, str]] = Counter()
        for ((high_left, high_right, low_left, low_right), high_color, low_color), count in (
            coupled_tensor_top.items()
        ):
            high_u, high_v = sorted((high_left, high_right))
            low_u, low_v = sorted((low_left, low_right))
            projected_from_coupled[(high_u, high_v, low_u, low_v, high_color, low_color)] += count
        assert projected_from_coupled == projected_tensor_top
        assert color_transition == color_formula(n)

        # Every low line has exactly 2^(n-2) high line lifts.  This is the
        # line-level commuting square, including its node endpoints and color.
        low_line_counts: Counter[tuple[int, int, str]] = Counter()
        for _, _, u, v, color in censuses[n - 1]["records"]:
            low_line_counts[(u, v, color)] += 1
        assert projected_low_lines == Counter(
            {key: lift_degree * count for key, count in low_line_counts.items()}
        )

        # A low blue line has blue lifts iff it is difference-striped, and then
        # exactly two.  Other low lines can also be images of high blue lines,
        # so blue is a transition label rather than a nested subobject.
        for mask, other, _, _, low_color in censuses[n - 1]["records"]:
            count = blue_lifts_by_low_line[min(mask, other)]
            if low_color == "B":
                assert count == (2 if is_difference_striped(mask, n - 1) else 0)

        formula = color_formula(n)
        colors = [formula[("B", "B")], formula[("B", "K")], formula[("K", "B")], formula[("K", "K")]]
        high_image_range = (min(map(len, high_edge_images.values())), max(map(len, high_edge_images.values())))
        low_source_range = (min(map(len, low_edge_sources.values())), max(map(len, low_edge_sources.values())))
        loop_compact = Counter()
        for (hc, lc, high_loop, low_loop), count in loop_transport.items():
            loop_compact[(hc + lc, "L" if high_loop else "E", "L" if low_loop else "E")] += count
        lines.append(
            f"{n} {len(transfer)} {sorted(parent_histogram.items())} {lift_degree} {colors} "
            f"[{len(coupled_tensor_top)},{len(projected_tensor_top)}] "
            f"[{high_image_range[0]}..{high_image_range[1]},"
            f"{low_source_range[0]}..{low_source_range[1]}] {sorted_counter(loop_compact)}"
        )

    lines.extend(
        [
            "",
            "MODE-B TWO-END LINE/DEFECT TOWER",
            "n->n-2 line-degree blue-lifts/blue-parent defects[high,low,kernel-dim] "
            "colors[BB,BK,KB,KK] loop[EE,EL,LE,LL] sheet-transport",
        ]
    )

    for n in range(5, 8):
        low_n = n - 2
        high_q = list(map(int, atlases[n]["q"]))
        low_q = list(map(int, atlases[low_n]["q"]))
        high_classes = list(map(int, atlases[n]["class_codes"]))
        low_classes = list(map(int, atlases[low_n]["class_codes"]))
        high_m = len(tile_schema(n))
        low_m = len(tile_schema(low_n))
        high_full = (1 << high_m) - 1
        low_full = (1 << low_m) - 1
        line_degree = 1 << (2 * n - 5)
        blue_degree = 1 << (n - 2)

        line_preimages: Counter[int] = Counter()
        blue_preimages: Counter[int] = Counter()
        colors: Counter[tuple[str, str]] = Counter()
        loop_status: Counter[tuple[bool, bool]] = Counter()
        sheet_transport: Counter[tuple[str, str, str]] = Counter()
        defect_images: dict[int, int] = {}

        for mask, other, high_u, high_v, high_color in censuses[n]["records"]:
            low_mask = mode_b_delete(mask, n)
            low_other = mode_b_delete(other, n)
            assert low_other == (low_mask ^ low_full)
            low_line = min(low_mask, low_other)
            low_u, low_v = sorted((low_q[low_mask], low_q[low_other]))
            low_color = "B" if is_grid_symmetric(low_mask, low_n) else "K"
            line_preimages[low_line] += 1
            if high_color == "B":
                blue_preimages[low_line] += 1
            colors[(high_color, low_color)] += 1
            high_loop = high_u == high_v
            low_loop = low_u == low_v
            loop_status[(high_loop, low_loop)] += 1
            high_sheet = line_sheet_status(high_q, high_classes, mask, other)
            low_sheet = line_sheet_status(low_q, low_classes, low_mask, low_other)
            sheet_transport[(high_color + low_color, high_sheet, low_sheet)] += 1

            high_defect = reflection_defect(mask, n)
            low_defect = reflection_defect(low_mask, low_n)
            assert mode_b_delete(high_defect, n) == low_defect
            previous = defect_images.setdefault(high_defect, low_defect)
            assert previous == low_defect

        assert set(line_preimages.values()) == {line_degree}
        for mask, other, _, _, low_color in censuses[low_n]["records"]:
            line = min(mask, other)
            expected = blue_degree if low_color == "B" else 0
            assert blue_preimages[line] == expected

        high_defect_fibres = Counter(
            reflection_defect(mask, n) for mask, _, _, _, _ in censuses[n]["records"]
        )
        low_defect_fibres = Counter(
            reflection_defect(mask, low_n)
            for mask, _, _, _, _ in censuses[low_n]["records"]
        )
        high_f = (n - 1) // 2
        low_f = (low_n - 1) // 2
        high_h = (high_m + high_f) // 2
        low_h = (low_m + low_f) // 2
        high_defect_dim = (high_m - high_f) // 2
        low_defect_dim = (low_m - low_f) // 2
        assert len(high_defect_fibres) == 1 << high_defect_dim
        assert len(low_defect_fibres) == 1 << low_defect_dim
        assert set(high_defect_fibres.values()) == {1 << (high_h - 1)}
        assert set(low_defect_fibres.values()) == {1 << (low_h - 1)}
        assert high_defect_dim - low_defect_dim == n - 3
        defect_image_fibres = Counter(defect_images.values())
        assert set(defect_image_fibres) == set(low_defect_fibres)
        assert set(defect_image_fibres.values()) == {1 << (n - 3)}

        high_blue = 1 << (high_h - 1)
        low_blue = 1 << (low_h - 1)
        high_total = 1 << (high_m - 1)
        low_total = 1 << (low_m - 1)
        expected_colors = Counter(
            {
                ("B", "B"): high_blue,
                ("B", "K"): 0,
                ("K", "B"): (line_degree - blue_degree) * low_blue,
                ("K", "K"): line_degree * (low_total - low_blue),
            }
        )
        assert colors == +expected_colors

        compact_loops = [
            loop_status[(False, False)],
            loop_status[(False, True)],
            loop_status[(True, False)],
            loop_status[(True, True)],
        ]
        compact_colors = [
            colors[("B", "B")],
            colors[("B", "K")],
            colors[("K", "B")],
            colors[("K", "K")],
        ]
        lines.append(
            f"{n}->{low_n} {line_degree} {blue_degree} "
            f"[{high_defect_dim},{low_defect_dim},{n - 3}] {compact_colors} "
            f"{compact_loops} {sorted_counter(sheet_transport)}"
        )

    lines.extend(
        [
            "",
            "NODE-QUOTIENT COMPOSITION TEST",
            "n unequal[entries,rows/total] max_error witness=(high,low,actual,Markov-product)",
        ]
    )

    for n in range(5, 8):
        high_q = list(map(int, atlases[n]["q"]))
        middle_q = list(map(int, atlases[n - 1]["q"]))
        low_q = list(map(int, atlases[n - 2]["q"]))
        high_fibres = Counter(high_q)
        middle_fibres = Counter(middle_q)

        first: Counter[tuple[int, int]] = Counter()
        second: Counter[tuple[int, int]] = Counter()
        direct: Counter[tuple[int, int]] = Counter()
        for mask, high_node in enumerate(high_q):
            middle_mask = endpoint_delete(mask, n, "top")
            low_mask = endpoint_delete(middle_mask, n - 1, "top")
            first[(high_node, middle_q[middle_mask])] += 1
            direct[(high_node, low_q[low_mask])] += 1
        for mask, middle_node in enumerate(middle_q):
            second[(middle_node, low_q[endpoint_delete(mask, n - 1, "top")])] += 1

        bad_entries = 0
        bad_rows: set[int] = set()
        maximum_error = Fraction(0)
        witness = None
        for high_node, high_size in high_fibres.items():
            for low_node in set(low_q):
                actual = Fraction(direct[(high_node, low_node)], high_size)
                markov_product = sum(
                    Fraction(first[(high_node, middle_node)], high_size)
                    * Fraction(second[(middle_node, low_node)], middle_size)
                    for middle_node, middle_size in middle_fibres.items()
                )
                if actual == markov_product:
                    continue
                bad_entries += 1
                bad_rows.add(high_node)
                error = abs(actual - markov_product)
                if error > maximum_error:
                    maximum_error = error
                    witness = (high_node, low_node, actual, markov_product)
        assert witness is not None
        high_node, low_node, actual, markov_product = witness
        witness_text = (
            f"({high_node},{low_node},{format_fraction(actual)},{format_fraction(markov_product)})"
        )
        lines.append(
            f"{n} [{bad_entries},{len(bad_rows)}/{len(high_fibres)}] "
            f"{format_fraction(maximum_error)} {witness_text}"
        )

    lines.extend(
        [
            "",
            "GENERAL IDENTITIES CHECKED FOR n=4..7",
            "d_top C_n = C_(n-1) d_top; d_bottom C_n = C_(n-1) d_bottom",
            "d_top S_n = S_(n-1) d_bottom",
            "d_top d_bottom = d_bottom d_top (where both composites are defined)",
            "P_top = P_bottom after converse merging",
            "each lower complement-line has exactly 2^(n-2) high line lifts",
            "BB masks = difference-striped masks; dim = n-2; BB lines = 2^(n-3)",
            "reflection-line orbits have size 1 exactly for blue and size 2 exactly for black",
            "every projected black pair has even multiplicity, including black self-loops",
            "black half-edge degree is even at every node (loops count twice)",
            "blue half-edge degree is odd exactly at self-converse nodes and zero otherwise",
            "merged loops retain class-sheet holonomy: same class or converse-sheet swap",
            "reflection defect gives 0 -> blue lines -> all lines -> defect space -> 0",
            "Mode-B defect fibres have dimension n-3; inherited/fresh blackness is distinguished",
            "Mode-B loop and class-sheet status are transported but are not hereditary",
            "",
            "PRESERVATION VERDICT",
            "tilings: honest endpoint-deletion maps and 2^(n-2)-fold extension fibres",
            "lines: honest regular set maps; multiplicity, color transition, and loop status retained by tensor",
            "nodes: no deletion function in general; retain weighted span M_n <- X_n -> M_(n-1)",
            "node span: not Markov-compositional; retain the reached lower tiling, not only its node",
            "blue/black: not preserved as a color; retain the two-time label (chi_n,chi_(n-1))",
            "simple node support: insufficient; it forgets parallel lines, self-loop mass, and cross-n coupling",
            "",
            "ALL ASSERTIONS PASSED",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-atlas", type=Path, default=DEFAULT_SMALL)
    parser.add_argument("--n7-atlas", type=Path, default=DEFAULT_N7)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--stdout", action="store_true")
    args = parser.parse_args()

    report = analyze(load_atlases(args.small_atlas, args.n7_atlas))
    args.output.write_text(report)
    if args.stdout:
        print(report, end="")


if __name__ == "__main__":
    main()
