#!/usr/bin/env python3
"""Finite-exact labelled structure of the endpoint-589 q=50 failure family.

The exactness of the family itself is inherited from the supplied complete
carrier-failure ledger.  Everything reported here is recomputed from that
ledger and the independently frozen q=50 active-carrier ledger.
"""

from __future__ import annotations

import argparse
import collections
import csv
import math
from pathlib import Path


POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
MASK30 = (1 << 30) - 1
EXPECTED_FAILURE_FNV = 0xFF421454F02D9099
EXPECTED_ACTIVE_FNV = 0xC075113890C7F5E1
HUB_BITS = (14, 24)  # labels 95 and 193
PETALS = frozenset((10, 15, 20, 30, 60, 63, 132, 170, 176, 190, 264, 290))
PETAL_MASK = sum(1 << bit for bit, label in enumerate(POOL) if label in PETALS)


def fail(message: str) -> "NoReturn":
    raise RuntimeError(message)


def require(condition: bool, message: str) -> None:
    if not condition:
        fail(message)


def fnv_add(state: int, word: int) -> int:
    for shift in range(0, 64, 8):
        state ^= (word >> shift) & 0xFF
        state = (state * 0x100000001B3) & ((1 << 64) - 1)
    return state


def mask_labels(mask: int) -> str:
    return "{" + ",".join(str(POOL[bit]) for bit in range(30)
                            if mask & (1 << bit)) + "}"


def read_failures(path: Path) -> list[int]:
    result: list[int] = []
    state = 0xCBF29CE484222325
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ["q", "r", "body_hex"],
                "failure header changed")
        for row in reader:
            q, r = int(row["q"]), int(row["r"])
            if (q, r) != (50, 589):
                continue
            body = int(row["body_hex"], 16)
            require(body <= MASK30 and body.bit_count() == 9,
                    "failure body escaped rank nine")
            result.append(body)
            state = fnv_add(state, body)
    require(len(result) == 20025 and len(set(result)) == len(result),
            "q50 failure count/distinctness changed")
    require(state == EXPECTED_FAILURE_FNV, "q50 failure FNV changed")
    return result


def read_active(path: Path) -> list[int]:
    result: list[int] = []
    state = 0xCBF29CE484222325
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames ==
                ["mask_hex", "rank", "joint", "margin_ticks"],
                "active header changed")
        for row in reader:
            mask = int(row["mask_hex"], 16)
            rank = int(row["rank"])
            require(mask <= MASK30 and mask.bit_count() == rank and
                    rank in (8, 9) and int(row["margin_ticks"]) > 0,
                    "active mask row changed")
            result.append(mask)
            state = fnv_add(state, mask)
    require(len(result) == 1398 and len(set(result)) == len(result),
            "active count/distinctness changed")
    require(state == EXPECTED_ACTIVE_FNV, "active FNV changed")
    rank_counts = collections.Counter(mask.bit_count() for mask in result)
    require(rank_counts == {8: 1347, 9: 51}, "active rank split changed")
    # The active deck is a clutter: no rank-nine member contains an active
    # rank-eight member.  Equal-rank members are distinct and incomparable.
    rank8 = set(mask for mask in result if mask.bit_count() == 8)
    for mask in result:
        if mask.bit_count() != 9:
            continue
        for bit in range(30):
            if mask & (1 << bit):
                require(mask ^ (1 << bit) not in rank8,
                        "active deck lost inclusion minimality")
    return result


def minimal_clutter(edges: set[int]) -> set[int]:
    ordered = sorted(edges, key=lambda mask: (mask.bit_count(), mask))
    minimal: list[int] = []
    for edge in ordered:
        if not any((prior & edge) == prior for prior in minimal):
            minimal.append(edge)
    return set(minimal)


def incidence_components(edges: set[int], universe_mask: int) -> list[int]:
    remaining = universe_mask
    components: list[int] = []
    while remaining:
        seed = remaining & -remaining
        component = seed
        changed = True
        while changed:
            changed = False
            for edge in edges:
                if edge & component and edge & ~component:
                    enlarged = component | edge
                    if enlarged != component:
                        component = enlarged
                        changed = True
        component &= universe_mask
        components.append(component)
        remaining &= ~component
    return sorted(components, key=lambda mask: (-mask.bit_count(), mask))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("failures", type=Path)
    parser.add_argument("active", type=Path)
    parser.add_argument("vertex_csv", type=Path)
    parser.add_argument("fibre_csv", type=Path)
    parser.add_argument("neither_csv", type=Path)
    args = parser.parse_args()

    bodies = read_failures(args.failures)
    active = read_active(args.active)
    body_set = set(bodies)

    support = 0
    common = MASK30
    for body in bodies:
        support |= body
        common &= body
        require(all(body & edge for edge in active),
                "listed failure is not an active-deck transversal")
    require(support == MASK30 and common == 0, "global support/common changed")

    body_degree = [sum(bool(body & (1 << bit)) for body in bodies)
                   for bit in range(30)]
    active_degree = [sum(bool(edge & (1 << bit)) for edge in active)
                     for bit in range(30)]
    require(len(set(body_degree)) == 30,
            "failure vertex degrees ceased to be distinct")
    top = tuple(sorted(range(30), key=lambda bit: body_degree[bit],
                       reverse=True)[:2])
    require(set(top) == set(HUB_BITS), "two dominant hubs changed")

    codegrees: list[tuple[int, int, int]] = []
    for left in range(30):
        for right in range(left + 1, 30):
            count = sum(bool(body & (1 << left)) and
                        bool(body & (1 << right))
                        for body in bodies)
            codegrees.append((count, left, right))
    codegrees.sort(reverse=True)
    require(codegrees[0] == (17454, 14, 24), "top codegree changed")

    private_hist: collections.Counter[int] = collections.Counter()
    private_examples: dict[int, int] = {}
    for body in bodies:
        private = 0
        for bit in range(30):
            singleton = 1 << bit
            if not body & singleton:
                continue
            if any((edge & body) == singleton for edge in active):
                private += 1
        private_hist[private] += 1
        private_examples.setdefault(private, body)
    require(private_hist == {4: 12, 5: 155, 6: 1315, 7: 5011,
                             8: 7881, 9: 5651},
            "private-witness histogram changed")
    petal_layers = collections.Counter((body & PETAL_MASK).bit_count()
                                       for body in bodies)
    require(petal_layers == {0: 870, 1: 4519, 2: 7803, 3: 5502,
                             4: 1293, 5: 38},
            "THM4234/4242 petal-layer split changed")

    with args.vertex_csv.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("bit", "label", "failure_degree", "active_degree",
                         "gcd_label_50", "gcd_label_589"))
        for bit, label in enumerate(POOL):
            writer.writerow((bit, label, body_degree[bit], active_degree[bit],
                             math.gcd(label, 50), math.gcd(label, 589)))

    hubs = sum(1 << bit for bit in HUB_BITS)
    outside = MASK30 ^ hubs
    fibre_rows: list[tuple[object, ...]] = []
    neither: list[int] = []
    for has95, has193 in ((0, 0), (0, 1), (1, 0), (1, 1)):
        selected = (has95 << 14) | (has193 << 24)
        fibre = [body for body in bodies if (body & hubs) == selected]
        residuals = [body & outside for body in fibre]
        target_rank = 9 - has95 - has193

        # Exact deletion/contraction quotient: active edges hit by selected
        # hubs disappear; both hub coordinates are deleted from every other
        # edge.  A fibre body is a transversal iff its residual is a
        # transversal of this reduced family.
        reduced = {edge & outside for edge in active if not edge & selected}
        require(0 not in reduced, "selected hubs alone hit a required edge")
        clutter = minimal_clutter(reduced)
        require(all(residual.bit_count() == target_rank and
                    all(residual & edge for edge in reduced)
                    for residual in residuals),
                "fibre quotient equivalence failed on supplied family")
        require(len(set(residuals)) == len(residuals),
                "fibre residual parameterization not injective")

        ranks = collections.Counter(edge.bit_count() for edge in reduced)
        clutter_ranks = collections.Counter(edge.bit_count() for edge in clutter)
        components = incidence_components(clutter, outside)
        fibre_common = outside
        fibre_support = 0
        for residual in residuals:
            fibre_common &= residual
            fibre_support |= residual
        fibre_rows.append((
            has95, has193, target_rank, len(fibre), len(reduced),
            ";".join(f"{rank}:{ranks[rank]}" for rank in sorted(ranks)),
            len(clutter),
            ";".join(f"{rank}:{clutter_ranks[rank]}"
                     for rank in sorted(clutter_ranks)),
            ";".join(str(component.bit_count()) for component in components),
            f"{fibre_common:08x}", f"{fibre_support:08x}",
        ))
        if not has95 and not has193:
            neither = fibre

    require([row[3] for row in fibre_rows] == [7, 1347, 1217, 17454],
            "hub-fibre counts changed")
    require(sum(int(row[3]) for row in fibre_rows) == 20025,
            "hub fibres do not partition family")
    require(len(neither) == 7, "neither-hub fibre changed")
    neither_common = MASK30
    neither_support = 0
    for body in neither:
        neither_common &= body
        neither_support |= body
    require(neither_common == 0x00100400 and neither_support == 0x26FC3707,
            "neither-hub common/support changed")

    with args.fibre_csv.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("has_95", "has_193", "residual_rank",
                         "transversal_count", "reduced_edges",
                         "reduced_rank_counts", "minimal_clutter_edges",
                         "minimal_clutter_rank_counts",
                         "clutter_component_sizes", "residual_common_hex",
                         "residual_support_hex"))
        writer.writerows(fibre_rows)

    with args.neither_csv.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("body_hex", "labels"))
        for body in neither:
            writer.writerow((f"{body:08x}", mask_labels(body)))

    # GCD coordinates are not sufficient: the two labels divisible by 19
    # already have radically different body degrees.
    bit95, bit190 = POOL.index(95), POOL.index(190)
    require(math.gcd(95, 589) == math.gcd(190, 589) == 19 and
            body_degree[bit95] == 18671 and body_degree[bit190] == 3465,
            "gcd hostile control changed")

    print("LRC14_ENDPOINT589_Q50_FAILURE_STRUCTURE_V1")
    print("FAMILY", len(bodies), "RANK 9 FAILURE_FNV",
          f"{EXPECTED_FAILURE_FNV:016x}", "SUPPORT", f"{support:08x}",
          "COMMON", f"{common:08x}")
    print("ACTIVE_CLUTTER", len(active), "RANK8 1347 RANK9 51 FNV",
          f"{EXPECTED_ACTIVE_FNV:016x}")
    print("HUBS bit14=95 degree=18671 bit24=193 degree=18801",
          "BOTH=17454 ONLY95=1217 ONLY193=1347 NEITHER=7")
    print("NEITHER_COMMON", f"{neither_common:08x}",
          mask_labels(neither_common), "SUPPORT", f"{neither_support:08x}")
    print("PRIVATE_LABEL_COUNTS",
          " ".join(f"{key}:{private_hist[key]}" for key in sorted(private_hist)))
    print("FIXED50_PETAL_LAYERS",
          " ".join(f"{key}:{petal_layers[key]}" for key in sorted(petal_layers)),
          "THM4234_ALREADY_COFINAL_K_LE_3", sum(petal_layers[k]
                                                for k in range(4)),
          "K_GE_4_REQUIRING_OTHER_MECHANISM",
          sum(value for key, value in petal_layers.items() if key >= 4))
    print("AUTOMORPHISM_CERTIFICATE 30_DISTINCT_VERTEX_DEGREES_IMPLY_"
          "TRIVIAL_LABEL_PERMUTATION_GROUP_AND_20025_SINGLETON_ORBITS")
    print("TOP_CODEGREES", " ".join(
        f"{POOL[left]}+{POOL[right]}:{count}"
        for count, left, right in codegrees[:10]))
    for row in fibre_rows:
        print("FIBRE", f"95={row[0]}", f"193={row[1]}",
              "RESIDUAL_RANK", row[2], "COUNT", row[3],
              "REDUCED", row[4], row[5], "CLUTTER", row[6], row[7],
              "COMPONENTS", row[8], "COMMON", row[9], "SUPPORT", row[10])
    print("GCD_HOSTILE gcd(95,589)=gcd(190,589)=19 BUT_DEGREES 18671 3465")
    print("SCOPE FINITE_EXACT_CONDITIONAL_ON_EXACT_INPUT_LEDGERS_"
          "NO_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
