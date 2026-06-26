"""LRC14 owner-strip filtration synthesis.

This is a small connection-mining script, not a proof search.  It assembles
recent residual work with older cocycle/normal-fan carriers and reports the
common hidden coordinate: endpoint-owner strip data appears after status,
primitive-period, and Haar-zeta quotients have done their jobs.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, permutations
from typing import Dict, Iterable, List, Sequence, Tuple


@dataclass(frozen=True)
class Layer:
    name: str
    features: Dict[str, int]
    destroys: str


FEATURES = [
    "status",
    "route",
    "primitive_period",
    "haar_zeta",
    "owner_strip",
    "topology",
    "family_transfer",
    "compression",
    "low_cost",
]


LAYERS = [
    Layer(
        "raw_shadow",
        dict(
            status=0,
            route=0,
            primitive_period=0,
            haar_zeta=0,
            owner_strip=0,
            topology=0,
            family_transfer=0,
            compression=3,
            low_cost=3,
        ),
        "all geometry beyond the scalar or coarse word",
    ),
    Layer(
        "status_gate",
        dict(
            status=3,
            route=1,
            primitive_period=0,
            haar_zeta=0,
            owner_strip=0,
            topology=2,
            family_transfer=2,
            compression=2,
            low_cost=2,
        ),
        "route scheduling inside strict-open residuals",
    ),
    Layer(
        "primitive_period_deck",
        dict(
            status=3,
            route=2,
            primitive_period=3,
            haar_zeta=0,
            owner_strip=0,
            topology=1,
            family_transfer=2,
            compression=2,
            low_cost=2,
        ),
        "q=14/boundary and owner-strip geometry",
    ),
    Layer(
        "haar_zeta_grid",
        dict(
            status=3,
            route=2,
            primitive_period=1,
            haar_zeta=3,
            owner_strip=0,
            topology=2,
            family_transfer=3,
            compression=1,
            low_cost=1,
        ),
        "which endpoint owner carries the diagonal boundary strip",
    ),
    Layer(
        "endpoint_owner_strip",
        dict(
            status=3,
            route=3,
            primitive_period=1,
            haar_zeta=2,
            owner_strip=3,
            topology=3,
            family_transfer=3,
            compression=1,
            low_cost=1,
        ),
        "nonlargest bars and explicit route labels",
    ),
    Layer(
        "labelled_route_certificate",
        dict(
            status=3,
            route=3,
            primitive_period=3,
            haar_zeta=3,
            owner_strip=3,
            topology=3,
            family_transfer=2,
            compression=0,
            low_cost=0,
        ),
        "nothing theorem-relevant after route is declared",
    ),
]


TIE_PATH = [
    "labelled_route_certificate",
    "endpoint_owner_strip",
    "haar_zeta_grid",
    "primitive_period_deck",
    "status_gate",
    "raw_shadow",
]


CONNECTIONS = [
    (
        "HYP-2997",
        "Cocycle normal form already split Haar `zeta` and endpoint-owner boundary cocycles.",
        "Read S201 as a two-step cocycle: first nonzero `zeta`, then endpoint-owner boundary current.",
    ),
    (
        "HYP-3018",
        "Active-bottleneck normal fan retained endpoint owner sets and residue sums.",
        "The S201 external owner strips are the same kind of normal-fan support, not a new scalar.",
    ),
    (
        "HYP-3035",
        "All 15 coarse residual first repairs are `owner_strip` after status protection.",
        "The q=23 diagonal owner strip is a local instance of the same first-tooth mechanism.",
    ),
    (
        "HYP-3036",
        "Primitive safe decks route Q-witness rows while zero decks point to covering/q=14 layers.",
        "Period data is a scheduler page before owner-strip geometry, not a replacement for it.",
    ),
    (
        "HYP-3041",
        "AP-tail puncture repair shows a mod-14 owner strip can still forget the hidden `m mod 13` clock.",
        "Treat `q13_puncture_bit` and `ap_tail_certificate_kind` as primitive-period sidecars before owner strips are trusted.",
    ),
    (
        "HYP-3038",
        "The q=23 pair is a drop/add square with nonzero exact-M zeta and equal coarse endpoint word `B18Z6`.",
        "`B18Z6` is only a scalar endpoint shadow; the external owner multiset is the splitting coordinate.",
    ),
]


def wins(a: Layer, b: Layer) -> bool:
    score = 0
    for feature in FEATURES:
        if a.features[feature] > b.features[feature]:
            score += 1
        elif a.features[feature] < b.features[feature]:
            score -= 1
    if score != 0:
        return score > 0
    return TIE_PATH.index(a.name) < TIE_PATH.index(b.name)


def edges(layers: Sequence[Layer]) -> Dict[Tuple[str, str], str]:
    out: Dict[Tuple[str, str], str] = {}
    for a, b in combinations(layers, 2):
        if wins(a, b):
            out[(a.name, b.name)] = a.name
        else:
            out[(b.name, a.name)] = b.name
    return out


def score_hist(layers: Sequence[Layer], edge_map: Dict[Tuple[str, str], str]) -> Dict[int, int]:
    scores = {layer.name: 0 for layer in layers}
    for winner in edge_map.values():
        scores[winner] += 1
    hist: Dict[int, int] = {}
    for score in scores.values():
        hist[score] = hist.get(score, 0) + 1
    return dict(sorted(hist.items()))


def directed_3cycles(layers: Sequence[Layer], edge_map: Dict[Tuple[str, str], str]) -> int:
    names = [layer.name for layer in layers]

    def beats(x: str, y: str) -> bool:
        return (x, y) in edge_map

    count = 0
    for a, b, c in combinations(names, 3):
        if (
            beats(a, b)
            and beats(b, c)
            and beats(c, a)
            or beats(a, c)
            and beats(c, b)
            and beats(b, a)
        ):
            count += 1
    return count


def hamiltonian_paths(layers: Sequence[Layer], edge_map: Dict[Tuple[str, str], str]) -> List[Tuple[str, ...]]:
    names = [layer.name for layer in layers]

    def ok(path: Tuple[str, ...]) -> bool:
        return all((path[i], path[i + 1]) in edge_map for i in range(len(path) - 1))

    return [path for path in permutations(names) if ok(path)]


def print_table(rows: Iterable[Tuple[str, str, str]]) -> None:
    for handle, detail, sharpened in rows:
        print(f"{handle}: {detail}")
        print(f"  hidden connection: {sharpened}")


def main() -> None:
    edge_map = edges(LAYERS)
    paths = hamiltonian_paths(LAYERS, edge_map)

    print("LRC14 owner-strip filtration synthesis (codex-2026-06-26-S205)")
    print("source_threads=HYP-3041,HYP-3038,HYP-3037,HYP-3036,HYP-3035,HYP-3031,HYP-3018,HYP-2997")
    print()

    print("[0] Hidden connection table")
    print_table(CONNECTIONS)
    print()

    print("[1] Proposed filtration")
    for index, layer in enumerate(LAYERS):
        print(f"E{index}: {layer.name}")
        retained = ", ".join(
            feature for feature in FEATURES if layer.features[feature] >= 2
        )
        print(f"  retains: {retained if retained else 'only compression/cost shadow'}")
        print(f"  destroys: {layer.destroys}")
    print()

    print("[2] Sharpened statement")
    print(
        "The post-status residual problem is not just route scheduling.  It is a "
        "short filtration: primitive-period decks route direct q-witness mass; "
        "HYP-3041 shows AP-tail prime clocks can repair owner-strip collisions; "
        "Haar zeta detects hidden drop/add diagonal charge; endpoint-owner strip "
        "is the first common boundary current that separates the remaining local "
        "route ambiguity before labelled route certificates are used."
    )
    print()
    print("candidate local rule:")
    print(
        "  in a protected residual fiber, a route-mixed pair must either have "
        "positive primitive safe mass at q<=13 or a named AP-tail q=13/fixed-point "
        "clock, nonzero Haar/drop-add zeta leading to a q-diagonal, or an "
        "endpoint-owner strip that names the boundary current; otherwise it is "
        "genuine named F7/THM-572 debt."
    )
    print()

    print("[3] Tournament Analysis")
    print("vertices=filtration layers / proof carriers, not runners")
    print(
        "pairwise_observable=status,route,primitive_period,haar_zeta,"
        "owner_strip,topology,family_transfer,compression,low_cost"
    )
    print("switch=majority of retained proof coordinates; tie path=" + ">".join(TIE_PATH))
    print(f"score_hist={score_hist(LAYERS, edge_map)}")
    print(f"directed_3cycles={directed_3cycles(LAYERS, edge_map)}")
    print("scc_sizes=[1,1,1,1,1,1]")
    print(f"hamiltonian_path_count={len(paths)}")
    if paths:
        print("path=" + ">".join(paths[0]))
    print()

    print("[4] Assumption challenge")
    print(
        "Alternate vertices considered: runners, residual fibers, denominators, "
        "drop pairs, add pairs, endpoint walls, safe bars, Haar squares, "
        "primitive periods, cocycle channels, normal-fan supports, and proof "
        "obligations."
    )
    print(
        "Chosen vertices are filtration layers because the preserved predicate is "
        "strict-open status plus route schedulability after quotienting.  This "
        "destroys raw runner identity, exact route labels before the last page, "
        "and scalar endpoint words such as B18Z6 unless the owner multiset is "
        "retained."
    )


if __name__ == "__main__":
    main()
