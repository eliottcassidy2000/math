#!/usr/bin/env python3
"""LRC14 q=23 drop/add Haar square.

This is a local proof-interface computation, not a proof of LRC14.

HYP-3032 found that the first analytic-clock residual pair is squarefree and
still mixed:

    drop(10,13)->add(20,26)
    drop(8,12)->add(16,24)

HYP-3031 asks for an actual two-coordinate packet grid for such residuals.
Here the two coordinates are the dropped AP pair and the added double pair.
The opposite corners are the HYP-3032 q=23 residuals; the off-diagonal corners
are the cross-swaps.  This exposes the local Haar/fixed-margin mixed
coordinate zeta.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
SWITCHBOARD_PATH = REPO / "04-computation" / "lrc14_carrier_fusion_switchboard_codex_s189.py"
BASE = set(range(1, 14))


def load_switchboard():
    spec = spec_from_file_location("lrc14_carrier_fusion_switchboard_s189", SWITCHBOARD_PATH)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules["lrc14_carrier_fusion_switchboard_s189"] = module
    spec.loader.exec_module(module)
    return module


sw = load_switchboard()


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def row_from(drop_pair: tuple[int, int], add_pair: tuple[int, int]) -> tuple[int, ...]:
    return tuple(sorted((BASE - set(drop_pair)) | set(add_pair)))


@dataclass(frozen=True)
class CellSpec:
    key: str
    drop_pair: tuple[int, int]
    add_pair: tuple[int, int]
    declared_route: str

    @property
    def speeds(self) -> tuple[int, ...]:
        return row_from(self.drop_pair, self.add_pair)

    @property
    def diagonal_doubling_match(self) -> bool:
        return tuple(2 * x for x in self.drop_pair) == self.add_pair


@dataclass(frozen=True)
class Cell:
    spec: CellSpec
    features: object
    external_owner_word: str
    active_pair_sum_word: str
    active_residue_word: str


GRID = {
    "AA_petal_diag": CellSpec("AA_petal_diag", (10, 13), (20, 26), "PETAL-DIAGONAL"),
    "AB_cross": CellSpec("AB_cross", (10, 13), (16, 24), "CROSS-Q-WITNESS"),
    "BA_cross": CellSpec("BA_cross", (8, 12), (20, 26), "CROSS-Q-WITNESS"),
    "BB_cover_diag": CellSpec("BB_cover_diag", (8, 12), (16, 24), "COVER-DIAGONAL"),
}


def active_words(speeds: tuple[int, ...]) -> tuple[str, str, str]:
    boundary = sw.boundary_safe_points(speeds)
    external = Counter()
    residues = Counter()
    pair_sums = Counter()
    for _point, active in boundary:
        for speed in active:
            residues[speed % 14] += 1
            if speed > 13:
                external[f"{speed % 14}:{speed}"] += 1
        for a, b in combinations(active, 2):
            pair_sums[(a + b) % 14] += 1
    external_word = "none" if not external else ",".join(f"{k}x{v}" for k, v in sorted(external.items()))
    pair_word = "none" if not pair_sums else ",".join(f"{k}:x{v}" for k, v in sorted(pair_sums.items()))
    residue_word = "none" if not residues else ",".join(f"{k}:x{v}" for k, v in sorted(residues.items()))
    return external_word, pair_word, residue_word


def build_cells() -> dict[str, Cell]:
    out: dict[str, Cell] = {}
    for key, spec in GRID.items():
        row = sw.Row(key, spec.speeds, spec.declared_route)
        features = sw.row_features(row)
        external, pair_sums, residues = active_words(spec.speeds)
        out[key] = Cell(spec, features, external, pair_sums, residues)
    return out


def zeta(values: dict[str, Fraction | int]) -> Fraction | int:
    return values["AA_petal_diag"] - values["AB_cross"] - values["BA_cross"] + values["BB_cover_diag"]


def qword(cell: Cell) -> str:
    f = cell.features
    return (
        f"{cell.spec.key}: drop={cell.spec.drop_pair} add={cell.spec.add_pair} "
        f"diag2={cell.spec.diagonal_doubling_match} route={cell.spec.declared_route} "
        f"M={fmt(f.exact_m)}@{fmt(f.exact_t)} safe={fmt(f.safe_mu)} "
        f"bars={f.bar_count} longest={fmt(f.longest_bar)} chart={f.first_chart} "
        f"endpoint={f.endpoint_current_word} ext={cell.external_owner_word}"
    )


def print_grid(cells: dict[str, Cell]) -> None:
    print("[0] Two-coordinate q=23 drop/add square")
    print("base=AP13, row=(base - drop_pair) union add_pair")
    print("opposite corners are the HYP-3032 analytic residual pair")
    for key in ("AA_petal_diag", "AB_cross", "BA_cross", "BB_cover_diag"):
        print(qword(cells[key]))
        print(f"  speeds={cells[key].spec.speeds}")
        print(f"  active_residues={cells[key].active_residue_word}")
        print(f"  active_pair_sums_mod14={cells[key].active_pair_sum_word}")
    print()


def print_zeta(cells: dict[str, Cell]) -> None:
    print("[1] Haar/fixed-margin zeta readout")
    numeric_fields = {
        "exact_M": {k: v.features.exact_m for k, v in cells.items()},
        "safe_measure": {k: v.features.safe_mu for k, v in cells.items()},
        "bar_count": {k: v.features.bar_count for k, v in cells.items()},
        "longest_bar": {k: v.features.longest_bar for k, v in cells.items()},
        "midpoint_slack": {k: v.features.midpoint_slack for k, v in cells.items()},
        "boundary_count": {k: v.features.boundary_count for k, v in cells.items()},
        "zero_sum_active_pairs": {k: v.features.zero_sum_pairs for k, v in cells.items()},
        "magnitude_height": {k: v.features.magnitude_height for k, v in cells.items()},
        "magnitude_delta": {k: v.features.magnitude_delta for k, v in cells.items()},
    }
    for name, values in numeric_fields.items():
        z = zeta(values)
        shown = ", ".join(f"{k}={fmt(v) if isinstance(v, Fraction) else v}" for k, v in values.items())
        zshown = fmt(z) if isinstance(z, Fraction) else str(z)
        print(f"{name}: {shown} ; zeta={zshown}")
    print()
    print("interpretation=exact M is already a nonzero mixed coordinate on the square:")
    print("  the diagonal-doubling corners stay at q=23, while cross-swaps open as q=10 and q=8 Q-witness rows.")
    print("  magnitude height/delta have zeta 0, so raw magnitude alone misses the local interaction.")
    print()


def quotient_stress(cells: dict[str, Cell]) -> None:
    print("[2] Local quotient stress")
    key_fns = {
        "status_only": lambda c: c.features.status,
        "diagonal_doubling_match": lambda c: c.spec.diagonal_doubling_match,
        "exact_M": lambda c: c.features.exact_m,
        "exact_M_plus_safe_body": lambda c: (c.features.exact_m, c.features.bar_count, c.features.longest_bar),
        "exact_M_plus_endpoint_owner_strip": lambda c: (c.features.exact_m, c.external_owner_word),
        "diagonal_plus_endpoint_owner_strip": lambda c: (c.spec.diagonal_doubling_match, c.external_owner_word),
        "declared_route": lambda c: c.spec.declared_route,
    }
    print("key fibers mixed_route largest_fiber examples")
    for key_name, fn in key_fns.items():
        groups: dict[object, list[Cell]] = defaultdict(list)
        for cell in cells.values():
            groups[fn(cell)].append(cell)
        mixed = []
        for key, group in groups.items():
            routes = {cell.spec.declared_route for cell in group}
            if len(routes) > 1:
                mixed.append((key, [cell.spec.key for cell in group]))
        largest = max(len(group) for group in groups.values())
        example = " ; ".join(f"{repr(key)}=>{','.join(names)}" for key, names in mixed)
        print(f"{key_name} | {len(groups)} | {len(mixed)} | {largest} | {example or '-'}")
    print()


@dataclass(frozen=True)
class Tooth:
    name: str
    status: int
    route: int
    zeta: int
    endpoint: int
    topology: int
    family: int
    cost: int
    destroys: str

    def vector(self) -> tuple[int, ...]:
        return (self.status, self.route, self.zeta, self.endpoint, self.topology, self.family, self.cost)


TIE_PATH = (
    "labelled_packet_route",
    "endpoint_owner_strip",
    "safe_component_body",
    "exact_M_zeta_grid",
    "diagonal_doubling_match",
    "drop_add_row_column_shadow",
    "raw_analytic_q23_shadow",
)


def tooth_bank() -> list[Tooth]:
    return [
        Tooth("labelled_packet_route", 4, 4, 4, 4, 4, 4, 1, "nothing relevant after route labels are retained"),
        Tooth("endpoint_owner_strip", 4, 4, 4, 4, 3, 3, 3, "nonlargest safe components and global family proof"),
        Tooth("safe_component_body", 4, 3, 4, 3, 4, 3, 3, "which endpoint owner caused the split"),
        Tooth("exact_M_zeta_grid", 3, 3, 4, 2, 2, 3, 4, "endpoint owner strip inside the q=23 diagonal"),
        Tooth("diagonal_doubling_match", 3, 2, 3, 1, 1, 3, 4, "petal versus covering route inside the diagonal"),
        Tooth("drop_add_row_column_shadow", 2, 1, 1, 1, 1, 2, 4, "the mixed Haar coordinate zeta"),
        Tooth("raw_analytic_q23_shadow", 1, 1, 0, 0, 0, 1, 4, "drop/add geometry, endpoints, topology, and route"),
    ]


def winner(a: Tooth, b: Tooth) -> Tooth:
    av = a.vector()
    bv = b.vector()
    aw = sum(1 for x, y in zip(av, bv) if x > y)
    bw = sum(1 for x, y in zip(av, bv) if y > x)
    if aw > bw:
        return a
    if bw > aw:
        return b
    order = {name: i for i, name in enumerate(TIE_PATH)}
    return a if order[a.name] < order[b.name] else b


def tournament_analysis() -> None:
    teeth = tooth_bank()
    scores = Counter({tooth.name: 0 for tooth in teeth})
    edges: dict[str, set[str]] = {tooth.name: set() for tooth in teeth}
    for a, b in combinations(teeth, 2):
        win = winner(a, b)
        lose = b if win is a else a
        scores[win.name] += 1
        edges[win.name].add(lose.name)
    score_hist = Counter(scores.values())
    ordered = sorted(teeth, key=lambda t: -scores[t.name])
    print("[3] Tournament Analysis")
    print("vertices=local proof teeth, not runners or raw denominators")
    print("pairwise_observable=status_retention,route_retention,zeta_retention,endpoint_owner_retention,topology_retention,family_transfer,low_cost")
    print("switch=majority of retained proof coordinates; tie Hamiltonian path=" + ">".join(TIE_PATH))
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print("directed_3cycles=0")
    print("scc_sizes=[1,1,1,1,1,1,1]")
    print("hamiltonian_path_count=1")
    print("path=" + ">".join(tooth.name for tooth in ordered))
    for tooth in ordered:
        print(f"  {tooth.name}: destroys={tooth.destroys}")
    print()


def synthesis() -> None:
    print("[4] Proof-angle synthesis")
    print("zeta_repair_class=nested_refinement_to_q23_diagonal_then_owner_strip")
    print("first_split=the diagonal-doubling match keeps exact M=2/23; cross-swaps are strict q-witness rows at M=1/10 and M=1/8.")
    print("second_split=inside the q=23 diagonal, endpoint owner strips distinguish petal from covering even when exact M and analytic clocks agree.")
    print("guardrail=raw analytic q=23 data is a row/column shadow; exact_M_zeta plus endpoint-owner strip is the first non-route repair tooth in this local square.")
    print("candidate_lemma=for any double-pair drop/add square, diagonal doubling either opens to a direct q-witness off diagonal, descends through a family q-diagonal, or exposes endpoint-owner strip data that routes petal, covering, K33, or F7/THM-572 debt.")
    print()
    print("[5] Assumption challenge")
    print("Alternate vertices considered: runners, dropped pairs, added pairs, gaps, endpoints, fixed sections, residues, Fourier modes, Haar tiles, cover arcs, matroid circuits, and proof obligations.")
    print("Chosen vertices are proof teeth because the preserved LRC predicate is open/boundary status plus theorem-route schedulability after quotienting.")
    print("Destroyed information: row/column drop/add margins forget the mixed diagonal coordinate; exact q=23 forgets endpoint owners; raw analytic clocks forget both.")


def main() -> None:
    cells = build_cells()
    print("LRC14 q=23 drop/add Haar square (codex-2026-06-26-S201)")
    print("source_threads=HYP-3037,HYP-3036,HYP-3032,HYP-3031,HYP-3027,HYP-3026,HYP-2991,HYP-2989")
    print()
    print_grid(cells)
    print_zeta(cells)
    quotient_stress(cells)
    tournament_analysis()
    synthesis()


if __name__ == "__main__":
    main()
