#!/usr/bin/env python3
"""LRC14 residual status-gate switchboard.

This is a lightweight proof-interface audit.  It deliberately reuses the
stored full-bank S188/S189 evidence instead of recomputing the HYP-2963 exact
packet bank.  The target is the remaining mixed-route residue after the coarse
Erdos-Turan plus Henselian-unit gate:

    mix route remains small, but mixed boundary/open status is already zero.

That suggests a status-first proof theorem: prove boundary/open purity before
trying to classify every open route.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations, permutations
from pathlib import Path
import re


REPO = Path(__file__).resolve().parents[1]
S188_OUT = REPO / "05-knowledge" / "results" / "lrc14_fiber_zipper_convergence_codex_s188.out"
S189_OUT = REPO / "05-knowledge" / "results" / "lrc14_carrier_fusion_switchboard_codex_s189.out"


@dataclass(frozen=True)
class GateEvidence:
    name: str
    fibers: int
    mixed_route: int
    mixed_status: int
    max_mixed: int


@dataclass(frozen=True)
class ResidualExample:
    name: str
    route: str
    m_value: str
    safe_mu: str
    q0: int
    exit_tooth: str
    reason: str


@dataclass(frozen=True)
class Carrier:
    name: str
    predicate: int
    compression: int
    route_split: int
    arithmetic: int
    topology: int
    proof_cost: int
    anti_address: int

    @property
    def vector(self) -> tuple[int, ...]:
        return (
            self.predicate,
            self.compression,
            self.route_split,
            self.arithmetic,
            self.topology,
            6 - self.proof_cost,
            self.anti_address,
        )


TIE_PATH = (
    "status_first_coarse_et_unit_gate",
    "magnitude_cocycle_fallback",
    "safe_stick_barcode_certificate",
    "henselian_unit_lift_rule",
    "exact_et_address_clock",
    "barcode_packet_zipper_address",
    "raw_automatic_word",
)


def read_text_auto(path: Path) -> str:
    data = path.read_bytes()
    for encoding in ("utf-8-sig", "utf-16"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            pass
    return data.decode("utf-8", errors="replace")


def read_s188() -> str:
    return read_text_auto(S188_OUT)


def read_s189() -> str:
    return read_text_auto(S189_OUT)


def parse_gate(text: str, name: str) -> GateEvidence:
    pattern = re.compile(
        rf"^\s+{re.escape(name)}\s+(\d+)\s+(\d+)\s+\d+\s+(\d+)\s+\d+\s+(\d+)",
        re.MULTILINE,
    )
    match = pattern.search(text)
    if not match:
        raise ValueError(f"missing gate row: {name}")
    fibers, mixed_route, mixed_status, max_mixed = (int(x) for x in match.groups())
    return GateEvidence(name, fibers, mixed_route, mixed_status, max_mixed)


def parse_target_gate(text: str) -> tuple[int, int, int]:
    match = re.search(
        r"coarse_et_unit_gate: fibers=(\d+) mixed_route=(\d+)\n"
        r"\s+largest_mixed_rows=(\d+)",
        text,
    )
    if not match:
        raise ValueError("missing target coarse gate summary")
    return tuple(int(x) for x in match.groups())


def parse_s189_barcode_stress(text: str) -> tuple[int, int, int]:
    match = re.search(r"automatic_plus_barcode_shape \| \d+ \| (\d+) \| (\d+) \| (\d+) \|", text)
    if not match:
        raise ValueError("missing S189 barcode stress row")
    fibers, mixed_status, mixed_route = (int(x) for x in match.groups())
    return fibers, mixed_status, mixed_route


def residual_examples() -> list[ResidualExample]:
    return [
        ResidualExample(
            "single swap 13->143",
            "COVERING-MOMENT",
            "11/144",
            "4307/194040",
            14,
            "strict-safe-stick / covering-moment",
            "q0=14 but safe_mu>0, so the LRC predicate is already strict-open",
        ),
        ResidualExample(
            "single swap 13->31",
            "Q-WITNESS",
            "1/13",
            "194039/6015240",
            13,
            "direct q-witness",
            "q0<=13 gives a denominator witness strictly above 1/14",
        ),
        ResidualExample(
            "single swap 13->59",
            "Q-WITNESS",
            "1/13",
            "53509/1635480",
            13,
            "direct q-witness",
            "q0<=13 gives a denominator witness strictly above 1/14",
        ),
        ResidualExample(
            "single swap 13->115",
            "Q-WITNESS",
            "1/13",
            "13469/405720",
            13,
            "direct q-witness",
            "q0<=13 gives a denominator witness strictly above 1/14",
        ),
    ]


def orient(a: Carrier, b: Carrier) -> str:
    av = a.vector
    bv = b.vector
    awins = sum(x > y for x, y in zip(av, bv))
    bwins = sum(x < y for x, y in zip(av, bv))
    if awins != bwins:
        return a.name if awins > bwins else b.name
    if av != bv:
        return a.name if av > bv else b.name
    return a.name if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else b.name


def tournament_fingerprint(carriers: list[Carrier]) -> tuple[Counter, int, list[int], int, list[str]]:
    names = [carrier.name for carrier in carriers]
    edges: dict[tuple[int, int], bool] = {}
    score = Counter()
    for i, j in combinations(range(len(carriers)), 2):
        winner = orient(carriers[i], carriers[j])
        i_wins = winner == names[i]
        edges[(i, j)] = i_wins
        score[i if i_wins else j] += 1

    score_hist = Counter(score[i] for i in range(len(carriers)))
    directed_3cycles = 0
    for a, b, c in combinations(range(len(carriers)), 3):
        ab = edges[(a, b)]
        ac = edges[(a, c)]
        bc = edges[(b, c)]
        if (ab and bc and not ac) or ((not ab) and (not bc) and ac):
            directed_3cycles += 1

    mat = [[False] * len(carriers) for _ in carriers]
    for i, j in combinations(range(len(carriers)), 2):
        if edges[(i, j)]:
            mat[i][j] = True
        else:
            mat[j][i] = True

    def reach(start: int, graph: list[list[bool]]) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            node = stack.pop()
            for nxt, ok in enumerate(graph[node]):
                if ok and nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        return seen

    reverse = [[mat[j][i] for j in range(len(mat))] for i in range(len(mat))]
    remaining = set(range(len(mat)))
    scc_sizes = []
    while remaining:
        seed = min(remaining)
        comp = reach(seed, mat) & reach(seed, reverse)
        scc_sizes.append(len(comp))
        remaining -= comp

    hp_count = 0
    for order in permutations(range(len(mat))):
        if all(mat[order[i]][order[i + 1]] for i in range(len(order) - 1)):
            hp_count += 1

    ranking = sorted(range(len(carriers)), key=lambda idx: (-score[idx], names[idx]))
    return score_hist, directed_3cycles, sorted(scc_sizes, reverse=True), hp_count, [names[idx] for idx in ranking]


def carriers() -> list[Carrier]:
    return [
        Carrier("status_first_coarse_et_unit_gate", 5, 5, 2, 4, 2, 2, 5),
        Carrier("magnitude_cocycle_fallback", 5, 2, 5, 5, 1, 3, 2),
        Carrier("safe_stick_barcode_certificate", 5, 3, 3, 1, 5, 3, 4),
        Carrier("henselian_unit_lift_rule", 4, 4, 2, 5, 1, 2, 5),
        Carrier("exact_et_address_clock", 5, 0, 5, 5, 0, 1, 0),
        Carrier("barcode_packet_zipper_address", 5, 0, 5, 4, 5, 1, 0),
        Carrier("raw_automatic_word", 1, 5, 0, 1, 0, 1, 5),
    ]


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, section boundaries, residue fibers, ET clocks,")
    print("    p-adic unit roots, zero-root scale debt, safe sticks, barcode")
    print("    bars, magnitude cocycles, route labels, and proof obligations.")
    print("  chosen vertices:")
    print("    residual proof obligations after the coarse ET+unit status gate.")
    print("  preserved LRC predicate:")
    print("    boundary/open status at threshold 1/14, not full theorem-route purity.")
    print("  destroyed information:")
    print("    exact route family, exact ET address, endpoint owners, and full")
    print("    magnitude/barcode identity until a certificate tooth is reattached.")
    print("  challenged assumption:")
    print("    route mixing after status convergence is not a counterexample signal;")
    print("    it is certificate scheduling debt.")
    print()


def main() -> None:
    s188 = read_s188()
    s189 = read_s189()
    coarse = parse_gate(s188, "coarse_et_unit_gate")
    magnitude = parse_gate(s188, "magnitude_cocycle")
    target_fibers, target_mixed, target_largest = parse_target_gate(s188)
    barcode_fibers, barcode_mixed_status, barcode_mixed_route = parse_s189_barcode_stress(s189)
    score_hist, c3, scc, hp, ranking = tournament_fingerprint(carriers())

    print("LRC14 residual status-gate switchboard (codex-2026-06-26-S192)")
    print("source_threads=HYP-3026,HYP-3025,HYP-3024,HYP-3023,HYP-3020,HYP-3017,HYP-3016,HYP-3015,HYP-2963")
    print()
    print_assumption_challenge()
    print("[1] Stored full-bank evidence")
    print(
        f"  full_bank_coarse_et_unit_gate fibers={coarse.fibers} "
        f"mixed_route={coarse.mixed_route} mixed_status={coarse.mixed_status} "
        f"max_mixed={coarse.max_mixed}"
    )
    print(
        f"  full_bank_magnitude_cocycle fibers={magnitude.fibers} "
        f"mixed_route={magnitude.mixed_route} mixed_status={magnitude.mixed_status} "
        f"max_mixed={magnitude.max_mixed}"
    )
    print(
        f"  target_word_coarse_gate fibers={target_fibers} "
        f"mixed_route={target_mixed} largest_mixed={target_largest}"
    )
    print(
        f"  s189_named_automatic_plus_barcode_shape fibers={barcode_fibers} "
        f"mixed_status={barcode_mixed_status} mixed_route={barcode_mixed_route}"
    )
    print("  interpretation=coarse ET+unit is already status-pure; magnitude is the route-pure fallback.")
    print()

    print("[2] Largest residual mixed-route sample from S188")
    print("  name | route | M | safe_mu | q0 | exit_tooth | reason")
    for row in residual_examples():
        print(
            "  "
            + " | ".join(
                (
                    row.name,
                    row.route,
                    row.m_value,
                    row.safe_mu,
                    str(row.q0),
                    row.exit_tooth,
                    row.reason,
                )
            )
        )
    print()

    print("[3] Status-first proof target")
    print("  The theorem should not first ask for route purity.  It should prove:")
    print("    in each automatic/residue fiber, coarse ET+Henselian-unit data cannot")
    print("    identify an AP/GW boundary equality atom with a strict-open packet.")
    print("  Once that is proved, the remaining route-mixed fibers are harmless for")
    print("  LRC14 itself and can be discharged by certificate teeth:")
    print("    q<=13 direct witness; AP/GW zero-bar equality; strict safe stick or")
    print("    Fejer/Haar certificate; unit-petal; K33/F7/THM-572 state lift;")
    print("    covering-moment; magnitude-cocycle family formula.")
    print()

    print("[4] Tournament Analysis")
    print("  vertices=proof obligations after the residual gate, not runners")
    print("  pairwise_observable=predicate, compression, route_split, arithmetic, topology, proof_cost, anti_address")
    print("  switch=majority comparison of carrier vectors")
    print("  tie_hamiltonian_path=" + ">".join(TIE_PATH))
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={c3}")
    print(f"  scc_sizes={scc}")
    print(f"  hamiltonian_path_count={hp}")
    print("  score_order=" + ">".join(ranking))
    print()

    print("[5] Next executable hook")
    print("  Add a cached packet-ledger mode to HYP-2963 so S188 residual fibers can be")
    print("  listed without recomputing exact M.  Then attach for each residual fiber")
    print("  the first successful tooth among q-witness, safe-stick/barcode, Fejer,")
    print("  magnitude formula, unit-petal, K33/F7, and covering-moment.")


if __name__ == "__main__":
    main()
