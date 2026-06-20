#!/usr/bin/env python3
"""HYP-2677 / T917: packet-sign tournament atlas for the LRC(14) far route.

This atlas starts from the exact one-missed-sector packet telescope used in
HYP-2674/HYP-2676, then tests the user's six-sign prompt in two ways:

1. The six sector signs can be read as the six edges of a K4 tournament.
2. The six sectors themselves can be read as tournament vertices, with arcs
   supplied by exact signed pair-pressure gauges.

The guardrail is important: signs alone are chamber data.  The script therefore
also reports packet magnitudes, opposite K4 edge-pair balances, additive
profiles, and the quotient information destroyed by each tournament model.
All LRC packet quantities are exact Fractions.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass, field
from itertools import combinations, permutations
from typing import Callable, Iterable

from lrc14_signed_packet_et_ruzsa_codex_s48 import (
    F,
    RowPacket,
    additive_profile,
    named_rows,
    normalized,
    packet_report,
    scan_b14_bank,
    sign_word,
)


LEX_EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EDGE_GAUGES = {
    "lex": LEX_EDGES,
    "opposite_pairs": ((0, 1), (2, 3), (0, 2), (1, 3), (0, 3), (1, 2)),
    "path_then_chords": ((0, 1), (1, 2), (2, 3), (0, 2), (1, 3), (0, 3)),
    "face_then_apex": ((0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)),
}
OPPOSITE_EDGE_PAIRS = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)


@dataclass
class AtlasRow:
    name: str
    Ep: tuple[int, ...]
    w: int
    delta: F
    wdelta: F
    packets: dict[int, F]
    packet_abs: dict[int, F]
    terms: list[F]
    runs: int
    profile: dict[str, object]
    aliases: list[str] = field(default_factory=list)

    @property
    def label(self) -> str:
        if not self.aliases:
            return self.name
        return self.name + "/" + "/".join(self.aliases)

    @property
    def sign(self) -> str:
        return sign_word(self.packets)

    @property
    def sector_abs(self) -> F:
        return sum(abs(v) for v in self.packets.values())

    @property
    def run_abs(self) -> F:
        return sum(abs(v) for v in self.terms)

    @property
    def signed_efficiency(self) -> F:
        return abs(self.wdelta) / self.run_abs if self.run_abs else F(0)

    @property
    def abel_bound(self) -> F:
        """Incoming THM-546 S2 signed Abel bound: |Delta_w| <= (6/49) V(E')/w."""
        return F(6 * self.runs, 49 * self.w)

    @property
    def abel_pressure(self) -> F:
        return abs(self.delta) / self.abel_bound if self.abel_bound else F(0)


def q(x: F | int | object) -> str:
    return str(x)


def sign_char(x: F) -> str:
    return "+" if x > 0 else "-" if x < 0 else "0"


def make_row(name: str, Ep: Iterable[int], w: int) -> AtlasRow:
    row = normalized(Ep)
    delta, packets, packet_abs, terms, runs = packet_report(row, w)
    return AtlasRow(
        name=name,
        Ep=row,
        w=w,
        delta=delta,
        wdelta=delta * w,
        packets=packets,
        packet_abs=packet_abs,
        terms=terms,
        runs=runs,
        profile=additive_profile(row),
    )


def add_alias(rows: dict[tuple[tuple[int, ...], int], AtlasRow], name: str, Ep: tuple[int, ...], w: int) -> None:
    key = (normalized(Ep), w)
    if key in rows:
        if name != rows[key].name and name not in rows[key].aliases:
            rows[key].aliases.append(name)
    else:
        rows[key] = make_row(name, Ep, w)


def selected_rows() -> list[AtlasRow]:
    rows: dict[tuple[tuple[int, ...], int], AtlasRow] = {}

    top_b14, _buckets = scan_b14_bank()
    for i, rep in enumerate(top_b14, 1):
        add_alias(rows, f"B14_top{i:02d}", rep.Ep, rep.w)

    for rep in named_rows():
        add_alias(rows, rep.name, rep.Ep, rep.w)

    extras = [
        ("B14_even_AP_boundary_probe", (0, 2, 4, 6, 8, 10, 12, 14), 16),
        ("B14_odd_near_AP_probe", (0, 1, 3, 5, 7, 9, 11, 13), 15),
        ("wide_two_cluster_probe", (0, 1, 2, 3, 30, 31, 32, 33), 35),
        ("mixed_squarefree_probe", (0, 2, 3, 5, 7, 11, 13, 17), 19),
    ]
    for name, Ep, w in extras:
        add_alias(rows, name, Ep, w)

    def key(row: AtlasRow) -> tuple[int, F, F, tuple[int, ...], int]:
        is_b14 = 0 if row.name.startswith("B14_top") else 1
        return (is_b14, -row.delta, -row.signed_efficiency, row.Ep, row.w)

    return sorted(rows.values(), key=key)


def empty_adj(n: int) -> list[list[bool]]:
    return [[False for _ in range(n)] for _ in range(n)]


def orient(adj: list[list[bool]], a: int, b: int) -> None:
    adj[a][b] = True
    adj[b][a] = False


def tournament_fingerprint(adj: list[list[bool]]) -> dict[str, object]:
    n = len(adj)
    scores = tuple(sum(1 for j in range(n) if i != j and adj[i][j]) for i in range(n))
    cycle3 = 0
    for tri in combinations(range(n), 3):
        local = [sum(1 for j in tri if i != j and adj[i][j]) for i in tri]
        if sorted(local) == [1, 1, 1]:
            cycle3 += 1

    def reach(start: int, reverse: bool = False) -> set[int]:
        seen = {start}
        todo = deque([start])
        while todo:
            u = todo.popleft()
            for v in range(n):
                edge = adj[v][u] if reverse else adj[u][v]
                if v not in seen and edge:
                    seen.add(v)
                    todo.append(v)
        return seen

    remaining = set(range(n))
    scc_sizes = []
    while remaining:
        start = min(remaining)
        comp = reach(start) & reach(start, reverse=True)
        scc_sizes.append(len(comp))
        remaining -= comp

    full = (1 << n) - 1
    dp = [[0 for _ in range(n)] for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += count

    return {
        "scores": scores,
        "score_hist": tuple(sorted(Counter(scores).items())),
        "score_sorted": tuple(sorted(scores, reverse=True)),
        "cycle3": cycle3,
        "scc": tuple(sorted(scc_sizes, reverse=True)),
        "hp": sum(dp[full]),
    }


def fp_key(fp: dict[str, object]) -> tuple[object, ...]:
    return (fp["score_sorted"], fp["cycle3"], fp["scc"], fp["hp"])


def fp_short(fp: dict[str, object]) -> str:
    return f"scores={fp['score_sorted']}, c3={fp['cycle3']}, scc={fp['scc']}, HP={fp['hp']}"


def k4_from_sign_word(signs: str, edge_order: tuple[tuple[int, int], ...]) -> list[list[bool]]:
    adj = empty_adj(4)
    for ch, (a, b) in zip(signs, edge_order):
        if ch == "-":
            orient(adj, b, a)
        else:
            orient(adj, a, b)
    return adj


def k4_all_permutation_types(signs: str) -> Counter[tuple[object, ...]]:
    out: Counter[tuple[object, ...]] = Counter()
    for permuted in set(permutations(signs)):
        fp = tournament_fingerprint(k4_from_sign_word("".join(permuted), LEX_EDGES))
        out[fp_key(fp)] += 1
    return out


def k4_type_summary(counts: Counter[tuple[object, ...]]) -> str:
    parts = []
    for key, count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
        scores, cycles, scc, hp = key
        parts.append(f"{count}x(scores={scores},c3={cycles},scc={scc},HP={hp})")
    return "; ".join(parts)


def edge_pair_balances(row: AtlasRow, edge_order: tuple[tuple[int, int], ...]) -> tuple[tuple[str, F], ...]:
    edge_to_sector = {edge: sector for sector, edge in enumerate(edge_order, 1)}
    out = []
    for e1, e2 in OPPOSITE_EDGE_PAIRS:
        s1 = edge_to_sector[e1]
        s2 = edge_to_sector[e2]
        total = row.packets[s1] + row.packets[s2]
        out.append((f"{s1}+{s2}", total))
    return tuple(out)


def sector_values(row: AtlasRow) -> tuple[F, ...]:
    return tuple(row.packets[s] for s in range(1, 7))


def opposite_sector(s: int) -> int:
    return 7 - s


def cyclic_base(s: int, t: int) -> int:
    step = (t - s) % 7
    if step in (1, 2, 3):
        return 1
    if step in (4, 5, 6):
        return -1
    raise ValueError((s, t))


PairGauge = Callable[[AtlasRow, int, int], F]


def scalar_signed_rank(row: AtlasRow, s: int, t: int) -> F:
    return row.packets[s] - row.packets[t]


def cyclic_pair_sum(row: AtlasRow, s: int, t: int) -> F:
    return cyclic_base(s, t) * (row.packets[s] + row.packets[t])


def opposite_compensated_cyclic(row: AtlasRow, s: int, t: int) -> F:
    paired = row.packets[s] + row.packets[t]
    opposite = row.packets[opposite_sector(s)] + row.packets[opposite_sector(t)]
    return cyclic_base(s, t) * (paired - opposite)


def absolute_pressure_cyclic(row: AtlasRow, s: int, t: int) -> F:
    paired = abs(row.packets[s]) + abs(row.packets[t])
    opposite = abs(row.packets[opposite_sector(s)]) + abs(row.packets[opposite_sector(t)])
    return cyclic_base(s, t) * (paired - opposite)


SECTOR_GAUGES: dict[str, PairGauge] = {
    "scalar_signed_rank": scalar_signed_rank,
    "cyclic_pair_sum": cyclic_pair_sum,
    "opposite_compensated": opposite_compensated_cyclic,
    "absolute_opposite": absolute_pressure_cyclic,
}


def sector_tournament(row: AtlasRow, gauge: PairGauge) -> tuple[list[list[bool]], int]:
    adj = empty_adj(6)
    ties = 0
    for s, t in combinations(range(1, 7), 2):
        obs = gauge(row, s, t)
        if obs > 0:
            orient(adj, s - 1, t - 1)
        elif obs < 0:
            orient(adj, t - 1, s - 1)
        else:
            ties += 1
            orient(adj, s - 1, t - 1)
    return adj, ties


def print_row_bank(rows: list[AtlasRow]) -> None:
    print("Rows in atlas")
    for row in rows:
        prof = row.profile
        print(
            f"  {row.label}: E'={row.Ep}, w={row.w}, "
            f"Delta={q(row.delta)}, wDelta={q(row.wdelta)}, sign={row.sign}, "
            f"eff={q(row.signed_efficiency)}, V={row.runs}, "
            f"AbelPressure={q(row.abel_pressure)}, exc={prof['sumset_excess']}, "
            f"K2={prof['K2']}, sqfree={prof['sqfree']}"
        )
    print()


def print_k4_edge_atlas(rows: list[AtlasRow]) -> None:
    print("K4 edge-sign quotient")
    print("  six sectors are assigned to K4 edges; + uses the gauge edge direction, - flips it")
    for row in rows:
        lex_fp = tournament_fingerprint(k4_from_sign_word(row.sign, EDGE_GAUGES["lex"]))
        all_types = k4_all_permutation_types(row.sign)
        balances = edge_pair_balances(row, EDGE_GAUGES["lex"])
        balance_text = ", ".join(f"{name}:{sign_char(value)}{q(abs(value))}" for name, value in balances)
        print(
            f"  {row.label}: sign={row.sign}, lex {fp_short(lex_fp)}, "
            f"all sector-edge permutations [{k4_type_summary(all_types)}], "
            f"opposite-pair balances {balance_text}"
        )
    print()


def print_gauge_examples(rows: list[AtlasRow]) -> None:
    print("K4 gauge examples")
    focus_names = {"KPS_third_pocket", "HYP2675_boundary_leader_base", "B14_top01"}
    for row in rows:
        if row.name not in focus_names and not any(alias in focus_names for alias in row.aliases):
            continue
        print(f"  {row.label}:")
        for gauge_name, edge_order in EDGE_GAUGES.items():
            fp = tournament_fingerprint(k4_from_sign_word(row.sign, edge_order))
            balances = edge_pair_balances(row, edge_order)
            balance_text = ", ".join(f"{name}:{sign_char(value)}{q(abs(value))}" for name, value in balances)
            print(f"    {gauge_name}: {fp_short(fp)}; pair balances {balance_text}")
    print()


def print_sector_vertex_atlas(rows: list[AtlasRow]) -> None:
    print("Six-sector vertex tournaments")
    print("  vertices are missed sectors 1..6; pairwise gauges use exact packet values")
    for row in rows:
        print(f"  {row.label}: sector packets={tuple(q(v) for v in sector_values(row))}")
        for gauge_name, gauge in SECTOR_GAUGES.items():
            adj, ties = sector_tournament(row, gauge)
            fp = tournament_fingerprint(adj)
            print(f"    {gauge_name}: {fp_short(fp)}, ties={ties}")
    print()


def row_family_summary(rows: list[AtlasRow]) -> None:
    print("Family-level contrasts")
    groups: dict[str, list[AtlasRow]] = defaultdict(list)
    for row in rows:
        if row.name.startswith("B14_top"):
            groups["B14_top12"].append(row)
        elif "third_pocket" in row.label:
            groups["KPS_third_pocket"].append(row)
        elif "HYP2675" in row.label:
            groups["HYP2675_wide_branch"].append(row)
        elif "dyadic" in row.label or "doubled" in row.label:
            groups["dyadic_doubled_resonances"].append(row)
        else:
            groups["other_named"].append(row)

    for group, members in groups.items():
        signs = Counter(row.sign for row in members)
        k4_types = Counter()
        sector_types = Counter()
        for row in members:
            k4_types[fp_key(tournament_fingerprint(k4_from_sign_word(row.sign, LEX_EDGES)))] += 1
            adj, _ties = sector_tournament(row, cyclic_pair_sum)
            sector_types[fp_key(tournament_fingerprint(adj))] += 1
        print(f"  {group}: count={len(members)}, signs={dict(signs)}")
        print(f"    K4 lex types={dict(k4_types)}")
        print(f"    cyclic-pair sector types={dict(sector_types)}")
    print()


def print_hypothesis_triage(rows: list[AtlasRow]) -> None:
    b14 = [row for row in rows if row.name.startswith("B14_top")]
    kps = [row for row in rows if "KPS_third_pocket" in row.label][0]
    boundary = [row for row in rows if "HYP2675_boundary" in row.label][0]
    true_wide = [row for row in rows if "HYP2675_true_wide" in row.label][0]

    b14_signs = Counter(row.sign for row in b14)
    b14_k4 = {fp_key(tournament_fingerprint(k4_from_sign_word(row.sign, LEX_EDGES))) for row in b14}
    kps_perm = k4_all_permutation_types(kps.sign)
    boundary_perm = k4_all_permutation_types(boundary.sign)

    print("Hypothesis triage")
    if b14_signs == Counter({"++++++": len(b14)}) and len(b14_k4) == 1:
        print("  PASS H1: B14 near-speed leaders collapse to one transitive K4 sign type.")
        print("    Reading: the K4 sign quotient locates the same-sign pocket but cannot separate its rows.")
    else:
        print(f"  FAIL H1: B14 signs/types vary: signs={dict(b14_signs)}, K4types={len(b14_k4)}")

    if len(kps_perm) > 1:
        print("  PASS H2: the KPS third-pocket sign word has multiple K4 types under sector-edge relabelling.")
        print(f"    Exact type distribution: {k4_type_summary(kps_perm)}")
    else:
        print("  FAIL H2: KPS third-pocket sign word is K4-rigid under relabelling.")

    if boundary.sign != true_wide.sign:
        print("  PASS H3: HYP-2675 boundary and true-wide rows split before magnitudes.")
        print(f"    boundary sign={boundary.sign}; true-wide sign={true_wide.sign}")
    else:
        print("  FAIL H3: HYP-2675 boundary and true-wide have same sign word.")

    kps_adj, _ = sector_tournament(kps, opposite_compensated_cyclic)
    kps_fp = tournament_fingerprint(kps_adj)
    b14_fps = {
        fp_key(tournament_fingerprint(sector_tournament(row, opposite_compensated_cyclic)[0]))
        for row in b14
    }
    if fp_key(kps_fp) not in b14_fps:
        print("  PASS H4: opposite-compensated sector tournament separates KPS from all B14 leaders.")
        print(f"    KPS opposite-compensated type: {fp_short(kps_fp)}")
    else:
        print("  WARN H4: KPS shares an opposite-compensated sector type with a B14 leader.")

    same_sign_rows = [row for row in rows if row.sign == "++++++"]
    same_sign_sector_types = {
        fp_key(tournament_fingerprint(sector_tournament(row, opposite_compensated_cyclic)[0]))
        for row in same_sign_rows
    }
    if len(same_sign_sector_types) > 1:
        print("  PASS H5: magnitudes inside the same `++++++` chamber create multiple sector-tournament types.")
        print(f"    same-sign opposite-compensated type count={len(same_sign_sector_types)}")
    else:
        print("  WARN H5: chosen sector gauges do not split the same-sign chamber.")

    print()
    print("Proof-route reading")
    print("  1. K4 edge signs are a finite chamber coordinate: useful for detecting the ++++++ pocket,")
    print("     but too lossy to prove cancellation by themselves.")
    print("  2. Opposite K4 edge-pair sums are a better typed address for the user's positive/negative")
    print("     prompt: they keep exact packet mass and expose cancellation before absolute values.")
    print("  3. The six-sector vertex gauges turn a negative packet-pair total into an arc flip in a")
    print("     fixed mod-7 cyclic tournament, so KPS-style cancellation becomes visible as topology.")
    print("  4. Incoming THM-546 S2 gives the rational signed Abel bound")
    print("     |Delta_w| <= (6/49) V(E')/w.  The tournament atlas is therefore a classifier")
    print("     for tight ungapped/same-sign packets, not a substitute for the analytic bound.")
    print("  5. The next sharp theorem should classify ++++++ rows by finite Ruzsa/Freiman model plus")
    print("     opposite-pair packet balances; the complement should prove that a cyclic-pair arc flip")
    print("     or small exact pair mass forces signed Erdos-Turan cancellation.")
    print("  6. No LRC(14) proof is claimed here.")


def main() -> None:
    rows = selected_rows()
    print("HYP-2677 / T917 LRC14 packet-sign tournament atlas")
    print("Arithmetic: exact Fractions for packet masses; no numeric approximation is used in the tables.")
    print("Pairwise observables are declared for each tournament quotient.")
    print()
    print_row_bank(rows)
    print_k4_edge_atlas(rows)
    print_gauge_examples(rows)
    print_sector_vertex_atlas(rows)
    row_family_summary(rows)
    print_hypothesis_triage(rows)
    print("PASS: packet tournament atlas complete.")


if __name__ == "__main__":
    main()
