#!/usr/bin/env python3
"""Residual certificate teeth for the LRC14 coarse status gate.

This lightweight pass parses the stored S194 status-topology output instead of
recomputing the full HYP-2963 packet bank.  It focuses only on the 38 packets
inside the 15 route-mixed fibers left by the coarse ET+Henselian-unit status
gate.

The tested angle is a small non-route scheduler:

    arc topology compact signature + unit-scale tooth

where the topology compact signature is

    closed beta0, closed beta1, open beta0, open beta1, safe topes, quotient defect

and the unit-scale tooth only records whether the displayed exact M is a unit
fraction with denominator at most 13, a larger unit fraction, or nonunit scale.

Tournament Analysis declaration:
  vertices: residual proof carriers, not runners;
  pairwise observable: status preservation, residual route purity, topology
    retention, scale retention, compression, noncircularity, and proof cost;
  switch/gauge: majority comparison of observable vectors;
  tie Hamiltonian path: topology+scale tooth > exact magnitude fallback >
    topology tooth > unit-scale tooth > coarse residual fiber.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path
import re


REPO = Path(__file__).resolve().parents[1]
S194_OUT = REPO / "05-knowledge" / "results" / "lrc14_status_topology_gate_codex_s194.out"


@dataclass(frozen=True)
class ResidualRow:
    fiber: int
    name: str
    status: str
    route: str
    M: Fraction
    mu: Fraction
    closed_b0: int
    closed_b1: int
    open_b0: int
    open_b1: int
    topes: int
    quotient_defect: int

    @property
    def topology(self) -> tuple[int, int, int, int, int, int]:
        return (
            self.closed_b0,
            self.closed_b1,
            self.open_b0,
            self.open_b1,
            self.topes,
            self.quotient_defect,
        )

    @property
    def unit_scale_tooth(self) -> str:
        if self.M.numerator == 1 and self.M.denominator <= 13:
            return "unit-q<=13"
        if self.M.numerator == 1:
            return "unit-q>=14"
        return "nonunit-scale"

    @property
    def topology_bucket(self) -> tuple[int, int, int]:
        return (self.topes, self.quotient_defect, self.open_b0 - self.closed_b0)


def parse_fraction(text: str) -> Fraction:
    return Fraction(text)


ROW_RE = re.compile(
    r"^\s{4}(?P<name>.*?)\s+"
    r"status=(?P<status>\S+)\s+"
    r"route=(?P<route>\S+)\s+"
    r"M=\s*(?P<M>\S+)\s+"
    r"mu=\s*(?P<mu>\S+)\s+"
    r"closed=\((?P<closed_b0>\d+),(?P<closed_b1>\d+)\)\s+"
    r"open=\((?P<open_b0>\d+),(?P<open_b1>\d+)\)\s+"
    r"topes=(?P<topes>\d+)\s+"
    r"qdef=(?P<qdef>\d+)\s+"
    r"apgw_cycle=(?P<apgw>\S+)"
)


def parse_residual_rows(path: Path = S194_OUT) -> list[ResidualRow]:
    rows: list[ResidualRow] = []
    active_fiber: int | None = None
    in_section = False

    for line in path.read_text().splitlines():
        if line.startswith("[2] Coarse ET+unit route-mixed fibers"):
            in_section = True
            continue
        if in_section and line.startswith("[3] "):
            break
        if not in_section:
            continue

        header = re.match(r"^\s{2}coarse_route_mixed_fiber\[(\d+)\]", line)
        if header:
            active_fiber = int(header.group(1))
            continue

        match = ROW_RE.match(line)
        if match and active_fiber is not None:
            data = match.groupdict()
            rows.append(
                ResidualRow(
                    fiber=active_fiber,
                    name=data["name"].strip(),
                    status=data["status"],
                    route=data["route"],
                    M=parse_fraction(data["M"]),
                    mu=parse_fraction(data["mu"]),
                    closed_b0=int(data["closed_b0"]),
                    closed_b1=int(data["closed_b1"]),
                    open_b0=int(data["open_b0"]),
                    open_b1=int(data["open_b1"]),
                    topes=int(data["topes"]),
                    quotient_defect=int(data["qdef"]),
                )
            )

    if not rows:
        raise ValueError(f"no residual rows parsed from {path}")
    return rows


@dataclass(frozen=True)
class CarrierAudit:
    name: str
    description: str
    status: int
    route: int
    topology: int
    scale: int
    compression: int
    noncircular: int
    proof_cost: int
    key_func: object

    @property
    def vector(self) -> tuple[int, ...]:
        return (
            self.status,
            self.route,
            self.topology,
            self.scale,
            self.compression,
            self.noncircular,
            8 - self.proof_cost,
        )


def carrier_definitions() -> tuple[CarrierAudit, ...]:
    return (
        CarrierAudit(
            "coarse_residual_fiber",
            "S194 coarse ET+unit residual fiber id",
            5,
            1,
            1,
            1,
            5,
            5,
            1,
            lambda row: (row.fiber,),
        ),
        CarrierAudit(
            "topology_tooth",
            "closed/open arc Betti, safe tope count, quotient defect",
            5,
            3,
            5,
            1,
            4,
            5,
            2,
            lambda row: row.topology,
        ),
        CarrierAudit(
            "unit_scale_tooth",
            "unit fraction q<=13 / unit q>=14 / nonunit scale",
            5,
            3,
            1,
            4,
            5,
            4,
            1,
            lambda row: (row.unit_scale_tooth,),
        ),
        CarrierAudit(
            "topology_bucket_plus_unit_scale",
            "safe-topes/quotient-defect topology bucket joined to unit-scale tooth",
            5,
            5,
            4,
            4,
            4,
            5,
            2,
            lambda row: (row.topology_bucket, row.unit_scale_tooth),
        ),
        CarrierAudit(
            "topology_plus_unit_scale",
            "full topology compact signature joined to unit-scale tooth",
            5,
            5,
            5,
            4,
            3,
            5,
            2,
            lambda row: (row.topology, row.unit_scale_tooth),
        ),
        CarrierAudit(
            "exact_M_fallback",
            "exact M value on the residual ledger",
            5,
            3,
            1,
            5,
            1,
            3,
            4,
            lambda row: (row.M,),
        ),
    )


def group_rows(rows: list[ResidualRow], carrier: CarrierAudit) -> dict[object, list[ResidualRow]]:
    groups: dict[object, list[ResidualRow]] = defaultdict(list)
    for row in rows:
        groups[carrier.key_func(row)].append(row)
    return groups


def audit_carrier(rows: list[ResidualRow], carrier: CarrierAudit) -> dict[str, object]:
    groups = group_rows(rows, carrier)
    mixed_route = [group for group in groups.values() if len({row.route for row in group}) > 1]
    mixed_status = [group for group in groups.values() if len({row.status for row in group}) > 1]
    return {
        "groups": groups,
        "fiber_count": len(groups),
        "mixed_route": len(mixed_route),
        "mixed_status": len(mixed_status),
        "max_mixed": max((len(group) for group in mixed_route), default=0),
        "mixed_groups": sorted(mixed_route, key=lambda group: (-len(group), group[0].fiber, group[0].name)),
    }


TIE_PATH = (
    "topology_plus_unit_scale",
    "topology_bucket_plus_unit_scale",
    "exact_M_fallback",
    "topology_tooth",
    "unit_scale_tooth",
    "coarse_residual_fiber",
)


def orient(a: CarrierAudit, b: CarrierAudit) -> str:
    av = a.vector
    bv = b.vector
    awins = sum(x > y for x, y in zip(av, bv))
    bwins = sum(x < y for x, y in zip(av, bv))
    if awins != bwins:
        return a.name if awins > bwins else b.name
    if av != bv:
        return a.name if av > bv else b.name
    return a.name if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else b.name


def tournament_fingerprint(carriers: tuple[CarrierAudit, ...]) -> dict[str, object]:
    names = [carrier.name for carrier in carriers]
    edges: dict[tuple[int, int], bool] = {}
    score = Counter()
    for i, j in combinations(range(len(carriers)), 2):
        winner = orient(carriers[i], carriers[j])
        i_wins = winner == names[i]
        edges[(i, j)] = i_wins
        score[i if i_wins else j] += 1

    c3 = 0
    for a, b, c in combinations(range(len(carriers)), 3):
        ab = edges[(a, b)]
        ac = edges[(a, c)]
        bc = edges[(b, c)]
        if (ab and bc and not ac) or ((not ab) and (not bc) and ac):
            c3 += 1

    matrix = [[False] * len(carriers) for _ in carriers]
    for i, j in combinations(range(len(carriers)), 2):
        if edges[(i, j)]:
            matrix[i][j] = True
        else:
            matrix[j][i] = True

    hp = 0
    for order in permutations(range(len(carriers))):
        if all(matrix[order[i]][order[i + 1]] for i in range(len(order) - 1)):
            hp += 1

    ranking = sorted(range(len(carriers)), key=lambda i: (-score[i], TIE_PATH.index(names[i])))
    return {
        "score_hist": dict(sorted(Counter(score[i] for i in range(len(carriers))).items())),
        "directed_3cycles": c3,
        "hamiltonian_path_count": hp,
        "score_order": [names[i] for i in ranking],
    }


def fmt(frac: Fraction) -> str:
    return str(frac.numerator) if frac.denominator == 1 else f"{frac.numerator}/{frac.denominator}"


def print_group(group: list[ResidualRow]) -> None:
    for row in sorted(group, key=lambda item: (item.route, item.name)):
        print(
            f"      fiber={row.fiber:02d} {row.name:<36} route={row.route:<16} "
            f"M={fmt(row.M):>7} mu={fmt(row.mu):>12} "
            f"topology={row.topology} unit={row.unit_scale_tooth}"
        )


def main() -> None:
    rows = parse_residual_rows()
    carriers = carrier_definitions()
    audits = {carrier.name: audit_carrier(rows, carrier) for carrier in carriers}

    print("=== LRC14 residual certificate teeth S197 ===")
    print(f"source={S194_OUT.relative_to(REPO)}")
    print(f"residual_rows={len(rows)} residual_fibers={len(set(row.fiber for row in rows))}")
    print(f"route_counts={dict(Counter(row.route for row in rows))}")
    print(f"status_counts={dict(Counter(row.status for row in rows))}")
    print()

    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, coarse residual fibers, exact routes, topology classes,")
    print("    scale teeth, safe-tope buckets, magnitude values, and proof teeth.")
    print("  chosen vertices:")
    print("    residual proof carriers, not runners or route labels.")
    print("  preserved LRC predicate:")
    print("    boundary/open status is already pure; the carrier schedules open")
    print("    residual route certificates without using route labels as a key.")
    print("  destroyed information:")
    print("    raw coarse ET address, exact magnitude except for a small unit/nonunit")
    print("    scale tooth, and row identity.")
    print()

    print("[1] Residual carrier audit")
    for carrier in carriers:
        audit = audits[carrier.name]
        print(
            f"  {carrier.name:<34} fibers={audit['fiber_count']:<2} "
            f"mixed_route={audit['mixed_route']:<2} max_mixed={audit['max_mixed']:<2} "
            f"mixed_status={audit['mixed_status']:<2} :: {carrier.description}"
        )
    print()

    print("[2] Remaining mixed pairs before the joined tooth")
    for name in ("topology_tooth", "unit_scale_tooth"):
        audit = audits[name]
        print(f"  {name}: mixed_route={audit['mixed_route']} max_mixed={audit['max_mixed']}")
        for group in audit["mixed_groups"]:
            print(f"    routes={dict(Counter(row.route for row in group))}")
            print_group(group)
    print()

    print("[3] Joined-tooth route scheduler")
    joined = audits["topology_plus_unit_scale"]
    print(
        f"  topology_plus_unit_scale fibers={joined['fiber_count']} "
        f"mixed_route={joined['mixed_route']} mixed_status={joined['mixed_status']}"
    )
    print("  compressed bucket version:")
    bucket = audits["topology_bucket_plus_unit_scale"]
    print(
        f"    topology_bucket_plus_unit_scale fibers={bucket['fiber_count']} "
        f"mixed_route={bucket['mixed_route']} mixed_status={bucket['mixed_status']}"
    )
    print("  pure route fibers by tooth:")
    groups = group_rows(rows, next(carrier for carrier in carriers if carrier.name == "topology_plus_unit_scale"))
    for key, group in sorted(groups.items(), key=lambda item: (str(item[0]), len(item[1]))):
        print(f"    tooth={key} size={len(group)} routes={dict(Counter(row.route for row in group))}")
    print()

    print("[4] Tournament Analysis over residual proof carriers")
    fp = tournament_fingerprint(carriers)
    print("  vertices_are=residual proof carriers, not runners")
    print("  score_hist=" + str(fp["score_hist"]))
    print("  directed_3cycles=" + str(fp["directed_3cycles"]))
    print("  hamiltonian_path_count=" + str(fp["hamiltonian_path_count"]))
    print("  score_order=" + " > ".join(fp["score_order"]))
    print()

    print("[5] Proof readout")
    print("  1. S194 leaves 15 route-mixed coarse ET+unit fibers, but all 38")
    print("     packets are strict-open; this pass only schedules certificates.")
    print("  2. Topology alone contracts the 15 mixed fibers to 3 mixed topology")
    print("     classes, but a common open-topology class still mixes routes.")
    print("  3. Unit-scale alone separates nonunit covering rows from most direct")
    print("     q-witness rows, but unit-scale covering rows still mix with")
    print("     unit q-witness rows.")
    print("  4. Joining topology compact signature with the unit-scale tooth gives")
    print("     21 residual fibers and 0 route-mixed fibers.")
    print("  5. The compressed topology bucket joined to the same unit-scale tooth")
    print("     gives the same purity on the stored S194 residual ledger.")


if __name__ == "__main__":
    main()
