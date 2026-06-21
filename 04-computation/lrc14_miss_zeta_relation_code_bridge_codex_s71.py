#!/usr/bin/env python3
"""HYP-2726/S71: bridge generated-word compatibility to relation-code packets.

HYP-2722 says cheap abstract q0-hiding atom moves must first pass the
HYP-2698 generated miss-zeta word cone.  HYP-2724 says the remaining carrier
error is organized by the low-support weight enumerator of the relation code
Lambda(E).  This scout checks whether those are one theorem or two sequential
filters.

For each sparse-tail challenger shape from the HYP-2722 frontier, aggregate the
exact generated-word metrics over coherent contexts, then attach the KPS
relation-code observables:

    dmin, A2, A3, A4 for primitive relations with |coef| <= 2.

The test asks whether relation-code features predict the generated q0/W12/U4
barriers directly, or only identify the finite low-support packets after the
generated-word compatibility filter has already fired.
"""

from __future__ import annotations

import importlib.util
import itertools
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
S71_PATH = ROOT / "04-computation" / "lrc14_miss_zeta_word_compatibility_codex_s71.py"
KPS_PATH = ROOT / "04-computation" / "lrc_q108_relation_code_mds_kps.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


s71 = load_module("s71_word_compat_bridge", S71_PATH)
kps = load_module("kps_relation_code_bridge", KPS_PATH)


def fmt_float(x: float | None) -> str:
    if x is None:
        return "None"
    return f"{x:+.3f}"


def fmt_frac(x: F | None) -> str:
    if x is None:
        return "None"
    return f"{x} ({float(x):.9f})"


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 3 or len(set(xs)) < 2 or len(set(ys)) < 2:
        return None
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den_x = sum((x - mx) ** 2 for x in xs)
    den_y = sum((y - my) ** 2 for y in ys)
    den = (den_x * den_y) ** 0.5
    if den == 0:
        return None
    return num / den


@dataclass(frozen=True)
class ShapeRow:
    size: int
    shape: tuple[int, ...]
    A2: int
    A3: int
    A4: int
    dmin: int
    dA3: int
    dA4: int
    min_q0: F
    min_w12: F
    min_u4: F
    min_b2: F
    min_tax: F
    min_d1: F

    def feature(self, name: str) -> float:
        return float(getattr(self, name))


def shape_rows(metrics: Iterable[object]) -> list[ShapeRow]:
    by_shape: dict[tuple[int, tuple[int, ...]], list[object]] = defaultdict(list)
    for metric in metrics:
        if metric.q0 > 0:
            by_shape[(metric.size, metric.shape)].append(metric)

    rows: list[ShapeRow] = []
    for (size, shape), group in sorted(by_shape.items()):
        counts, dmin = kps.support_spectrum(shape, B=2, max_support=min(4, size))
        kcounts, _ = kps.support_spectrum(tuple(range(size)), B=2, max_support=min(4, size))
        rows.append(
            ShapeRow(
                size=size,
                shape=shape,
                A2=counts.get(2, 0),
                A3=counts.get(3, 0),
                A4=counts.get(4, 0),
                dmin=dmin if dmin is not None else 99,
                dA3=kcounts.get(3, 0) - counts.get(3, 0),
                dA4=kcounts.get(4, 0) - counts.get(4, 0),
                min_q0=min(m.q0 for m in group),
                min_w12=min(m.low_leak_2 for m in group if m.low_leak_2 is not None),
                min_u4=min(m.u4_norm for m in group if m.u4_norm is not None),
                min_b2=min(m.b2_norm for m in group if m.b2_norm is not None),
                min_tax=min(m.normalized_tax for m in group if m.normalized_tax is not None),
                min_d1=min(m.cheap_dist_1 for m in group if m.cheap_dist_1 is not None),
            )
        )
    return rows


def print_correlation_block(title: str, rows: list[ShapeRow]) -> None:
    features = ["A2", "A3", "A4", "dmin", "dA3", "dA4"]
    targets = ["min_q0", "min_w12", "min_u4", "min_b2", "min_tax", "min_d1"]
    print(title)
    for feature in features:
        entries: list[str] = []
        xs = [row.feature(feature) for row in rows]
        for target in targets:
            ys = [row.feature(target) for row in rows]
            r = pearson(xs, ys)
            if r is not None:
                entries.append(f"{target}={fmt_float(r)}")
        if entries:
            print(f"  {feature:<5} " + " ".join(entries))


def print_extreme_rows(rows: list[ShapeRow]) -> None:
    print("\nEXTREME SHAPE WITNESSES")
    specs = [
        ("smallest q0", lambda r: r.min_q0),
        ("smallest W1+W2 leakage", lambda r: r.min_w12),
        ("smallest U4/q0", lambda r: r.min_u4),
        ("smallest B2/q0", lambda r: r.min_b2),
        ("closest cheap r=1 direction", lambda r: r.min_d1),
    ]
    for title, key in specs:
        row = min(rows, key=key)
        print(f"  {title}:")
        print(
            f"    size={row.size} shape={row.shape} "
            f"A2/A3/A4={row.A2}/{row.A3}/{row.A4} dmin={row.dmin} "
            f"dA3={row.dA3} dA4={row.dA4}"
        )
        print(
            f"    q0={fmt_frac(row.min_q0)} W12={fmt_frac(row.min_w12)} "
            f"U4/q0={fmt_frac(row.min_u4)} B2/q0={fmt_frac(row.min_b2)} "
            f"tax={fmt_frac(row.min_tax)} d1={fmt_frac(row.min_d1)}"
        )


def tournament_bridge(rows: list[ShapeRow]) -> None:
    print("\nTOURNAMENT ANALYSIS")
    print("  vertices: sparse-tail challenger shapes, not runners or arcs")
    print("  generated-risk edge: smaller (q0, W12, U4, cheap-r1 distance)")
    print("  relation-proxy edge: smaller dmin, then larger A3, then larger A4")
    print("  tie Hamiltonian path: lexicographic shape order")

    for size in sorted({row.size for row in rows}):
        group = [row for row in rows if row.size == size]
        if len(group) < 2:
            continue

        def risk_key(row: ShapeRow):
            return (row.min_q0, row.min_w12, row.min_u4, row.min_d1, row.shape)

        def relation_key(row: ShapeRow):
            return (row.dmin, -row.A3, -row.A4, row.shape)

        risk_edges = set()
        relation_edges = set()
        risk_score = Counter()
        relation_score = Counter()
        flips = 0
        total = 0
        for a, b in itertools.combinations(group, 2):
            if risk_key(a) <= risk_key(b):
                risk_edges.add((a.shape, b.shape))
                risk_score[a.shape] += 1
                risk_winner = a.shape
            else:
                risk_edges.add((b.shape, a.shape))
                risk_score[b.shape] += 1
                risk_winner = b.shape
            if relation_key(a) <= relation_key(b):
                relation_edges.add((a.shape, b.shape))
                relation_score[a.shape] += 1
                relation_winner = a.shape
            else:
                relation_edges.add((b.shape, a.shape))
                relation_score[b.shape] += 1
                relation_winner = b.shape
            total += 1
            if risk_winner != relation_winner:
                flips += 1

        cycles = 0
        for a, b, c in itertools.combinations([row.shape for row in group], 3):
            if (a, b) in risk_edges and (b, c) in risk_edges and (c, a) in risk_edges:
                cycles += 1
            if (a, c) in risk_edges and (c, b) in risk_edges and (b, a) in risk_edges:
                cycles += 1

        print(f"  size={size}: vertices={len(group)} pairs={total} proxy_flips={flips}")
        print(
            "    risk_score_hist="
            + str(dict(sorted(Counter(risk_score[row.shape] for row in group).items())))
        )
        print(
            "    relation_score_hist="
            + str(dict(sorted(Counter(relation_score[row.shape] for row in group).items())))
        )
        print(f"    risk_directed_3cycles={cycles}")
        risk_path = " > ".join(str(row.shape) for row in sorted(group, key=risk_key)[:8])
        rel_path = " > ".join(str(row.shape) for row in sorted(group, key=relation_key)[:8])
        print(f"    risk path head: {risk_path}")
        print(f"    relation path head: {rel_path}")


def print_synthesis(rows: list[ShapeRow]) -> None:
    print("\nSYNTHESIS")
    print("  Relation-code features are real signal but not a replacement for")
    print("  generated-word compatibility.  Globally A3/A4/dmin correlate strongly")
    print("  with q0, U4, and B2, but that is partly a size effect.  Within fixed")
    print("  shape size, the signs are mixed, and size=3 has the worst q0/U4")
    print("  witnesses while the simple |coef|<=2 relation spectrum is flat.")
    print("  Therefore the proof stack should stay sequential:")
    print("    generated miss-zeta death-chain/context-merge compatibility first;")
    print("    relation-code A3/MDS packet selection second;")
    print("    factorial q0/Vitali atom boundary evaluation last.")
    print("  The useful theorem target is not 'A3 implies compatibility'.  It is")
    print("  'compatibility reduces surviving packets to a low-support relation")
    print("  ledger whose leading finite cases are organized by A3/dmin/MDS'.")


def main() -> None:
    print("HYP-2726/S71 miss-zeta compatibility / relation-code bridge")
    print("Exact generated-word metrics from HYP-2722; relation spectra from HYP-2724.\n")
    metrics = s71.frontier_metrics()
    rows = shape_rows(metrics)
    print(f"\nAGGREGATE: tests={len(metrics)} unique_shapes={len(rows)}")
    print_correlation_block("\nGLOBAL PEARSON (unique shapes; size-confounded)", rows)
    print("\nPER-SIZE PEARSON (unique shapes)")
    for size in sorted({row.size for row in rows}):
        group = [row for row in rows if row.size == size]
        print(f"\nsize={size} n={len(group)}")
        print_correlation_block("", group)
    print_extreme_rows(rows)
    tournament_bridge(rows)
    print_synthesis(rows)


if __name__ == "__main__":
    main()
