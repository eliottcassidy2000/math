#!/usr/bin/env python3
"""NORK pinch-template audit for the LRC14 labelled-packet proof route.

NORK means "No Open Residual Kernel": the counterexample-shaped bucket in
HYP-2956 where qdiv >= 14, the strict safe set has no open interval, and the
row is not the AP/Goddyn-Wong boundary atom.

This script extends the packet-migration gauntlet in two ways:

* it mines endpoint-owner pinch templates from positive-open rows, so the
  proof target becomes a finite atlas of how near-F6 rows escape; and
* it optionally adds a first bounded 4-swap AP-neighborhood scan.

Tournament Analysis vertices are proof carriers and pinch templates, not raw
runners.  The pair observable is which carrier preserves the NORK predicate
and the labels needed to discharge it before scalarizing.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path
import argparse
import os
import sys


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s146 = load_module(
    "nork_s146_taut_boundary",
    REPO / "04-computation" / "lrc14_haar_baire_taut_boundary_s146.py",
)


KNOWN_REPLACEMENTS = {
    (12, 24): "GW",
    (12, 36): "K33",
    (10, 20): "P10",
    (13, 26): "P13",
}


@dataclass(frozen=True)
class RowSpec:
    bank: str
    name: str
    speeds: tuple[int, ...]
    holes: tuple[int, ...]
    adds: tuple[int, ...]
    qdiv: int
    qroute: str


@dataclass(frozen=True)
class Classified:
    bank: str
    name: str
    speeds: tuple[int, ...]
    holes: tuple[int, ...]
    adds: tuple[int, ...]
    qdiv: int
    qroute: str
    status: str
    family: str
    mass: Fraction
    closed_points: int
    skeleton: tuple[tuple[str, tuple[str, ...]], ...]
    atom_keys: tuple[str, ...]
    transfer: str
    pinch: tuple[str, str, tuple[str, ...], tuple[str, ...], str, str] | None
    template: str


@dataclass(frozen=True)
class Carrier:
    name: str
    score: tuple[int, int, int, int, int, int]
    note: str


def qdiv(speeds: tuple[int, ...], cap: int = 160) -> int:
    for d in range(2, cap + 1):
        if all(v % d for v in speeds):
            return d
    return cap + 1


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def row_name(holes: tuple[int, ...], adds: tuple[int, ...]) -> str:
    if not holes and not adds:
        return "AP"
    if len(holes) == 1 and len(adds) == 1:
        return f"{holes[0]}->{adds[0]}"
    return f"drop({','.join(map(str, holes))})->add({','.join(map(str, adds))})"


def atom_keys(holes: tuple[int, ...], adds: tuple[int, ...]) -> tuple[str, ...]:
    remaining = list(adds)
    keys: list[str] = []
    for h in holes:
        for a in list(remaining):
            key = KNOWN_REPLACEMENTS.get((h, a))
            if key is not None:
                keys.append(key)
                remaining.remove(a)
                break
    return tuple(keys)


def generate_bank(k: int, add_max: int) -> tuple[list[RowSpec], Counter[str]]:
    bank = f"{k}-swap add<={add_max}"
    rows: list[RowSpec] = []
    qroutes: Counter[str] = Counter()
    if k == 0:
        q = qdiv(AP)
        rows.append(RowSpec(bank, "AP", AP, (), (), q, "qdiv>=14 exact"))
        qroutes["qdiv>=14 exact"] += 1
        return rows, qroutes
    for holes in combinations(AP, k):
        base = set(AP) - set(holes)
        for adds in combinations(range(14, add_max + 1), k):
            speeds = tuple(sorted(base | set(adds)))
            if len(speeds) != 13 or not primitive(speeds):
                continue
            q = qdiv(speeds)
            if q < 14:
                qroute = "qdiv<14 strict witness"
            elif q == 14:
                qroute = "qdiv=14 exact"
                rows.append(RowSpec(bank, row_name(holes, adds), speeds, holes, adds, q, qroute))
            else:
                qroute = "qdiv>14 exact"
                rows.append(RowSpec(bank, row_name(holes, adds), speeds, holes, adds, q, qroute))
            qroutes[qroute] += 1
    return rows, qroutes


def boundary_skeleton(speeds: tuple[int, ...], pts: tuple[Fraction, ...]) -> tuple[tuple[str, tuple[str, ...]], ...]:
    return tuple((str(t), s146.active_owners(speeds, t)) for t in pts)


APGW_SKELETON = boundary_skeleton(AP, s146.threshold_safe_points(AP))


def owner_norm(owners: tuple[str, ...]) -> tuple[str, ...]:
    out: list[str] = []
    for owner in owners:
        side = owner[-1]
        speed = int(owner[:-1])
        shell = speed % 27
        shell = min(shell, 27 - shell)
        out.append(f"{side}:v{speed}:m14={speed % 14}:c27={shell}:g{gcd(shell, 27)}")
    return tuple(sorted(out))


def atom_class(keys: tuple[str, ...]) -> str:
    if "K33" in keys:
        return "K33"
    if keys and set(keys) <= {"P10", "P13", "GW"}:
        return "+".join(keys)
    return "unlabelled"


def q_class(q: int) -> str:
    if q < 14:
        return "q<14"
    if q == 14:
        return "q=14"
    return "q>14"


def pinch_component(
    speeds: tuple[int, ...],
    comps: list[tuple[Fraction, Fraction]],
) -> tuple[str, str, tuple[str, ...], tuple[str, ...], str, str] | None:
    if not comps:
        return None
    a, b = min(comps, key=lambda interval: (interval[1] - interval[0], interval[0]))
    left = s146.active_owners(speeds, a)
    right = s146.active_owners(speeds, b)
    length = b - a
    slack = s146.midpoint_slack(speeds, a, b)
    return (str(a), str(b), left, right, str(length), str(slack))


def classify_family(
    status: str,
    name: str,
    q: int,
    keys: tuple[str, ...],
    skeleton: tuple[tuple[str, tuple[str, ...]], ...],
) -> str:
    if status == "covered":
        return "F6-covered-zero-open"
    if status == "boundary_only":
        if skeleton == APGW_SKELETON and (name == "AP" or keys == ("GW",)):
            return "F1-AP/GW-boundary"
        return "F6-new-boundary-kernel"
    if "K33" in keys:
        return "F3-positive-K33"
    if keys and set(keys) <= {"P10", "P13", "GW"}:
        return "F2-positive-unit-petal"
    if q > 14:
        return "F5-positive-covering"
    return "F4-positive-q14-front"


def classify_row(row: RowSpec) -> Classified:
    comps = s146.safe_open_components(row.speeds)
    mass = s146.interval_measure(comps)
    keys = atom_keys(row.holes, row.adds)
    transfer = s146.transfer_label(row.speeds)
    if mass > 0:
        status = "positive_open"
        pts: tuple[Fraction, ...] = ()
        skeleton: tuple[tuple[str, tuple[str, ...]], ...] = ()
        pinch = pinch_component(row.speeds, comps)
    else:
        pinch = None
        pts = s146.threshold_safe_points(row.speeds)
        skeleton = boundary_skeleton(row.speeds, pts)
        status = "boundary_only" if pts else "covered"
    family = classify_family(status, row.name, row.qdiv, keys, skeleton)
    if pinch is None:
        template = family
    else:
        _, _, left, right, _, _ = pinch
        template = "|".join(
            [
                family,
                q_class(row.qdiv),
                atom_class(keys),
                ",".join(owner_norm(left)) or "-",
                ",".join(owner_norm(right)) or "-",
            ]
        )
    return Classified(
        row.bank,
        row.name,
        row.speeds,
        row.holes,
        row.adds,
        row.qdiv,
        row.qroute,
        status,
        family,
        mass,
        len(pts),
        skeleton,
        keys,
        transfer,
        pinch,
        template,
    )


def classify_bank(rows: list[RowSpec], workers: int, progress_every: int) -> list[Classified]:
    if workers <= 1:
        out: list[Classified] = []
        for i, row in enumerate(rows, 1):
            out.append(classify_row(row))
            if progress_every and i % progress_every == 0:
                print(f"    exact classified {i}/{len(rows)}", file=sys.stderr, flush=True)
        return out
    out = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(classify_row, row) for row in rows]
        for i, fut in enumerate(as_completed(futures), 1):
            out.append(fut.result())
            if progress_every and i % progress_every == 0:
                print(f"    exact classified {i}/{len(rows)}", file=sys.stderr, flush=True)
    return out


def fraction_text(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge / NORK quotient")
    print("  alternate vertices considered:")
    print("    runners, gaps, residues, boundary endpoints, safe components,")
    print("    fixed-margin packet classes, C27 shells, K33 incidence, moment signs,")
    print("    and proof obligations.")
    print("  chosen vertices:")
    print("    qdiv>=14 AP-mutation rows plus their tightest positive Haar/Baire")
    print("    pinch front, or their zero-open boundary skeleton if no front exists.")
    print("  preserved LRC predicate:")
    print("    whether the strict safe set U(S) has positive Haar interior; if not,")
    print("    whether the threshold support is exactly the AP/GW skeleton or a new")
    print("    F6 NORK obstruction.")
    print("  destroyed information:")
    print("    exact M away from threshold, global wide-decorrelation estimates, and")
    print("    row identity inside already-positive families.")
    print("  challenged assumption:")
    print("    positive-open rows are not all alike; the proof should keep their")
    print("    endpoint-owner pinch template before replacing them by a scalar mass.")


def summarize_bank(bank: str, generated: int, qroutes: Counter[str], rows: list[Classified]) -> Counter[str]:
    status_counts = Counter(row.status for row in rows)
    family_counts = Counter(row.family for row in rows)
    threats = [row for row in rows if row.family.startswith("F6")]
    positives = [row for row in rows if row.status == "positive_open"]
    by_mass = sorted(positives, key=lambda row: (row.mass, row.bank, row.name))
    by_width = sorted(
        positives,
        key=lambda row: (Fraction(row.pinch[4]) if row.pinch else Fraction(1), row.bank, row.name),
    )
    template_counts = Counter(row.template for row in positives)
    template_best: dict[str, Classified] = {}
    for row in positives:
        old = template_best.get(row.template)
        if old is None:
            template_best[row.template] = row
            continue
        old_width = Fraction(old.pinch[4]) if old.pinch else Fraction(1)
        new_width = Fraction(row.pinch[4]) if row.pinch else Fraction(1)
        if (new_width, row.mass, row.name) < (old_width, old.mass, old.name):
            template_best[row.template] = row

    print(f"[1] {bank}")
    print(f"  generated_rows={generated}")
    print(f"  qroute_counts={dict(sorted(qroutes.items()))}")
    print(f"  exact_qdiv>=14_rows={len(rows)}")
    print(f"  exact_status_counts={dict(sorted(status_counts.items()))}")
    print(f"  family_counts={dict(sorted(family_counts.items()))}")
    print(f"  pinch_template_count={len(template_counts)}")
    if threats:
        print("  NORK/F6 threats:")
        for row in threats[:20]:
            sk = "APGW" if row.skeleton == APGW_SKELETON else "OTHER"
            print(
                f"    {row.name:36s} status={row.status:13s} qdiv={row.qdiv:<3d} "
                f"closed={row.closed_points:<2d} skeleton={sk:5s} transfer={row.transfer}"
            )
        if len(threats) > 20:
            print(f"    ... {len(threats) - 20} more")
    print("  smallest positive masses:")
    for row in by_mass[:8]:
        print_front(row)
    print("  tightest positive pinch fronts:")
    for row in by_width[:8]:
        print_front(row)
    print("  most common pinch templates:")
    for template, count in template_counts.most_common(6):
        rep = template_best[template]
        width = rep.pinch[4] if rep.pinch else "-"
        print(f"    count={count:<6d} width_min={width:>10s} rep={rep.name} template={template[:160]}")
    print()
    return family_counts


def print_front(row: Classified) -> None:
    if row.pinch is None:
        print(f"    {row.name:36s} mass={fraction_text(row.mass):>10s} no-front")
        return
    a, b, left, right, width, slack = row.pinch
    print(
        f"    {row.name:36s} family={row.family:24s} mass={fraction_text(row.mass):>10s} "
        f"width={width:>10s} qdiv={row.qdiv:<3d} keys={','.join(row.atom_keys) or '-':8s} "
        f"({a},{b}) {left}->{right} slack={slack}"
    )


def carrier_tournament() -> tuple[list[Carrier], dict[str, set[str]]]:
    carriers = [
        Carrier("F6_NORK_sink", (6, 6, 6, 5, 6, 6), "only counterexample-shaped packet bucket"),
        Carrier("boundary_skeleton", (6, 5, 6, 5, 5, 6), "zero-open endpoint-owner code"),
        Carrier("pinch_template", (5, 6, 5, 6, 5, 6), "positive-open escape with owners retained"),
        Carrier("fixed_margin_packet", (5, 5, 5, 6, 6, 6), "family identity before scalarization"),
        Carrier("C27_K33_labels", (4, 4, 5, 6, 6, 6), "petal/K33 discharge and state-lift address"),
        Carrier("qdiv_gate", (5, 3, 3, 2, 5, 4), "direct denominator strict witness"),
        Carrier("raw_scalar_mass", (2, 4, 2, 1, 2, 2), "safe after discharge, weak before labels"),
        Carrier("raw_runner_set", (1, 1, 1, 1, 1, 1), "too lossy"),
    ]
    out = {c.name: set() for c in carriers}
    tie = {c.name: i for i, c in enumerate(carriers)}
    for a, b in combinations(carriers, 2):
        awins = sum(x > y for x, y in zip(a.score, b.score))
        bwins = sum(y > x for x, y in zip(a.score, b.score))
        if awins > bwins or (awins == bwins and tie[a.name] < tie[b.name]):
            out[a.name].add(b.name)
        else:
            out[b.name].add(a.name)
    return carriers, out


def hamiltonian_paths(vertices: list[str], out: dict[str, set[str]]) -> int:
    idx = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    dp: dict[tuple[int, int], int] = {}
    for v in vertices:
        dp[(1 << idx[v], idx[v])] = 1
    for mask in range(1 << n):
        for last in range(n):
            cur = dp.get((mask, last), 0)
            if not cur:
                continue
            for w in out[vertices[last]]:
                j = idx[w]
                if mask & (1 << j):
                    continue
                dp[(mask | (1 << j), j)] = dp.get((mask | (1 << j), j), 0) + cur
    full = (1 << n) - 1
    return sum(dp.get((full, j), 0) for j in range(n))


def sccs(vertices: list[str], out: dict[str, set[str]]) -> list[list[str]]:
    rev = {v: set() for v in vertices}
    for v, ws in out.items():
        for w in ws:
            rev[w].add(v)
    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in out[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)
    seen.clear()
    comps: list[list[str]] = []

    def rdfs(v: str, comp: list[str]) -> None:
        seen.add(v)
        comp.append(v)
        for w in rev[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[str] = []
            rdfs(v, comp)
            comps.append(sorted(comp))
    return sorted(comps, key=lambda c: (-len(c), c))


def print_tournament_analysis() -> None:
    carriers, out = carrier_tournament()
    vertices = [c.name for c in carriers]
    scores = {v: len(out[v]) for v in vertices}
    c3 = 0
    for a, b, c in combinations(vertices, 3):
        if (b in out[a] and c in out[b] and a in out[c]) or (
            c in out[a] and b in out[c] and a in out[b]
        ):
            c3 += 1
    print("[2] Tournament Analysis")
    print("  vertices are proof carriers/pinch templates, not runners.")
    print("  pair observable:")
    print("    which carrier preserves NORK status, boundary owners, packet labels,")
    print("    state-lift visibility, and anti-scalarization before discharge.")
    print("  switch:")
    print("    majority over the six retention coordinates; ties follow the listed order.")
    for carrier in carriers:
        print(f"    {carrier.name:22s} vector={carrier.score} note={carrier.note}")
    print(
        "  fingerprint="
        + str(
            {
                "score_hist": dict(sorted(Counter(scores.values()).items())),
                "directed_3_cycles": c3,
                "hamiltonian_paths": hamiltonian_paths(vertices, out),
                "sccs": sccs(vertices, out),
            }
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--one-max", type=int, default=420)
    parser.add_argument("--two-max", type=int, default=60)
    parser.add_argument("--three-max", type=int, default=34)
    parser.add_argument("--four-max", type=int, default=24)
    parser.add_argument("--workers", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--progress-every", type=int, default=10000)
    args = parser.parse_args()

    print("=" * 78)
    print("LRC14 NORK pinch-template audit")
    print("=" * 78)
    print_assumption_challenge()
    print()

    specs = [(0, 13), (1, args.one_max), (2, args.two_max), (3, args.three_max)]
    if args.four_max >= 17:
        specs.append((4, args.four_max))

    total_families: Counter[str] = Counter()
    any_f6 = False
    for k, add_max in specs:
        rows, qroutes = generate_bank(k, add_max)
        generated = sum(qroutes.values())
        classified = classify_bank(rows, args.workers, args.progress_every)
        total_families.update(summarize_bank(rows[0].bank if rows else f"{k}-swap add<={add_max}", generated, qroutes, classified))
        any_f6 = any_f6 or any(row.family.startswith("F6") and row.family != "F1-AP/GW-boundary" for row in classified)

    print_tournament_analysis()
    print()
    print("[3] NORK readout")
    if any_f6:
        print("  A non-AP/GW F6 NORK packet appeared.  This is a live proof obstruction.")
    else:
        print("  No non-AP/GW F6 NORK packet appeared in this audit.")
    print("  Total family counts:")
    for family, count in sorted(total_families.items()):
        print(f"    {family:28s} {count}")
    print("  Proof-shaped takeaway:")
    print("    Replace the scalar question 'is the safe measure small?' by the labelled")
    print("    question 'which endpoint-owner pinch template makes it positive?'.")
    print("    The LRC route should prove that every qdiv>=14 non-AP/GW packet either")
    print("    enters one of these positive pinch templates or exposes a genuine F6/F7")
    print("    obstruction.")


if __name__ == "__main__":
    main()
