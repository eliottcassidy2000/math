#!/usr/bin/env python3
"""S150: packet-migration gauntlet for the LRC14 AP/GW boundary source core.

This script continues the HYP-2947/HYP-2951 packet route.  HYP-2951 checked
that in AP one-swap rows through add<=160 and AP two-swap rows through add<=40,
the only rows with threshold support but no strict Haar-open witness are AP and
Goddyn-Wong.  Here we stress that boundary-only claim farther:

  * one-swap AP rows through add<=420;
  * two-swap AP rows through add<=60;
  * three-swap AP rows through add<=30.

The script uses the q-divisibility witness first.  If qdiv(S)<14, then t=1/q
already gives a strict open witness, so exact interval work is unnecessary.
Only qdiv>=14 rows are sent to the exact S146 Haar/Baire interval classifier.

Tournament Analysis vertices are packet states and proof carriers, not runners.
The preserved predicate is whether the row stays in the zero-regular-open
boundary source core after qdiv, boundary-owner, C27, and state-lift labels are
attached.
"""

from __future__ import annotations

from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))
DELTA = Fraction(1, 14)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s146 = load_module(
    "s150_s146_taut_boundary",
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
class ExactClass:
    bank: str
    name: str
    speeds: tuple[int, ...]
    holes: tuple[int, ...]
    adds: tuple[int, ...]
    qdiv: int
    status: str
    mass: Fraction
    closed_points: int
    skeleton: tuple[tuple[str, tuple[str, ...]], ...]
    first_front: tuple[str, str, tuple[str, ...], tuple[str, ...], str] | None
    transfer: str
    atom_keys: tuple[str, ...]
    packet_state: str


def qdiv(speeds: tuple[int, ...], cap: int = 120) -> int:
    for d in range(2, cap + 1):
        if all(v % d for v in speeds):
            return d
    return cap + 1


def row_name(holes: tuple[int, ...], adds: tuple[int, ...]) -> str:
    if not holes and not adds:
        return "AP"
    if len(holes) == 1 and len(adds) == 1:
        return f"{holes[0]}->{adds[0]}"
    h = ",".join(map(str, holes))
    a = ",".join(map(str, adds))
    return f"drop({h})->add({a})"


def atom_keys(holes: tuple[int, ...], adds: tuple[int, ...]) -> tuple[str, ...]:
    remaining = list(adds)
    keys: list[str] = []
    for h in holes:
        for a in list(remaining):
            key = KNOWN_REPLACEMENTS.get((h, a))
            if key is None:
                continue
            keys.append(key)
            remaining.remove(a)
            break
    return tuple(keys)


def generate_bank(k: int, add_max: int) -> tuple[list[RowSpec], Counter[str]]:
    bank = f"{k}-swap add<={add_max}"
    rows: list[RowSpec] = []
    qroute_counts: Counter[str] = Counter()
    if k == 0:
        q = qdiv(AP)
        rows.append(RowSpec(bank, "AP", AP, (), (), q, "qdiv>=14 exact"))
        qroute_counts["qdiv>=14 exact"] += 1
        return rows, qroute_counts

    for holes in combinations(AP, k):
        base = set(AP) - set(holes)
        for adds in combinations(range(14, add_max + 1), k):
            speeds = tuple(sorted(base | set(adds)))
            if len(speeds) != 13 or gcd(*speeds) != 1:
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
            qroute_counts[qroute] += 1
    return rows, qroute_counts


def boundary_skeleton(speeds: tuple[int, ...], pts: tuple[Fraction, ...]) -> tuple[tuple[str, tuple[str, ...]], ...]:
    return tuple((str(t), s146.active_owners(speeds, t)) for t in pts)


def classify_exact(row: RowSpec) -> ExactClass:
    comps = s146.safe_open_components(row.speeds)
    mass = s146.interval_measure(comps)
    keys = atom_keys(row.holes, row.adds)
    transfer = s146.transfer_label(row.speeds)
    if mass > 0:
        a, b = comps[0]
        first_front = (
            str(a),
            str(b),
            s146.active_owners(row.speeds, a),
            s146.active_owners(row.speeds, b),
            str(s146.midpoint_slack(row.speeds, a, b)),
        )
        status = "positive_open"
        pts: tuple[Fraction, ...] = ()
        skeleton: tuple[tuple[str, tuple[str, ...]], ...] = ()
    else:
        first_front = None
        pts = s146.threshold_safe_points(row.speeds)
        skeleton = boundary_skeleton(row.speeds, pts)
        status = "boundary_only" if pts else "covered"

    packet_state = packet_state_for(row, status, keys, skeleton)
    return ExactClass(
        row.bank,
        row.name,
        row.speeds,
        row.holes,
        row.adds,
        row.qdiv,
        status,
        mass,
        len(pts),
        skeleton,
        first_front,
        transfer,
        keys,
        packet_state,
    )


def apgw_skeleton() -> tuple[tuple[str, tuple[str, ...]], ...]:
    pts = s146.threshold_safe_points(AP)
    return boundary_skeleton(AP, pts)


APGW_SKELETON = apgw_skeleton()


def packet_state_for(
    row: RowSpec,
    status: str,
    keys: tuple[str, ...],
    skeleton: tuple[tuple[str, tuple[str, ...]], ...],
) -> str:
    if status == "covered":
        return "counterexample-shaped-covered"
    if status == "boundary_only":
        if skeleton == APGW_SKELETON and (row.name == "AP" or keys == ("GW",)):
            return "AP/GW-boundary-source"
        if skeleton == APGW_SKELETON:
            return "AP/GW-skeleton-impostor"
        return "unknown-boundary-source"
    if "K33" in keys:
        return "positive-K33-state-lift"
    if set(keys) and set(keys) <= {"P10", "P13", "GW"}:
        return "positive-unit-petal-or-GW-strip"
    if row.qdiv > 14:
        return "positive-covering-off-apex"
    return "positive-unlabelled-q14-front"


def classify_bank(rows: list[RowSpec], workers: int, progress_every: int) -> list[ExactClass]:
    if not rows:
        return []
    out: list[ExactClass] = []
    if workers <= 1:
        for i, row in enumerate(rows, 1):
            out.append(classify_exact(row))
            if progress_every and i % progress_every == 0:
                print(f"    exact classified {i}/{len(rows)}", file=sys.stderr, flush=True)
        return out
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(classify_exact, row) for row in rows]
        for i, fut in enumerate(as_completed(futs), 1):
            out.append(fut.result())
            if progress_every and i % progress_every == 0:
                print(f"    exact classified {i}/{len(rows)}", file=sys.stderr, flush=True)
    return out


def fmt_fraction(x: Fraction) -> str:
    return str(x)


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    runners, residues mod 14, C27 shell pairs, boundary endpoints,")
    print("    open Haar fronts, wall-crossing events, Kpq/K33 owners, exact")
    print("    denominator packets, and proof obligations.")
    print("  chosen exact objects:")
    print("    qdiv-gated AP mutation rows, then S146 regular-open safe intervals")
    print("    and finite boundary-owner skeletons for qdiv>=14 rows.")
    print("  quotient preserves:")
    print("    whether U(S)={t:min||vt||>1/14} has positive Haar interior,")
    print("    whether zero-interior threshold support has the AP/GW owner skeleton,")
    print("    and the visible C27/K33 replacement labels.")
    print("  quotient destroys:")
    print("    exact M(S) away from the threshold, unbounded row identity, and")
    print("    analytic wide/decorrelation data.  Those must be attached later.")
    print("  challenged assumption:")
    print("    a row with AP-like residues is AP/GW-hard.  The test is stricter:")
    print("    it must keep zero regular-open witness and the labelled boundary")
    print("    source skeleton after qdiv and C27/K33 labels are retained.")


def print_bank_summary(
    bank: str,
    qroute_counts: Counter[str],
    exact: list[ExactClass],
    sample_count: int,
) -> None:
    status_counts = Counter(row.status for row in exact)
    packet_counts = Counter(row.packet_state for row in exact)
    print(f"[1] {bank}")
    print(f"  generated_rows={sample_count}")
    print(f"  qroute_counts={dict(sorted(qroute_counts.items()))}")
    print(f"  exact_qdiv>=14_rows={len(exact)}")
    print(f"  exact_status_counts={dict(sorted(status_counts.items()))}")
    print(f"  packet_state_counts={dict(sorted(packet_counts.items()))}")

    boundary = [row for row in exact if row.status == "boundary_only"]
    covered = [row for row in exact if row.status == "covered"]
    positives = sorted(
        (row for row in exact if row.status == "positive_open"),
        key=lambda r: (r.mass, r.bank, r.name),
    )
    if boundary:
        print("  boundary-only rows:")
        for row in boundary[:20]:
            sk = "APGW" if row.skeleton == APGW_SKELETON else "OTHER"
            print(
                f"    {row.name:34s} qdiv={row.qdiv:<3d} closed_pts={row.closed_points:<2d} "
                f"skeleton={sk:5s} keys={','.join(row.atom_keys) or '-':8s} transfer={row.transfer}"
            )
        if len(boundary) > 20:
            print(f"    ... {len(boundary) - 20} more")
    if covered:
        print("  COVERED rows, possible counterexample-shaped support:")
        for row in covered[:20]:
            print(f"    {row.name:34s} qdiv={row.qdiv:<3d} transfer={row.transfer}")
    print("  smallest exact positive Haar fronts among qdiv>=14 rows:")
    for row in positives[:10]:
        front = row.first_front
        if front is None:
            front_text = "-"
        else:
            a, b, left, right, slack = front
            front_text = f"({a},{b}) {left}->{right} slack={slack}"
        print(
            f"    {row.name:34s} mass={fmt_fraction(row.mass):>10s} "
            f"qdiv={row.qdiv:<3d} keys={','.join(row.atom_keys) or '-':8s} {front_text}"
        )


@dataclass(frozen=True)
class Carrier:
    name: str
    score: tuple[int, int, int, int, int, int]
    note: str


def carrier_tournament() -> tuple[list[Carrier], dict[str, set[str]]]:
    carriers = [
        Carrier("qdiv_gate", (6, 5, 4, 3, 4, 6), "strict q-witness before exact interval work"),
        Carrier("Haar_open_front", (5, 6, 5, 4, 4, 6), "regular-open positive witness/discharge"),
        Carrier("boundary_owner_skeleton", (5, 5, 6, 5, 5, 6), "zero-interior endpoint debt with owners"),
        Carrier("C27_transfer_labels", (4, 4, 5, 6, 5, 6), "unit/nonunit shell ownership"),
        Carrier("unital_affine_packet", (4, 4, 4, 6, 6, 6), "branch-local chart plus order/depth"),
        Carrier("K33_state_lift_flag", (3, 3, 4, 5, 6, 6), "HYP-2908/THM-572 endpoint address"),
        Carrier("wide_decorrelation_tail", (3, 5, 3, 3, 4, 5), "unbounded nonlocal escape control"),
        Carrier("raw_runner_residue", (1, 1, 1, 1, 1, 1), "too lossy alone"),
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


def tournament_fingerprint(vertices: list[str], out: dict[str, set[str]]) -> dict[str, object]:
    scores = {v: len(out[v]) for v in vertices}
    c3 = 0
    for a, b, c in combinations(vertices, 3):
        if (b in out[a] and c in out[b] and a in out[c]) or (
            c in out[a] and b in out[c] and a in out[b]
        ):
            c3 += 1
    return {
        "score_hist": dict(sorted(Counter(scores.values()).items())),
        "directed_3_cycles": c3,
        "hamiltonian_paths": hamiltonian_paths(vertices, out),
        "sccs": strongly_connected_components(vertices, out),
    }


def strongly_connected_components(vertices: list[str], out: dict[str, set[str]]) -> list[list[str]]:
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


def print_tournament_analysis(packet_counts: Counter[str]) -> None:
    print("[2] Tournament Analysis")
    print("  vertices:")
    print("    proof carriers and observed packet states, not runners or arcs.")
    print("  pairwise observable:")
    print("    which carrier preserves the zero-regular-open LRC predicate while")
    print("    destroying the least boundary-owner/C27/K33/state-lift data.")
    print("  switch/gauge:")
    print("    majority over qdiv retention, Haar interior, boundary code, C27")
    print("    owner data, state-lift fit, and anti-scalar resistance; ties follow")
    print("    the printed Hamiltonian path.")
    carriers, out = carrier_tournament()
    names = [c.name for c in carriers]
    fp = tournament_fingerprint(names, out)
    print("  carrier scores:")
    for c in carriers:
        print(f"    {c.name:24s} vector={c.score} note={c.note}")
    order = sorted(names, key=lambda n: -len(out[n]))
    print(f"  hamiltonian_path={' > '.join(order)}")
    print(f"  fingerprint={fp}")
    print("  observed packet states:")
    for state, count in sorted(packet_counts.items()):
        print(f"    {state:36s} {count}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--one-max", type=int, default=420)
    parser.add_argument("--two-max", type=int, default=60)
    parser.add_argument("--three-max", type=int, default=30)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=5000)
    args = parser.parse_args()

    print("=" * 78)
    print("S150 LRC14 packet-migration gauntlet")
    print("=" * 78)
    print_assumption_challenge()
    print()

    bank_specs = [
        (0, 13),
        (1, args.one_max),
        (2, args.two_max),
        (3, args.three_max),
    ]
    all_packet_counts: Counter[str] = Counter()
    any_boundary_impostor = False
    any_covered = False
    for k, add_max in bank_specs:
        rows, qcounts = generate_bank(k, add_max)
        sample_count = sum(qcounts.values())
        exact = classify_bank(rows, workers=args.workers, progress_every=args.progress_every)
        all_packet_counts.update(row.packet_state for row in exact)
        any_boundary_impostor = any_boundary_impostor or any(
            row.packet_state in {"AP/GW-skeleton-impostor", "unknown-boundary-source"}
            for row in exact
        )
        any_covered = any_covered or any(row.status == "covered" for row in exact)
        print_bank_summary(rows[0].bank if rows else f"{k}-swap add<={add_max}", qcounts, exact, sample_count)
        print()

    print_tournament_analysis(all_packet_counts)
    print()
    print("[3] Readout")
    if any_covered:
        print("  A covered qdiv>=14 row appeared.  Treat this as a possible")
        print("  counterexample-shaped support and audit immediately.")
    else:
        print("  No covered qdiv>=14 row appeared in this gauntlet.")
    if any_boundary_impostor:
        print("  A non-AP/GW boundary-only packet appeared.  This is a new")
        print("  boundary-source obligation, not a Haar-positive discharge.")
    else:
        print("  No non-AP/GW boundary-only packet appeared in this gauntlet.")
    print("  Evidence added by S150:")
    print(f"    one-swap AP rows through add<={args.one_max};")
    print(f"    two-swap AP rows through add<={args.two_max};")
    print(f"    three-swap AP rows through add<={args.three_max}.")
    print("  This supports the packet-source-core picture: outside AP/GW, rows")
    print("  seen here either have a strict qdiv/Haar witness, migrate to an")
    print("  open front, or carry labels for the C27/K33/state-lift branch.")


if __name__ == "__main__":
    main()
