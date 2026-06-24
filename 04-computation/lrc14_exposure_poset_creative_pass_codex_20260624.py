#!/usr/bin/env python3
"""Creative LRC14 exposure-poset proof pass.

The recent repo frontier says that a strict LRC14 counterexample is no longer
just a 13-speed row.  It would have to be an unexposed labelled source packet:
no q-witness, no AP/GW boundary atom, no positive Haar bridge, no K33/state-lift
address, and no finite Fejer/Toeplitz dual.

This script audits that slogan as a finite proof object.  Rows are tested
against certificate channels, then the channels are compared by a tournament.

Tournament Analysis declaration:

  vertices:
    certificate/exposure channels, not runners or arcs.

  pairwise observable:
    severity-weighted exposure of audited rows, plus retention of the LRC
    predicate fields: q-threshold, open-vs-boundary status, exact safe
    component, Fejer degree/margin, packet-family label, state-lift/petal
    address, and anti-scalarization guard.

  switch/gauge:
    channel A beats B if it exposes more weighted hard mass not exposed by B;
    ties are broken by a declared retention vector and then the Hamiltonian
    path shown in CHANNELS.

  challenged assumption:
    a proof channel should not be judged by how many easy rows it sees.  It is
    judged by how much of the low-margin, high-degree, q>=14 frontier it
    exposes while retaining the labels needed for a formal handoff.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, log10
from pathlib import Path
import argparse
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


s157 = load_module(
    "exp_s157_fejer",
    REPO / "04-computation" / "lrc14_fourier_toeplitz_fejer_fullbank_codex_s157.py",
)


@dataclass(frozen=True)
class ExposureAudit:
    name: str
    family: str
    speeds: tuple[int, ...]
    qdiv: int
    safe_mu: Fraction
    component_count: int
    largest_width: Fraction
    fejer_degree: int | None
    fejer_qform: float | None
    invisible_speeds: tuple[int, ...]
    channels: frozenset[str]

    @property
    def severity(self) -> float:
        """Heuristic weight for rows close to the current proof frontier."""
        weight = 1.0
        if self.qdiv >= 14:
            weight += 3.0
        if self.qdiv > 14:
            weight += 2.0
        if self.safe_mu == 0:
            weight += 4.0
        else:
            weight += min(4.0, 0.02 / float(self.safe_mu))
        if self.fejer_degree is not None:
            weight += min(3.0, self.fejer_degree / 100.0)
        if self.fejer_qform is not None and self.fejer_qform < 0:
            weight += min(3.0, -log10(max(-self.fejer_qform, 1e-18)) / 4.0)
        if "K33_STATE_LIFT" in self.channels:
            weight += 2.0
        if "HARD_FEJER_MARGIN" in self.channels:
            weight += 2.0
        return weight


CHANNELS = (
    "Q_WITNESS",
    "AP_GW_TAUT_BOUNDARY",
    "OPEN_HAAR_BRIDGE",
    "FEJER_PSD_DUAL",
    "K33_STATE_LIFT",
    "C27_PETAL_EXIT",
    "LATE_COVERING_PRESSURE",
    "HARD_FEJER_MARGIN",
    "UNEXPOSED_SOURCE_KERNEL",
)


RETENTION: dict[str, tuple[int, int, int, int, int, int, int]] = {
    # q, boundary, open interval, harmonic dual, packet label, state label,
    # anti-scalarization.
    "Q_WITNESS": (10, 2, 3, 1, 2, 1, 6),
    "AP_GW_TAUT_BOUNDARY": (9, 10, 1, 3, 8, 3, 8),
    "OPEN_HAAR_BRIDGE": (7, 7, 10, 4, 5, 3, 7),
    "FEJER_PSD_DUAL": (6, 5, 8, 10, 6, 3, 9),
    "K33_STATE_LIFT": (7, 6, 6, 5, 9, 10, 9),
    "C27_PETAL_EXIT": (8, 7, 7, 4, 9, 4, 8),
    "LATE_COVERING_PRESSURE": (10, 5, 5, 6, 7, 5, 8),
    "HARD_FEJER_MARGIN": (6, 5, 7, 10, 8, 4, 10),
    "UNEXPOSED_SOURCE_KERNEL": (1, 1, 1, 1, 1, 1, 1),
}


def fmt_fraction(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def fmt_float(x: float | None) -> str:
    if x is None:
        return "-"
    if abs(x) < 1e-3:
        return f"{x:.6e}"
    return f"{x:.9f}"


def primitive(row: tuple[int, ...]) -> bool:
    g = 0
    for v in row:
        g = gcd(g, v)
    return g == 1


def qdiv(row: tuple[int, ...], cap: int = 240) -> int:
    for d in range(2, cap + 1):
        if all(v % d != 0 for v in row):
            return d
    return cap + 1


def unique_rows(rows: list) -> list:
    seen: set[tuple[int, ...]] = set()
    out = []
    for row in rows:
        speeds = tuple(sorted(set(row.speeds)))
        if len(speeds) != 13 or not primitive(speeds) or speeds in seen:
            continue
        seen.add(speeds)
        out.append(s157.Row(row.name, row.family, speeds))
    return out


def build_rows(
    single_limit: int,
    two_swap_limit: int,
    alias_depth: int,
    lcm_tail_max: int,
    include_hard_tail: bool,
) -> list:
    rows = []
    rows.extend(s157.named_rows())
    rows.extend(s157.build_one_swap_rows(single_limit))
    rows.extend(
        s157.build_hyp2963_rows(
            single_limit=single_limit,
            two_swap_limit=two_swap_limit,
            alias_depth=alias_depth,
            lcm_tail_max=lcm_tail_max,
        )
    )
    if include_hard_tail:
        # Pull in the S162 hard rows without importing the whole interval-budget
        # script. These are the rows where Fejer formalization is most relevant.
        def repl(holes: tuple[int, ...], adds: tuple[int, ...]) -> tuple[int, ...]:
            return tuple(sorted((set(AP) - set(holes)) | set(adds)))

        Row = s157.Row
        rows.extend(
            [
                Row("hard single 6->63", "hard Fejer one-swap", repl((6,), (63,))),
                Row("hard single 6->86", "hard Fejer one-swap", repl((6,), (86,))),
                Row(
                    "hard two drop(12,13)->add(14,29)",
                    "small margin two-swap",
                    repl((12, 13), (14, 29)),
                ),
                Row(
                    "hard two drop(6,10)->add(30,34)",
                    "hard Fejer two-swap",
                    repl((6, 10), (30, 34)),
                ),
                Row(
                    "hard two drop(10,12)->add(24,30)",
                    "hard Fejer two-swap",
                    repl((10, 12), (24, 30)),
                ),
            ]
        )
    return unique_rows(rows)


def classify_channels(row, audit, q: int) -> frozenset[str]:
    channels: set[str] = set()
    safe_mu = audit.safe_measure
    cert = audit.certificate
    text = f"{row.name} {row.family}".lower()
    if q <= 13:
        channels.add("Q_WITNESS")
    if safe_mu == 0 and (
        "ap" in text
        or "gw" in text
        or row.speeds in {AP, tuple(list(range(1, 12)) + [13, 24])}
    ):
        channels.add("AP_GW_TAUT_BOUNDARY")
    if safe_mu > 0:
        channels.add("OPEN_HAAR_BRIDGE")
    if cert is not None:
        channels.add("FEJER_PSD_DUAL")
        if cert.degree >= 120 or cert.qform > -1e-5:
            channels.add("HARD_FEJER_MARGIN")
    if "k33" in text or "state" in text:
        channels.add("K33_STATE_LIFT")
    if "petal" in text or "p10" in text or "c27" in text:
        channels.add("C27_PETAL_EXIT")
    if q > 14 or "covering" in text:
        channels.add("LATE_COVERING_PRESSURE")
    if not channels or channels <= {"LATE_COVERING_PRESSURE"}:
        channels.add("UNEXPOSED_SOURCE_KERNEL")
    return frozenset(channels)


def audit_rows(rows: list, max_degree: int, negative_tol: float) -> list[ExposureAudit]:
    ghat = s157.ghat_table(max_degree)
    out: list[ExposureAudit] = []
    for row in rows:
        q = qdiv(row.speeds)
        row_audit = s157.first_fejer_certificate(row, max_degree, ghat, negative_tol)
        cert = row_audit.certificate
        channels = classify_channels(row, row_audit, q)
        out.append(
            ExposureAudit(
                name=row.name,
                family=row.family,
                speeds=row.speeds,
                qdiv=q,
                safe_mu=row_audit.safe_measure,
                component_count=row_audit.component_count,
                largest_width=row_audit.largest_width,
                fejer_degree=cert.degree if cert else None,
                fejer_qform=cert.qform if cert else None,
                invisible_speeds=row_audit.invisible_speeds_at_first,
                channels=channels,
            )
        )
    return out


def channel_exposure(audits: list[ExposureAudit]) -> dict[str, dict[str, float]]:
    stats: dict[str, dict[str, float]] = {
        ch: {"rows": 0.0, "weight": 0.0, "unique_weight": 0.0, "frontier": 0.0}
        for ch in CHANNELS
    }
    for audit in audits:
        severity = audit.severity
        for ch in audit.channels:
            stats[ch]["rows"] += 1.0
            stats[ch]["weight"] += severity
            if len(audit.channels) == 1:
                stats[ch]["unique_weight"] += severity
            if audit.qdiv >= 14 and (audit.safe_mu == 0 or audit.safe_mu <= Fraction(1, 1000)):
                stats[ch]["frontier"] += severity
    return stats


def tournament(audits: list[ExposureAudit]):
    rank = {name: idx for idx, name in enumerate(CHANNELS)}
    stats = channel_exposure(audits)
    n = len(CHANNELS)

    def retention_cmp(a: str, b: str) -> int:
        av = RETENTION[a]
        bv = RETENTION[b]
        awins = sum(x > y for x, y in zip(av, bv))
        bwins = sum(y > x for x, y in zip(av, bv))
        if awins != bwins:
            return 1 if awins > bwins else -1
        return 1 if rank[a] < rank[b] else -1

    matrix: dict[tuple[str, str], bool] = {}
    for a, b in combinations(CHANNELS, 2):
        only_a = sum(x.severity for x in audits if a in x.channels and b not in x.channels)
        only_b = sum(x.severity for x in audits if b in x.channels and a not in x.channels)
        if abs(only_a - only_b) > 1e-9:
            awins = only_a > only_b
        else:
            awins = retention_cmp(a, b) > 0
        matrix[(a, b)] = awins

    def edge(i: int, j: int) -> bool:
        if i == j:
            raise ValueError("loop")
        a = CHANNELS[i]
        b = CHANNELS[j]
        if i < j:
            return matrix[(a, b)]
        return not matrix[(b, a)]

    outdeg = [sum(1 for j in range(n) if i != j and edge(i, j)) for i in range(n)]
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if edge(i, j) and edge(j, k) and edge(k, i):
            c3 += 1
        if edge(i, k) and edge(k, j) and edge(j, i):
            c3 += 1

    reach = [[False] * n for _ in range(n)]
    for i in range(n):
        dq = deque([i])
        reach[i][i] = True
        while dq:
            u = dq.popleft()
            for v in range(n):
                if u != v and edge(u, v) and not reach[i][v]:
                    reach[i][v] = True
                    dq.append(v)
    seen = [False] * n
    sccs = []
    for i in range(n):
        if seen[i]:
            continue
        comp = [j for j in range(n) if reach[i][j] and reach[j][i]]
        for j in comp:
            seen[j] = True
        sccs.append(tuple(CHANNELS[j] for j in comp))

    return stats, {
        "score_hist": dict(sorted(Counter(outdeg).items())),
        "scores": tuple(sorted(zip(outdeg, CHANNELS), reverse=True)),
        "c3": c3,
        "sccs": tuple(sccs),
    }


def print_header(args, rows: list, audits: list[ExposureAudit]) -> None:
    print("LRC14 EXPOSURE-POSET CREATIVE PROOF PASS")
    print("=" * 78)
    print(
        f"rows={len(rows)}, max_degree={args.max_degree}, negative_tol={args.negative_tol}, "
        f"single_limit={args.single_limit}, two_swap_limit={args.two_swap_limit}, "
        f"alias_depth={args.alias_depth}, lcm_tail_max={args.lcm_tail_max}"
    )
    print()
    print("[0] Quotient declaration")
    print("  considered vertices:")
    print("    runners, gaps, fixed circle sections, wall crossings, residues,")
    print("    safe components, denominator events, Fourier modes, packet families,")
    print("    and proof obligations.")
    print("  chosen vertices:")
    print("    certificate/exposure channels.")
    print("  preserved LRC predicate:")
    print("    q-threshold, open-vs-boundary status, safe-component size,")
    print("    divisor-curried Fejer certificate, packet labels, and K33/petal")
    print("    side channels.")
    print("  destroyed information:")
    print("    full endpoint-owner incidence and exact family quantifiers; this is")
    print("    a theorem-target finder, not a formal proof.")
    print("  challenged assumption:")
    print("    counting easy exposed rows is not enough; proof channels must expose")
    print("    the low-margin q>=14 frontier without erasing labels.")
    print()
    zero = [x for x in audits if x.safe_mu == 0]
    positive = [x for x in audits if x.safe_mu > 0]
    fejer_hits = [x for x in positive if "FEJER_PSD_DUAL" in x.channels]
    print("[1] Exposure census")
    print(f"  zero-safe rows={len(zero)}: {[x.name for x in zero[:8]]}")
    print(f"  positive-safe rows={len(positive)}")
    print(f"  positive rows with Fejer PSD dual={len(fejer_hits)}")
    print(
        "  unexposed rows="
        f"{sum(1 for x in audits if x.channels == frozenset({'UNEXPOSED_SOURCE_KERNEL'}))}"
    )
    print(
        "  rows carrying UNEXPOSED_SOURCE_KERNEL="
        f"{sum(1 for x in audits if 'UNEXPOSED_SOURCE_KERNEL' in x.channels)}"
    )
    print()


def print_channel_table(stats: dict[str, dict[str, float]]) -> None:
    print("[2] Channel exposure table")
    print(f"  {'channel':28s} {'rows':>6s} {'weight':>10s} {'frontier':>10s} {'unique':>10s}")
    for ch in CHANNELS:
        s = stats[ch]
        print(
            f"  {ch:28s} {int(s['rows']):6d} {s['weight']:10.3f} "
            f"{s['frontier']:10.3f} {s['unique_weight']:10.3f}"
        )
    print()


def print_frontier(audits: list[ExposureAudit], limit: int) -> None:
    frontier = sorted(
        [x for x in audits if x.qdiv >= 14],
        key=lambda x: (-x.severity, x.safe_mu, -(x.fejer_degree or 0), x.name),
    )
    print(f"[3] Highest-severity q>=14 rows (top {limit})")
    print(
        f"  {'row':42s} {'q':>3s} {'mu':>11s} {'deg':>5s} "
        f"{'qform':>13s} channels"
    )
    for x in frontier[:limit]:
        print(
            f"  {x.name[:42]:42s} {x.qdiv:3d} {fmt_fraction(x.safe_mu):>11s} "
            f"{str(x.fejer_degree) if x.fejer_degree is not None else '-':>5s} "
            f"{fmt_float(x.fejer_qform):>13s} {','.join(sorted(x.channels))}"
        )
    print()


def print_tournament(fp: dict[str, object]) -> None:
    print("[4] Exposure-channel tournament")
    print("  pair observable:")
    print("    severity-weighted rows exposed by one channel and not the other;")
    print("    ties use retention of q, boundary, open interval, harmonic dual,")
    print("    packet label, state label, and anti-scalarization.")
    print("  tie Hamiltonian path:")
    print("    " + " -> ".join(CHANNELS))
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3_cycles={fp['c3']}")
    print(f"  sccs={fp['sccs']}")
    print("  scores:")
    for score, name in fp["scores"]:
        print(f"    {name:28s} {score}")
    print()


def print_theorem_target(audits: list[ExposureAudit]) -> None:
    zero = [x for x in audits if x.safe_mu == 0]
    non_ap_zero = [x for x in zero if "AP_GW_TAUT_BOUNDARY" not in x.channels]
    hard = [x for x in audits if "HARD_FEJER_MARGIN" in x.channels]
    print("[5] Creative proof target")
    print("  Proposed exposure theorem:")
    print("    Every primitive LRC14 source packet either has a q-witness, is one of")
    print("    the AP/GW taut boundary atoms, exposes a positive Haar bridge with a")
    print("    familywise Fejer interval certificate, or carries a labelled K33/petal")
    print("    handoff whose state-lift theorem discharges it.")
    print()
    print("  In this audit:")
    print(f"    non-AP/GW zero-safe rows={len(non_ap_zero)}")
    print(f"    hard Fejer-margin rows={len(hard)}")
    print("  The missing global lemma is not another scalar invariant.  It is a")
    print("  no-hidden-kernel lemma: once the packet labels are retained, the only")
    print("  zero-open kernel is AP/GW, and every positive packet has a uniform")
    print("  interval-certificate assembly path.")
    print()
    print("  Suggested next formal sublemma:")
    print("    q>=14 and non-AP/GW source packets have either a positive endpoint")
    print("    bridge whose largest component supports a Fejer certificate before")
    print("    the first invisible speed, or a named K33/C27 state-lift label.")
    print("    The hard rows below show where the quantified family proof should start:")
    for x in sorted(hard, key=lambda y: (-(y.fejer_degree or 0), y.safe_mu, y.name))[:8]:
        print(
            f"      {x.name}: q={x.qdiv}, mu={fmt_fraction(x.safe_mu)}, "
            f"degree={x.fejer_degree}, qform={fmt_float(x.fejer_qform)}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=120)
    parser.add_argument("--two-swap-limit", type=int, default=30)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--max-degree", type=int, default=320)
    parser.add_argument("--negative-tol", type=float, default=1e-8)
    parser.add_argument("--frontier-limit", type=int, default=24)
    parser.add_argument("--no-hard-tail", action="store_true")
    args = parser.parse_args()

    rows = build_rows(
        single_limit=args.single_limit,
        two_swap_limit=args.two_swap_limit,
        alias_depth=args.alias_depth,
        lcm_tail_max=args.lcm_tail_max,
        include_hard_tail=not args.no_hard_tail,
    )
    audits = audit_rows(rows, args.max_degree, args.negative_tol)
    stats, fp = tournament(audits)

    print_header(args, rows, audits)
    print_channel_table(stats)
    print_frontier(audits, args.frontier_limit)
    print_tournament(fp)
    print_theorem_target(audits)


if __name__ == "__main__":
    main()
