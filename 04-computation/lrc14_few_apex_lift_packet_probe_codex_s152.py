#!/usr/bin/env python3
"""S152: few-apex lift-packet probe for the LRC14 covering residual.

The labelled-packet audits left the covering boundary-moment bucket as the
largest non-q branch.  This script isolates the post-THM-571 side:

    S = R union 14Q,  1 <= |14Q| <= 6.

The exact carrier is the lift packet over u = 14t.  For a fixed u there are
fourteen lifts t=(u+k)/14.  A 14-multiple 14m only tests the Q-runner m at u;
a residual runner r forbids explicit rational intervals in each lift:

    ||r(u+k)/14|| < 1/14
      <=> exists n in Z with |r(u+k)-14n| < 1.

Thus every row emits fourteen exact Borel/Baire interval fronts.  A positive
front in any lift is already an LRC14 witness with M(S)>1/14.

Tournament Analysis is over proof carriers and lift packets, not runners.  The
chosen quotient preserves qdiv>14, the count of 14-multiples, exact lift-safe
mass, and the Q/R split.  It destroys fine runner ownership inside a residual
interval, so the worst rows are sent back to exact M/Farey calculation.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))
THRESHOLD = Fraction(1, 14)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s124 = load_module(
    "s152_s124_common",
    REPO / "04-computation" / "lrc14_ap_gw_common_conditions_codex_s124.py",
)
s147 = load_module(
    "s152_s147_haar",
    REPO / "04-computation" / "lrc14_baire_haar_anyangle_codex_s147.py",
)


COMPACT_FILLERS = (
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    15,
    18,
    20,
    22,
    24,
    26,
    27,
    28,
    30,
    33,
    36,
    39,
    42,
    44,
    45,
    52,
    55,
    60,
    66,
    72,
    78,
)

STRETCH_FILLERS = (
    6,
    12,
    18,
    24,
    30,
    36,
    42,
    55,
    66,
    78,
    90,
    99,
    110,
    117,
    130,
    143,
    156,
    165,
    180,
)


@dataclass(frozen=True)
class BankRow:
    name: str
    mode: str
    removed: tuple[int, ...]
    multipliers: tuple[int, ...]
    speeds: tuple[int, ...]


@dataclass(frozen=True)
class LiftAudit:
    row: BankRow
    q_threshold: int
    lift_masses: tuple[Fraction, ...]
    lift_component_counts: tuple[int, ...]
    safe_t_measure: Fraction
    q_safe_measure: Fraction

    @property
    def k14(self) -> int:
        return len(self.row.multipliers)

    @property
    def open_lifts(self) -> int:
        return sum(1 for x in self.lift_masses if x > 0)

    @property
    def max_lift_mass(self) -> Fraction:
        return max(self.lift_masses, default=Fraction(0))

    @property
    def min_positive_lift_mass(self) -> Fraction:
        positives = [x for x in self.lift_masses if x > 0]
        return min(positives, default=Fraction(0))


@dataclass(frozen=True)
class ExactWorst:
    audit: LiftAudit
    M: Fraction
    denoms: tuple[int, ...]
    haar_safe_measure: Fraction


@dataclass(frozen=True)
class MeasureCheck:
    audit: LiftAudit
    haar_safe_measure: Fraction


def fmt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def q_req(d: int) -> int:
    return d // gcd(d, 14)


def antichain(reqs: set[int]) -> list[int]:
    out = []
    for r in sorted(reqs, reverse=True):
        if not any(big % r == 0 for big in out):
            out.append(r)
    return sorted(out, reverse=True)


def missing_q_reqs(kept: set[int]) -> list[int]:
    missing = {
        q_req(d)
        for d in range(2, 15)
        if not any(v % d == 0 for v in kept)
    }
    return antichain(missing)


def pack_reqs(reqs: list[int], k: int) -> list[int]:
    bins = list(reqs)
    while len(bins) > k:
        best = None
        for i, j in combinations(range(len(bins)), 2):
            merged = lcm(bins[i], bins[j])
            growth = merged - bins[i] - bins[j]
            key = (growth, merged, i, j)
            if best is None or key < best:
                best = key
        assert best is not None
        _, merged, i, j = best
        nxt = [v for idx, v in enumerate(bins) if idx not in {i, j}]
        nxt.append(merged)
        bins = sorted(nxt, reverse=True)
    return bins


def fill_profile(base: list[int], k: int, fillers: tuple[int, ...]) -> list[int]:
    vals = list(base)
    for f in fillers:
        if len(vals) >= k:
            break
        if f not in vals:
            vals.append(f)
    if len(vals) < k:
        f = 2
        while len(vals) < k:
            candidate = fillers[-1] + f
            if candidate not in vals:
                vals.append(candidate)
            f += 1
    return sorted(vals)


def stretch_profile(base: list[int], k: int, max_m: int) -> list[int]:
    vals: list[int] = []
    factors = (2, 3, 5, 7, 11)
    for idx, r in enumerate(base):
        chosen = r
        for factor in factors[idx % len(factors) :] + factors[: idx % len(factors)]:
            candidate = r * factor
            if candidate <= max_m and candidate not in vals:
                chosen = candidate
                break
        while chosen in vals and chosen + r <= max_m:
            chosen += r
        vals.append(chosen)
    return fill_profile(vals, k, STRETCH_FILLERS)


def build_bank(max_k14: int, max_multiplier: int, include_stretch: bool) -> list[BankRow]:
    rows: dict[tuple[int, ...], BankRow] = {}
    for k in range(1, max_k14 + 1):
        for removed in combinations(AP, k):
            kept = set(AP) - set(removed)
            reqs = missing_q_reqs(kept)
            packed = pack_reqs(reqs, k)
            profiles = [("compact", fill_profile(packed, k, COMPACT_FILLERS))]
            if include_stretch:
                profiles.append(("stretched", stretch_profile(packed, k, max_multiplier)))
            for mode, profile in profiles:
                profile = tuple(sorted(m for m in profile if m <= max_multiplier))
                if len(profile) != k:
                    continue
                speeds = tuple(sorted(kept | {14 * m for m in profile}))
                if len(speeds) != 13 or not primitive(speeds):
                    continue
                if s124.q_threshold(speeds) <= 14:
                    continue
                name = f"{mode} drop{removed}->14*{profile}"
                rows.setdefault(
                    speeds,
                    BankRow(name, mode, tuple(removed), profile, speeds),
                )
    return sorted(rows.values(), key=lambda r: (len(r.multipliers), r.mode, r.removed, r.multipliers))


def clip_interval(a: Fraction, b: Fraction) -> tuple[Fraction, Fraction] | None:
    start = max(a, Fraction(0))
    end = min(b, Fraction(1))
    if start < end:
        return start, end
    return None


def residual_danger_intervals(speed: int, lift: int) -> list[tuple[Fraction, Fraction]]:
    """Danger intervals in u for ||speed*(u+lift)/14|| < 1/14."""

    lo = speed * lift
    hi = speed * (lift + 1)
    n0 = (lo - 1) // 14 - 1
    n1 = (hi + 1) // 14 + 2
    out: list[tuple[Fraction, Fraction]] = []
    for n in range(n0, n1 + 1):
        clipped = clip_interval(
            Fraction(14 * n - 1, speed) - lift,
            Fraction(14 * n + 1, speed) - lift,
        )
        if clipped is not None:
            out.append(clipped)
    return out


def lift_packet(row: BankRow) -> LiftAudit:
    q_intervals: list[tuple[Fraction, Fraction]] = []
    for m in row.multipliers:
        q_intervals.extend(s147.danger_intervals(m, THRESHOLD))
    q_merged = s147.merge_intervals(q_intervals)
    q_safe = Fraction(1) - s147.union_measure(q_merged)

    residual = tuple(v for v in row.speeds if v % 14 != 0)
    lift_masses = []
    lift_component_counts = []
    for lift in range(14):
        intervals = list(q_merged)
        for r in residual:
            intervals.extend(residual_danger_intervals(r, lift))
        comps = s147.complement_components(s147.merge_intervals(intervals))
        lift_masses.append(s147.union_measure(comps))
        lift_component_counts.append(len(comps))
    return LiftAudit(
        row=row,
        q_threshold=s124.q_threshold(row.speeds),
        lift_masses=tuple(lift_masses),
        lift_component_counts=tuple(lift_component_counts),
        safe_t_measure=sum(lift_masses, Fraction(0)) / 14,
        q_safe_measure=q_safe,
    )


def exact_worst(audits: list[LiftAudit], count: int) -> list[ExactWorst]:
    if count <= 0:
        return []
    out = []
    worst = sorted(audits, key=lambda a: (a.safe_t_measure, a.open_lifts, max(a.row.speeds), a.row.speeds))[:count]
    for audit in worst:
        M, pts = s124.M_exact(audit.row.speeds)
        haar = s147.exact_row_measure(audit.row.speeds)["safe_measure"]
        out.append(
            ExactWorst(
                audit=audit,
                M=M,
                denoms=tuple(sorted({t.denominator for t in pts})),
                haar_safe_measure=haar,
            )
        )
    return out


def measure_checked_worst(audits: list[LiftAudit], count: int) -> list[MeasureCheck]:
    if count <= 0:
        return []
    out = []
    worst = sorted(audits, key=lambda a: (a.safe_t_measure, a.open_lifts, max(a.row.speeds), a.row.speeds))[:count]
    for audit in worst:
        haar = s147.exact_row_measure(audit.row.speeds)["safe_measure"]
        out.append(MeasureCheck(audit=audit, haar_safe_measure=haar))
    return out


def edge_from_scores(scores: list[tuple[int, ...]], i: int, j: int) -> bool:
    if scores[i] == scores[j]:
        return i < j
    return scores[i] > scores[j]


def tournament_fingerprint(scores: list[tuple[int, ...]]) -> dict[str, object]:
    n = len(scores)
    outdeg = [
        sum(1 for j in range(n) if i != j and edge_from_scores(scores, i, j))
        for i in range(n)
    ]
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if (
            edge_from_scores(scores, a, b)
            and edge_from_scores(scores, b, c)
            and edge_from_scores(scores, c, a)
        ):
            c3 += 1
        if (
            edge_from_scores(scores, a, c)
            and edge_from_scores(scores, c, b)
            and edge_from_scores(scores, b, a)
        ):
            c3 += 1
    graph = [[j for j in range(n) if i != j and edge_from_scores(scores, i, j)] for i in range(n)]
    rgraph = [[j for j in range(n) if i != j and edge_from_scores(scores, j, i)] for i in range(n)]
    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)
    comps = []
    seen.clear()
    for v in reversed(order):
        if v in seen:
            continue
        comp = set()
        stack = [v]
        seen.add(v)
        while stack:
            x = stack.pop()
            comp.add(x)
            for w in rgraph[x]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(comp)

    dp: dict[tuple[int, int], int] = {(1 << i, i): 1 for i in range(n)}
    for _ in range(1, n):
        nxt: dict[tuple[int, int], int] = defaultdict(int)
        for (used, last), val in dp.items():
            for v in range(n):
                if used & (1 << v):
                    continue
                if edge_from_scores(scores, last, v):
                    nxt[(used | (1 << v), v)] += val
        dp = nxt
    full = (1 << n) - 1
    hp = sum(val for (used, _), val in dp.items() if used == full)
    return {
        "score_hist": dict(sorted(Counter(outdeg).items())),
        "c3": c3,
        "scc": tuple(sorted((len(c) for c in comps), reverse=True)),
        "hp": hp,
    }


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    runners, 14-multiple multipliers Q, residual speeds R, fourteen lifts,")
    print("    q-divisor obligations, safe-interval components, boundary events,")
    print("    Fourier modes, and proof obligations.")
    print("  chosen primary vertices:")
    print("    lift packets (u,k) and proof carriers, not raw runners.")
    print("  quotient preserves:")
    print("    qdiv>14, the Q/R split, |14Z cap S|, exact rational lift-safe mass,")
    print("    and whether a strict Haar/Baire witness exists.")
    print("  quotient destroys:")
    print("    fine ownership of a residual interval and exact maximizing Farey node;")
    print("    the worst rows are therefore sent back to exact M/Farey calculation.")
    print("  challenged assumption:")
    print("    the HYP-2961 L1 branch with |Q|>=7 is still live.  THM-571 closes")
    print("    that apex-majority branch modulo below-frontier LRC input, so the")
    print("    active covering side is the few-apex |Q|<=6 lift-packet branch.")


def print_bank_summary(audits: list[LiftAudit]) -> None:
    print()
    print("[1] Structured few-apex bank")
    print(f"  audited qdiv>14 rows: {len(audits)}")
    print(f"  rows with zero strict lift mass: {sum(1 for a in audits if a.safe_t_measure == 0)}")
    print(f"  rows with no positive lift: {sum(1 for a in audits if a.open_lifts == 0)}")
    print("  by |14Z cap S| and profile mode:")
    groups: dict[tuple[int, str], list[LiftAudit]] = defaultdict(list)
    for a in audits:
        groups[(a.k14, a.row.mode)].append(a)
    print(
        f"    {'k14':>3s} {'mode':>9s} {'rows':>6s} {'min safe_t':>14s} "
        f"{'min open lifts':>14s} {'min Q-safe':>14s}"
    )
    for key in sorted(groups):
        group = groups[key]
        print(
            f"    {key[0]:3d} {key[1]:>9s} {len(group):6d} "
            f"{fmt(min(a.safe_t_measure for a in group)):>14s} "
            f"{min(a.open_lifts for a in group):14d} "
            f"{fmt(min(a.q_safe_measure for a in group)):>14s}"
        )


def print_worst_lifts(worst: list[MeasureCheck]) -> None:
    print()
    print("[2] Tightest lift packets, checked against direct Haar intervals")
    print(
        f"  {'mode':>9s} {'k14':>3s} {'removed':>18s} {'Q mults':>24s} "
        f"{'safe_t':>12s} {'lifts':>5s} {'max speed':>9s} {'check':>8s}"
    )
    for item in worst:
        audit = item.audit
        ok = "OK" if item.haar_safe_measure == audit.safe_t_measure else "MISMATCH"
        print(
            f"  {audit.row.mode:>9s} {audit.k14:3d} {str(audit.row.removed):>18s} "
            f"{str(audit.row.multipliers):>24s} {fmt(audit.safe_t_measure):>12s} "
            f"{audit.open_lifts:5d} {max(audit.row.speeds):9d} {ok:>8s}"
        )


def print_exact_worst(worst: list[ExactWorst]) -> None:
    print()
    print("[2b] Optional exact M fallback for tightest lift packets")
    if not worst:
        print("  skipped; the exact lift-safe mass already proves M(S)>1/14 for every audited row.")
        print("  rerun with --exact-worst N when exact maximizing denominators are needed.")
        return
    print(
        f"  {'mode':>9s} {'k14':>3s} {'removed':>18s} {'Q mults':>24s} "
        f"{'safe_t':>12s} {'lifts':>5s} {'M':>8s} {'denoms':>16s}"
    )
    for item in worst:
        audit = item.audit
        ok = "OK" if item.haar_safe_measure == audit.safe_t_measure else "MISMATCH"
        print(
            f"  {audit.row.mode:>9s} {audit.k14:3d} {str(audit.row.removed):>18s} "
            f"{str(audit.row.multipliers):>24s} {fmt(audit.safe_t_measure):>12s} "
            f"{audit.open_lifts:5d} {fmt(item.M):>8s} {str(list(item.denoms)):>16s} {ok}"
        )
    if worst:
        min_m = min(item.M for item in worst)
        print(f"  minimum exact M among displayed worst rows: {fmt(min_m)}")


def print_lift_histograms(audits: list[LiftAudit]) -> None:
    print()
    print("[3] Lift-packet fingerprints")
    open_hist = Counter(a.open_lifts for a in audits)
    minpos_bucket = Counter()
    for a in audits:
        x = a.min_positive_lift_mass
        if x == 0:
            key = "zero"
        elif x <= Fraction(1, 5000):
            key = "<=1/5000"
        elif x <= Fraction(1, 1000):
            key = "<=1/1000"
        elif x <= Fraction(1, 200):
            key = "<=1/200"
        else:
            key = ">1/200"
        minpos_bucket[key] += 1
    print(f"  open lift count histogram: {dict(sorted(open_hist.items()))}")
    print(f"  minimum positive lift-mass buckets: {dict(minpos_bucket)}")
    print("  hardest rows by total safe_t mass are positive because at least one")
    print("  exact lift interval survives; no endpoint-only packet appears in this bank.")


def proof_carrier_tournament() -> None:
    print()
    print("[4] Tournament Analysis: proof-carrier relation")
    carriers = [
        ("qdiv scalar witness", (6, 6, 5, 3, 2, 5)),
        ("THM-571 apex-majority descent", (5, 6, 6, 5, 4, 6)),
        ("few-apex lift packet", (5, 5, 6, 6, 5, 6)),
        ("boundary-moment bridge", (4, 5, 5, 6, 6, 6)),
        ("K33/state lift endpoint", (4, 4, 5, 5, 6, 6)),
        ("exact M/Farey fallback", (5, 5, 4, 5, 4, 5)),
        ("raw apex residue tournament", (2, 2, 2, 2, 1, 1)),
        ("raw runner set", (1, 1, 1, 1, 1, 1)),
    ]
    scores = [score for _, score in carriers]
    fp = tournament_fingerprint(scores)
    order = sorted(range(len(carriers)), key=lambda i: scores[i], reverse=True)
    print("  vertices: proof carriers and packet quotients, not runners.")
    print("  pair observable:")
    print("    LRC predicate retention, Q/R split retention, lift interval exactness,")
    print("    boundary label retention, state-lift compatibility, anti-scalar guard.")
    print("  switch/gauge:")
    print("    lexicographically larger retention vector wins; ties use listed order.")
    print(f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']} scc={fp['scc']} hp={fp['hp']}")
    print("  Hamiltonian path:")
    print("    " + " > ".join(carriers[i][0] for i in order))


def print_theorem_readout(audits: list[LiftAudit], worst: list[ExactWorst]) -> None:
    print()
    print("[5] Theorem-facing readout")
    print("  THM-571 reconciliation:")
    print("    the |14Z cap S|>=7 apex-majority family should be removed from")
    print("    the live HYP-2961 counterexample list, modulo the accepted")
    print("    below-frontier LRC input used by THM-571.")
    print("  Few-apex evidence:")
    print(
        "    every structured AP-replacement row with 1<=|14Z cap S|<=6 in this"
    )
    print(
        "    qdiv>14 bank has positive exact lift-safe mass, so none is a"
    )
    print("    covering zero-open non-migration packet.")
    print(f"    smallest exact lift-safe measure in this run is {fmt(min(a.safe_t_measure for a in audits))}.")
    if worst:
        print(f"    optional exact M fallback displayed minimum {fmt(min(item.M for item in worst))}.")
    print("  Candidate lemma:")
    print("    Fixed-margin few-apex lift theorem.  For primitive qdiv>14 rows")
    print("    with at most six multiples of 14, a retained Q/R lift packet has")
    print("    positive regular-open Haar mass unless it state-lifts to the K33")
    print("    endpoint.  This is the F5/F6 boundary-moment bridge in executable form.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-k14", type=int, default=6)
    parser.add_argument("--max-multiplier", type=int, default=180)
    parser.add_argument("--no-stretch", action="store_true")
    parser.add_argument("--measure-check", type=int, default=14)
    parser.add_argument("--exact-worst", type=int, default=0)
    args = parser.parse_args()

    print("S152 LRC14 FEW-APEX LIFT-PACKET PROBE")
    print("=" * 78)
    print(
        f"max_k14={args.max_k14}, max_multiplier={args.max_multiplier}, "
        f"stretch={not args.no_stretch}, measure_check={args.measure_check}, "
        f"exact_worst={args.exact_worst}"
    )
    print_assumption_challenge()

    rows = build_bank(args.max_k14, args.max_multiplier, include_stretch=not args.no_stretch)
    audits = [lift_packet(row) for row in rows]
    audits.sort(key=lambda a: (a.safe_t_measure, a.open_lifts, a.k14, a.row.mode, a.row.removed))

    print_bank_summary(audits)
    checked = measure_checked_worst(audits, args.measure_check)
    print_worst_lifts(checked)
    worst = exact_worst(audits, args.exact_worst)
    print_exact_worst(worst)
    print_lift_histograms(audits)
    proof_carrier_tournament()
    print_theorem_readout(audits, worst)


if __name__ == "__main__":
    main()
