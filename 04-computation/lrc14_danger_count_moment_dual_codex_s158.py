#!/usr/bin/env python3
"""S158: danger-count moment dual route for LRC14.

This is intentionally different from the recent few-apex lift-packet and
apex-aperture work.  Instead of choosing a denominator-14 apex or lifting by
u=14t, look at the integer-valued danger count

    N_S(t) = #{s in S : ||s t|| < 1/14}.

If S were a strict counterexample with no safe time, then N_S(t) >= 1 for
every t.  Therefore every polynomial P on {0,...,13} with

    P(0)=1,    P(n)<=0 for n=1,...,13

is a valid dual certificate:

    safe_mu(S) = mu(N=0) >= E[P(N)].

The script computes the exact distribution of N by sweeping rational danger
arc endpoints and searches exact factorial-moment polynomial duals.  Positive
dual value proves a positive safe set using only the retained danger-count
moments, not an interval/lift witness.

Tournament Analysis uses moment-dual carriers as vertices.  The pairwise
observable is the proof information retained: cover predicate, count
distribution, low-degree dual strength, AP/GW equality visibility, and
resistance to scalar-only collapse.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb, gcd
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


s147 = load_module(
    "s158_haar_events",
    REPO / "04-computation" / "lrc14_baire_haar_anyangle_codex_s147.py",
)


@dataclass(frozen=True)
class Row:
    name: str
    speeds: tuple[int, ...]
    family: str


@dataclass(frozen=True)
class DualCandidate:
    degree: int
    active: tuple[int, ...]
    coeffs: tuple[Fraction, ...]
    values: tuple[Fraction, ...]

    def expectation(self, moments: tuple[Fraction, ...]) -> Fraction:
        return Fraction(1) + sum(
            self.coeffs[k - 1] * moments[k] for k in range(1, self.degree + 1)
        )


@dataclass(frozen=True)
class RowAudit:
    row: Row
    distribution: tuple[tuple[int, Fraction], ...]
    moments: tuple[Fraction, ...]
    best_by_degree: tuple[tuple[int, Fraction, tuple[int, ...]], ...]
    first_positive_degree: int | None

    @property
    def safe_measure(self) -> Fraction:
        return dict(self.distribution).get(0, Fraction(0))

    @property
    def mean_danger(self) -> Fraction:
        return self.moments[1]

    @property
    def max_count(self) -> int:
        return max((n for n, mass in self.distribution if mass), default=0)


def fmt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def danger_count_distribution(speeds: tuple[int, ...]) -> tuple[tuple[int, Fraction], ...]:
    """Exact distribution of N_S(t) by a rational endpoint sweep."""

    events: dict[Fraction, int] = {Fraction(0): 0, Fraction(1): 0}
    for speed in speeds:
        for lo, hi in s147.danger_intervals(speed, THRESHOLD):
            if hi <= lo:
                continue
            events[lo] = events.get(lo, 0) + 1
            events[hi] = events.get(hi, 0) - 1

    pts = sorted(events)
    active = 0
    dist: dict[int, Fraction] = {}
    for idx, point in enumerate(pts):
        active += events[point]
        nxt = pts[idx + 1] if idx + 1 < len(pts) else Fraction(1)
        if nxt > point:
            dist[active] = dist.get(active, Fraction(0)) + (nxt - point)
    return tuple(sorted((n, mass) for n, mass in dist.items() if mass))


def binom(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    return comb(n, k)


def factorial_moments(
    distribution: tuple[tuple[int, Fraction], ...], max_degree: int
) -> tuple[Fraction, ...]:
    return tuple(
        sum(Fraction(binom(n, k)) * mass for n, mass in distribution)
        for k in range(max_degree + 1)
    )


def solve_linear(
    matrix: list[list[Fraction]], rhs: list[Fraction]
) -> tuple[Fraction, ...] | None:
    n = len(rhs)
    aug = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    pivot = 0
    for col in range(n):
        row = next((r for r in range(pivot, n) if aug[r][col]), None)
        if row is None:
            return None
        aug[pivot], aug[row] = aug[row], aug[pivot]
        scale = aug[pivot][col]
        aug[pivot] = [x / scale for x in aug[pivot]]
        for r in range(n):
            if r == pivot or not aug[r][col]:
                continue
            factor = aug[r][col]
            aug[r] = [aug[r][j] - factor * aug[pivot][j] for j in range(n + 1)]
        pivot += 1
    return tuple(aug[i][-1] for i in range(n))


def polynomial_values(coeffs: tuple[Fraction, ...], max_count: int) -> tuple[Fraction, ...]:
    values = []
    for n in range(max_count + 1):
        values.append(
            Fraction(1)
            + sum(coeffs[k - 1] * binom(n, k) for k in range(1, len(coeffs) + 1))
        )
    return tuple(values)


def build_dual_library(max_degree: int, max_count: int) -> dict[int, tuple[DualCandidate, ...]]:
    """Enumerate exact LP vertices for each degree.

    The feasible set is independent of S.  A degree-R vertex has R active
    constraints P(n)=0 among n=1,...,max_count.
    """

    library: dict[int, list[DualCandidate]] = {}
    for degree in range(1, max_degree + 1):
        candidates: list[DualCandidate] = []
        for active in combinations(range(1, max_count + 1), degree):
            matrix = [
                [Fraction(binom(n, k)) for k in range(1, degree + 1)]
                for n in active
            ]
            coeffs = solve_linear(matrix, [Fraction(-1)] * degree)
            if coeffs is None:
                continue
            values = polynomial_values(coeffs, max_count)
            if all(values[n] <= 0 for n in range(1, max_count + 1)):
                candidates.append(DualCandidate(degree, active, coeffs, values))
        library[degree] = candidates
    return {degree: tuple(candidates) for degree, candidates in library.items()}


def best_dual(
    library: dict[int, tuple[DualCandidate, ...]],
    degree: int,
    moments: tuple[Fraction, ...],
) -> tuple[Fraction, tuple[int, ...], DualCandidate] | None:
    best: tuple[Fraction, tuple[int, ...], DualCandidate] | None = None
    for candidate in library[degree]:
        value = candidate.expectation(moments)
        if best is None or value > best[0]:
            best = (value, candidate.active, candidate)
    return best


def named_rows() -> list[Row]:
    def replace(holes: tuple[int, ...], adds: tuple[int, ...], name: str, family: str) -> Row:
        speeds = tuple(sorted((set(AP) - set(holes)) | set(adds)))
        return Row(name, speeds, family)

    return [
        Row("AP", AP, "AP/GW equality"),
        Row("GW 12->24", tuple(list(range(1, 12)) + [13, 24]), "AP/GW equality"),
        replace((12,), (36,), "near/K33 12->36", "K33 state-lift"),
        replace((10,), (20,), "petal 10->20", "unit petal"),
        replace((13,), (26,), "petal 13->26", "unit petal"),
        replace((10, 12), (20, 24), "P10+GW splice", "two-block splice"),
        replace((10, 12), (20, 36), "P10+K33 splice", "two-block K33 splice"),
        replace((12,), (84,), "covering 12->84", "covering comb"),
        replace((12,), (168,), "covering 12->168", "covering comb"),
        replace((6,), (14,), "few-apex 6->14", "few-apex comb"),
        replace((6,), (28,), "few-apex 6->28", "few-apex comb"),
        replace((12,), (14 * 12,), "few-apex 12->168", "few-apex comb"),
    ]


def one_swap_rows(limit: int) -> list[Row]:
    rows: list[Row] = []
    for hole in AP:
        for add in range(14, limit + 1):
            if add in AP or add == hole:
                continue
            speeds = tuple(sorted((set(AP) - {hole}) | {add}))
            if not primitive(speeds):
                continue
            rows.append(Row(f"drop({hole})->{add}", speeds, "one-swap AP bank"))
    return rows


def audit_row(
    row: Row,
    library: dict[int, tuple[DualCandidate, ...]],
    max_degree: int,
) -> RowAudit:
    distribution = danger_count_distribution(row.speeds)
    moments = factorial_moments(distribution, max_degree)
    best_rows = []
    first_positive: int | None = None
    for degree in range(1, max_degree + 1):
        best = best_dual(library, degree, moments)
        if best is None:
            continue
        value, active, _ = best
        best_rows.append((degree, value, active))
        if first_positive is None and value > 0:
            first_positive = degree
    return RowAudit(row, distribution, moments, tuple(best_rows), first_positive)


def distribution_text(distribution: tuple[tuple[int, Fraction], ...], limit: int = 6) -> str:
    shown = sorted(distribution)[:limit]
    chunks = [f"{n}:{fmt(mass)}" for n, mass in shown]
    if len(distribution) > limit:
        chunks.append("...")
    return "{" + ", ".join(chunks) + "}"


def print_scope() -> None:
    print("[0] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    runners, danger arcs, exact endpoints, safe intervals, lift packets,")
    print("    danger-count states N=0..13, factorial moments, polynomial duals,")
    print("    Fourier modes, boundary owners, proof obligations.")
    print("  chosen primary vertices:")
    print("    danger-count states and moment-dual certificates, not runners or lifts.")
    print("  quotient preserves:")
    print("    the cover predicate N(t)>=1, exact Haar distribution of N, AP/GW")
    print("    equality visibility, and low-degree dual proof obligations.")
    print("  quotient destroys:")
    print("    where the safe interval sits, which endpoint owns a wall, and Q/R")
    print("    lift geometry; these must be recovered only after the moment gate.")
    print("  challenged assumption:")
    print("    the next proof has to be an interval/lift construction.  A cover has")
    print("    a purely count-theoretic shadow: N(t)>=1 everywhere.")
    print()


def print_dual_theorem() -> None:
    print("[1] Exact danger-count dual theorem")
    print("  Let N(t)=#{s in S: ||s t||<1/14}.  Then")
    print("    safe_mu(S) = mu(N=0).")
    print("  For any polynomial P on {0,...,13} with")
    print("    P(0)=1 and P(n)<=0 for every n>=1,")
    print("  the pointwise inequality 1_{N=0} >= P(N) gives")
    print("    safe_mu(S) >= E[P(N)].")
    print("  Using the factorial basis P(n)=1+sum y_k*C(n,k),")
    print("  E[P(N)] depends only on the moments E[C(N,k)].")
    print()


def print_named(audits: list[RowAudit], max_degree: int) -> None:
    tail_degree = min(13, max_degree)
    print("[2] Named-row danger-count dual audit")
    print(
        f"  {'row':24s} {'safe_mu':>12s} {'first +deg':>10s} "
        f"{'best d8':>12s} {'best d9':>12s} "
        f"{('best d' + str(tail_degree)):>12s} {'dist prefix'}"
    )
    for audit in audits:
        best = {degree: value for degree, value, _ in audit.best_by_degree}
        first = "-" if audit.first_positive_degree is None else str(audit.first_positive_degree)
        print(
            f"  {audit.row.name:24s} {fmt(audit.safe_measure):>12s} {first:>10s} "
            f"{fmt(best.get(8, Fraction(0))):>12s} "
            f"{fmt(best.get(9, Fraction(0))):>12s} "
            f"{fmt(best.get(tail_degree, Fraction(0))):>12s} "
            f"{distribution_text(audit.distribution)}"
        )
    print()
    print("  readout:")
    print("    AP and GW stay at nonpositive dual value until the full degree-13")
    print("    interpolant, where the best bound is exactly 0.")
    print("    Positive rows often need degree 8 or 9 before the count-only dual")
    print("    sees the safe mass; this is a high-degree count obstruction, not a")
    print("    low-degree second-moment proof.")
    print()


def print_bank_summary(audits: list[RowAudit], max_degree: int) -> None:
    if not audits:
        return
    positive = [a for a in audits if a.safe_measure > 0]
    zero = [a for a in audits if a.safe_measure == 0]
    first_hist: dict[str, int] = {}
    for audit in audits:
        key = "none" if audit.first_positive_degree is None else str(audit.first_positive_degree)
        first_hist[key] = first_hist.get(key, 0) + 1
    d9_cert = 0
    worst_d9: tuple[Fraction, RowAudit] | None = None
    for audit in audits:
        best = {degree: value for degree, value, _ in audit.best_by_degree}
        value = best.get(9, Fraction(0))
        if value > 0:
            d9_cert += 1
        if audit.safe_measure > 0 and (worst_d9 is None or value < worst_d9[0]):
            worst_d9 = (value, audit)

    print("[3] One-swap AP-bank moment gate")
    print(f"  audited rows: {len(audits)}")
    print(f"  zero safe-measure rows: {len(zero)}")
    print(f"  positive safe-measure rows: {len(positive)}")
    print(f"  first positive dual-degree histogram: {dict(sorted(first_hist.items()))}")
    print(f"  rows certified by degree <=9 dual: {d9_cert}")
    if worst_d9 is not None:
        value, audit = worst_d9
        print(
            "  worst degree-9 lower bound among positive rows: "
            f"{fmt(value)} at {audit.row.name}, safe_mu={fmt(audit.safe_measure)}"
        )
    tight = sorted(positive, key=lambda a: (a.safe_measure, a.row.name))[:8]
    print("  smallest positive one-swap rows:")
    for audit in tight:
        best = {degree: value for degree, value, _ in audit.best_by_degree}
        first = "-" if audit.first_positive_degree is None else str(audit.first_positive_degree)
        print(
            f"    {audit.row.name:18s} safe={fmt(audit.safe_measure):>10s} "
            f"first={first:>2s} d9={fmt(best.get(9, Fraction(0))):>12s} "
            f"d13={fmt(best.get(min(13, max_degree), Fraction(0))):>10s}"
        )
    print()


def print_active_patterns(audits: list[RowAudit]) -> None:
    print("[4] Active-set fingerprints")
    selected = [a for a in audits if a.row.name in {"AP", "GW 12->24", "near/K33 12->36", "covering 12->84"}]
    for audit in selected:
        best = {degree: (value, active) for degree, value, active in audit.best_by_degree}
        print(f"  {audit.row.name}:")
        for degree in (8, 9, 13):
            if degree in best:
                value, active = best[degree]
                print(f"    degree {degree:2d}: value={fmt(value):>12s} active={active}")
    print()
    print("  The active zero sets are not endpoint labels.  They are count-level")
    print("  atoms: which danger multiplicities the dual polynomial must touch.")
    print("  AP/GW are invisible to low-degree count duals for the same reason they")
    print("  are equality atoms: their N-distribution keeps mass in many count levels.")
    print()


def tournament_analysis() -> None:
    print("[5] Tournament Analysis: proof-carrier relation")
    carriers = [
        ("full danger-count distribution", (6, 6, 6, 4, 5, 6)),
        ("degree-13 interpolating dual", (6, 5, 6, 3, 4, 5)),
        ("degree-9 moment dual", (5, 4, 5, 5, 5, 5)),
        ("gK8/Delsarte sector moment", (4, 4, 4, 6, 5, 5)),
        ("Bonferroni pair moment", (3, 3, 3, 5, 4, 3)),
        ("raw Haar safe measure", (3, 2, 5, 2, 3, 2)),
        ("endpoint/lift packet", (4, 5, 4, 4, 4, 4)),
        ("raw runner set", (1, 1, 1, 1, 1, 1)),
    ]
    wins = {name: 0 for name, _ in carriers}
    for i, (left_name, left_score) in enumerate(carriers):
        for right_name, right_score in carriers[i + 1 :]:
            if left_score > right_score:
                wins[left_name] += 1
            elif right_score > left_score:
                wins[right_name] += 1
            elif left_name < right_name:
                wins[left_name] += 1
            else:
                wins[right_name] += 1
    hist: dict[int, int] = {}
    for score in wins.values():
        hist[score] = hist.get(score, 0) + 1
    order = sorted(carriers, key=lambda item: (-wins[item[0]], item[0]))
    print("  vertices: proof carriers, not runners.")
    print("  pair observable:")
    print("    cover predicate retention, exact count distribution, low-degree")
    print("    certificate strength, AP/GW equality visibility, finite-atlas fit,")
    print("    and anti-scalarization.")
    print("  switch/gauge:")
    print("    lexicographic retention vector; ties use carrier name order.")
    print(f"  fingerprint: score_hist={dict(sorted(hist.items()))} c3=0 hp=1")
    print("  Hamiltonian path:")
    print("    " + " > ".join(name for name, _ in order))
    print()


def print_readout(max_degree: int) -> None:
    print("[6] Theorem-facing readout")
    print("  New lemma target:")
    print("    Danger-count moment gate.  A strict LRC14 counterexample must have")
    print("    safe_mu=0, hence every feasible count-dual polynomial has nonpositive")
    print("    expectation.  Therefore any packet with a positive degree-R dual bound")
    print("    is discharged without locating the witness interval.")
    print("  What this pass found:")
    print("    The route is not a cheap pair/second-moment proof.  Named hard rows")
    print("    remain invisible through degree 6.  Degree 8/9 begins to separate")
    print("    positive rows from AP/GW; degree 13 recovers exact safe mass.")
    print("  How this differs from HYP-2968:")
    print("    HYP-2968 proves positivity by exact Q/R lift intervals.  This route")
    print("    forgets all lift geometry and asks whether the integer multiplicity")
    print("    distribution alone already forces a safe set.")
    print(f"  live next step: search for a universal degree <= {max_degree} dual family")
    print("    on fixed-margin labelled packets, or prove why AP/GW are the only")
    print("    rows whose danger-count distribution defeats every low-degree dual.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-degree", type=int, default=13)
    parser.add_argument("--one-swap-limit", type=int, default=80)
    parser.add_argument(
        "--skip-bank",
        action="store_true",
        help="only audit named rows; skip the one-swap AP bank",
    )
    args = parser.parse_args()
    if args.max_degree > 13 or args.max_degree < 1:
        raise SystemExit("--max-degree must be between 1 and 13")

    print("S158 LRC14 DANGER-COUNT MOMENT-DUAL PROBE")
    print("=" * 78)
    print(f"max_degree={args.max_degree}, one_swap_limit={args.one_swap_limit}")
    print_scope()
    print_dual_theorem()
    library = build_dual_library(args.max_degree, 13)
    print("[library]")
    print(
        "  feasible dual vertices by degree: "
        + str({degree: len(library[degree]) for degree in range(1, args.max_degree + 1)})
    )
    print()

    named_audits = [audit_row(row, library, args.max_degree) for row in named_rows()]
    print_named(named_audits, args.max_degree)
    if not args.skip_bank:
        bank = one_swap_rows(args.one_swap_limit)
        bank_audits = [audit_row(row, library, args.max_degree) for row in bank]
        print_bank_summary(bank_audits, args.max_degree)
    print_active_patterns(named_audits)
    tournament_analysis()
    print_readout(args.max_degree)


if __name__ == "__main__":
    main()
