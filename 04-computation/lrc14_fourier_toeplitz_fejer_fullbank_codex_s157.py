#!/usr/bin/env python3
"""S157: Fourier-Toeplitz PSD dual route for LRC14.

Let

    C_S(t) = #{v in S : ||v t|| < 1/14},      F_S(t)=C_S(t)-1.

If S were a strict LRC14 counterexample, then the open danger arcs would cover
the circle and F_S(t) would be nonnegative almost everywhere.  Therefore every
finite Toeplitz moment matrix

    T_d(S) = [hat F_S(i-j)]_{0<=i,j<=d}

must be positive semidefinite.  A single vector p with p^*T_d p < 0 is a
dual certificate for a strict safe interval.

This script uses the phase-sensitive Fejer vector

    p_j = exp(-2*pi*i*j*x) / sqrt(d+1),  j=0,...,d,

centered at an exact safe component midpoint x.  Its quadratic form is

    Q_d(x) = c_0 + 2*sum_{k=1}^d (1-k/(d+1))*c_k*cos(2*pi*k*x),

where c_k=hat F_S(k).  Negative Q_d(x) proves T_d is not PSD.  This is not
claimed to be the true first negative eigenvalue; it is an explicit
Toeplitz-vector certificate.

The Fourier coefficients are sparse and "curried" by divisibility:

    c_0 = |S|/7 - 1,
    c_k = sum_{v in S, v|k} sin(pi*(k/v)/7)/(pi*(k/v)),  k>0.

Thus low modes only see the divisor fibers of the row.  That is the concrete
Fourier analogue of the repo's recurring n*2/n+2 recursion and Farey
denominator packet language.

Tournament Analysis uses proof carriers/Fourier packets as vertices, not
runners.  The chosen quotient preserves the nonnegativity necessary condition
and an explicit PSD-failure vector, while destroying exact endpoint ownership.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import cos, gcd, pi, sin
from pathlib import Path
import argparse
import sys
import time


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))
N = 14


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s147 = load_module(
    "s157_haar_baire",
    REPO / "04-computation" / "lrc14_baire_haar_anyangle_codex_s147.py",
)
s2963 = load_module(
    "s157_hyp2963_packets",
    REPO / "04-computation" / "lrc14_labelled_packet_counterexample_classifier_codex_20260624.py",
)


@dataclass(frozen=True)
class Row:
    name: str
    family: str
    speeds: tuple[int, ...]


@dataclass(frozen=True)
class SafeComponent:
    lo: Fraction
    hi: Fraction

    @property
    def width(self) -> Fraction:
        return self.hi - self.lo

    @property
    def midpoint(self) -> Fraction:
        return (self.lo + self.hi) / 2


@dataclass(frozen=True)
class FejerCertificate:
    degree: int
    qform: float
    center: Fraction
    width: Fraction
    safe_measure: Fraction
    component_count: int


@dataclass(frozen=True)
class RowAudit:
    row: Row
    safe_measure: Fraction
    component_count: int
    largest_width: Fraction
    certificate: FejerCertificate | None
    best_degree: int
    best_qform: float
    invisible_speeds_at_first: tuple[int, ...]


PROOF_VERTICES = (
    "exact_nonnegative_function",
    "Toeplitz_PSD_cone",
    "Fejer_center_certificate",
    "divisor_curried_Fourier_fibers",
    "HYP2973_count_distribution",
    "HYP2972_twist_ladder",
    "endpoint_lift_packets",
    "raw_safe_interval",
    "raw_runner_set",
)

PROOF_SCORES: dict[str, tuple[int, int, int, int, int, int]] = {
    # cover predicate, phase sensitivity, certificate concreteness,
    # packet compatibility, exactness path, anti-scalarization
    "exact_nonnegative_function": (9, 9, 4, 5, 9, 8),
    "Toeplitz_PSD_cone": (8, 9, 6, 6, 8, 8),
    "Fejer_center_certificate": (8, 8, 9, 6, 6, 7),
    "divisor_curried_Fourier_fibers": (7, 8, 7, 8, 6, 9),
    "HYP2973_count_distribution": (8, 3, 8, 7, 8, 4),
    "HYP2972_twist_ladder": (8, 6, 9, 7, 7, 7),
    "endpoint_lift_packets": (7, 7, 8, 9, 8, 8),
    "raw_safe_interval": (8, 5, 8, 4, 7, 4),
    "raw_runner_set": (2, 1, 1, 1, 1, 1),
}

TIE_PATH = PROOF_VERTICES


def fmt_fraction(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def fmt_float(x: float) -> str:
    return f"{x:.9g}"


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def row_replace(holes: tuple[int, ...], adds: tuple[int, ...], name: str, family: str) -> Row:
    return Row(name, family, tuple(sorted((set(AP) - set(holes)) | set(adds))))


def named_rows() -> list[Row]:
    return [
        Row("AP", "AP/GW equality", AP),
        Row("GW 12->24", "AP/GW equality", tuple(list(range(1, 12)) + [13, 24])),
        row_replace((12,), (36,), "near/K33 12->36", "K33 state-lift"),
        row_replace((10,), (20,), "petal 10->20", "unit petal"),
        row_replace((13,), (26,), "petal 13->26", "unit petal"),
        row_replace((10, 12), (20, 24), "P10+GW", "two-block splice"),
        row_replace((10, 12), (20, 36), "P10+K33", "two-block K33 splice"),
        row_replace((12,), (84,), "covering 12->84", "covering comb"),
        row_replace((12,), (168,), "covering 12->168", "covering comb"),
        row_replace((6,), (14,), "few-apex 6->14", "few-apex comb"),
        row_replace((6,), (28,), "few-apex 6->28", "few-apex comb"),
    ]


def build_hyp2963_rows(
    single_limit: int,
    two_swap_limit: int,
    alias_depth: int,
    lcm_tail_max: int,
) -> list[Row]:
    rows = s2963.build_bank(single_limit, two_swap_limit, alias_depth, lcm_tail_max)
    return [Row(name, family, tuple(speeds)) for name, family, speeds in rows]


def build_one_swap_rows(limit: int) -> list[Row]:
    out: list[Row] = []
    for hole in AP:
        for add in range(14, limit + 1):
            speeds = tuple(sorted((set(AP) - {hole}) | {add}))
            if len(speeds) == 13 and primitive(speeds):
                out.append(Row(f"drop({hole})->{add}", "one-swap AP bank", speeds))
    return out


def ghat_table(max_degree: int) -> list[float]:
    return [0.0] + [sin(pi * m / 7.0) / (pi * m) for m in range(1, max_degree + 1)]


def fourier_coefficients(
    speeds: tuple[int, ...],
    max_degree: int,
    ghat: list[float],
) -> list[float]:
    coeffs = [0.0 for _ in range(max_degree + 1)]
    coeffs[0] = len(speeds) / 7.0 - 1.0
    for v in speeds:
        for m in range(1, max_degree // v + 1):
            coeffs[m * v] += ghat[m]
    return coeffs


def safe_components(speeds: tuple[int, ...]) -> tuple[Fraction, tuple[SafeComponent, ...]]:
    exact = s147.exact_row_measure(speeds)
    comps = tuple(SafeComponent(lo, hi) for lo, hi in exact["safe_components"])
    return exact["safe_measure"], comps


def first_fejer_certificate(
    row: Row,
    max_degree: int,
    ghat: list[float],
    negative_tol: float,
) -> RowAudit:
    safe_measure, comps = safe_components(row.speeds)
    if not comps:
        return RowAudit(row, safe_measure, 0, Fraction(0), None, 0, len(row.speeds) / 7.0 - 1.0, ())

    largest = max(comps, key=lambda comp: comp.width)
    center = largest.midpoint
    x = float(center)
    coeffs = fourier_coefficients(row.speeds, max_degree, ghat)

    running_sum = 0.0
    running_weighted = 0.0
    best_degree = 0
    best_qform = coeffs[0]
    first: FejerCertificate | None = None

    for degree in range(1, max_degree + 1):
        term = coeffs[degree] * cos(2.0 * pi * degree * x)
        running_sum += term
        running_weighted += degree * term
        qform = coeffs[0] + 2.0 * (running_sum - running_weighted / (degree + 1))
        if qform < best_qform:
            best_degree = degree
            best_qform = qform
        if qform < -negative_tol:
            first = FejerCertificate(
                degree,
                qform,
                center,
                largest.width,
                safe_measure,
                len(comps),
            )
            break

    degree_for_visibility = first.degree if first else max_degree
    invisible = tuple(v for v in row.speeds if v > degree_for_visibility)
    return RowAudit(
        row,
        safe_measure,
        len(comps),
        largest.width,
        first,
        best_degree,
        best_qform,
        invisible,
    )


def audit_rows(
    rows: list[Row],
    max_degree: int,
    negative_tol: float,
    progress_every: int = 0,
) -> list[RowAudit]:
    ghat = ghat_table(max_degree)
    audits: list[RowAudit] = []
    start = time.time()
    for idx, row in enumerate(rows, 1):
        audits.append(first_fejer_certificate(row, max_degree, ghat, negative_tol))
        if progress_every and idx % progress_every == 0:
            hits = sum(a.certificate is not None for a in audits)
            misses = sum(a.safe_measure > 0 and a.certificate is None for a in audits)
            max_degree_seen = max((a.certificate.degree for a in audits if a.certificate), default=0)
            print(
                f"    progress {idx}/{len(rows)} elapsed={time.time()-start:.1f}s "
                f"hits={hits} misses={misses} max_first_degree={max_degree_seen}"
            )
    return audits


def summarize_audits(label: str, audits: list[RowAudit]) -> None:
    zero = [a for a in audits if a.safe_measure == 0]
    positive = [a for a in audits if a.safe_measure > 0]
    hits = [a for a in positive if a.certificate is not None]
    misses = [a for a in positive if a.certificate is None]
    degree_hist = Counter(a.certificate.degree for a in hits if a.certificate)
    bucket_hist = Counter()
    for audit in hits:
        assert audit.certificate
        degree = audit.certificate.degree
        lo = 10 * (degree // 10)
        bucket_hist[f"{lo:03d}-{lo+9:03d}"] += 1

    print(f"[summary] {label}")
    print(f"  rows={len(audits)} zero_safe={len(zero)} positive_safe={len(positive)}")
    print(f"  Fejer PSD-vector hits={len(hits)} misses_at_cap={len(misses)}")
    if hits:
        max_hit = max(hits, key=lambda a: a.certificate.degree if a.certificate else -1)
        min_margin = min(hits, key=lambda a: abs(a.certificate.qform) if a.certificate else 10**9)
        assert max_hit.certificate and min_margin.certificate
        print(
            f"  max_first_degree={max_hit.certificate.degree} "
            f"at {max_hit.row.name} ({max_hit.row.family})"
        )
        print(
            f"  smallest_negative_margin={fmt_float(min_margin.certificate.qform)} "
            f"at {min_margin.row.name}, degree={min_margin.certificate.degree}"
        )
        print(f"  first-degree histogram top={degree_hist.most_common(16)}")
        print(f"  first-degree buckets={dict(sorted(bucket_hist.items()))}")
    if zero:
        print("  zero-safe rows:")
        for audit in zero[:8]:
            print(f"    {audit.row.name} ({audit.row.family})")
    if misses:
        print("  positive rows with no Fejer certificate at cap:")
        for audit in sorted(misses, key=lambda a: (a.safe_measure, a.row.name))[:12]:
            print(
                f"    {audit.row.name:36s} family={audit.row.family:18s} "
                f"safe={fmt_fraction(audit.safe_measure):>10s} "
                f"largest_width={fmt_fraction(audit.largest_width):>10s} "
                f"best=({audit.best_degree},{fmt_float(audit.best_qform)})"
            )
    print()


def print_named_table(audits: list[RowAudit]) -> None:
    print("[2] Named-row Fourier-Toeplitz Fejer certificates")
    print(
        f"  {'row':24s} {'safe':>11s} {'comp':>4s} {'width':>11s} "
        f"{'degree':>6s} {'qform':>14s} {'center':>12s} invisible speeds"
    )
    for audit in audits:
        cert = audit.certificate
        if cert is None:
            degree = "-"
            qform = "-"
            center = "-"
            width = fmt_fraction(audit.largest_width)
        else:
            degree = str(cert.degree)
            qform = fmt_float(cert.qform)
            center = fmt_fraction(cert.center)
            width = fmt_fraction(cert.width)
        print(
            f"  {audit.row.name[:24]:24s} {fmt_fraction(audit.safe_measure):>11s} "
            f"{audit.component_count:4d} {width:>11s} {degree:>6s} "
            f"{qform:>14s} {center:>12s} {audit.invisible_speeds_at_first}"
        )
    print()


def print_hardest(audits: list[RowAudit], limit: int) -> None:
    hits = [a for a in audits if a.certificate is not None]
    hardest = sorted(hits, key=lambda a: (a.certificate.degree, a.row.name), reverse=True)
    print(f"[4] Hardest Fejer certificates (top {limit})")
    print(
        f"  {'row':42s} {'family':16s} {'degree':>6s} {'safe':>11s} "
        f"{'width':>11s} {'qform':>14s} invisible"
    )
    for audit in hardest[:limit]:
        cert = audit.certificate
        assert cert is not None
        print(
            f"  {audit.row.name[:42]:42s} {audit.row.family[:16]:16s} "
            f"{cert.degree:6d} {fmt_fraction(audit.safe_measure):>11s} "
            f"{fmt_fraction(cert.width):>11s} {fmt_float(cert.qform):>14s} "
            f"{audit.invisible_speeds_at_first}"
        )
    print()


def print_fourier_formula() -> None:
    print("[1] Fourier-Toeplitz necessary condition")
    print("  Define C_S(t)=#{v in S : ||v*t||<1/14} and F_S=C_S-1.")
    print("  A strict counterexample would have F_S(t)>=0 almost everywhere.")
    print("  Hence every Toeplitz moment matrix [hat F_S(i-j)] must be PSD.")
    print("  A single vector p with p^*T_d p<0 is a dual certificate for a")
    print("  strict safe interval.")
    print()
    print("  Divisor-curried Fourier coefficients:")
    print("    c_0 = |S|/7 - 1 = 6/7 for 13-speed LRC14 rows")
    print("    c_k = sum_{v|k, v in S} sin(pi*(k/v)/7)/(pi*(k/v)), k>0")
    print("  Low modes only see the speed divisor fibers.  This is the")
    print("  phase-sensitive cousin of HYP-2973's count distribution.")
    print()
    print("  Fejer vector used here:")
    print("    p_j=exp(-2*pi*i*j*x)/sqrt(d+1), centered at a safe-component midpoint.")
    print("    Q_d(x)=c_0+2*sum_{k<=d}(1-k/(d+1))*c_k*cos(2*pi*k*x).")
    print("    Q_d(x)<0 proves T_d is not PSD, but may occur after the true first")
    print("    negative eigenvalue.  Floating margins are reported; exact interval")
    print("    enclosures are the next formalization step.")
    print()


def tournament_fingerprint() -> dict[str, object]:
    tie_rank = {name: idx for idx, name in enumerate(TIE_PATH)}
    n = len(PROOF_VERTICES)
    mask = 0
    bit = 0
    for i, j in combinations(range(n), 2):
        a = PROOF_VERTICES[i]
        b = PROOF_VERTICES[j]
        av = PROOF_SCORES[a]
        bv = PROOF_SCORES[b]
        awins = sum(x > y for x, y in zip(av, bv))
        bwins = sum(y > x for x, y in zip(av, bv))
        if awins > bwins or (awins == bwins and tie_rank[a] < tie_rank[b]):
            mask |= 1 << bit
        bit += 1

    def edge(i: int, j: int) -> bool:
        if i == j:
            raise ValueError("loop")
        if i > j:
            return not edge(j, i)
        bit_idx = 0
        for a, b in combinations(range(n), 2):
            if a == i and b == j:
                return bool((mask >> bit_idx) & 1)
            bit_idx += 1
        raise AssertionError("unreachable")

    outdeg = [sum(1 for j in range(n) if i != j and edge(i, j)) for i in range(n)]
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if edge(a, b) and edge(b, c) and edge(c, a):
            c3 += 1
        if edge(a, c) and edge(c, b) and edge(b, a):
            c3 += 1
    return {
        "score_hist": dict(sorted(Counter(outdeg).items())),
        "c3": c3,
        "leaders": tuple(
            name for _, name in sorted(zip(outdeg, PROOF_VERTICES), reverse=True)[:4]
        ),
    }


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    runners, danger arcs, exact endpoints, safe components, count states,")
    print("    Fourier modes, Toeplitz rows, Fejer centers, denominator ladders,")
    print("    endpoint owners, packet labels, and proof obligations.")
    print("  chosen vertices:")
    print("    Fourier modes and explicit Toeplitz PSD-test vectors.")
    print("  preserved LRC predicate:")
    print("    a strict counterexample requires F_S=C_S-1>=0 a.e.; any negative")
    print("    Toeplitz quadratic form contradicts that and proves a strict safe set.")
    print("  destroyed information:")
    print("    endpoint-owner labels and exact interval atlas, except for the midpoint")
    print("    used to place the Fejer vector.")
    print("  challenged assumption:")
    print("    HYP-2973's danger-count distribution may not be the only moment")
    print("    object; phase-sensitive Fourier fibers can certify rows while retaining")
    print("    divisor/Farey structure.")
    print()


def print_tournament_analysis() -> None:
    fp = tournament_fingerprint()
    print("[5] Tournament Analysis")
    print("  vertices are proof carriers/Fourier packets, not runners.")
    print("  pair observable:")
    print("    retention of cover predicate, phase information, certificate")
    print("    concreteness, packet compatibility, exactness path, and")
    print("    anti-scalarization.")
    print("  switch/gauge:")
    print("    majority over the six retention coordinates; ties follow the displayed")
    print("    Hamiltonian path.")
    print("  Hamiltonian path:")
    print("    " + " > ".join(TIE_PATH))
    print(
        "  fingerprint: "
        f"score_hist={fp['score_hist']} c3={fp['c3']} leaders={fp['leaders']}"
    )
    print()


def print_theorem_readout(max_degree: int, bank_label: str, bank_audits: list[RowAudit]) -> None:
    positive_misses = [a for a in bank_audits if a.safe_measure > 0 and a.certificate is None]
    positive_hits = [a for a in bank_audits if a.safe_measure > 0 and a.certificate is not None]
    print("[6] Theorem-facing readout")
    print("  Necessary condition:")
    print("    Every strict counterexample must pass all finite Toeplitz PSD tests for")
    print("    F_S=C_S-1.  A negative Fejer quadratic form is a concrete dual witness.")
    print("  Current finite-bank evidence:")
    if positive_hits:
        max_hit = max(a.certificate.degree for a in positive_hits if a.certificate)
        print(
            f"    all {len(positive_hits)} positive {bank_label} rows audited here have a "
            f"Fejer PSD violation by degree <= {max_hit}."
        )
    print(f"    positive misses at cap {max_degree}: {len(positive_misses)}")
    print("    AP/GW remain zero-safe equality atoms, as expected.")
    print("  Proof target:")
    print("    Fourier-Toeplitz packet gate.  Outside AP/GW, a packet should expose")
    print("    either a bounded-degree divisor-curried Toeplitz violation, a")
    print("    HYP-2973 count-dual, or a labelled C27/K33/state-lift obstruction.")
    print("  Important caveat:")
    print("    This script gives floating Fejer-vector certificates, not exact")
    print("    interval-arithmetic proof objects and not true first eigenvalue degrees.")
    print("    The next formal step is to interval-enclose the reported negative")
    print("    trigonometric sums row-family-wise.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-degree", type=int, default=512)
    parser.add_argument("--negative-tol", type=float, default=1e-8)
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--skip-full-bank", action="store_true")
    parser.add_argument("--progress-every", type=int, default=5000)
    parser.add_argument("--hardest", type=int, default=20)
    args = parser.parse_args()

    print("S157 LRC14 FOURIER-TOEPLITZ PSD DUAL")
    print("=" * 78)
    print(
        f"max_degree={args.max_degree}, negative_tol={args.negative_tol}, "
        f"single_limit={args.single_limit}, two_swap_limit={args.two_swap_limit}, "
        f"alias_depth={args.alias_depth}, lcm_tail_max={args.lcm_tail_max}"
    )
    print_assumption_challenge()
    print_fourier_formula()

    named = audit_rows(named_rows(), args.max_degree, args.negative_tol)
    print_named_table(named)
    summarize_audits("named rows", named)

    one_swap = audit_rows(build_one_swap_rows(args.single_limit), args.max_degree, args.negative_tol)
    summarize_audits(f"one-swap AP bank through add<={args.single_limit}", one_swap)

    bank_audits: list[RowAudit] = []
    if not args.skip_full_bank:
        rows = build_hyp2963_rows(
            args.single_limit,
            args.two_swap_limit,
            args.alias_depth,
            args.lcm_tail_max,
        )
        print(f"[3] Full HYP-2963-bank Fejer scan rows={len(rows)}")
        bank_audits = audit_rows(rows, args.max_degree, args.negative_tol, args.progress_every)
        bank_label = "HYP-2963 packet-bank"
        summarize_audits("HYP-2963 packet bank", bank_audits)
        print_hardest(bank_audits, args.hardest)
    else:
        print("[3] Full HYP-2963-bank Fejer scan skipped by flag.")
        bank_audits = named + one_swap
        bank_label = "named/one-swap sample"

    print_tournament_analysis()
    print_theorem_readout(args.max_degree, bank_label, bank_audits)


if __name__ == "__main__":
    main()
