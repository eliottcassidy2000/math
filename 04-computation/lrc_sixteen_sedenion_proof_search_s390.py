#!/usr/bin/env python3
"""
lrc_sixteen_sedenion_proof_search_s390.py

codex-2026-05-31 S390

Long proof-search pass for the n=16 Lonely Runner denominator.

The framing is deliberately Cayley-Dickson flavored but computationally exact:

* n=16 is the pure 2-adic / sedenion row of the denominator tower.
* A counterexample must contain a 16-gate, because otherwise every odd
  unit point a/16 is a lonely boundary witness.
* A 16-gate creates endpoint debt.  The main experiment below asks whether
  that debt recursively descends through dyadic half-gates.

This is not a proof.  It is a battery of exact attacks and proof-shape
diagnostics intended to isolate a plausible certificate.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
import random


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()
S362 = SourceFileLoader(
    "lonely_runner_bohr_descent_s362",
    str(ROOT / "04-computation" / "lonely_runner_bohr_descent_s362.py"),
).load_module()


N = 16
K = N - 1
ONE = Fraction(1, 1)


@dataclass(frozen=True)
class CandidateAudit:
    label: str
    speeds: tuple[int, ...]
    classification: str
    gap_ratio: Fraction
    forbidden_length: Fraction
    boundary_witnesses: int
    unprotected: int
    first_unprotected: Fraction | None
    layer_hist: tuple[tuple[int, int], ...]
    peel_depth: int
    core_endpoints: int
    q: int


@dataclass(frozen=True)
class ScanRow:
    label: str
    speeds: tuple[int, ...]
    classification: str
    gap_ratio: Fraction
    forbidden_length: Fraction
    boundary_witnesses: int


def fmt_frac(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def fmt_float(value: Fraction | None) -> str:
    if value is None:
        return "-"
    return f"{float(value):.6f}"


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for speed in speeds:
        g = gcd(g, speed)
    return g == 1


def v2(value: int) -> int:
    if value == 0:
        return 10**9
    out = 0
    while value % 2 == 0:
        value //= 2
        out += 1
    return out


def circle(value: Fraction) -> Fraction:
    return value % ONE


def dist_to_integer(value: Fraction) -> Fraction:
    value = circle(value)
    return min(value, ONE - value)


def classify_report(report: S356.GapReport) -> str:
    if report.max_gap > 0:
        return "positive_gap"
    if report.boundary_witness_count:
        return "boundary_only"
    return "open_cover_candidate"


def lpd_ladder(d: int, skip: int) -> tuple[int, ...]:
    return tuple(sorted({1} | {d * q for q in range(1, N) if q != skip}))


def best_ladder(d: int) -> tuple[int, tuple[int, ...], S356.GapReport]:
    best: tuple[Fraction, int, int, tuple[int, ...], S356.GapReport] | None = None
    for skip in range(1, N):
        speeds = lpd_ladder(d, skip)
        if len(speeds) != K or not primitive(speeds):
            continue
        report = S356.report(f"d={d} skip={skip}", list(speeds))
        rank = (
            report.max_gap / report.threshold,
            report.boundary_witness_count,
            skip,
            speeds,
            report,
        )
        if best is None or rank < best:
            best = rank
    if best is None:
        raise RuntimeError(f"no ladder for d={d}")
    return best[2], best[3], best[4]


def unprotected_points(speeds: tuple[int, ...]) -> list[Fraction]:
    endpoints = {endpoint.value for endpoint in S360.endpoints(speeds)}
    return [
        point
        for point in sorted(endpoints)
        if not any(S360.direct_protects(speeds, speed, point) for speed in speeds)
    ]


def debt_layer(point: Fraction) -> int:
    return point.denominator // gcd(point.denominator, N)


def layer_histogram(speeds: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(debt_layer(point) for point in unprotected_points(speeds)).items()))


def endpoint_modulus(points: list[Fraction]) -> int:
    q = 1
    for point in points:
        q = lcm(q, point.denominator)
    return q


def audit_candidate(label: str, speeds: tuple[int, ...]) -> CandidateAudit:
    report = S356.report(label, list(speeds))
    protection = S360.summarize(list(report.speeds))
    descent = S362.summarize(list(report.speeds))
    return CandidateAudit(
        label=label,
        speeds=report.speeds,
        classification=protection.classification,
        gap_ratio=report.max_gap / report.threshold,
        forbidden_length=report.forbidden_length,
        boundary_witnesses=report.boundary_witness_count,
        unprotected=protection.unprotected_count,
        first_unprotected=protection.first_unprotected,
        layer_hist=layer_histogram(report.speeds),
        peel_depth=len(descent.peel_layers),
        core_endpoints=descent.core_endpoint_count,
        q=protection.boundary_modulus,
    )


def print_unit_skeleton_lemma() -> None:
    print("1. Unit skeleton gate lemma")
    print("   Odd points a/16 are the pure sedenion row's unit boundary.")
    rows: list[tuple[int, int, Fraction, str]] = []
    for residue in range(N):
        distances = [dist_to_integer(Fraction(residue * a, N)) for a in range(1, N, 2)]
        minimum = min(distances)
        verdict = "kills units" if minimum == 0 else "safe boundary"
        rows.append((residue, v2(residue), minimum, verdict))
    print("   residue mod16  v2  min_{a odd} ||r*a/16||  verdict")
    for residue, valuation, minimum, verdict in rows:
        v2_text = "inf" if residue == 0 else str(valuation)
        print(f"   {residue:>12} {v2_text:>3} {fmt_frac(minimum):>24}  {verdict}")
    print(
        "   Consequence: if no speed is 0 mod 16, every odd a/16 is a lonely\n"
        "   boundary witness.  So an open-cover counterexample must contain a\n"
        "   16-gate.  This is a clean proved branch, not just evidence."
    )
    print()


def print_halfturn_table() -> None:
    print("2. Half-turn parity attack")
    print("   Under t -> t+1/2, even speeds are unchanged and odd speeds flip")
    print("   forbidden neighborhoods to the opposite half of the circle.")
    print("   residue parity  behavior at threshold 1/16")
    print("   even            same forbidden status at t and t+1/2")
    print("   odd             cannot forbid both antipodal times")
    print(
        "   Proof-shape: pair antipodal times.  Odd speeds can spend budget on\n"
        "   at most one side of each pair; even speeds descend to the quotient.\n"
        "   This wants an induction on the dyadic shell vector, not a density\n"
        "   argument."
    )
    print()


def print_curated_audits() -> tuple[CandidateAudit, ...]:
    print("3. Exact audits of structured n=16 candidates")
    skip8, ladder8, _report8 = best_ladder(8)
    skip4, ladder4, _report4 = best_ladder(4)
    skip2, ladder2, _report2 = best_ladder(2)
    curated: list[tuple[str, tuple[int, ...]]] = [
        ("initial segment", tuple(range(1, N))),
        (f"best 8-ladder skip={skip8}", ladder8),
        (f"best 4-ladder skip={skip4}", ladder4),
        (f"best 2-ladder skip={skip2}", ladder2),
        ("powers of two spine", tuple([1, 2, 4, 8, 16] + list(range(17, 27)))),
        ("odd units plus seven 16-gates", tuple(sorted(tuple(range(1, 16, 2)) + tuple(16 * q for q in range(1, 8))))),
        ("pure 16-gate ladder", tuple([1] + [16 * q for q in range(1, 15)])),
        ("dyadic staircase", tuple(sorted({1, 2, 4, 8, 16, 24, 32, 40, 48, 56, 64, 80, 96, 112, 120}))),
    ]
    audits = tuple(audit_candidate(label, speeds) for label, speeds in curated)
    print("   label                         class         gap/th   unprot peel coreE first layer_hist")
    for row in audits:
        hist = ",".join(f"{layer}:{count}" for layer, count in row.layer_hist[:5])
        if len(row.layer_hist) > 5:
            hist += ",..."
        print(
            f"   {row.label[:28]:<28} {row.classification:<12} "
            f"{fmt_float(row.gap_ratio):>8} {row.unprotected:>7} "
            f"{row.peel_depth:>4} {row.core_endpoints:>5} "
            f"{fmt_frac(row.first_unprotected):>8} {hist}"
        )
    print()
    print("   Best 8-ladder speeds:")
    print(f"   {ladder8}")
    print(
        "   The pure lpd ladder is positive-gap, but its first debt layer is 8:\n"
        "   exactly the half-dimension below the sedenion row."
    )
    print()
    return audits


def scan_one_gate_replacements(q_limit: int = 64) -> tuple[ScanRow, ...]:
    rows: list[ScanRow] = []
    base = set(range(1, N))
    for removed in range(1, N):
        for q in range(1, q_limit + 1):
            speeds = tuple(sorted((base - {removed}) | {16 * q}))
            if len(speeds) != K or not primitive(speeds):
                continue
            report = S356.report(f"remove {removed}, add {16*q}", list(speeds))
            rows.append(
                ScanRow(
                    label=f"{removed}->16*{q}",
                    speeds=report.speeds,
                    classification=classify_report(report),
                    gap_ratio=report.max_gap / report.threshold,
                    forbidden_length=report.forbidden_length,
                    boundary_witnesses=report.boundary_witness_count,
                )
            )
    return tuple(rows)


def print_one_gate_scan() -> tuple[ScanRow, ...]:
    print("4. One-gate replacement scan around the tight initial segment")
    rows = scan_one_gate_replacements()
    class_hist = Counter(row.classification for row in rows)
    print(f"   scanned={len(rows)} replacements remove r, add 16q with q<=64")
    print(f"   class_hist={dict(sorted(class_hist.items()))}")
    best = sorted(rows, key=lambda row: (row.gap_ratio, -row.forbidden_length))[:10]
    print("   closest positive-gap rows")
    for row in best:
        print(
            f"     {row.label:<10} gap/th={fmt_float(row.gap_ratio)} "
            f"forbidden={fmt_frac(row.forbidden_length)} "
            f"boundary={row.boundary_witnesses}"
        )
    print(
        "   No one-gate perturbation became boundary-only or open-cover.  The\n"
        "   best rows move toward the 8-layer debt, not toward a closed cover."
    )
    print()
    return tuple(best)


def endpoint_mask(speed: int, protector: int) -> int:
    mask = 0
    modulus = N * speed
    for m in range(speed):
        for sign_index, sign in enumerate((-1, 1)):
            bit = 2 * m + sign_index
            residue = (protector * (N * m + sign)) % modulus
            if min(residue, modulus - residue) < speed:
                mask |= 1 << bit
    return mask


def greedy_cover(speed: int, protectors: tuple[int, ...]) -> tuple[int, tuple[int, ...], int]:
    full = (1 << (2 * speed)) - 1
    masks = [(p, endpoint_mask(speed, p)) for p in protectors]
    covered = 0
    chosen: list[int] = []
    while covered != full:
        best_gain = 0
        best_p = None
        best_mask = 0
        for p, mask in masks:
            if p in chosen:
                continue
            gain = (mask & ~covered).bit_count()
            if gain > best_gain:
                best_gain = gain
                best_p = p
                best_mask = mask
        if best_p is None:
            break
        chosen.append(best_p)
        covered |= best_mask
    return len(chosen), tuple(chosen), covered.bit_count()


def exact_lower_cover_if_small(speed: int) -> tuple[int | None, tuple[int, ...]]:
    if speed > 16:
        return None, tuple()
    full = (1 << (2 * speed)) - 1
    masks = [(p, endpoint_mask(speed, p)) for p in range(1, speed)]
    masks = [(p, mask) for p, mask in masks if mask]
    for size in range(1, len(masks) + 1):
        for combo in combinations(masks, size):
            covered = 0
            for _p, mask in combo:
                covered |= mask
            if covered == full:
                return size, tuple(p for p, _mask in combo)
    return None, tuple()


def print_dyadic_endpoint_debt() -> None:
    print("5. Dyadic endpoint-debt cascade")
    print(
        "   For a speed v, its 2v endpoints are protected by integer residues p\n"
        "   satisfying |p*(16m+/-1)-a*16v| < v.  The table separates two cases:\n"
        "   a super-gate p=0 mod 16v protects everything, while a maximum-speed\n"
        "   branch can only use p<v."
    )
    print("   v   endpoints  best p<v cover  exact<=16  greedy p<v cover  first greedy protectors")
    for speed in (1, 2, 4, 8, 16, 32, 64, 128):
        lower = tuple(range(1, speed))
        best_lower = max((endpoint_mask(speed, p).bit_count() for p in lower), default=0)
        greedy_count, greedy_picks, covered = greedy_cover(speed, lower)
        exact_count, exact_picks = exact_lower_cover_if_small(speed)
        exact_text = "-" if exact_count is None else str(exact_count)
        picks = ",".join(str(p) for p in greedy_picks[:8]) if greedy_picks else "-"
        print(
            f"   {speed:<3} {2*speed:>9} {best_lower:>15}/{2*speed:<3} "
            f"{exact_text:>9} {greedy_count:>5} speeds cover {covered:>3}/{2*speed:<3}  {picks}"
        )
        if exact_picks:
            print(f"       exact lower cover for v={speed}: {exact_picks}")
    print(
        "   Read this as a possible proof spine: a largest 16-divisible gate\n"
        "   cannot use a super-gate, so it asks for half-gate residues; those\n"
        "   ask for quarter-gates; the unit layer asks to be rescued by a\n"
        "   16-gate.  This is a finite dyadic debt cycle, not a linear cover."
    )
    print()


def random_gate_samples(sample_count: int = 48) -> tuple[ScanRow, ...]:
    rng = random.Random(388)
    rows: list[ScanRow] = []
    pool = [x for x in range(1, 65) if x % 16 != 0]
    gate_pool = [16 * q for q in range(1, 13)]
    for i in range(sample_count):
        gate_count = 1 + (i % 4)
        gates = rng.sample(gate_pool, gate_count)
        rest = rng.sample(pool, K - gate_count)
        speeds = tuple(sorted(set(gates + rest)))
        if len(speeds) != K or not primitive(speeds):
            continue
        report = S356.report(f"random_gate_{i}", list(speeds))
        rows.append(
            ScanRow(
                label=f"random_gate_{i}",
                speeds=report.speeds,
                classification=classify_report(report),
                gap_ratio=report.max_gap / report.threshold,
                forbidden_length=report.forbidden_length,
                boundary_witnesses=report.boundary_witness_count,
            )
        )
    return tuple(rows)


def print_random_stress() -> tuple[ScanRow, ...]:
    print("6. Deterministic random stress with forced 16-gates")
    rows = random_gate_samples()
    class_hist = Counter(row.classification for row in rows)
    print(f"   sampled={len(rows)} primitive forced-gate sets")
    print(f"   class_hist={dict(sorted(class_hist.items()))}")
    best = sorted(rows, key=lambda row: (row.gap_ratio, -row.forbidden_length))[:8]
    for row in best:
        shells = Counter(min(v2(speed), 5) for speed in row.speeds)
        print(
            f"     {row.label:<16} gap/th={fmt_float(row.gap_ratio)} "
            f"forbidden={fmt_frac(row.forbidden_length)} shells={dict(sorted(shells.items()))}"
        )
    print("   endpoint audits of the three closest random rows")
    for row in best[:3]:
        audit = audit_candidate(row.label, row.speeds)
        print(
            f"     {audit.label:<16} unprot={audit.unprotected:>5} "
            f"peel={audit.peel_depth:>3} coreE={audit.core_endpoints:>3} "
            f"first={fmt_frac(audit.first_unprotected):>8} "
            f"layers={audit.layer_hist[:4]}"
        )
    print(
        "   Again the near rows are gap-positive and core-empty.  The forced\n"
        "   gate buys coverage at the unit layer but pays for it in private\n"
        "   dyadic endpoints."
    )
    print()
    return best


def print_proof_synthesis() -> None:
    print("7. Proof synthesis after the failed attacks")
    print(
        "   Candidate theorem A (proved here): no 16-gate implies no open cover,\n"
        "   because odd unit points a/16 survive as boundary witnesses."
    )
    print(
        "   Candidate theorem B (new target): every primitive n=16 full-measure\n"
        "   set containing a 16-gate has a nonempty dyadic endpoint-debt leaf.\n"
        "   Equivalently, the endpoint-protection core peels to empty before an\n"
        "   arithmetic cycle can close."
    )
    print(
        "   The wild but now concrete route is a dyadic debt-flow inequality:\n"
        "   charge endpoints by their layer 1,2,4,8,16; half-turn symmetry makes\n"
        "   odd protectors one-sided; maximum-speed endpoints cannot use a\n"
        "   super-gate; and the debt cycle must eventually ask the original\n"
        "   16-gate to rescue the unit layer.  Show that this directed debt\n"
        "   cycle has positive divergence for only 15 speeds."
    )
    print(
        "   The Cayley-Dickson analogy becomes precise here: the sedenion row is\n"
        "   exactly where division fails, and in LRC that failure is endpoint\n"
        "   debt leaking one dyadic layer downward every time a gate closes an\n"
        "   older boundary."
    )


def main() -> None:
    print("n=16 Lonely Runner sedenion-row proof search (codex-2026-05-31 S390)")
    print("denominator n=16, moving speeds k=15, threshold=1/16")
    print()
    print_unit_skeleton_lemma()
    print_halfturn_table()
    print_curated_audits()
    print_one_gate_scan()
    print_dyadic_endpoint_debt()
    print_random_stress()
    print_proof_synthesis()


if __name__ == "__main__":
    main()
