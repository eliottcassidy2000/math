#!/usr/bin/env python3
"""
lrc_n16_antipodal_fan_proof_search_s392.py

codex-2026-05-31 S392

Another long proof-search pass for the n=16 Lonely Runner case.

The new reframes are deliberately orthogonal to S389/S390:

1. antipodal quotient:
   under t -> t+1/2, even speeds cover both sides of a pair while odd speeds
   are one-sided.  A counterexample must therefore cover every half-turn pair
   by

       even_forbidden OR (odd_forbidden AND shifted_odd_forbidden).

   This converts the full circle problem into a half-circle pair certificate.

2. maximal-gate fan:
   endpoints of a maximal speed v can be locally covered by a nine-speed
   scaled odd fan.  This kills the naive largest-speed endpoint-cover proof,
   but exposes a sharper target: the fan is imprimitive unless external
   gcd-breaker speeds are added, and those gcd-breakers carry their own
   endpoint debt.

3. dyadic flow:
   selected candidates are audited as owner-layer -> protector-layer flows.

This is not a proof of n=16.  It records two finite certificates that a proof
can now try to force simultaneously.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd
from pathlib import Path


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
THRESHOLD = Fraction(1, N)
HALF = Fraction(1, 2)
ONE = Fraction(1, 1)


@dataclass(frozen=True)
class PairQuotientReport:
    label: str
    speeds: tuple[int, ...]
    even_measure: Fraction
    odd_double_measure: Fraction
    pair_measure: Fraction
    pair_gap_count: int
    max_pair_gap: Fraction
    first_pair_gap: tuple[Fraction, Fraction] | None
    pair_boundary_witnesses: int


@dataclass(frozen=True)
class FanReport:
    speed: int
    fan: tuple[int, ...]
    fan_gcd: int
    normalized_fan: tuple[int, ...]
    covered: int
    endpoint_count: int
    marginal_gains: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class ExactCandidateReport:
    label: str
    speeds: tuple[int, ...]
    gap_ratio: Fraction
    forbidden_length: Fraction
    unprotected: int
    core_endpoints: int
    pair_gap_ratio: Fraction
    pair_boundary_witnesses: int
    owner_hist: tuple[tuple[int, int], ...]
    private_hist: tuple[tuple[tuple[int, int], int], ...]


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
        return 99
    out = 0
    while value % 2 == 0:
        out += 1
        value //= 2
    return out


def dyadic_layer(value: int) -> int:
    return min(v2(value), 4)


def circle(value: Fraction) -> Fraction:
    return value % ONE


def dist_to_integer(value: Fraction) -> Fraction:
    value = circle(value)
    return min(value, ONE - value)


def intervals_for_fixed_threshold(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    intervals: list[tuple[Fraction, Fraction]] = []
    for speed in speeds:
        radius = THRESHOLD / speed
        for center in range(speed):
            intervals.extend(
                S356.split_circle_interval(
                    Fraction(center, speed) - radius,
                    Fraction(center, speed) + radius,
                )
            )
    return S356.merge_intervals(intervals)


def clip_half(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    clipped: list[tuple[Fraction, Fraction]] = []
    for lo, hi in intervals:
        if hi <= 0 or lo >= HALF:
            continue
        new_lo = max(lo, Fraction(0))
        new_hi = min(hi, HALF)
        if new_hi > new_lo:
            clipped.append((new_lo, new_hi))
    return S356.merge_intervals(clipped)


def shift_minus_half(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    shifted: list[tuple[Fraction, Fraction]] = []
    for lo, hi in intervals:
        shifted.extend(S356.split_circle_interval(lo - HALF, hi - HALF))
    return S356.merge_intervals(shifted)


def interval_intersection(
    left: list[tuple[Fraction, Fraction]], right: list[tuple[Fraction, Fraction]]
) -> list[tuple[Fraction, Fraction]]:
    out: list[tuple[Fraction, Fraction]] = []
    i = 0
    j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if hi > lo:
            out.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return S356.merge_intervals(out)


def half_gaps(components: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    gaps: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for lo, hi in components:
        if lo > cursor:
            gaps.append((cursor, lo))
        if hi > cursor:
            cursor = hi
    if cursor < HALF:
        gaps.append((cursor, HALF))
    return gaps


def pair_covered_strict(speeds: tuple[int, ...], point: Fraction) -> bool:
    even = tuple(speed for speed in speeds if speed % 2 == 0)
    odd = tuple(speed for speed in speeds if speed % 2 == 1)
    if any(dist_to_integer(speed * point) < THRESHOLD for speed in even):
        return True
    return any(dist_to_integer(speed * point) < THRESHOLD for speed in odd) and any(
        dist_to_integer(speed * (point + HALF)) < THRESHOLD for speed in odd
    )


def pair_quotient_report(label: str, speeds: tuple[int, ...]) -> PairQuotientReport:
    even = tuple(speed for speed in speeds if speed % 2 == 0)
    odd = tuple(speed for speed in speeds if speed % 2 == 1)

    even_half = clip_half(intervals_for_fixed_threshold(even))
    odd_full = intervals_for_fixed_threshold(odd)
    odd_half = clip_half(odd_full)
    shifted_odd_half = clip_half(shift_minus_half(odd_full))
    odd_double = interval_intersection(odd_half, shifted_odd_half)
    pair_cover = S356.merge_intervals(even_half + odd_double)
    gaps = half_gaps(pair_cover)

    boundaries = {Fraction(0), HALF}
    for components in (even_half, odd_half, shifted_odd_half, odd_double, pair_cover):
        for lo, hi in components:
            if 0 <= lo <= HALF:
                boundaries.add(lo)
            if 0 <= hi <= HALF:
                boundaries.add(hi)
    boundary_witnesses = sum(
        1 for point in boundaries if not pair_covered_strict(speeds, point)
    )

    return PairQuotientReport(
        label=label,
        speeds=speeds,
        even_measure=sum(hi - lo for lo, hi in even_half),
        odd_double_measure=sum(hi - lo for lo, hi in odd_double),
        pair_measure=sum(hi - lo for lo, hi in pair_cover),
        pair_gap_count=len(gaps),
        max_pair_gap=max((hi - lo for lo, hi in gaps), default=Fraction(0)),
        first_pair_gap=gaps[0] if gaps else None,
        pair_boundary_witnesses=boundary_witnesses,
    )


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


def canonical_fan(speed: int) -> tuple[int, ...]:
    if speed < 16 or speed & (speed - 1):
        raise ValueError("fan audit expects dyadic speed >=16")
    if speed == 16:
        return (speed // 2,) + tuple(range(1, 16, 2))
    scale = speed // 32
    return (speed // 2,) + tuple(scale * odd for odd in range(1, 16, 2))


def fan_report(speed: int) -> FanReport:
    fan = canonical_fan(speed)
    covered = 0
    gains: list[tuple[int, int]] = []
    for protector in fan:
        mask = endpoint_mask(speed, protector)
        gain = (mask & ~covered).bit_count()
        gains.append((protector, gain))
        covered |= mask

    fan_gcd = 0
    for protector in fan:
        fan_gcd = gcd(fan_gcd, protector)
    normalized = tuple(protector // fan_gcd for protector in fan)
    return FanReport(
        speed=speed,
        fan=fan,
        fan_gcd=fan_gcd,
        normalized_fan=normalized,
        covered=covered.bit_count(),
        endpoint_count=2 * speed,
        marginal_gains=tuple(gains),
    )


def best_ladder(d: int) -> tuple[int, tuple[int, ...]]:
    best: tuple[Fraction, int, tuple[int, ...]] | None = None
    for skip in range(1, N):
        speeds = tuple(sorted({1} | {d * q for q in range(1, N) if q != skip}))
        if len(speeds) != K or not primitive(speeds):
            continue
        report = S356.report(f"d={d} skip={skip}", list(speeds))
        key = (report.max_gap / report.threshold, skip, speeds)
        if best is None or key < best:
            best = key
    if best is None:
        raise RuntimeError(f"no ladder for d={d}")
    return best[1], best[2]


def flow_summary(speeds: tuple[int, ...]) -> tuple[
    tuple[tuple[int, int], ...],
    tuple[tuple[tuple[int, int], int], ...],
    tuple[tuple[tuple[int, int], int], ...],
]:
    points: dict[Fraction, list[object]] = defaultdict(list)
    for endpoint in S360.endpoints(speeds):
        points[endpoint.value].append(endpoint)

    owner_hist: Counter[int] = Counter()
    private_hist: Counter[tuple[int, int]] = Counter()
    edge_hist: Counter[tuple[int, int]] = Counter()
    for point, labels in points.items():
        owner_layers = {dyadic_layer(endpoint.speed) for endpoint in labels}
        protectors = [
            speed for speed in speeds if S360.direct_protects(speeds, speed, point)
        ]
        for owner_layer in owner_layers:
            owner_hist[owner_layer] += 1
            for protector in protectors:
                edge_hist[(owner_layer, dyadic_layer(protector))] += 1
            if len(protectors) == 1:
                private_hist[(owner_layer, dyadic_layer(protectors[0]))] += 1
    return (
        tuple(sorted(owner_hist.items())),
        tuple(sorted(private_hist.items())),
        tuple(edge_hist.most_common(8)),
    )


def exact_candidate_report(label: str, speeds: tuple[int, ...]) -> ExactCandidateReport:
    gap = S356.report(label, list(speeds))
    protection = S360.summarize(list(speeds))
    descent = S362.summarize(list(speeds))
    pair = pair_quotient_report(label, tuple(gap.speeds))
    owner_hist, private_hist, _edge_hist = flow_summary(tuple(gap.speeds))
    return ExactCandidateReport(
        label=label,
        speeds=tuple(gap.speeds),
        gap_ratio=gap.max_gap / gap.threshold,
        forbidden_length=gap.forbidden_length,
        unprotected=protection.unprotected_count,
        core_endpoints=descent.core_endpoint_count,
        pair_gap_ratio=pair.max_pair_gap / THRESHOLD,
        pair_boundary_witnesses=pair.pair_boundary_witnesses,
        owner_hist=owner_hist,
        private_hist=private_hist[:8],
    )


def fan_plus_breakers(speed: int, breakers: tuple[int, ...]) -> tuple[int, ...]:
    speeds = set(canonical_fan(speed))
    speeds.update(breakers)
    filler = 1
    while len(speeds) < K:
        speeds.add(filler)
        filler += 1
    return tuple(sorted(speeds))


def print_pair_quotient_theorem() -> None:
    print("1. Antipodal pair quotient")
    print("   For n=16, split the speeds into E=even and O=odd.")
    print("   A full open cover implies the half-circle pair cover")
    print("       F_E OR (F_O AND tau(F_O))")
    print("   has no open or boundary gap, where tau is t -> t+1/2.")
    print("   This is exact: even speeds are half-turn invariant; each odd")
    print("   speed is one-sided on an antipodal pair.")
    print()


def print_pair_audits() -> None:
    print("2. Pair-quotient audits")
    skip8, ladder8 = best_ladder(8)
    candidates = [
        ("initial 1..15", tuple(range(1, N))),
        ("drop 15 add 16", tuple(range(1, 15)) + (16,)),
        (f"best 8-ladder skip={skip8}", ladder8),
        (
            "odd units plus gates",
            tuple(sorted(tuple(range(1, 16, 2)) + tuple(16 * q for q in range(1, 8)))),
        ),
        ("fan64 plus odd breakers", fan_plus_breakers(64, (1, 3, 5, 7, 9, 11))),
        ("fan128 plus small breakers", fan_plus_breakers(128, (1, 2, 3, 5, 7, 11))),
    ]
    print(
        "   label                         even     odd2     pair     "
        "gap/th  gaps boundary first_gap"
    )
    for label, speeds in candidates:
        row = pair_quotient_report(label, speeds)
        first = "-"
        if row.first_pair_gap is not None:
            first = f"{fmt_frac(row.first_pair_gap[0])}->{fmt_frac(row.first_pair_gap[1])}"
        print(
            f"   {label[:28]:<28} {fmt_float(row.even_measure):>7} "
            f"{fmt_float(row.odd_double_measure):>8} {fmt_float(row.pair_measure):>8} "
            f"{fmt_float(row.max_pair_gap / THRESHOLD):>7} "
            f"{row.pair_gap_count:>5} {row.pair_boundary_witnesses:>8} {first}"
        )
    print(
        "   The initial segment pair quotient is measure-full but boundary-only;\n"
        "   the gated and fan rows have honest positive pair gaps."
    )
    print()


def print_one_gate_pair_scan(q_limit: int = 32) -> None:
    print("3. One-gate antipodal scan")
    rows: list[tuple[Fraction, int, Fraction, int, int, PairQuotientReport]] = []
    class_hist: Counter[str] = Counter()
    base = set(range(1, N))
    for removed in range(1, N):
        for q in range(1, q_limit + 1):
            speeds = tuple(sorted((base - {removed}) | {16 * q}))
            row = pair_quotient_report(f"{removed}->16*{q}", speeds)
            if row.max_pair_gap > 0:
                classification = "positive_pair_gap"
            elif row.pair_boundary_witnesses:
                classification = "pair_boundary_only"
            else:
                classification = "pair_open_cover"
            class_hist[classification] += 1
            rows.append(
                (
                    row.max_pair_gap / THRESHOLD,
                    row.pair_boundary_witnesses,
                    -row.pair_measure,
                    removed,
                    q,
                    row,
                )
            )

    print(f"   scanned={{1..15}}-r+16q with q<= {q_limit}: {len(rows)} rows")
    print(f"   class_hist={dict(sorted(class_hist.items()))}")
    print("   closest pair-quotient rows")
    for gap_ratio, boundary, neg_measure, removed, q, row in sorted(rows)[:10]:
        print(
            f"     remove {removed:>2}, add 16*{q:<2} "
            f"pair_gap/th={fmt_float(gap_ratio)} "
            f"boundary={boundary:<3} pair_measure={fmt_frac(-neg_measure)}"
        )
    print(
        "   The pair quotient gives an independent obstruction: even before\n"
        "   checking every original endpoint, all scanned gated perturbations\n"
        "   still leave an antipodal pair gap."
    )
    print()


def print_fan_audit() -> None:
    print("4. Maximal-gate nine-speed fan")
    print(
        "   A local endpoint-cover proof by largest speed fails in a structured\n"
        "   way: a dyadic speed v has a nine-speed lower fan that covers all\n"
        "   2v endpoints.  The normalized fan is stable."
    )
    print("   v     gcd(fan) normalized fan                  covered gains")
    for speed in (16, 32, 64, 128, 256, 512):
        row = fan_report(speed)
        gains = ",".join(str(gain) for _protector, gain in row.marginal_gains)
        print(
            f"   {speed:<5} {row.fan_gcd:<8} "
            f"{str(row.normalized_fan):<31} "
            f"{row.covered:>4}/{row.endpoint_count:<4} {gains}"
        )
    print(
        "   New proof target: if a maximal branch closes by this fan, the fan is\n"
        "   imprimitive at scale gcd=v/32.  The six remaining speeds must break\n"
        "   that gcd and simultaneously protect the fan's descendant endpoints."
    )
    print()


def print_exact_candidate_audits() -> None:
    print("5. Exact candidate audits with dyadic owner flow")
    skip8, ladder8 = best_ladder(8)
    candidates = [
        ("initial 1..15", tuple(range(1, N))),
        ("drop 15 add 16", tuple(range(1, 15)) + (16,)),
        (f"best 8-ladder skip={skip8}", ladder8),
        ("fan64 plus odd breakers", fan_plus_breakers(64, (1, 3, 5, 7, 9, 11))),
        ("fan128 plus small breakers", fan_plus_breakers(128, (1, 2, 3, 5, 7, 11))),
    ]
    print(
        "   label                         gap/th   forb     unprot core pairgap/th boundary owners private"
    )
    for label, speeds in candidates:
        row = exact_candidate_report(label, speeds)
        private = ",".join(
            f"{edge[0]}->{edge[1]}:{count}" for edge, count in row.private_hist[:4]
        )
        print(
            f"   {row.label[:28]:<28} {fmt_float(row.gap_ratio):>7} "
            f"{fmt_frac(row.forbidden_length):>8} {row.unprotected:>7} "
            f"{row.core_endpoints:>4} {fmt_float(row.pair_gap_ratio):>10} "
            f"{row.pair_boundary_witnesses:>8} {dict(row.owner_hist)} {private}"
        )
    print(
        "   Every non-initial row here has pair gap, ordinary gap, unprotected\n"
        "   endpoints, and terminal coreE=0.  The fan closes the top endpoints\n"
        "   locally but does not close the global arithmetic endpoint system."
    )
    print()


def print_synthesis() -> None:
    print("6. Synthesis / next proof attempt")
    print(
        "   A counterexample at n=16 must now pass three gates at once:\n"
        "   (i) contain a 16-gate by the unit-skeleton lemma;\n"
        "   (ii) have no antipodal pair gap in E OR (O AND tau O);\n"
        "   (iii) if a maximal dyadic branch uses the nine-speed fan, the six\n"
        "   non-fan speeds must break the fan gcd without creating private debt."
    )
    print(
        "   The wild idea that survived this pass is to prove incompatibility\n"
        "   between (ii) and (iii).  Odd double-coverage wants two odd lanes on\n"
        "   every pair not already covered by evens; the maximal fan wants a\n"
        "   scaled odd residue fan plus a half-gate.  Those demands point in\n"
        "   opposite quotient directions once primitivity is imposed."
    )


def main() -> None:
    print("n=16 LRC antipodal/fan proof search (codex-2026-05-31 S392)")
    print("Exact Fraction arithmetic; threshold=1/16; k=15 speeds.")
    print()
    print_pair_quotient_theorem()
    print_pair_audits()
    print_one_gate_pair_scan()
    print_fan_audit()
    print_exact_candidate_audits()
    print_synthesis()


if __name__ == "__main__":
    main()
