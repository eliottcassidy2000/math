#!/usr/bin/env python3
"""
lrc_n14_creative_reframes_s440.py

codex-2026-05-31 S440

Creative reframes for the n=14 Lonely Runner frontier.

This pass deliberately does not search harder for speed sets.  It extracts
proof-shaped ledgers from the existing near-misses:

1. small-denominator rows as a modulus-cover market;
2. the antipodal quotient t ~ t+1/2 as a parity fold;
3. the local endpoint-cover "fan tax" of a 14-gate;
4. owner-wise endpoint debt in gate-heavy rows;
5. product-depth potentials on the 2 x 7 CRT denominator grid.

The target is a proof grammar for n=14: no 14-gate leaves the unit boundary,
while a 14-gate must pay a heptagonal local fan and then exports debt into the
product of the 2- and 7-adic depth directions.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd, prod
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
S420 = SourceFileLoader(
    "lrc_integer_programming_modes_s420",
    str(ROOT / "04-computation" / "lrc_integer_programming_modes_s420.py"),
).load_module()


N = 14
K = 13
ONE = Fraction(1, 1)
HALF = Fraction(1, 2)
PRIMES = (2, 7)


@dataclass(frozen=True)
class ReframeRow:
    label: str
    speeds: tuple[int, ...]
    classification: str
    gap_ratio: Fraction
    unprotected: int
    peel_depth: int
    core_endpoints: int
    missing_moduli: tuple[int, ...]
    fragile_moduli: tuple[int, ...]
    pair_gap_ratio: Fraction
    pair_components: int
    depth_hist: tuple[tuple[tuple[int, ...], int], ...]
    frontier_mass: Fraction
    denominator_pressure: int


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def fmt_float(value: Fraction) -> str:
    return f"{float(value):.6f}"


def circle(value: Fraction) -> Fraction:
    return value % ONE


def initial() -> tuple[int, ...]:
    return tuple(range(1, N))


def drop_add(dropped: int, added: int) -> tuple[int, ...]:
    speeds = set(initial())
    speeds.remove(dropped)
    speeds.add(added)
    return tuple(sorted(speeds))


def seven_ladder() -> tuple[int, ...]:
    return (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91)


def s380_gate_ladder() -> tuple[int, ...]:
    return (1, 14, 28, 42, 56, 70, 98, 112, 126, 140, 154, 168, 182)


def seven_with_42() -> tuple[int, ...]:
    return tuple(sorted((set(seven_ladder()) - {84}) | {42}))


def boundary_swap_12_24() -> tuple[int, ...]:
    return drop_add(12, 24)


def canonical(raw: tuple[int, ...]) -> tuple[int, ...]:
    return S356.normalize_speed_set(list(raw))


def vp(value: int, prime: int) -> int:
    value = abs(value)
    if value == 0:
        raise ValueError("v_p(0) is not finite here")
    out = 0
    while value % prime == 0:
        out += 1
        value //= prime
    return out


def extra_depth(point: Fraction) -> tuple[int, ...]:
    return tuple(max(0, vp(point.denominator, prime) - vp(N, prime)) for prime in PRIMES)


def depth_scale(depth: tuple[int, ...]) -> int:
    return prod(prime**height for prime, height in zip(PRIMES, depth))


def unprotected_values(speeds: tuple[int, ...]) -> set[Fraction]:
    values = {endpoint.value for endpoint in S360.endpoints(speeds)}
    return {
        value
        for value in values
        if not any(S360.direct_protects(speeds, speed, value) for speed in speeds)
    }


def moduli_ledger(speeds: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...], dict[int, int]]:
    payers = {
        modulus: sum(1 for speed in speeds if speed % modulus == 0)
        for modulus in range(2, N + 1)
    }
    missing = tuple(modulus for modulus, count in payers.items() if count == 0)
    fragile = tuple(modulus for modulus, count in payers.items() if count == 1)
    return missing, fragile, payers


def clip_interval(
    interval: tuple[Fraction, Fraction],
    lo: Fraction,
    hi: Fraction,
) -> tuple[Fraction, Fraction] | None:
    start = max(interval[0], lo)
    end = min(interval[1], hi)
    if start < end:
        return (start, end)
    return None


def intersect_interval_lists(
    left: list[tuple[Fraction, Fraction]],
    right: list[tuple[Fraction, Fraction]],
) -> list[tuple[Fraction, Fraction]]:
    out: list[tuple[Fraction, Fraction]] = []
    for a_lo, a_hi in left:
        for b_lo, b_hi in right:
            lo = max(a_lo, b_lo)
            hi = min(a_hi, b_hi)
            if lo < hi:
                out.append((lo, hi))
    return S356.merge_intervals(out) if out else []


def antipodal_pair_cover(
    speeds: tuple[int, ...]
) -> tuple[list[tuple[Fraction, Fraction]], list[tuple[Fraction, Fraction]]]:
    """Return quotient pair-covered intervals and quotient gaps in t-units.

    A quotient point q in [0,1/2) is pair-covered when q and q+1/2 are both
    forbidden.  If the quotient has any gap, at least one of each corresponding
    antipodal pair is a lonely witness.
    """

    components = S356.merge_intervals(S356.forbidden_intervals(speeds))
    lower: list[tuple[Fraction, Fraction]] = []
    upper_shifted: list[tuple[Fraction, Fraction]] = []
    for component in components:
        clipped_lower = clip_interval(component, Fraction(0), HALF)
        if clipped_lower is not None:
            lower.append(clipped_lower)
        clipped_upper = clip_interval(component, HALF, ONE)
        if clipped_upper is not None:
            upper_shifted.append((clipped_upper[0] - HALF, clipped_upper[1] - HALF))

    pair = intersect_interval_lists(
        S356.merge_intervals(lower) if lower else [],
        S356.merge_intervals(upper_shifted) if upper_shifted else [],
    )
    scaled = [(2 * lo, 2 * hi) for lo, hi in pair]
    scaled_components = S356.merge_intervals(scaled) if scaled else []
    scaled_gaps = S356.circular_gaps(scaled_components)
    gaps = [(lo / 2, hi / 2) for lo, hi in scaled_gaps]
    return pair, gaps


def pair_gap_ratio(speeds: tuple[int, ...]) -> tuple[Fraction, int]:
    pair, gaps = antipodal_pair_cover(speeds)
    max_gap = max((hi - lo for lo, hi in gaps), default=Fraction(0))
    return max_gap / Fraction(1, N), len(pair)


def row(label: str, raw_speeds: tuple[int, ...]) -> ReframeRow:
    speeds = canonical(raw_speeds)
    protection = S360.summarize(list(speeds))
    descent = S362.summarize(list(speeds))
    missing, fragile, _payers = moduli_ledger(speeds)
    pair_ratio, pair_components = pair_gap_ratio(speeds)
    unprotected = unprotected_values(speeds)
    depth_hist = Counter(extra_depth(point) for point in unprotected)
    frontier_mass = sum(Fraction(1, depth_scale(extra_depth(point))) for point in unprotected)
    denominator_pressure = sum(depth_scale(extra_depth(point)) for point in unprotected)
    return ReframeRow(
        label=label,
        speeds=speeds,
        classification=protection.classification,
        gap_ratio=protection.max_gap / protection.threshold,
        unprotected=protection.unprotected_count,
        peel_depth=len(descent.peel_layers),
        core_endpoints=descent.core_endpoint_count,
        missing_moduli=missing,
        fragile_moduli=fragile,
        pair_gap_ratio=pair_ratio,
        pair_components=pair_components,
        depth_hist=tuple(sorted(depth_hist.items())),
        frontier_mass=frontier_mass,
        denominator_pressure=denominator_pressure,
    )


def fmt_depth(depth: tuple[int, ...]) -> str:
    return "{" + ",".join(f"{prime}:+{height}" for prime, height in zip(PRIMES, depth)) + "}"


def fmt_depth_hist(hist: tuple[tuple[tuple[int, ...], int], ...]) -> str:
    if not hist:
        return "-"
    return " ".join(f"{fmt_depth(depth)}:{count}" for depth, count in hist)


def print_header(title: str) -> None:
    print("=" * 92)
    print(title)
    print("=" * 92)


def print_reframe_table(rows: tuple[ReframeRow, ...]) -> None:
    print_header("N=14 CANDIDATE ROWS THROUGH FOUR PROOF LENSES")
    print(
        "gap/th is the ordinary LRC complement gap.  pair-gap/th is the "
        "antipodal quotient gap: a gap there means some t or t+1/2 is safe."
    )
    print()
    header = (
        "label                    class          gap/th   pair/th  "
        "unprot peel core missing fragile"
    )
    print(header)
    print("-" * len(header))
    for item in rows:
        missing = ",".join(map(str, item.missing_moduli)) or "-"
        fragile = ",".join(map(str, item.fragile_moduli)) or "-"
        print(
            f"{item.label:<24} {item.classification:<14} "
            f"{fmt_float(item.gap_ratio):>7} {fmt_float(item.pair_gap_ratio):>8} "
            f"{item.unprotected:>6} {item.peel_depth:>4} {item.core_endpoints:>4} "
            f"{missing:<8} {fragile}"
        )
    print()


def print_product_depth_table(rows: tuple[ReframeRow, ...]) -> None:
    print_header("PRODUCT-DEPTH ENDPOINT DEBT")
    print(
        "Depths are extra denominator exponents beyond the base n=14 lattice. "
        "frontier_mass=sum 1/(2^a 7^b); denominator_pressure=sum 2^a 7^b."
    )
    print()
    header = "label                    unprotected depth histogram                  frontier  pressure"
    print(header)
    print("-" * len(header))
    for item in rows:
        print(
            f"{item.label:<24} {item.unprotected:>11} "
            f"{fmt_depth_hist(item.depth_hist):<32} "
            f"{fmt(item.frontier_mass):>8} {item.denominator_pressure:>9}"
        )
    print()


def owner_debt_ledger(label: str, raw_speeds: tuple[int, ...], limit: int = 10) -> None:
    speeds = canonical(raw_speeds)
    bad = unprotected_values(speeds)
    by_owner_labels: Counter[int] = Counter()
    by_owner_unique: dict[int, set[Fraction]] = defaultdict(set)
    for endpoint in S360.endpoints(speeds):
        if endpoint.value in bad:
            by_owner_labels[endpoint.speed] += 1
            by_owner_unique[endpoint.speed].add(endpoint.value)

    print(f"[{label}] owner-wise exposed endpoint labels")
    print("  owner  speed_mode        labels  unique")
    for owner, count in by_owner_labels.most_common(limit):
        mode = S420.speed_mode(owner)
        print(
            f"  {owner:>5}  2^{mode.dyadic_height}*{mode.odd_core:<7} "
            f"{count:>6} {len(by_owner_unique[owner]):>7}"
        )
    print()


def local_owner14_fan_tax() -> None:
    print_header("LOCAL 14-GATE FAN TAX")
    targets = S420.endpoint_values_for_owner(N, N)
    candidates = tuple(range(1, N))
    result = S420.set_cover_result(
        "owner 14 endpoints, lower 1..13",
        N,
        targets,
        candidates,
    )
    print(
        "A 14-gate has 28 endpoints.  If it is maximal, its endpoints must be "
        "covered by lower columns.  This local set cover is the n=14 analogue "
        "of the n=16 nine-cover, but with a 2 x 7 CRT flavor."
    )
    print()
    print(f"exact lower-cover size={result.exact_size} columns={result.exact_columns}")
    print(f"forced columns by private endpoint rows={result.forced_columns}")
    print()

    columns = S420.build_cover_columns(N, targets, candidates)
    union_without: dict[int, int] = {}
    for column in columns:
        mask = 0
        for other in columns:
            if other.speed != column.speed:
                mask |= other.mask
        union_without[column.speed] = mask

    print("  p  mode         covers private private_depths")
    for column in columns:
        private_mask = column.mask & ~union_without[column.speed]
        private_depths: Counter[tuple[int, ...]] = Counter()
        work = private_mask
        while work:
            bit = work & -work
            target_index = bit.bit_length() - 1
            private_depths[extra_depth(targets[target_index])] += 1
            work -= bit
        mode = S420.speed_mode(column.speed)
        private_text = fmt_depth_hist(tuple(sorted(private_depths.items())))
        print(
            f"  {column.speed:>2}  2^{mode.dyadic_height}*{mode.odd_core:<5} "
            f"{column.size:>6} {private_mask.bit_count():>7} {private_text}"
        )
    print()
    print(
        "Reframe: a single 14-gate locally demands the six unit residue columns, "
        "the half-gate 7, and at least one even bridge.  So a global proof can "
        "charge every 14-gate against a heptagonal fan before counting the "
        "deeper product-depth endpoints it creates."
    )
    print()


def print_owner_ledgers() -> None:
    print_header("OWNER-WISE ENDPOINT DEBT")
    print(
        "Unprotected endpoint values are counted back against the speed whose "
        "forbidden interval owns that endpoint.  This identifies which gates "
        "export the debt."
    )
    print()
    owner_debt_ledger("seven-ladder", seven_ladder())
    owner_debt_ledger("S380 gate ladder", s380_gate_ladder())
    owner_debt_ledger("drop13 add182", drop_add(13, 182))


def print_antipodal_interpretation(rows: tuple[ReframeRow, ...]) -> None:
    print_header("ANTIPODAL / CRT FOLD INTERPRETATION")
    print(
        "The quotient t~t+1/2 is not a new theorem by itself, but it forces "
        "an even/odd split.  Even speeds cover both sides of a pair in the "
        "same way.  Odd speeds must cover one side through a low band and the "
        "other through a high band.  At n=14 this pairs the 2-adic fold with "
        "the mod-7 unit cycle."
    )
    print()
    for item in rows:
        even = sum(1 for speed in item.speeds if speed % 2 == 0)
        odd = len(item.speeds) - even
        gates = sum(1 for speed in item.speeds if speed % N == 0)
        print(
            f"{item.label:<24} odd={odd:>2} even={even:>2} "
            f"14-gates={gates:>2} pair_components={item.pair_components:>4} "
            f"pair_gap/th={fmt_float(item.pair_gap_ratio)}"
        )
    print()
    print(
        "Creative target: prove that any n=14 row with enough even/14-gate "
        "mass to kill unit points leaves an antipodal high-band deficit unless "
        "it imports odd unit columns; those odd columns are exactly the local "
        "fan-tax columns in the 14-gate cover."
    )
    print()


def print_new_models() -> None:
    print_header("NEW PROOF MODELS GENERATED")
    print("1. Gate fan tax.")
    print(
        "   Before a 14-gate can help globally, its own endpoints demand a "
        "local cover of size eight: six unit residues, 7, and an even bridge."
    )
    print("2. Product-depth debt.")
    print(
        "   The S380 ladder pays every small-denominator row and uses only "
        "forced columns, but leaves debt at {2:+1,7:+1} and {2:+3,7:+1}."
    )
    print("3. Antipodal CRT fold.")
    print(
        "   The even fold t~t+1/2 separates even columns from odd high/low "
        "columns.  For n=14 the missing proof should couple this with the "
        "mod-7 unit cycle."
    )
    print("4. Owner-charge descent.")
    print(
        "   Endpoint debt should be charged to owner speeds, not just to the "
        "set as a whole.  Gate-heavy candidates expose which owners export "
        "which denominator depths."
    )
    print()


def main() -> None:
    examples = (
        row("initial", initial()),
        row("boundary 12->24", boundary_swap_12_24()),
        row("drop13 add14", drop_add(13, 14)),
        row("drop13 add182", drop_add(13, 182)),
        row("seven-ladder", seven_ladder()),
        row("seven with 42", seven_with_42()),
        row("S380 gate ladder", s380_gate_ladder()),
    )

    print("n=14 LRC creative proof reframes (codex-2026-05-31 S440)")
    print("All endpoint and interval computations are exact over Fraction.\n")
    print_reframe_table(examples)
    print_product_depth_table(examples)
    local_owner14_fan_tax()
    print_owner_ledgers()
    print_antipodal_interpretation(examples)
    print_new_models()


if __name__ == "__main__":
    main()
