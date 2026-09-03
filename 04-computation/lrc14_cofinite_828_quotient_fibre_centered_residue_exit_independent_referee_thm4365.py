#!/usr/bin/env python3
"""Import-free exact referee for the proposed cofinite THM-4363 fibre.

The script reconstructs the fixed open-tooth cover from the definitions in
tracked THM-4363.  It does not import any project module or scout artifact.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")

H = 420
ANCHOR = 840
FIXED = (3, 39, 11, 1691, 3371, 5051, 6731, 8411, 10091, 525, 945)
COLLARS = (19, 20, 259, 260, 299, 300, 539, 540, 579, 580, 819, 820)
RESIDUAL = tuple(k for k in range(840) if k not in set(COLLARS))

MOD = 47194
X1 = Fraction(1303, MOD)
NEXT_LEFT = Fraction(2603, 94234)  # left wall of 6731@186

EXPECTED_STATUS_HASH = "4d34447b9eca8c8a9302a0f799a56300b4b96135cd0c3a245a97960975f9a347"
EXPECTED_COMPLETED_HASH = "64704c712a0e3e70ca0ecb3264834b3f61c4d473ddc20ea8c44ccc5c2c616d11"


@dataclass(frozen=True, order=True)
class Tooth:
    speed: int
    address: int
    left: Fraction
    right: Fraction


@dataclass(frozen=True)
class Trace:
    status: str
    chain: tuple[Tooth, ...]
    exit: Fraction | None
    extensions: tuple[Fraction, ...]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def floor_fraction(x: Fraction) -> int:
    return x.numerator // x.denominator


def ceil_fraction(x: Fraction) -> int:
    return -((-x.numerator) // x.denominator)


def component(k: int) -> tuple[Fraction, Fraction]:
    return Fraction(14 * k + 1, 11760), Fraction(14 * k + 13, 11760)


def tooth(speed: int, address: int) -> Tooth:
    return Tooth(
        speed,
        address,
        Fraction(14 * address - 1, 14 * speed),
        Fraction(14 * address + 1, 14 * speed),
    )


def teeth_meeting(speed: int, left: Fraction, right: Fraction) -> tuple[Tooth, ...]:
    lo = floor_fraction(speed * left) - 2
    hi = ceil_fraction(speed * right) + 2
    return tuple(
        z
        for n in range(lo, hi + 1)
        for z in (tooth(speed, n),)
        if z.left < right and left < z.right
    )


def fixed_bank(k: int) -> tuple[Tooth, ...]:
    left, right = component(k)
    return tuple(z for speed in FIXED for z in teeth_meeting(speed, left, right))


def trace_from_bank(bank: tuple[Tooth, ...], k: int) -> Trace:
    left, right = component(k)
    cursor = left
    chain: list[Tooth] = []
    extensions: list[Fraction] = []
    while cursor <= right:
        active = tuple(z for z in bank if z.left < cursor < z.right)
        if not active:
            return Trace("missing", tuple(chain), cursor, tuple(extensions))
        farthest = max(z.right for z in active)
        chosen = min(
            (z for z in active if z.right == farthest),
            key=lambda z: (z.speed, z.address),
        )
        require(chosen.right > cursor, "selected tooth failed to advance")
        extensions.append(chosen.right - cursor)
        chain.append(chosen)
        cursor = chosen.right
    return Trace("span" if len(chain) == 1 else "renew", tuple(chain), None, tuple(extensions))


def trace_with_parameter(k: int, p: int, cached_bank: tuple[Tooth, ...]) -> Trace:
    left, right = component(k)
    return trace_from_bank(cached_bank + teeth_meeting(p, left, right), k)


def positive_fixed_gaps(bank: tuple[Tooth, ...], k: int) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    """Positive-length components of I_k minus the union of fixed teeth.

    Each returned tuple is (length,left,right).  Zero-length open-endpoint
    holes are irrelevant to the no-completion certificates below: finding
    one positive gap not coverable by P already forces the component missing.
    """
    left, right = component(k)
    intervals = sorted((max(left, z.left), min(right, z.right)) for z in bank)
    cursor = left
    ans: list[tuple[Fraction, Fraction, Fraction]] = []
    for a, b in intervals:
        if a > cursor:
            ans.append((a - cursor, cursor, a))
        if b > cursor:
            cursor = b
    if cursor < right:
        ans.append((right - cursor, cursor, right))
    return tuple(ans)


def one_p_tooth_strictly_contains_gap(p: int, a: Fraction, b: Fraction) -> bool:
    """Whether one open P-tooth contains both endpoints of [a,b]."""
    # (14n-1)/(14p) < a and b < (14n+1)/(14p), equivalently
    # p*b-1/14 < n < p*a+1/14.
    lower = p * b - Fraction(1, 14)
    upper = p * a + Fraction(1, 14)
    least = floor_fraction(lower) + 1
    greatest = ceil_fraction(upper) - 1
    return least <= greatest


def physical_chain(trace: Trace) -> tuple[tuple[int, int], ...]:
    return tuple((z.speed, z.address) for z in trace.chain)


def centered_rho(p: int) -> int:
    raw = (1303 * p) % MOD
    # Declared convention: -MOD/2 < rho <= MOD/2.  At the ambiguous
    # half-modulus class this chooses +23597.
    return raw if raw <= MOD // 2 else raw - MOD


def predicted_k23_exit(p: int) -> tuple[int, int, Fraction]:
    rho = centered_rho(p)
    n = (1303 * p - rho) // MOD
    if abs(rho) >= 3371:
        return rho, n, X1
    return rho, n, X1 + Fraction(3371 - rho, MOD * p)


def circle_distance(speed: int, x: Fraction) -> Fraction:
    y = speed * x
    r = y - floor_fraction(y)
    return min(r, 1 - r)


def binding_set(p: int, x: Fraction) -> tuple[Fraction, tuple[int, ...]]:
    vals = tuple((w, circle_distance(w, x)) for w in (ANCHOR,) + FIXED + (p,))
    level = min(v for _, v in vals)
    return level, tuple(w for w, v in vals if v == level)


def status_bytes(traces: dict[int, Trace]) -> bytes:
    code = {"missing": "M", "span": "S", "renew": "R"}
    return "".join(f"{k}|{code[traces[k].status]}\n" for k in RESIDUAL).encode("ascii")


def completed_bytes(traces: dict[int, Trace]) -> bytes:
    lines = []
    for k in RESIDUAL:
        tr = traces[k]
        chain = "" if tr.status == "missing" else ",".join(
            f"{z.speed}@{z.address}" for z in tr.chain
        )
        lines.append(f"{k}|{tr.status}|{chain}\n")
    return "".join(lines).encode("ascii")


def main() -> None:
    require(len(RESIDUAL) == 828, "residual address count changed")
    require(gcd(1303, MOD) == 1, "1303 is not a unit modulo 47194")

    banks = {k: fixed_bank(k) for k in RESIDUAL}
    fixed = {k: trace_from_bank(banks[k], k) for k in RESIDUAL}
    census = Counter(tr.status for tr in fixed.values())
    require(census == Counter({"missing": 546, "span": 276, "renew": 6}), "fixed census changed")

    status_stream = status_bytes(fixed)
    completed_stream = completed_bytes(fixed)
    status_hash = sha256(status_stream).hexdigest()
    completed_hash = sha256(completed_stream).hexdigest()
    require((len(status_stream), status_hash) == (4860, EXPECTED_STATUS_HASH), "status quotient mismatch")
    require((len(completed_stream), completed_hash) == (10940, EXPECTED_COMPLETED_HASH), "completed quotient mismatch")
    require(next(k for k in RESIDUAL if fixed[k].status == "missing") == 23, "fixed first missing address changed")

    extension_records = tuple(
        (extension, k, j, fixed[k].chain[j].speed, fixed[k].chain[j].address)
        for k in RESIDUAL
        if fixed[k].status != "missing"
        for j, extension in enumerate(fixed[k].extensions)
    )
    min_extension = min(x[0] for x in extension_records)
    extension_minimizers = tuple(x for x in extension_records if x[0] == min_extension)
    require(min_extension == Fraction(1, 28665), "completed reach margin changed")
    require(
        extension_minimizers
        == (
            (Fraction(1, 28665), 216, 1, 945, 244),
            (Fraction(1, 28665), 496, 1, 945, 559),
            (Fraction(1, 28665), 776, 1, 945, 874),
        ),
        "completed reach minimizers changed",
    )

    gaps = {
        k: tuple(sorted(positive_fixed_gaps(banks[k], k), reverse=True))
        for k in RESIDUAL
        if fixed[k].status == "missing"
    }
    require(all(gaps.values()), "a fixed-missing component has no positive gap")
    largest = {k: rows[0] for k, rows in gaps.items()}
    min_largest_length = min(row[0] for row in largest.values())
    min_largest_records = tuple((k,) + largest[k] for k in sorted(largest) if largest[k][0] == min_largest_length)
    require(min_largest_length == Fraction(7003, 594127807), "largest-gap margin changed")
    require(
        min_largest_records
        == (
            (410, Fraction(7003, 594127807), Fraction(57597, 117754), Fraction(69103, 141274)),
            (429, Fraction(7003, 594127807), Fraction(72171, 141274), Fraction(60157, 117754)),
        ),
        "largest-gap minimizers changed",
    )
    diameter_threshold = Fraction(1, 7 * min_largest_length)
    require(diameter_threshold == Fraction(84875401, 7003), "diameter threshold changed")
    require(12119 < diameter_threshold < 12120, "integer diameter threshold changed")
    require(Fraction(1, 7 * 12121) < min_largest_length, "P=12121 tail inequality failed")

    # Completed chains are automatically stable already for P>=4095: an
    # active P-tooth can advance a cursor by strictly less than its diameter.
    require(Fraction(1, 7 * 4095) == min_extension, "completed-chain threshold changed")

    # Alignment-free tail: the largest fixed gap in each missing component is
    # longer than every P-tooth for P>=12121, so no missing component closes.
    # Exact finite bridge 11019..12119: for every row/component exhibit at
    # least one positive fixed gap that no open P-tooth strictly contains.
    finite_safe_count = 0
    finite_min_uncertified = None
    for p in range(11019, 12120, 2):
        require(p > max(FIXED), "parameter is not a new distinct speed")
        for k, rows in gaps.items():
            witness = next(
                ((length, a, b) for length, a, b in rows if not one_p_tooth_strictly_contains_gap(p, a, b)),
                None,
            )
            require(witness is not None, f"finite bridge closes fixed-missing component P={p}, k={k}")
        finite_safe_count += 1
        if finite_min_uncertified is None:
            finite_min_uncertified = p
    require(finite_safe_count == 551 and finite_min_uncertified == 11019, "finite bridge range changed")

    # The immediate predecessor is a genuine hostile row: exactly two fixed-
    # missing components become renew components, so 11019 is the least odd
    # cofinite threshold, while 12121 is only the alignment-free threshold.
    p_bad = 11017
    bad_traces = {k: trace_with_parameter(k, p_bad, banks[k]) for k in RESIDUAL}
    bad_differences = tuple(k for k in RESIDUAL if (
        bad_traces[k].status != fixed[k].status
        or (
            fixed[k].status != "missing"
            and physical_chain(bad_traces[k]) != physical_chain(fixed[k])
        )
    ))
    require(bad_differences == (147, 692), "P=11017 hostile address set changed")
    require(
        physical_chain(bad_traces[147])
        == ((1691, 296), (11017, 1929), (10091, 1767), (525, 92), (11, 2)),
        "P=11017 k=147 hostile chain changed",
    )
    require(
        physical_chain(bad_traces[692])
        == ((11, 9), (525, 433), (10091, 8324), (11017, 9088), (1691, 1395)),
        "P=11017 k=692 hostile chain changed",
    )
    for p in (11019, 12119, 12121):
        # Direct recurrence controls at the two junctions and analytic-tail start.
        traces = {k: trace_with_parameter(k, p, banks[k]) for k in RESIDUAL}
        require(status_bytes(traces) == status_stream, f"direct status control failed P={p}")
        require(completed_bytes(traces) == completed_stream, f"direct completed control failed P={p}")

    # First missing trace.  The first two fixed reaches dominate all P-tooth
    # advances for P>=10141.
    left23, _ = component(23)
    x0 = tooth(945, 26).right
    prefix_margins = (x0 - left23, X1 - x0)
    require(prefix_margins == (Fraction(13, 105840), Fraction(92, 4459833)), "k=23 prefix margins changed")
    require(Fraction(1, 7 * 10141) < min(prefix_margins), "k=23 prefix stability failed")
    require(tooth(3371, 93).right == X1, "x1 fixed endpoint changed")
    require(tooth(6731, 186).left == NEXT_LEFT, "next fixed left wall changed")
    future_fixed_lefts = tuple(sorted(
        z.left for z in banks[23] if z.left >= X1 and z.right > X1
    ))
    require(future_fixed_lefts and future_fixed_lefts[0] == NEXT_LEFT > X1, "next fixed start is not exact")
    next_gap = NEXT_LEFT - X1
    require(next_gap == Fraction(2110, 158831407), "post-x1 fixed gap changed")
    universal_exit_threshold = Fraction(1, 7 * next_gap)
    require(universal_exit_threshold == Fraction(22690201, 2110), "universal exit threshold changed")
    require(10753 < universal_exit_threshold < 10754, "universal exit integer bracket changed")

    finite_exit_active = 0
    finite_exit_min_margin = None
    finite_exit_min_record = None
    for p in range(10141, 10754, 2):
        rho, n, predicted = predicted_k23_exit(p)
        tr = trace_with_parameter(23, p, banks[23])
        require(tr.status == "missing" and tr.exit == predicted, f"finite k23 law failed P={p}")
        require(physical_chain(tr)[:2] == ((945, 26), (3371, 93)), "fixed k23 prefix changed")
        if abs(rho) < 3371:
            finite_exit_active += 1
            margin = NEXT_LEFT - predicted
            require(margin > 0, f"active P tooth reaches next fixed tooth P={p}")
            if finite_exit_min_margin is None or margin < finite_exit_min_margin:
                finite_exit_min_margin = margin
                finite_exit_min_record = (p, rho, n, margin)
            require(physical_chain(tr)[-1] == (p, n), "active terminal tooth changed")
        else:
            require(predicted == X1 and physical_chain(tr) == ((945, 26), (3371, 93)), "inactive trace changed")
        level, binders = binding_set(p, predicted)
        require(level == Fraction(1, 14), f"exit clearance changed P={p}")
        if abs(rho) < 3371:
            require(binders == (p,), f"active exit binder changed P={p}")
        elif abs(rho) == 3371:
            require(set(binders) == {3371, p}, f"boundary binder set changed P={p}")
        else:
            require(binders == (3371,), f"inactive exit binder changed P={p}")
    require(
        finite_exit_min_record == (10465, -3171, 289, Fraction(19, 493079405)),
        "finite active post-wall margin changed",
    )

    # P>=10755 is automatic because a P-tooth's whole diameter is below the
    # next fixed gap.  P=10139 shows the consecutive-odd lower bound is sharp.
    require(Fraction(1, 7 * 10755) < next_gap, "analytic k23 tail inequality failed")
    p_hostile = 10139
    hostile_rho, hostile_n, hostile_prediction = predicted_k23_exit(p_hostile)
    hostile_trace = trace_with_parameter(23, p_hostile, banks[23])
    hostile_cross = hostile_prediction - NEXT_LEFT
    require((hostile_rho, hostile_n) == (-3203, 280), "P=10139 residue data changed")
    require(hostile_cross == Fraction(31, 68245609) > 0, "P=10139 crossing margin changed")
    require(
        physical_chain(hostile_trace)
        == ((945, 26), (3371, 93), (10139, 280), (6731, 186), (10091, 279)),
        "P=10139 hostile chain changed",
    )
    require(hostile_trace.exit == Fraction(3907, 141274) != hostile_prediction, "P=10139 is not hostile")

    # Residue count, including the centered half-modulus convention.
    odd_classes = tuple(range(1, MOD, 2))
    rho_values = tuple(centered_rho(p) for p in odd_classes)
    require(len(set(rho_values)) == MOD // 2, "odd residue map is not bijective")
    require(23597 in rho_values and -23597 not in rho_values, "half-modulus convention changed")
    active_classes = tuple(p for p in odd_classes if abs(centered_rho(p)) < 3371)
    boundary_classes = tuple(p for p in odd_classes if abs(centered_rho(p)) == 3371)
    require(len(active_classes) == 3370, "active odd residue count changed")
    require(len(boundary_classes) == 2, "activity boundary class count changed")
    require(boundary_classes == (3371, 43823), "boundary residue classes changed")
    boundary_controls = []
    for p in (43823, 50565):
        rho, n, x = predicted_k23_exit(p)
        tr = trace_with_parameter(23, p, banks[23])
        level, binders = binding_set(p, x)
        require(abs(rho) == 3371 and x == X1 and tr.exit == X1, "boundary endpoint trace changed")
        require(level == Fraction(1, 14) and set(binders) == {3371, p}, "boundary endpoint binders changed")
        boundary_controls.append((p, rho, n, binders))
    inverse_1303 = pow(1303, -1, MOD)
    require(inverse_1303 == 1485 and centered_rho(inverse_1303) == 1, "rho=1 class changed")

    # Every P=1485+47194*t for t>=1 is odd, lies in the cofinite quotient
    # fibre, and has the distinct exit x1+3370/(47194 P).
    infinite_controls = []
    for t in (1, 2, 17):
        p = inverse_1303 + MOD * t
        rho, n, x = predicted_k23_exit(p)
        require(p >= 11019 and rho == 1, "infinite active class control failed")
        require(x == X1 + Fraction(3370, MOD * p), "infinite exit formula failed")
        infinite_controls.append((p, n, x))
    require(len({x for _, _, x in infinite_controls}) == len(infinite_controls), "sample exits collide")

    # Primitivity/distinctness in both theorem ranges is automatic.
    require(gcd(ANCHOR, *FIXED) == 1, "fixed row is not primitive")
    require(max(FIXED) == 10091 < 10139, "parameter separation changed")
    require(len(FIXED) == len(set(FIXED)), "fixed tail has a duplicate")

    print("THM4365 COFINITE FIBRE CLEANROOM REFEREE")
    print("fixed_census", dict(sorted(census.items())))
    print("status_stream", len(status_stream), status_hash)
    print("completed_stream", len(completed_stream), completed_hash)
    print("completed_min_reach", min_extension, extension_minimizers)
    print("missing_min_largest_gap", min_largest_length, min_largest_records)
    print("alignment_free_threshold", diameter_threshold, "least_integer", 12120, "least_odd", 12121)
    print("exact_odd_cofinite_threshold", 11019, "finite_rows_checked", finite_safe_count)
    print("hostile_q_predecessor", p_bad, bad_differences)
    print("hostile_q_chain_147", physical_chain(bad_traces[147]))
    print("hostile_q_chain_692", physical_chain(bad_traces[692]))
    print("k23_prefix_margins", prefix_margins)
    print("k23_next_fixed_gap", next_gap)
    print("k23_analytic_tail_threshold", universal_exit_threshold, "least_odd", 10755)
    print("k23_finite_band", "10141..10753", "active", finite_exit_active)
    print("k23_finite_min_margin", finite_exit_min_record)
    print("k23_hostile_10139", hostile_rho, hostile_n, hostile_cross, physical_chain(hostile_trace), hostile_trace.exit)
    print("centered_convention", "-23597<rho<=23597", "half_class", 23597)
    print("odd_classes", len(odd_classes), "active", len(active_classes), "boundary", len(boundary_classes))
    print("boundary_controls", tuple(boundary_controls))
    print("rho_one_class", inverse_1303, "infinite_controls", tuple(infinite_controls))
    print("endpoint_rule", "open; |rho|=3371 inactive and safe with binders {3371,P}")
    print("scope", "Q stable for every odd P>=11019; 12121 is sufficient but not minimal")
    print("PASS")


if __name__ == "__main__":
    main()
