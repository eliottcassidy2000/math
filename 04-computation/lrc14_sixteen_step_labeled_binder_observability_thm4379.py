#!/usr/bin/env python3
"""Primary exact verifier for THM-4379's labeled first-exit binder observer.

The emitted record is deliberately physical and labelled.  It keeps the
component/copy data, the selected tooth chain with Euclidean addresses and
owner parities, the first uncovered point, the clearance, and every binding
tooth with its wall side.  In particular, a boundary tie really contains the
numeric varying speed; it is not replaced by the abstract role name ``P``.

The verifier has two independent paths for the local record:

* the centered-residue formula inherited from THM-4365; and
* a direct exact open-tooth renewal trace followed by a distance calculation
  against all thirteen relative speeds.

It then exhausts the 23,597 half-phases (at two distinct scales), proves the
closed marked-interval cover under P -> P+2, audits all THM-4367 structural
active scale-fibre shapes, and checks sharp hostile and boundary rows.
Nothing here is an LRC(14) decrement.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from math import ceil, floor, gcd
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")


A = 3371
S = 1303
M = 14 * A
N = M // 2
TAIL = 11019
H = 420
ANCHOR = 840
FIXED = (3, 39, 11, 1691, A, 5051, 6731, 8411, 10091, 525, 945)
SPEEDS_WITHOUT_PARAMETER = (ANCHOR,) + FIXED
UNITS14 = (1, 3, 5, 9, 11, 13)
BASE_EXIT = Fraction(S, M)
CLEARANCE = Fraction(1, 14)
COMPONENT_K = 23
COPY_BIT = 0
INV_S = pow(S, -1, M)

CHECKS = 0


def require(condition: bool, message: str) -> None:
    """Optimization-stable assertion."""
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("check failed: " + message)


def centered_rho(parameter: int) -> int:
    """The representative of S*parameter modulo M in (-M/2,M/2]."""
    value = (S * parameter) % M
    if 2 * value > M:
        value -= M
    return value


def phase(parameter: int) -> int:
    """Half-coordinate z=(A-rho)/2 in Z/NZ."""
    rho = centered_rho(parameter)
    require((A - rho) % 2 == 0, "half-phase parity")
    return ((A - rho) // 2) % N


def euclidean_address(parameter: int) -> int:
    rho = centered_rho(parameter)
    numerator = S * parameter - rho
    require(numerator % M == 0, "Euclidean address integrality")
    return numerator // M


def owner_parity(speed: int, address: int) -> int:
    """THM-4363 owner parity (2n-epsilon*w) mod 2 at epsilon=0."""
    return (2 * address - COPY_BIT * speed) % 2


@dataclass(frozen=True)
class Tooth:
    speed: int
    address: int
    left: Fraction
    right: Fraction


@dataclass(frozen=True)
class Trace:
    status: str
    chain: tuple[Tooth, ...]
    frontiers: tuple[Fraction, ...]
    exit_point: Fraction


@dataclass(frozen=True)
class Record:
    """The declared complete physical first-exit record."""
    component: int
    component_mod_h: int
    copy_bit: int
    status: str
    # Each selected tooth is (numeric speed, Euclidean address, owner parity).
    chain: tuple[tuple[int, int, int], ...]
    frontiers: tuple[Fraction, ...]
    first_uncovered: Fraction
    clearance: Fraction
    # Each binder is (numeric speed, address, L/R wall side, owner parity).
    binders: tuple[tuple[int, int, str, int], ...]


def tooth(speed: int, address: int) -> Tooth:
    return Tooth(
        speed,
        address,
        Fraction(14 * address - 1, 14 * speed),
        Fraction(14 * address + 1, 14 * speed),
    )


def component(k: int) -> tuple[Fraction, Fraction]:
    return Fraction(14 * k + 1, 28 * H), Fraction(14 * k + 13, 28 * H)


def meeting_teeth(speed: int, left: Fraction, right: Fraction) -> tuple[Tooth, ...]:
    lo = floor(speed * left - CLEARANCE)
    hi = ceil(speed * right + CLEARANCE)
    answer = []
    for address in range(lo - 1, hi + 2):
        item = tooth(speed, address)
        if item.left < right and left < item.right:
            answer.append(item)
    return tuple(answer)


COMPONENT_LEFT, COMPONENT_RIGHT = component(COMPONENT_K)
FIXED_TEETH = tuple(
    item
    for speed in FIXED
    for item in meeting_teeth(speed, COMPONENT_LEFT, COMPONENT_RIGHT)
)


def selection_key(item: Tooth) -> tuple[Fraction, Fraction, int, int]:
    # THM-4365's farthest-right rule and its exact deterministic tie order.
    return item.right, -item.left, -item.speed, -item.address


def active_parameter_tooth(parameter: int, frontier: Fraction) -> Tooth | None:
    numerator = parameter * frontier.numerator
    denominator = frontier.denominator
    quotient = numerator // denominator
    active = []
    for address in (quotient, quotient + 1):
        if 14 * abs(numerator - address * denominator) < denominator:
            active.append(tooth(parameter, address))
    require(len(active) <= 1, "one speed has two active teeth")
    return active[0] if active else None


def direct_trace(parameter: int) -> Trace:
    """Independent component-23 renewal trace with exact open endpoints."""
    frontier = COMPONENT_LEFT
    chain: list[Tooth] = []
    frontiers: list[Fraction] = []
    while True:
        active = [item for item in FIXED_TEETH if item.left < frontier < item.right]
        extra = active_parameter_tooth(parameter, frontier)
        if extra is not None:
            active.append(extra)
        if not active:
            return Trace("missing", tuple(chain), tuple(frontiers), frontier)
        winner = max(active, key=selection_key)
        frontiers.append(frontier)
        chain.append(winner)
        frontier = winner.right
        if frontier > COMPONENT_RIGHT:
            status = "span" if len(chain) == 1 else "renew"
            return Trace(status, tuple(chain), tuple(frontiers), frontier)


def circle_distance(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def binding_teeth(parameter: int, exit_point: Fraction) -> tuple[Fraction, tuple[tuple[int, int, str, int], ...]]:
    """Compute clearance and all physical binders directly from distances."""
    values = tuple(
        (speed, circle_distance(speed * exit_point))
        for speed in SPEEDS_WITHOUT_PARAMETER + (parameter,)
    )
    clearance = min(value for _, value in values)
    binders = []
    for speed, value in values:
        if value != clearance:
            continue
        position = speed * exit_point
        lower = position.numerator // position.denominator
        residue = position - lower
        if residue == CLEARANCE:
            address, side = lower, "R"
        elif residue == 1 - CLEARANCE:
            address, side = lower + 1, "L"
        else:
            raise RuntimeError(
                f"binding point for speed {speed} is not a 1/14 wall: {position}"
            )
        binders.append((speed, address, side, owner_parity(speed, address)))
    return clearance, tuple(sorted(binders))


def record_from_direct_geometry(parameter: int) -> Record:
    row = direct_trace(parameter)
    clearance, binders = binding_teeth(parameter, row.exit_point)
    return Record(
        COMPONENT_K,
        COMPONENT_K % H,
        COPY_BIT,
        row.status,
        tuple(
            (item.speed, item.address, owner_parity(item.speed, item.address))
            for item in row.chain
        ),
        row.frontiers,
        row.exit_point,
        clearance,
        binders,
    )


FIXED_CHAIN = ((945, 26, 0), (A, 93, 0))
FIXED_FRONTIERS = (COMPONENT_LEFT, Fraction(73, 2646))
FIXED_BINDER = (A, 93, "R", 0)


def record_from_centered_law(parameter: int) -> Record:
    """Closed formula obtained from THM-4365, including boundary ties."""
    rho = centered_rho(parameter)
    address = euclidean_address(parameter)
    z = phase(parameter)
    strict_active = 1 <= z <= A - 1

    chain = FIXED_CHAIN
    frontiers = FIXED_FRONTIERS
    if strict_active:
        chain += ((parameter, address, owner_parity(parameter, address)),)
        frontiers += (BASE_EXIT,)
        exit_point = Fraction(14 * address + 1, 14 * parameter)
        binders = ((parameter, address, "R", owner_parity(parameter, address)),)
    elif z == 0:
        require(rho == A, "positive boundary phase")
        exit_point = BASE_EXIT
        binders = tuple(
            sorted(
                (
                    FIXED_BINDER,
                    (parameter, address, "R", owner_parity(parameter, address)),
                )
            )
        )
    elif z == A:
        require(rho == -A, "negative boundary phase")
        exit_point = BASE_EXIT
        binders = tuple(
            sorted(
                (
                    FIXED_BINDER,
                    (parameter, address, "L", owner_parity(parameter, address)),
                )
            )
        )
    else:
        require(abs(rho) > A, "strict inactive phase")
        exit_point = BASE_EXIT
        binders = (FIXED_BINDER,)

    return Record(
        COMPONENT_K,
        COMPONENT_K % H,
        COPY_BIT,
        "missing",
        chain,
        frontiers,
        exit_point,
        CLEARANCE,
        binders,
    )


def parameter_is_visible(record: Record, parameter: int) -> bool:
    return any(speed == parameter for speed, _, _, _ in record.binders)


def marked_phase(z: int) -> bool:
    """The varying speed binds iff z is in the closed interval J=[0,A]."""
    return 0 <= z % N <= A


def first_marked_hit(z: int) -> int:
    for time in range(17):
        if marked_phase(z - S * time):
            return time
    raise RuntimeError(f"phase {z} has no marked hit through time 16")


def phase_representative(z: int, lift: int = 0) -> int:
    """Least odd tail parameter in phase z, followed by an optional M-lift."""
    rho_mod_m = (A - 2 * (z % N)) % M
    parameter = (INV_S * rho_mod_m) % M
    if parameter == 0:
        parameter = M
    while parameter < TAIL:
        parameter += M
    parameter += lift * M
    require(parameter % 2 == 1 and parameter >= TAIL, "tail representative")
    require(phase(parameter) == z % N, "phase inverse")
    return parameter


def output_word(parameter: int, horizon: int) -> tuple[Record, ...]:
    return tuple(record_from_centered_law(parameter + 2 * time) for time in range(horizon + 1))


def decode_visible_parameter(word: tuple[Record, ...]) -> int | None:
    """Decode P from the first numeric nonfixed binder q=P+2t."""
    for time, record in enumerate(word):
        varying = tuple(
            speed
            for speed, _, _, _ in record.binders
            if speed != A
        )
        if varying:
            require(len(varying) == 1, "unique varying binder in one record")
            return varying[0] - 2 * time
    return None


def representative_metric_class(a: int, g0: int) -> tuple[int, int]:
    """One literal admissible (b,kappa) for each THM-4367 scale shape."""
    kappa_mod_14 = pow(g0, -1, 14)
    b = (INV_S * (A * kappa_mod_14 - a)) % M
    if b == 0:
        b = M
    while b < TAIL or gcd(a, b) != 1:
        b += M
    numerator = a + S * b
    require(numerator % A == 0, "metric class divisibility")
    kappa = numerator // A
    require(b % 2 == 1 and gcd(a, b) == 1, "metric class primitive parity")
    require((g0 * kappa) % 14 == 1, "metric class scale residue")
    return b, kappa


def full_record_role_projection(parameter: int) -> tuple[object, ...]:
    """Replace every varying chain/binder tooth by its role, erasing q and n."""
    record = record_from_centered_law(parameter)
    projected_chain = tuple(
        (("P", owner) if speed == parameter else (speed, address, owner))
        for speed, address, owner in record.chain
    )
    projected_binders = tuple(
        (
            ("P",) if speed == parameter else (speed, address),
            side,
            owner,
        )
        for speed, address, side, owner in record.binders
    )
    return (
        record.component,
        record.component_mod_h,
        record.copy_bit,
        record.status,
        projected_chain,
        record.frontiers,
        record.first_uncovered,
        record.clearance,
        projected_binders,
    )


def main() -> None:
    require((A, S, M, N, TAIL, INV_S) == (3371, 1303, 47194, 23597, 11019, 1485), "published constants")
    require(gcd(S, M) == 1 and N == 7 * A, "primitive rotation")
    require(COMPONENT_LEFT == Fraction(323, 11760), "component-23 left endpoint")
    require(FIXED_FRONTIERS[1] == tooth(945, 26).right, "first fixed frontier")
    require(BASE_EXIT == tooth(A, 93).right, "fixed D2 exit")

    # Direct local-geometry replay at every half-phase and two genuinely
    # distinct scale lifts.  This attacks the formula, rather than merely
    # re-evaluating its own cases.
    geometry_rows = 0
    owner_values = set()
    for z in range(N):
        for lift in (0, 1):
            parameter = phase_representative(z, lift)
            predicted = record_from_centered_law(parameter)
            direct = record_from_direct_geometry(parameter)
            require(direct == predicted, f"direct record mismatch P={parameter}")
            require(direct.status == "missing", "component 23 must remain missing")
            require(direct.clearance == CLEARANCE, "safe boundary clearance")
            owner_values.update(owner for _, _, owner in direct.chain)
            owner_values.update(owner for _, _, _, owner in direct.binders)
            geometry_rows += 1
    require(geometry_rows == 2 * N, "two-scale phase geometry universe")
    require(owner_values == {0}, "component-23 owner parity is constant zero")

    # Exact closed marked-interval return clock.  J includes both open-tooth
    # boundary ties z=0,A, unlike THM-4374's strict metric-active interval.
    hit_counts = Counter()
    unresolved_counts = {}
    baseline = record_from_centered_law(11123)
    require(not parameter_is_visible(baseline, 11123), "baseline hostile is fixed-only")
    for z in range(N):
        hit = first_marked_hit(z)
        hit_counts[hit] += 1
        require(marked_phase(z - S * hit), "reported marked hit")
        for time in range(hit):
            require(not marked_phase(z - S * time), "marked hit is first")
        parameter = phase_representative(z)
        for time in range(hit):
            current = parameter + 2 * time
            record = record_from_centered_law(current)
            require(record == baseline, "pre-hit physical record must be constant")
            require(not parameter_is_visible(record, current), "pre-hit parameter leaked")
        current = parameter + 2 * hit
        record = record_from_centered_law(current)
        require(parameter_is_visible(record, current), "hit must expose numeric parameter")

    expected_hits = {0: A + 1, 16: 680}
    expected_hits.update({time: S for time in range(1, 16)})
    require(dict(sorted(hit_counts.items())) == expected_hits, "closed marked first-hit distribution")

    for horizon in range(17):
        unresolved = tuple(
            z
            for z in range(N)
            if all(not marked_phase(z - S * time) for time in range(horizon + 1))
        )
        unresolved_counts[horizon] = len(unresolved)
        if horizon <= 15:
            expected = N - (A + 1) - horizon * S
            require(len(unresolved) == expected, f"unresolved count H={horizon}")
            require(unresolved == tuple(range(A + horizon * S + 1, N)), f"unresolved interval H={horizon}")
        else:
            require(not unresolved, "horizon 16 must cover every phase")

    # The first visible numeric binder is a literal inverse decoder for the
    # whole infinite tail.  Two lifts per phase make scale blindness hostile.
    decoded_rows = 0
    for z in range(N):
        for lift in (0, 1):
            parameter = phase_representative(z, lift)
            word = output_word(parameter, 16)
            require(decode_visible_parameter(word) == parameter, "horizon-16 decoder")
            decoded_rows += 1
    require(decoded_rows == 2 * N, "decoder universe")

    # Sharp global lower bound: the same unresolved phase at two scales has
    # sixteen identical observations (times 0,...,15), then different labelled
    # P-teeth at time 16.
    hostile_p = 11123
    hostile_p_lift = hostile_p + M
    require(phase(hostile_p) == phase(hostile_p_lift) == 22927, "sharp hostile phase")
    hostile_word_1 = output_word(hostile_p, 16)
    hostile_word_2 = output_word(hostile_p_lift, 16)
    require(hostile_word_1[:16] == hostile_word_2[:16], "horizon 15 collision")
    require(all(record == baseline for record in hostile_word_1[:16]), "hostile baseline block 1")
    require(all(record == baseline for record in hostile_word_2[:16]), "hostile baseline block 2")
    require(hostile_word_1[16] != hostile_word_2[16], "hostile split at time 16")
    require(hostile_word_1[16].binders == ((11155, 308, "R", 0),), "hostile binder 1")
    require(hostile_word_2[16].binders == ((58349, 1611, "R", 0),), "hostile binder 2")
    require(hostile_word_1[16].first_uncovered == Fraction(4313, 156170), "hostile exit 1")
    require(hostile_word_2[16].first_uncovered == Fraction(22555, 816886), "hostile exit 2")

    # Boundary hostiles: open-tooth equality does not erase the varying binder.
    boundary_rows = []
    for target_rho in (-A, A):
        parameter = next(
            p
            for p in range(TAIL, TAIL + M, 2)
            if centered_rho(p) == target_rho
        )
        record = record_from_centered_law(parameter)
        boundary_rows.append((parameter, target_rho, euclidean_address(parameter), record.binders))
        require(record.first_uncovered == BASE_EXIT, "boundary metric is baseline")
        require(len(record.binders) == 2, "boundary must retain a two-tooth tie")
        require(parameter_is_visible(record, parameter), "boundary must expose numeric parameter")
    require(boundary_rows[0] == (43823, -3371, 1210, ((3371, 93, "R", 0), (43823, 1210, "L", 0))), "negative boundary record")
    require(boundary_rows[1] == (50565, 3371, 1396, ((3371, 93, "R", 0), (50565, 1396, "R", 0))), "positive boundary record")

    # Exhaust every THM-4367 structural active scale-fibre shape.  Since the
    # physical binder is the numeric speed bg, the full record is already
    # injective at horizon zero, even when all current metrics coincide.
    structural_types = 0
    structural_points = 0
    structural_pairs = 0
    maximum_scale_fibre = 0
    for a in range(2, 6741, 2):
        scale_cap = 6740 // a
        for g0 in UNITS14:
            if g0 > scale_cap:
                continue
            scales = tuple(range(g0, scale_cap + 1, 14))
            if not scales:
                continue
            structural_types += 1
            b, kappa = representative_metric_class(a, g0)
            current_records = set()
            current_metrics = set()
            for g in scales:
                parameter = b * g
                require(centered_rho(parameter) == A - a * g, "active structural cell")
                record = record_from_centered_law(parameter)
                require(record.binders == ((parameter, euclidean_address(parameter), "R", 0),), "active unique numeric binder")
                require(record.first_uncovered == Fraction(kappa, 14 * b), "active common metric")
                current_records.add(record)
                current_metrics.add(record.first_uncovered)
            require(len(current_metrics) == 1, "one current metric fibre")
            require(len(current_records) == len(scales), "full record resolves every scale")
            structural_points += len(scales)
            structural_pairs += len(scales) * (len(scales) - 1) // 2
            maximum_scale_fibre = max(maximum_scale_fibre, len(scales))
    require(structural_types == 6106, "structural type count")
    require(structural_points == 13427, "structural scale point count")
    require(structural_pairs == 281073, "structural scale pair count")
    require(maximum_scale_fibre == 241, "sharp metric scale-fibre size")

    # Loss ledger: normalize every occurrence of the varying chain/binder
    # tooth to the abstract role P.  Boundary/tie information alone does not
    # buy the shorter horizon: THM-4374's sharp pair survives through H=16.
    metric_hostile_1, metric_hostile_2 = 253031, 258645
    role_word_1 = tuple(full_record_role_projection(metric_hostile_1 + 2 * t) for t in range(18))
    role_word_2 = tuple(full_record_role_projection(metric_hostile_2 + 2 * t) for t in range(18))
    require(role_word_1[:17] == role_word_2[:17], "role-normalized horizon-16 hostile")
    require(role_word_1[17] != role_word_2[17], "role-normalized split at 17")
    require(record_from_centered_law(metric_hostile_1) != record_from_centered_law(metric_hostile_2), "full numeric binders split metric hostile immediately")

    print("THM-4379 labeled first-exit binder observer primary: PASS")
    print(f"checks={CHECKS}")
    print(f"constants=A:{A},S:{S},M:{M},N:{N},tail:{TAIL}")
    print(f"direct_geometry_rows={geometry_rows}")
    print("first_marked_hit=" + ",".join(f"{time}:{hit_counts[time]}" for time in range(17)))
    print("unresolved_phase_cells=" + ",".join(f"{horizon}:{unresolved_counts[horizon]}" for horizon in range(17)))
    print("global_full_record_horizon=16_shifts_17_observations;sharp")
    print(f"sharp_h15_same_phase_pair={hostile_p},{hostile_p_lift};phase={phase(hostile_p)}")
    print(f"sharp_h16_binders={hostile_word_1[16].binders},{hostile_word_2[16].binders}")
    print("boundaries=" + ";".join(f"P={p},rho={rho},n={n},binders={binders}" for p, rho, n, binders in boundary_rows))
    print(f"active_structural_types={structural_types};points={structural_points};pairs={structural_pairs}")
    print("active_full_record_max_fibre=" + ",".join(f"{h}:1" for h in range(17)))
    print("metric_comparison=241,94x15,49,1_at_H17;full_physical=1_from_H0_on_active_fibres")
    print(f"role_normalized_metric_hostile={metric_hostile_1},{metric_hostile_2};equal_through_H16")
    print("scope=THM-4365 fixed h=420 odd P>=11019 tail;all rows safe;LRC(14) OPEN")


if __name__ == "__main__":
    main()
