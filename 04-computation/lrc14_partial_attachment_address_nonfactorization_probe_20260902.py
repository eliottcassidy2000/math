#!/usr/bin/env python3
"""Exact attachment-address probe on the h=420 anchored LRC(14) bank.

This is exploratory scratch work, not a canon theorem and not an LRC(14)
closure claim.  It reuses the exact open-tooth/farthest-reach semantics of
THM-4335/4348, but compares a deliberately coarse owner/support quotient
with the full endpoint-to-component attachment record.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relative: str):
    spec = spec_from_file_location(name, ROOT / relative)
    assert spec is not None and spec.loader is not None
    module = module_from_spec(spec)
    # dataclasses consult sys.modules while resolving postponed annotations.
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


renewal = load_module(
    "renewal4348",
    "04-computation/lrc14_third_tooth_competition_probe_third_tooth_20260902.py",
)
current = load_module(
    "current4345",
    "04-computation/lrc14_anchor_current_strip_probe_anchor_current_20260902.py",
)


H = 420
COMMON = tuple(11 + 1680 * k for k in range(7)) + (525, 945, 1365, 1575)
BASE_TEETH = tuple(
    tuple(
        t
        for w in COMMON
        for t in renewal.meeting_teeth(w, *renewal.anchor_component(H, k))
    )
    for k in range(2 * H)
)


def role_map(extra: int) -> dict[int, str]:
    labels = [f"D{i}" for i in range(7)] + [f"C{i}" for i in range(4)]
    ans = dict(zip(COMMON, labels, strict=True))
    ans[extra] = "P"
    return ans


@dataclass(frozen=True)
class RowRecord:
    extra: int
    # Old quotient: owner sheet and unaddressed chain support, as a multiset.
    coarse: tuple[tuple[tuple[int, str, tuple[str, ...]], int], ...]
    # Full endpoint attachment address, including physical component and tooth
    # address at every selected occurrence.
    attachments: tuple[tuple[object, ...], ...]
    statuses: tuple[tuple[str, int], ...]
    first_exit: tuple[object, ...] | None


def active_trace_cached(extra: int, k: int):
    """Same exact greedy trace as THM-4348, with common teeth cached."""
    left, right = renewal.anchor_component(H, k)
    teeth = BASE_TEETH[k] + renewal.meeting_teeth(extra, left, right)
    cursor = left
    chain = []
    frontiers = []
    while True:
        eligible = [t for t in teeth if t.left < cursor < t.right]
        if not eligible:
            return "missing", tuple(chain), tuple(frontiers)
        winner = max(eligible, key=renewal.selection_key)
        frontiers.append(cursor)
        chain.append(winner)
        cursor = winner.right
        if cursor > right:
            return ("span" if len(chain) == 1 else "renew"), tuple(chain), tuple(frontiers)


def row_record(extra: int) -> RowRecord:
    assert extra > 0 and extra % 2 == 1 and extra not in COMMON
    labels = role_map(extra)
    speeds = COMMON + (extra,)
    coarse_counter: Counter[tuple[int, str, tuple[str, ...]]] = Counter()
    attachments: list[tuple[object, ...]] = []
    statuses: Counter[str] = Counter()
    first_exit: tuple[object, ...] | None = None
    occurrence_counter: Counter[tuple[int, str, str]] = Counter()

    for k in range(2 * H):
        status, chain, frontiers = active_trace_cached(extra, k)
        statuses[status] += 1
        epsilon = int(k >= H)
        roles = tuple(labels[t.w] for t in chain)
        # This is the old completed-cover quotient: it knows whether a
        # component is missing/span/renew, the owner sheet, and the support
        # word of every completed cover.  On a missing component it retains
        # only the failure occurrence, not the partial endpoint itinerary.
        # That is precisely the candidate quotient under audit.
        coarse_counter[(epsilon, status, roles if status != "missing" else ())] += 1

        if status == "missing" and first_exit is None:
            cursor = (
                frontiers[-1]
                if frontiers
                else renewal.anchor_component(H, k)[0]
            )
            # active_trace stores the frontier before each selected tooth.  On
            # failure after a partial chain the missing point is the last
            # selected right endpoint.
            if chain:
                cursor = chain[-1].right
            full = (2 * H,) + speeds
            distances = tuple(current.circle_distance(v * cursor) for v in full)
            clearance = min(distances)
            binding = tuple(v for v, d in zip(full, distances) if d == clearance)
            assert clearance >= Fraction(1, 14)
            first_exit = (k, k % H, epsilon, cursor, clearance, binding)

        # A selected chain endpoint is an occurrence, even when the component
        # later fails.  Record each outgoing right side and its next consumer.
        for i, tooth in enumerate(chain):
            source_role = labels[tooth.w]
            key = (epsilon, source_role, "R")
            ordinal = occurrence_counter[key]
            occurrence_counter[key] += 1
            next_tooth = chain[i + 1] if i + 1 < len(chain) else None
            attachments.append(
                (
                    H,
                    epsilon,
                    source_role,
                    "R",
                    ordinal,
                    k,
                    k % H,
                    tooth.n,
                    None if next_tooth is None else labels[next_tooth.w],
                    None if next_tooth is None else next_tooth.n,
                    status,
                )
            )

    return RowRecord(
        extra=extra,
        coarse=tuple(sorted(coarse_counter.items(), key=repr)),
        attachments=tuple(attachments),
        statuses=tuple(sorted(statuses.items())),
        first_exit=first_exit,
    )


def live_filters(extra: int) -> bool:
    """Cheap inherited filters for the anchored 420-wall search universe."""
    if extra <= 0 or extra % 2 == 0 or extra in COMMON:
        return False
    speeds = COMMON + (extra,)
    # The fixed half-turn endpoints from THM-4330 must both fail; their
    # distance vectors agree, so test t=1/2+1/(28h).
    t = Fraction(1, 2) + Fraction(1, 28 * H)
    if min(current.circle_distance(v * t) for v in (2 * H,) + speeds) >= Fraction(1, 14):
        return False
    # THM-366 completeness is already supplied by COMMON, but assert it.
    if any(not any(v % d == 0 for v in (2 * H,) + speeds) for d in range(2, 15)):
        return False
    return True


def search(extras: Iterable[int]) -> tuple[RowRecord, RowRecord] | None:
    buckets: defaultdict[object, list[RowRecord]] = defaultdict(list)
    for extra in extras:
        if not live_filters(extra):
            continue
        record = row_record(extra)
        for earlier in buckets[record.coarse]:
            if earlier.attachments != record.attachments:
                return earlier, record
        buckets[record.coarse].append(record)
    return None


def reached_next_blocker_change(a: int, b: int) -> tuple[object, ...] | None:
    """First common reached right endpoint whose next winner changes."""
    labels_a = role_map(a)
    labels_b = role_map(b)
    speeds_a = COMMON + (a,)
    speeds_b = COMMON + (b,)
    first_address_only = None
    for k in range(2 * H):
        status_a, chain_a, _ = active_trace_cached(a, k)
        status_b, chain_b, _ = active_trace_cached(b, k)
        for i in range(min(len(chain_a), len(chain_b)) - 1):
            source_a, source_b = chain_a[i], chain_b[i]
            if source_a != source_b:
                break
            next_a, next_b = chain_a[i + 1], chain_b[i + 1]
            if next_a != next_b:
                assert renewal.residue_winner(source_a.w, source_a.n, speeds_a) == (
                    next_a.w,
                    next_a.n,
                )
                assert renewal.residue_winner(source_b.w, source_b.n, speeds_b) == (
                    next_b.w,
                    next_b.n,
                )
                epsilon = int(k >= H)
                record = (
                    "role" if labels_a[next_a.w] != labels_b[next_b.w] else "address",
                    H,
                    epsilon,
                    labels_a[source_a.w],
                    "R",
                    i,
                    k,
                    k % H,
                    source_a.n,
                    labels_a[next_a.w],
                    next_a.n,
                    labels_b[next_b.w],
                    next_b.n,
                    status_a,
                    status_b,
                )
                if record[0] == "role":
                    return record
                if first_address_only is None:
                    first_address_only = record
    return first_address_only


def changed_attachment_rows(a: RowRecord, b: RowRecord) -> tuple[object, ...]:
    aset = set(a.attachments)
    bset = set(b.attachments)
    return tuple(sorted(aset ^ bset, key=repr))


def trace_map(extra: int) -> tuple[tuple[object, ...], ...]:
    labels = role_map(extra)
    answer = []
    for k in range(2 * H):
        status, chain, _ = active_trace_cached(extra, k)
        answer.append((status, tuple((labels[t.w], t.n) for t in chain)))
    return tuple(answer)


def local_signed_current(t: Fraction, speeds: tuple[int, ...]) -> tuple[int, int, int]:
    lower = sum(current.circle_distance(w * t) < Fraction(1, 14) for w in speeds)
    upper = sum(
        current.circle_distance(w * (t + Fraction(1, 2))) < Fraction(1, 14)
        for w in speeds
    )
    return lower, upper, lower - upper


def one_sided_current_at(
    x: Fraction, speeds_a: tuple[int, ...], speeds_b: tuple[int, ...]
) -> tuple[object, ...]:
    walls = set(current.tail_events(speeds_a)) | set(current.tail_events(speeds_b))
    previous = max(w for w in walls if w < x)
    following = min(w for w in walls if w > x)
    left_test = (previous + x) / 2
    right_test = (x + following) / 2
    return (
        x,
        previous,
        following,
        local_signed_current(left_test, speeds_a),
        local_signed_current(left_test, speeds_b),
        local_signed_current(right_test, speeds_a),
        local_signed_current(right_test, speeds_b),
    )


def component_exit(extra: int, k: int) -> tuple[object, ...]:
    status, chain, _ = active_trace_cached(extra, k)
    assert status == "missing"
    cursor = chain[-1].right if chain else renewal.anchor_component(H, k)[0]
    full = (2 * H,) + COMMON + (extra,)
    distances = tuple(current.circle_distance(v * cursor) for v in full)
    clearance = min(distances)
    binding = tuple(v for v, d in zip(full, distances) if d == clearance)
    assert clearance >= Fraction(1, 14)
    return cursor, clearance, binding


def audit_resonant_side_blowup() -> tuple[int, tuple[object, ...]]:
    """Audit the exact 7u|h nearest-boundary obstruction.

    At a signed u-wall the anchor is zero.  Its two nearest anchor-safe
    boundary points split: the side outside the u-tooth is killed by 13u,
    while the other side is killed by u.  Hence merely resolving a resonant
    wall into its two adjacent component attachments cannot produce a witness.
    """
    checks = 0
    witness = None
    for L in range(1, 8):
        for u in range(1, 16, 2):
            h = 7 * L * u
            delta = Fraction(1, 28 * h)
            for n in range(u):
                for sigma in (-1, 1):
                    wall = Fraction(14 * n + sigma, 14 * u)
                    assert current.circle_distance(2 * h * wall) == 0
                    assert current.circle_distance(u * wall) == Fraction(1, 14)
                    assert current.circle_distance(13 * u * wall) == Fraction(1, 14)
                    inward = wall - sigma * delta
                    outward = wall + sigma * delta
                    assert current.circle_distance(2 * h * inward) == Fraction(1, 14)
                    assert current.circle_distance(2 * h * outward) == Fraction(1, 14)
                    assert current.circle_distance(u * inward) < Fraction(1, 14)
                    assert current.circle_distance(u * outward) > Fraction(1, 14)
                    assert current.circle_distance(13 * u * inward) > Fraction(1, 14)
                    assert current.circle_distance(13 * u * outward) < Fraction(1, 14)
                    checks += 1
                    if (L, u, n, sigma) == (20, 3, 0, 1):
                        witness = (h, u, wall, inward, outward)
    # The universal audit above intentionally uses L<=7.  Freeze the live
    # h=420,u=3 instance separately (L=20).
    L, u, n, sigma = 20, 3, 0, 1
    h = 7 * L * u
    delta = Fraction(1, 28 * h)
    wall = Fraction(14 * n + sigma, 14 * u)
    inward, outward = wall - sigma * delta, wall + sigma * delta
    assert current.circle_distance(2 * h * inward) == Fraction(1, 14)
    assert current.circle_distance(u * inward) < Fraction(1, 14)
    assert current.circle_distance(13 * u * outward) < Fraction(1, 14)
    witness = (
        h,
        u,
        wall,
        inward,
        outward,
        current.circle_distance(u * inward),
        current.circle_distance(13 * u * outward),
    )
    return checks, witness


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")
    print("LRC14_PARTIAL_ATTACHMENT_ADDRESS_NONFACTORIZATION_PROBE=FINITE_EXACT")
    print(f"gate=h={H};420|h;row={{2h}} union COMMON union {{P}}")
    print("common=" + ",".join(map(str, COMMON)))

    # The two canonical THM-4330/4335 controls already form the minimum
    # cardinality hostile bank: no comparison is possible with one row.
    controls = [row_record(1287), row_record(9009)]
    print("controls")
    for record in controls:
        profile = current.exact_profile(H, COMMON + (record.extra,))
        print(
            f" P={record.extra};statuses={record.statuses};"
            f"coarse_equal_to_first={record.coarse == controls[0].coarse};"
            f"attachment_equal_to_first={record.attachments == controls[0].attachments};"
            f"first_exit={record.first_exit};core_energy={profile['core_energy']};"
            f"core_q_cubic={profile['core_q_cubic']}"
        )

    a, b = controls
    assert a.coarse == b.coarse
    assert a.attachments != b.attachments
    pa = current.exact_profile(H, COMMON + (a.extra,))
    pb = current.exact_profile(H, COMMON + (b.extra,))
    delta = changed_attachment_rows(a, b)
    next_change = reached_next_blocker_change(a.extra, b.extra)
    assert next_change is not None
    traces_a = trace_map(a.extra)
    traces_b = trace_map(b.extra)
    status_map_equal = tuple(row[0] for row in traces_a) == tuple(
        row[0] for row in traces_b
    )
    changed_components = tuple(
        k for k, (left, right) in enumerate(zip(traces_a, traces_b)) if left != right
    )
    role_changed_components = tuple(
        k
        for k, (left, right) in enumerate(zip(traces_a, traces_b))
        if tuple(role for role, _ in left[1]) != tuple(role for role, _ in right[1])
    )
    changed_nonmissing = tuple(
        k for k in changed_components if traces_a[k][0] != "missing"
    )
    assert not changed_nonmissing
    hostile_k = int(next_change[6])
    hostile_source_n = int(next_change[8])
    hostile_source_speed = 1575
    hostile_x = renewal.tooth(hostile_source_speed, hostile_source_n).right
    local_current = one_sided_current_at(
        hostile_x, COMMON + (a.extra,), COMMON + (b.extra,)
    )
    print("minimum_two_row_safe_control_hostile")
    print(f" P_pair=({a.extra},{b.extra})")
    print(f" statuses={a.statuses}")
    print(f" coarse={a.coarse}")
    print(f" coarse_sha256={sha256(repr(a.coarse).encode()).hexdigest()}")
    print(f" coarse_equal={a.coarse == b.coarse}")
    print(f" attachment_equal={a.attachments == b.attachments}")
    print(f" attachment_symmetric_difference={len(delta)}")
    print(f" status_map_equal={status_map_equal}")
    print(f" changed_component_traces={len(changed_components)}")
    print(f" role_changed_component_traces={len(role_changed_components)}")
    print(f" changed_nonmissing_component_traces={len(changed_nonmissing)}")
    print(f" first_changed_components={changed_components[:12]}")
    print(f" first_role_changed_components={role_changed_components[:12]}")
    print(f" first_attachment_difference={delta[0] if delta else None}")
    print(f" reached_next_blocker_change={next_change}")
    print(f" local_trace_A_k{hostile_k}={traces_a[hostile_k]}")
    print(f" local_trace_B_k{hostile_k}={traces_b[hostile_k]}")
    print(f" local_exit_A_k{hostile_k}={component_exit(a.extra, hostile_k)}")
    print(f" local_exit_B_k{hostile_k}={component_exit(b.extra, hostile_k)}")
    print(f" one_sided_signed_current={local_current}")
    print(f" first_exit_A={a.first_exit}")
    print(f" first_exit_B={b.first_exit}")
    print(f" first_exit_changed={a.first_exit != b.first_exit}")
    print(f" core_energy_A={pa['core_energy']}")
    print(f" core_energy_B={pb['core_energy']}")
    print(f" core_energy_delta={pa['core_energy'] - pb['core_energy']}")
    print(f" signed_current_changed={pa['core_energy'] != pb['core_energy']}")
    print(f" core_q_cubic_A={pa['core_q_cubic']}")
    print(f" core_q_cubic_B={pb['core_q_cubic']}")
    print(f" core_q_cubic_delta={pa['core_q_cubic'] - pb['core_q_cubic']}")
    resonance_checks, resonance_witness = audit_resonant_side_blowup()
    print("seven_u_divides_h_resonant_side_blowup")
    print(f" symbolic_grid_checks={resonance_checks}")
    print(f" h420_u3_witness={resonance_witness}")
    print(" split=inner_attachment_killed_by_u;outer_attachment_killed_by_13u")
    print("VERDICT=PASS")


if __name__ == "__main__":
    main()
