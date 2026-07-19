#!/usr/bin/env python3
"""Exact A_12 chip-walk carrier for the shallow n=12 LRC frontier.

Domain
------
The input is a residue-labelled packet

    W=(w_r)_(r in F_13^*),       w_r > 0,       w_r == r (mod 13).

This is exactly the shallow/full-residue branch of the denominator-13
tightness problem.  It is *not* the deep branch and it does not by itself
establish the covering hypotheses needed by the LRC14 inverse statement.

Put t=(j+u)/13, with j in F_13 and 0 <= u < 1.  Away from an event
u=k/w_r, runner r is 1/13-dangerous on exactly the Cayley edge

    E_r(u)={-floor(w_r u) r^-1, -(floor(w_r u)+1) r^-1}.

Initially the degree excess over one is 11 delta_0.  At u=k/w_r the edge
slides by the A_12 root

    delta_(-(k+1)r^-1) - delta_(-(k-1)r^-1).

All events at the same rational time must be applied simultaneously.  Hence
the original continuous covering predicate is equivalent, on the whole
unbounded domain above, to the ballot condition

    11 delta_0 + C(prefix) >= 0

for every grouped prefix of the mechanical root word.  At an event itself
the closed danger set contains the union of the adjacent edges, so checking
the generic chambers is complete.  The state space is the 11-chip simplex
of size C(23,11)=1,352,078, independent of the heights.

The still-open proof target isolated by this carrier is:

  A12 BALLOT RIGIDITY.  Every accepted mechanical word is a dilation ray
  W=c{1,...,12}, with 13 not dividing c.

Equivalently, every accepted packet which is a sink under legal coordinate
lowerings w_r -> w_r-13 is a dilation ray.  This script does not assert that
lemma.  It tests the exact carrier independently, exhausts heights 0..2,
and searches for counterexamples to weaker descent potentials.

Dilation rays have an exact Farey-comb/toothpick self-similarity.  If
r=ca (mod 13), u=(m+v)/c, and phi_(c,m)(x)=c^-1(x-m), then

    E_(ca)((m+v)/c) = phi_(c,m)(E_a(v)).

Thus c{1,...,12} consists of c affine copies of the Farey-12 base word.  A
proof of the converse can be cleanly targeted as a regeneration theorem:
an accepted primitive mechanical word has just the one base comb.

Tournament / carrier challenge
------------------------------
The predicate-preserving vertices are the 13 sheet cuts, not necessarily the
12 runners.  Orient a pair of cuts by lexicographic current deficit; root
negation plus u -> 1-u is the switch/gauge and 0,1,...,12 is the tie
Hamiltonian path.  This telemetry tournament is transitive and therefore
forgets the essential chronology.  The alternative event tournament (one
vertex per labelled wall crossing, ordered by k/w_r, ties contracted) retains
the Christoffel interleavings but is unbounded.  The grouped A_12 word is the
finite-state quotient that preserves coverage.  It destroys runner identity
unless root letters retain (r,k), and it destroys the arithmetic-realizability
language if event order is forgotten.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from heapq import heapify, heappop, heappush
from itertools import product
from math import comb, gcd
from pathlib import Path
from random import Random
from typing import Iterable


P = 13
RESIDUES = tuple(range(1, P))
INITIAL_EXCESS = (11,) + (0,) * 12
CHIP_SIMPLEX_SIZE = comb(23, 11)


def fmt(value: F | None) -> str:
    if value is None:
        return "survive"
    return f"{value.numerator}/{value.denominator}"


def validate_packet(packet: tuple[int, ...]) -> None:
    assert len(packet) == 12
    assert all(speed > 0 for speed in packet)
    assert all(speed % P == residue for residue, speed in zip(RESIDUES, packet))


def heights(packet: tuple[int, ...]) -> tuple[int, ...]:
    validate_packet(packet)
    return tuple((speed - residue) // P for residue, speed in zip(RESIDUES, packet))


def packet_from_heights(height_vector: tuple[int, ...]) -> tuple[int, ...]:
    assert len(height_vector) == 12
    assert all(height >= 0 for height in height_vector)
    return tuple(residue + P * height for residue, height in zip(RESIDUES, height_vector))


def dilation_packet(scale: int) -> tuple[int, ...]:
    assert scale > 0 and scale % P
    packet = [0] * 12
    for multiplier in RESIDUES:
        speed = scale * multiplier
        packet[speed % P - 1] = speed
    result = tuple(packet)
    validate_packet(result)
    return result


def dilation_scale(packet: tuple[int, ...]) -> int | None:
    common = 0
    for speed in packet:
        common = gcd(common, speed)
    if sorted(packet) == [common * multiplier for multiplier in RESIDUES]:
        return common
    return None


def primitive_normalization(packet: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
    """Divide by gcd and re-label the resulting full residue packet."""
    validate_packet(packet)
    common = 0
    for speed in packet:
        common = gcd(common, speed)
    assert common > 0 and common % P
    primitive = [0] * 12
    for speed in packet:
        reduced = speed // common
        assert reduced % P
        assert primitive[reduced % P - 1] == 0
        primitive[reduced % P - 1] = reduced
    result = tuple(primitive)
    validate_packet(result)
    assert gcd(*result) == 1
    return common, result


def affine_carry_heights(scale: int) -> tuple[int, ...]:
    """The twelve explicit height rays for scale=13q+s."""
    assert scale > 0 and scale % P
    quotient, residue_scale = divmod(scale, P)
    inverse_scale = pow(residue_scale, -1, P)
    return tuple(
        quotient * ((inverse_scale * residue) % P)
        + (residue_scale * ((inverse_scale * residue) % P)) // P
        for residue in RESIDUES
    )


def distance_mod_one(value: F) -> F:
    phase = value % 1
    return min(phase, 1 - phase)


def is_lonely(packet: tuple[int, ...], time: F) -> bool:
    return all(distance_mod_one(speed * time) > F(1, P) for speed in packet)


@dataclass(frozen=True)
class WalkAudit:
    accepted: bool
    first_tear: F | None
    first_tear_right: F | None
    deficient_sheets: tuple[int, ...]
    event_groups: int
    final_excess: tuple[int, ...] | None
    final_current: tuple[int, ...] | None
    prefix_minimum: tuple[int, ...] | None

    def lonely_witness(self, packet: tuple[int, ...]) -> tuple[int, F] | None:
        if self.first_tear is None:
            return None
        assert self.first_tear_right is not None
        sheet = self.deficient_sheets[0]
        u = (self.first_tear + self.first_tear_right) / 2
        time = F(sheet, P) + u / P
        assert is_lonely(packet, time)
        return sheet, time


def chipwalk(packet: tuple[int, ...], *, stop_at_tear: bool = False) -> WalkAudit:
    """Run the exact grouped mechanical A_12 root word."""
    validate_packet(packet)
    heap: list[tuple[F, int, int, int, int]] = []
    for index, (residue, speed) in enumerate(zip(RESIDUES, packet)):
        if speed > 1:
            heap.append((F(1, speed), index, 1, speed, pow(residue, -1, P)))
    heapify(heap)

    excess = list(INITIAL_EXCESS)
    current = [0] * P
    prefix_minimum = [0] * P
    first_tear: F | None = None
    first_tear_right: F | None = None
    deficient: tuple[int, ...] = ()
    groups = 0

    while heap:
        event_time = heap[0][0]
        updates: list[tuple[int, int]] = []
        while heap and heap[0][0] == event_time:
            _, index, integer, speed, inverse = heappop(heap)
            departure = (-(integer - 1) * inverse) % P
            arrival = (-(integer + 1) * inverse) % P
            updates.append((departure, arrival))
            if integer + 1 < speed:
                heappush(
                    heap,
                    (F(integer + 1, speed), index, integer + 1, speed, inverse),
                )

        increment = [0] * P
        for departure, arrival in updates:
            increment[departure] -= 1
            increment[arrival] += 1
        for sheet in range(P):
            excess[sheet] += increment[sheet]
            current[sheet] += increment[sheet]
            prefix_minimum[sheet] = min(prefix_minimum[sheet], current[sheet])
            assert excess[sheet] == INITIAL_EXCESS[sheet] + current[sheet]
        assert sum(excess) == 11 and sum(current) == 0
        groups += 1

        if first_tear is None and min(excess) < 0:
            first_tear = event_time
            first_tear_right = heap[0][0] if heap else F(1)
            deficient = tuple(sheet for sheet, value in enumerate(excess) if value < 0)
            if stop_at_tear:
                return WalkAudit(
                    False,
                    first_tear,
                    first_tear_right,
                    deficient,
                    groups,
                    None,
                    None,
                    None,
                )

    expected_final = (0,) * 12 + (11,)
    assert tuple(excess) == expected_final
    assert tuple(current) == tuple(
        expected_final[sheet] - INITIAL_EXCESS[sheet] for sheet in range(P)
    )
    accepted = first_tear is None
    if accepted:
        assert tuple(prefix_minimum) == (-11,) + (0,) * 12
    return WalkAudit(
        accepted,
        first_tear,
        first_tear_right,
        deficient,
        groups,
        tuple(excess),
        tuple(current),
        tuple(prefix_minimum),
    )


def original_circle_cover(packet: tuple[int, ...]) -> tuple[bool, F | None]:
    """Independent exact check in t-space at every danger boundary/chamber."""
    validate_packet(packet)
    boundaries = {F(0), F(1)}
    for speed in packet:
        for integer in range(speed):
            for sign in (-1, 1):
                point = F(P * integer + sign, P * speed) % 1
                boundaries.add(point)
    ordered = sorted(boundaries)
    probes = set(ordered)
    probes.update((left + right) / 2 for left, right in zip(ordered, ordered[1:]))
    for time in sorted(probes):
        if is_lonely(packet, time):
            return False, time
    return True, None


def edge_at(packet: tuple[int, ...], residue: int, u: F) -> tuple[int, int]:
    speed = packet[residue - 1]
    integer = (speed * u.numerator) // u.denominator
    inverse = pow(residue, -1, P)
    return ((-integer * inverse) % P, (-(integer + 1) * inverse) % P)


def graph_topology(packet: tuple[int, ...], u: F) -> tuple[int, int, int]:
    """Return (components, cycle rank, maximum degree) in a generic chamber."""
    parent = list(range(P))
    degrees = [0] * P

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    cycle_rank = 0
    for residue in RESIDUES:
        left, right = edge_at(packet, residue, u)
        degrees[left] += 1
        degrees[right] += 1
        root_left, root_right = find(left), find(right)
        if root_left == root_right:
            cycle_rank += 1
        else:
            parent[root_left] = root_right
    components = len({find(vertex) for vertex in range(P)})
    return components, cycle_rank, max(degrees)


def chamber_midpoints(packet: tuple[int, ...]) -> tuple[F, ...]:
    events = sorted(
        {F(integer, speed) for speed in packet for integer in range(1, speed)}
    )
    boundaries = (F(0), *events, F(1))
    return tuple((left + right) / 2 for left, right in zip(boundaries, boundaries[1:]))


def audit_toothpick_self_similarity(scales: tuple[int, ...]) -> dict[int, Counter]:
    base = dilation_packet(1)
    base_chambers = chamber_midpoints(base)
    assert len(base_chambers) == 46
    records: dict[int, Counter] = {}
    for scale in scales:
        packet = dilation_packet(scale)
        inverse_scale = pow(scale % P, -1, P)
        for block in range(scale):
            for v in base_chambers:
                u = F(block, scale) + v / scale
                for multiplier in RESIDUES:
                    residue = scale * multiplier % P
                    actual = set(edge_at(packet, residue, u))
                    base_edge = edge_at(base, multiplier, v)
                    transported = {
                        (inverse_scale * (sheet - block)) % P for sheet in base_edge
                    }
                    assert actual == transported

        profile: Counter[tuple[int, int]] = Counter()
        stars = 0
        for u in chamber_midpoints(packet):
            components, cycle_rank, max_degree = graph_topology(packet, u)
            profile[(components, cycle_rank)] += 1
            stars += max_degree == 12
        assert profile == Counter({(1, 0): 24 * scale, (2, 1): 22 * scale})
        assert stars == 2 * scale
        records[scale] = profile
    return records


def direct_option_masks(max_height: int = 2) -> tuple[int, dict[tuple[int, int], int], int]:
    """Direct modular masks, deliberately independent of the root updates."""
    options = tuple(
        residue + P * height
        for residue in RESIDUES
        for height in range(max_height + 1)
    )
    events = sorted(
        {F(integer, speed) for speed in options for integer in range(1, speed)}
    )
    boundaries = (F(0), *events, F(1))
    midpoints = tuple((left + right) / 2 for left, right in zip(boundaries, boundaries[1:]))
    masks: dict[tuple[int, int], int] = {}
    for residue in RESIDUES:
        for height in range(max_height + 1):
            speed = residue + P * height
            mask = 0
            for chamber, u in enumerate(midpoints):
                dangerous_sheets = []
                for sheet in range(P):
                    time = F(sheet, P) + u / P
                    if distance_mod_one(speed * time) < F(1, P):
                        dangerous_sheets.append(sheet)
                assert len(dangerous_sheets) == 2
                for sheet in dangerous_sheets:
                    mask |= 1 << (P * chamber + sheet)
            masks[(residue, height)] = mask
    total_bits = P * len(midpoints)
    return len(midpoints), masks, (1 << total_bits) - 1


def exhaustive_height_box(max_height: int = 2) -> tuple[list[tuple[int, ...]], str, int]:
    chambers, masks, all_mask = direct_option_masks(max_height)
    half = 6
    left_records = []
    right_records = []
    for choice in product(range(max_height + 1), repeat=half):
        mask = 0
        for residue, height in zip(range(1, half + 1), choice):
            mask |= masks[(residue, height)]
        left_records.append((choice, mask))
    for choice in product(range(max_height + 1), repeat=12 - half):
        mask = 0
        for residue, height in zip(range(half + 1, P), choice):
            mask |= masks[(residue, height)]
        right_records.append((choice, mask))

    accepted = []
    decision_digest = sha256()
    for left_choice, left_mask in left_records:
        need = all_mask ^ left_mask
        for right_choice, right_mask in right_records:
            choice = left_choice + right_choice
            survives = (need & right_mask) == need
            decision_digest.update(bytes(choice) + bytes((survives,)))
            if survives:
                accepted.append(choice)
    assert len(left_records) * len(right_records) == (max_height + 1) ** 12
    return accepted, decision_digest.hexdigest(), chambers


def lowering(packet: tuple[int, ...], coordinate: int) -> tuple[int, ...]:
    assert packet[coordinate] > P
    result = list(packet)
    result[coordinate] -= P
    lowered = tuple(result)
    validate_packet(lowered)
    return lowered


def exhaustive_h01_descent_counterexamples() -> tuple[int, list[tuple]]:
    cache: dict[tuple[int, ...], WalkAudit] = {}
    for choice in product((0, 1), repeat=12):
        cache[choice] = chipwalk(packet_from_heights(choice), stop_at_tear=True)

    regressions = []
    for choice, audit in cache.items():
        if audit.accepted:
            continue
        legal = [coordinate for coordinate, height in enumerate(choice) if height]
        if not legal:
            continue
        lowered = []
        for coordinate in legal:
            next_choice = choice[:coordinate] + (0,) + choice[coordinate + 1 :]
            next_audit = cache[next_choice]
            lowered.append((coordinate + 1, next_audit.first_tear))
        if all(
            next_tear is not None and next_tear < audit.first_tear
            for _, next_tear in lowered
        ):
            regressions.append((choice, audit.first_tear, tuple(lowered)))

    accepted_count = sum(audit.accepted for audit in cache.values())
    assert accepted_count == 2
    assert len(regressions) == 2
    return accepted_count, regressions


def audit_perturbations(scales: Iterable[int]) -> tuple[int, int]:
    upward_failures = 0
    downward_failures = 0
    for scale in scales:
        packet = dilation_packet(scale)
        assert chipwalk(packet).accepted
        for coordinate in range(12):
            lifted = list(packet)
            lifted[coordinate] += P
            lifted_packet = tuple(lifted)
            audit = chipwalk(lifted_packet, stop_at_tear=True)
            assert not audit.accepted
            assert audit.lonely_witness(lifted_packet) is not None
            upward_failures += 1
            if packet[coordinate] > P:
                lowered_packet = lowering(packet, coordinate)
                audit = chipwalk(lowered_packet, stop_at_tear=True)
                assert not audit.accepted
                assert audit.lonely_witness(lowered_packet) is not None
                downward_failures += 1
    return upward_failures, downward_failures


def audit_random_high_packets(count: int, seed: int) -> tuple[int, Counter[str]]:
    rng = Random(seed)
    accepted = 0
    tear_denominators: Counter[str] = Counter()
    for _ in range(count):
        choice = tuple(rng.randrange(0, 81) for _ in RESIDUES)
        packet = packet_from_heights(choice)
        audit = chipwalk(packet, stop_at_tear=True)
        accepted += audit.accepted
        if audit.accepted:
            assert dilation_scale(packet) is not None
        else:
            assert audit.lonely_witness(packet) is not None
            assert audit.first_tear is not None
            tear_denominators[str(audit.first_tear.denominator)] += 1
    return accepted, tear_denominators


def main() -> None:
    assert CHIP_SIMPLEX_SIZE == 1_352_078

    controls = {
        "AP": dilation_packet(1),
        "2AP": dilation_packet(2),
        "3AP": dilation_packet(3),
        "single_lift_r1": packet_from_heights((1,) + (0,) * 11),
        "all_plus_13": packet_from_heights((1,) * 12),
        "mixed_tail": (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 24, 25),
    }
    control_records = {}
    for name, packet in controls.items():
        validate_packet(packet)
        walk = chipwalk(packet)
        direct, direct_witness = original_circle_cover(packet)
        assert walk.accepted == direct
        if not walk.accepted:
            assert walk.lonely_witness(packet) is not None
            assert direct_witness is not None
        control_records[name] = (packet, walk, direct_witness)

    toothpick_scales = (1, 2, 3, 5, 14)
    toothpick = audit_toothpick_self_similarity(toothpick_scales)
    for scale in toothpick_scales + (17, 25, 27, 97):
        packet = dilation_packet(scale)
        assert heights(packet) == affine_carry_heights(scale)
        common, primitive = primitive_normalization(packet)
        assert common == scale and primitive == dilation_packet(1)
        assert chipwalk(packet).accepted == chipwalk(primitive).accepted

    accepted_h2, h2_digest, chambers_h2 = exhaustive_height_box(2)
    expected_h2 = [
        heights(dilation_packet(1)),
        heights(dilation_packet(2)),
        heights(dilation_packet(3)),
    ]
    assert accepted_h2 == expected_h2
    for choice in accepted_h2:
        packet = packet_from_heights(choice)
        assert chipwalk(packet).accepted
        assert original_circle_cover(packet)[0]
        assert dilation_scale(packet) in (1, 2, 3)

    all_high = set(product((1, 2), repeat=12))
    assert not all_high.intersection(accepted_h2)
    endpoint_only_false_positives = len(all_high)

    accepted_h01, descent_regressions = exhaustive_h01_descent_counterexamples()

    perturbation_scales = (1, 2, 3, 4, 5, 7, 12, 14, 17, 25, 27, 97)
    upward_failures, downward_failures = audit_perturbations(perturbation_scales)
    random_count = 1_000
    random_seed = 20260718
    random_accepted, random_tear_denominators = audit_random_high_packets(
        random_count, random_seed
    )
    assert random_accepted == 0

    source_digest = sha256(Path(__file__).read_bytes()).hexdigest()
    print("A12 SHALLOW CHIP-WALK / DESCENT AUDIT")
    print(f"source_sha256={source_digest}")
    print("domain=all positive residue-labelled w_r == r mod 13")
    print("scope=shallow full-residue chi_13=0 only; not deep and not INVcov")
    print("certificate_claim=exhaustive direct search only for h_r in {0,1,2}; arbitrary winding rigidity remains open")
    print(f"finite_state=11-chip simplex size C(23,11)={CHIP_SIMPLEX_SIZE}")
    print("exact_equivalence=chi_13(W)=0 iff 11*delta_0+C(prefix)>=0 for every grouped prefix")
    print("events_at_equal_fractions=grouped before testing; closed event points are covered by adjacent-edge union")

    print("\nCONTROL CROSS-CHECK: ROOT WALK VS ORIGINAL t-CIRCLE")
    for name, (packet, audit, direct_witness) in control_records.items():
        witness = audit.lonely_witness(packet)
        print(
            f"{name}: W={packet} accepted={audit.accepted} "
            f"first_tear={fmt(audit.first_tear)} deficient={audit.deficient_sheets} "
            f"root_witness={witness} direct_witness={direct_witness}"
        )

    print("\nDILATION / FAREY-COMB TOOTHPICK SELF-SIMILARITY")
    print("identity: E_ca((m+v)/c)=phi_(c,m)(E_a(v)), phi(x)=c^-1(x-m) mod 13")
    print("carry_ray: for c=13q+s and pi_s(r)=<s^-1 r>_13, h_r=q*pi_s(r)+floor(s*pi_s(r)/13)")
    print("gcd_reduction: divide W by gcd(W), re-label residues; circle covering is unchanged by the surjective time map")
    print("base_word=46 Farey-12 chambers: 24 trees + 22 (two-component, one-cycle) chambers")
    for scale in toothpick_scales:
        print(
            f"c={scale}: chambers={46*scale} profile={dict(toothpick[scale])} "
            f"star_chambers={2*scale}"
        )

    print("\nEXHAUSTIVE DIRECT-MASK BOX h_r IN {0,1,2}")
    print(f"global_chambers={chambers_h2} bits={P*chambers_h2}")
    print(f"packets={3**12} accepted={len(accepted_h2)} decision_sha256={h2_digest}")
    for choice in accepted_h2:
        packet = packet_from_heights(choice)
        print(f"accepted_h={choice} W={packet} dilation_scale={dilation_scale(packet)}")
    print("descent_status=all accepted nodes are dilation sinks; non-dilation descent implication is vacuous in this box")

    print("\nCOUNTEREXAMPLES TO WEAKER DESCENT SUMMARIES")
    print(
        f"endpoint_current_only: all {endpoint_only_false_positives} packets with "
        "h_r in {1,2} have the universal correct endpoint but none passes the full word"
    )
    print(
        f"h_r in {{0,1}}: accepted={accepted_h01}; loose packets for which every legal "
        f"lowering makes first-tear time strictly earlier={len(descent_regressions)}"
    )
    for choice, tear, lowered in descent_regressions:
        packet = packet_from_heights(choice)
        print(f"regression W={packet} h={choice} tear={fmt(tear)} lowerings={[(r,fmt(t)) for r,t in lowered]}")
    print("conclusion=first-tear time, endpoint current, and scalar demand are not descent energies")

    print("\nADVERSARIAL TESTS (EVIDENCE ONLY)")
    print(
        f"dilation_scales={perturbation_scales}: accepted_all=True; "
        f"one-coordinate +13 failures={upward_failures}; legal one-coordinate -13 failures={downward_failures}"
    )
    print(
        f"deterministic_random heights 0..80: seed={random_seed} packets={random_count} "
        f"accepted={random_accepted} first_tear_denominator_top={random_tear_denominators.most_common(12)}"
    )

    print("\nPROOF-READY OPEN TARGET")
    print("A12 ballot rigidity: every arithmetic-realizable accepted root word is c*[12].")
    print("primitive form: gcd(W)=1 and accepted implies W=[12].")
    print("regeneration form: every accepted primitive mechanical word is the single Farey-12 comb.")
    print("sink form: every accepted coordinate-lowering sink is a dilation ray.")
    print("warning=the tested first-tear descent potential is false; a proof must retain the full prefix-minimum stalk/event word")

    print("\nTOURNAMENT / ASSUMPTION CHALLENGE")
    print("cut_vertices=13 sheets; observable=lexicographic current deficit; switch=root negation plus time reversal")
    print("tie_hamiltonian_path=0>1>...>12; fingerprint=transitive scores 0..12, zero cycles, 13 singleton SCCs")
    print("event_vertices=labelled crossings (r,k), pair observable=order of k/w_r, simultaneous ties contracted")
    print("preserves=mechanical chronology and root ballot predicate when root labels are retained")
    print("destroys=cut tournament loses chronology; unlabelled roots lose runner identity; arbitrary root words lose arithmetic realizability")


if __name__ == "__main__":
    main()
