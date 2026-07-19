#!/usr/bin/env python3
"""Exact finite stress test for the A12 primitive essential-height bridge.

This audit connects three already separated facts without conflating them:

* the all-height shallow carrier: one positive speed w_r == r (mod 13) for
  each r in F_13^*, represented by its grouped A_12 root-current word;
* THM-770: the full-residue classification in the height box 0 <= h_r <= 12;
* THM-1142 (formerly the colliding THM-1125): the exact essential-region
  identity for a one-coordinate replacement.

For closed danger teeth

    D_w={t : ||wt|| <= 1/13}

and W=(w_1,...,w_12), define the open core-essential region

    E_r(W)=(R/Z) \ union_(s != r) D_(w_s).

If W covers the circle, then E_r(W) is contained in D_(w_r).  The exact
De Morgan identity, with no endpoint or measure loss, is

    (R/Z) \ (union_(s != r) D_(w_s) union D_(w_r-13))
       = E_r(W) \ D_(w_r-13).

Thus a legal lowering w_r -> w_r-13 preserves chi_13=0 if and only if

    E_r(W) subset D_(w_r-13).                              (E)

The smallest sufficient missing shallow lemma is:

  PRIMITIVE ESSENTIAL-HEIGHT DESCENT (PEHD13).
  If W is an accepted primitive full-residue packet and max h_r >= 13,
  then some coordinate with h_r >= 13 satisfies (E).

Iteration enters THM-770's height-12 box.  If at least one lowering occurred,
the final reverse step would be an accepted proper one-coordinate lift of a
dilated AP, contradicting THM-795.  PEHD13 therefore proves uniform *shallow*
rigidity.  It says nothing about THM-769's deep s>=2 branch and hence is not,
alone, n=12 sporadic emptiness or the compact LRC14 theorem.

The primitive clause is necessary.  Every dilation c[12] is accepted and is
a lowering-sink at arbitrarily large height by THM-795.  This script produces
literal strict witnesses in E_r(c[12]) \ D_(w_r-13) for all coordinates of
three high rays.

Finite exact tests
------------------
1. A direct 758-chamber meet-in-the-middle scan of all 4^12 height vectors
   h_r in {0,1,2,3}.  Exactly [12],2[12],3[12],4[12] survive.
2. Four sparse high shells h_r in {0,L}, L=13,25,50,100.  Each of the 4096
   rows is checked exactly; only [12] survives.
3. Full ternary coordinate cubes w_r in {w_r^0-13,w_r^0,w_r^0+13} around
   c[12] for c=14,17,97.  Each 3^12-row cube has only its centre ray.
4. Exact essential-region counterwitnesses for every legal lowering of those
   three accepted high rays.

All searches use the sheet decomposition t=(j+u)/13.  At a generic u,
runner r covers exactly

    {-floor(w_r u)r^-1, -(floor(w_r u)+1)r^-1}.

The truth value changes only at u=k/w_r.  The union of all candidate event
times gives a common atomic chamber complex; covering all thirteen sheets in
every chamber is equivalent to the original continuous closed-danger cover.
At an event the closed tooth contains the union of the adjacent generic
edges, so no endpoint case is dropped.

The exhaustive boxes are evidence and carrier validation, not PEHD13: they
contain no accepted non-dilation, so its implication is vacuous there.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from heapq import heapify, heappop, heappush
from itertools import product
from math import gcd
from pathlib import Path
from typing import Iterable


P = 13
RESIDUES = tuple(range(1, P))
DELTA = F(1, P)


def fmt(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def validate_packet(packet: tuple[int, ...]) -> None:
    assert len(packet) == 12
    assert all(speed > 0 for speed in packet)
    assert all(speed % P == residue for residue, speed in zip(RESIDUES, packet))


def packet_from_heights(height_vector: tuple[int, ...]) -> tuple[int, ...]:
    assert len(height_vector) == 12 and min(height_vector) >= 0
    packet = tuple(
        residue + P * height for residue, height in zip(RESIDUES, height_vector)
    )
    validate_packet(packet)
    return packet


def heights(packet: tuple[int, ...]) -> tuple[int, ...]:
    validate_packet(packet)
    return tuple((speed - residue) // P for residue, speed in zip(RESIDUES, packet))


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
    return common if sorted(packet) == [common * a for a in RESIDUES] else None


def packet_gcd(packet: tuple[int, ...]) -> int:
    common = 0
    for speed in packet:
        common = gcd(common, speed)
    return common


def distance_mod_one(value: F) -> F:
    phase = value % 1
    return min(phase, 1 - phase)


def is_strictly_core_safe(packet: tuple[int, ...], omitted: int, time: F) -> bool:
    return all(
        distance_mod_one(speed * time) > DELTA
        for coordinate, speed in enumerate(packet)
        if coordinate != omitted
    )


@dataclass(frozen=True)
class ScanResult:
    label: str
    rows: int
    chambers: int
    bits: int
    accepted: tuple[tuple[int, ...], ...]
    digest: str


def atomic_masks(
    option_groups: tuple[tuple[int, ...], ...]
) -> tuple[int, dict[tuple[int, int], int], int, str]:
    assert len(option_groups) == 12
    for residue, group in zip(RESIDUES, option_groups):
        assert group and len(set(group)) == len(group)
        assert all(speed > 0 and speed % P == residue for speed in group)

    all_options = tuple(speed for group in option_groups for speed in group)
    events = sorted(
        {F(integer, speed) for speed in all_options for integer in range(1, speed)}
    )
    boundaries = (F(0), *events, F(1))
    midpoints = tuple((left + right) / 2 for left, right in zip(boundaries, boundaries[1:]))
    masks: dict[tuple[int, int], int] = {}
    atlas_digest = sha256()
    for event in events:
        atlas_digest.update(f"{event.numerator}/{event.denominator}\n".encode())

    for residue, group in zip(RESIDUES, option_groups):
        inverse = pow(residue, -1, P)
        for option_index, speed in enumerate(group):
            mask = 0
            for chamber, u in enumerate(midpoints):
                integer = (speed * u.numerator) // u.denominator
                left = (-integer * inverse) % P
                right = (-(integer + 1) * inverse) % P
                assert left != right
                mask |= 1 << (P * chamber + left)
                mask |= 1 << (P * chamber + right)
            masks[(residue, option_index)] = mask
            atlas_digest.update(
                f"O|{residue}|{option_index}|{speed}|{mask.bit_count()}\n".encode()
            )

    bit_count = P * len(midpoints)
    return len(midpoints), masks, (1 << bit_count) - 1, atlas_digest.hexdigest()


def scan_options(label: str, option_groups: tuple[tuple[int, ...], ...]) -> ScanResult:
    chambers, masks, all_mask, atlas_digest = atomic_masks(option_groups)
    left_groups = option_groups[:6]
    right_groups = option_groups[6:]
    left_records = []
    right_records = []

    for choices in product(*(range(len(group)) for group in left_groups)):
        mask = 0
        speeds = []
        for residue, option_index, group in zip(range(1, 7), choices, left_groups):
            mask |= masks[(residue, option_index)]
            speeds.append(group[option_index])
        left_records.append((tuple(speeds), mask))

    for choices in product(*(range(len(group)) for group in right_groups)):
        mask = 0
        speeds = []
        for residue, option_index, group in zip(range(7, 13), choices, right_groups):
            mask |= masks[(residue, option_index)]
            speeds.append(group[option_index])
        right_records.append((tuple(speeds), mask))

    accepted = []
    for left_speeds, left_mask in left_records:
        needed = all_mask ^ left_mask
        for right_speeds, right_mask in right_records:
            if (needed & right_mask) == needed:
                accepted.append(left_speeds + right_speeds)

    rows = len(left_records) * len(right_records)
    assert rows == product_size(len(group) for group in option_groups)
    decision_digest = sha256()
    decision_digest.update(f"{label}|{rows}|{chambers}|{atlas_digest}\n".encode())
    for packet in accepted:
        decision_digest.update((",".join(map(str, packet)) + "\n").encode())
    return ScanResult(
        label,
        rows,
        chambers,
        P * chambers,
        tuple(accepted),
        decision_digest.hexdigest(),
    )


def product_size(values: Iterable[int]) -> int:
    result = 1
    for value in values:
        result *= value
    return result


def h_box_groups(max_height: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(residue + P * height for height in range(max_height + 1))
        for residue in RESIDUES
    )


def sparse_shell_groups(high_height: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        (residue, residue + P * high_height) for residue in RESIDUES
    )


def ray_cube_groups(scale: int) -> tuple[tuple[int, ...], ...]:
    centre = dilation_packet(scale)
    return tuple((speed - P, speed, speed + P) for speed in centre)


@dataclass(frozen=True)
class EssentialWitness:
    scale: int
    residue: int
    original: int
    lowered: int
    tear_u: F
    next_u: F
    sheet: int
    time: F


def first_tear(packet: tuple[int, ...]) -> tuple[F | None, tuple[int, ...]]:
    """First deficient generic chamber in the exact grouped root walk."""
    validate_packet(packet)
    heap: list[tuple[F, int, int, int, int]] = []
    for index, (residue, speed) in enumerate(zip(RESIDUES, packet)):
        if speed > 1:
            heap.append((F(1, speed), index, 1, speed, pow(residue, -1, P)))
    heapify(heap)
    excess = [11] + [0] * 12
    while heap:
        event_time = heap[0][0]
        updates = []
        while heap and heap[0][0] == event_time:
            _, index, integer, speed, inverse = heappop(heap)
            updates.append(
                ((-(integer - 1) * inverse) % P, (-(integer + 1) * inverse) % P)
            )
            if integer + 1 < speed:
                heappush(
                    heap,
                    (F(integer + 1, speed), index, integer + 1, speed, inverse),
                )
        for departure, arrival in updates:
            excess[departure] -= 1
            excess[arrival] += 1
        deficient = tuple(sheet for sheet, value in enumerate(excess) if value < 0)
        if deficient:
            return event_time, deficient
    return None, ()


def loose_proxy_regressions() -> tuple[int, tuple[tuple, ...]]:
    """Exhaust h in {0,1}; seek failure of first-tear descent before proof use."""
    cache = {}
    for height_vector in product((0, 1), repeat=12):
        cache[height_vector] = first_tear(packet_from_heights(height_vector))[0]

    regressions = []
    for height_vector, tear in cache.items():
        if tear is None:
            continue
        legal = [coordinate for coordinate, height in enumerate(height_vector) if height]
        if not legal:
            continue
        lowered_records = []
        for coordinate in legal:
            lowered = (
                height_vector[:coordinate]
                + (0,)
                + height_vector[coordinate + 1 :]
            )
            lowered_records.append((coordinate + 1, cache[lowered]))
        if all(lowered_tear is not None and lowered_tear < tear
               for _, lowered_tear in lowered_records):
            regressions.append((height_vector, tear, tuple(lowered_records)))

    accepted = sum(tear is None for tear in cache.values())
    assert accepted == 2
    assert len(regressions) == 2
    return accepted, tuple(regressions)


def lowering_tear(packet: tuple[int, ...], coordinate: int) -> EssentialWitness:
    """Find a strict point of E_r(W) \ D_(w_r-13) by the root walk."""
    validate_packet(packet)
    assert packet[coordinate] > P
    lowered = list(packet)
    lowered[coordinate] -= P
    lowered_packet = tuple(lowered)
    validate_packet(lowered_packet)

    heap: list[tuple[F, int, int, int, int]] = []
    for index, (residue, speed) in enumerate(zip(RESIDUES, lowered_packet)):
        if speed > 1:
            heap.append((F(1, speed), index, 1, speed, pow(residue, -1, P)))
    heapify(heap)
    excess = [11] + [0] * 12

    while heap:
        event_time = heap[0][0]
        updates = []
        while heap and heap[0][0] == event_time:
            _, index, integer, speed, inverse = heappop(heap)
            updates.append(
                ((-(integer - 1) * inverse) % P, (-(integer + 1) * inverse) % P)
            )
            if integer + 1 < speed:
                heappush(
                    heap,
                    (F(integer + 1, speed), index, integer + 1, speed, inverse),
                )
        for departure, arrival in updates:
            excess[departure] -= 1
            excess[arrival] += 1
        assert sum(excess) == 11
        deficient = [sheet for sheet, value in enumerate(excess) if value < 0]
        if deficient:
            next_time = heap[0][0] if heap else F(1)
            u = (event_time + next_time) / 2
            sheet = deficient[0]
            time = F(sheet, P) + u / P
            assert is_strictly_core_safe(packet, coordinate, time)
            assert distance_mod_one(lowered_packet[coordinate] * time) > DELTA
            assert distance_mod_one(packet[coordinate] * time) <= DELTA
            scale = dilation_scale(packet)
            assert scale is not None
            return EssentialWitness(
                scale,
                coordinate + 1,
                packet[coordinate],
                lowered_packet[coordinate],
                event_time,
                next_time,
                sheet,
                time,
            )
    raise AssertionError("lowered dilation coordinate unexpectedly remained accepted")


def essential_sink_audit(scales: tuple[int, ...]) -> tuple[EssentialWitness, ...]:
    witnesses = []
    for scale in scales:
        packet = dilation_packet(scale)
        for coordinate in range(12):
            witnesses.append(lowering_tear(packet, coordinate))
    return tuple(witnesses)


def print_scan(result: ScanResult) -> None:
    print(
        f"{result.label}: rows={result.rows} chambers={result.chambers} "
        f"bits={result.bits} accepted={len(result.accepted)} digest={result.digest}"
    )
    for packet in result.accepted:
        print(
            f"  W={packet} h={heights(packet)} gcd={packet_gcd(packet)} "
            f"dilation={dilation_scale(packet)}"
        )


def main() -> None:
    h3 = scan_options("height_box_0_3", h_box_groups(3))
    expected_h3 = tuple(dilation_packet(scale) for scale in (1, 2, 3, 4))
    assert h3.rows == 4**12
    assert h3.accepted == expected_h3

    sparse_results = []
    for high_height in (13, 25, 50, 100):
        result = scan_options(
            f"sparse_heights_0_{high_height}", sparse_shell_groups(high_height)
        )
        assert result.rows == 2**12
        assert result.accepted == (dilation_packet(1),)
        sparse_results.append(result)

    cube_results = []
    cube_scales = (14, 17, 97)
    for scale in cube_scales:
        result = scan_options(f"ray_cube_pm13_c{scale}", ray_cube_groups(scale))
        assert result.rows == 3**12
        assert result.accepted == (dilation_packet(scale),)
        cube_results.append(result)

    essential_witnesses = essential_sink_audit(cube_scales)
    assert len(essential_witnesses) == 36
    witness_digest = sha256()
    for witness in essential_witnesses:
        witness_digest.update(
            (
                f"{witness.scale}|{witness.residue}|{witness.original}|"
                f"{witness.lowered}|{fmt(witness.tear_u)}|{fmt(witness.next_u)}|"
                f"{witness.sheet}|{fmt(witness.time)}\n"
            ).encode()
        )

    proxy_accepted, proxy_regressions = loose_proxy_regressions()

    source_digest = sha256(Path(__file__).read_bytes()).hexdigest()
    print("A12 H3 / ESSENTIAL-HEIGHT REGENERATION BRIDGE AUDIT")
    print(f"source_sha256={source_digest}")
    print("closed_danger=D_w={t: ||wt|| <= 1/13}; essential_region E_r is open")
    print("exact_lowerability=E_r(W) subset D_(w_r-13); no a.e./endpoint relaxation")
    print("scope=shallow full-residue only; deep s>=2 and compact LRC14 residual are not proved")
    print("finite_claims=the explicit boxes below; PEHD13 remains an open all-height lemma")

    print("\nEXHAUSTIVE HEIGHT-THREE BOX")
    print_scan(h3)
    print("status=only four dilation rays; accepted non-dilation antecedent count is zero, so PEHD13 is vacuous here")

    print("\nEXACT SPARSE HIGH SHELLS")
    for result in sparse_results:
        print_scan(result)
    print("status=no primitive or nonprimitive high-shell survivor beyond the AP centre in these four slices")

    print("\nEXACT TERNARY +/-13 CUBES AROUND HIGH DILATIONS")
    for result in cube_results:
        print_scan(result)
    print("status=the full 3^12 coordinate cube has no accepted neighbour at any tested scale")

    print("\nESSENTIAL-REGION COUNTERWITNESSES TO UNQUALIFIED DESCENT")
    print(
        f"accepted_high_rays={cube_scales} legal_lowerings={len(essential_witnesses)} "
        f"lowerable=0 witness_sha256={witness_digest.hexdigest()}"
    )
    for scale in cube_scales:
        sample = next(
            witness for witness in essential_witnesses
            if witness.scale == scale and witness.residue == 1
        )
        print(
            f"c={scale} representative r=1: {sample.original}->{sample.lowered}, "
            f"tear_chamber_u=({fmt(sample.tear_u)},{fmt(sample.next_u)}), "
            f"sheet={sample.sheet}, strict_t={fmt(sample.time)} in E_1\\D_lower"
        )
    print("conclusion=primitive-or-dilation alternative is mandatory; accepted dilation rays are high sinks")

    print("\nNONACCEPTED DESCENT-POTENTIAL COUNTEREXAMPLES")
    print(
        f"height_box_0_1 accepted={proxy_accepted}; loose rows whose every legal lowering "
        f"tears strictly earlier={len(proxy_regressions)}"
    )
    for height_vector, tear, lowered_records in proxy_regressions:
        print(
            f"  W={packet_from_heights(height_vector)} tear={fmt(tear)} "
            f"lowerings={[(residue, fmt(lowered_tear)) for residue, lowered_tear in lowered_records]}"
        )
    print("proxy_verdict=first-tear time cannot orient an essential-height descent, even in the smallest box")

    print("\nSMALLEST EXACT SHALLOW PROOF OBLIGATION")
    print("PEHD13: accepted + primitive + full residues + max(h)>=13 => a coordinate r with h_r>=13 and E_r subset D_(w_r-13).")
    print("iteration: PEHD13 -> height<=12 -> THM-770 dilation; THM-795 forbids a nonempty reverse lift chain.")
    print("equivalent regeneration target: an accepted primitive mechanical geodesic has only the single 46-chamber Farey-12 comb.")
    print("negative=first-tear time and THM-1010 max-speed descent are not valid potentials; containment must retain the whole essential stalk.")

    print("\nBRIDGE TO UNIFORM n=12 / COMPACT LRC14")
    print("positive=PEHD13 would close arbitrary winding in THM-769's shallow full-residue branch.")
    print("missing_after_PEHD13=THM-769 deep binding scales s>=2 (at least two off-sheet runners) remain separate.")
    print("compact=rho<13 plus the auxiliary 1/13 cover is a stronger sufficient compact target, not the full sharp LRC14 residual; PEHD13 feeds it only when a shallow 12-core is first forced.")
    print("therefore=no valid direct inference from the h<=3 certificate to sporadic emptiness or INVcov.")

    print("\nCARRIER / TOURNAMENT AUDIT")
    print("predicate_vertices=13 sheet cuts; state=11-chip excess plus grouped labelled A12 root word")
    print("pair_observable=wall-crossing order k/w_r; gauge=time reversal/root negation; ties=contracted before testing")
    print("tie_hamiltonian_path=chronological event order, then sheet labels 0..12")
    print("preserves=closed danger coverage and strict essential witnesses")
    print("destroys=cut-deficit tournament loses chronology; unlabelled current loses runner identity; scalar descent loses E_r position")


if __name__ == "__main__":
    main()
