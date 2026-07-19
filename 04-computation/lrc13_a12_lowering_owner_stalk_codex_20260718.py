#!/usr/bin/env python3
"""All-height local owner-stalk criterion for lowering an A12 coordinate.

Let W=(w_r), w_r == r (mod 13), be a shallow full-residue packet and let

    D_w={t in R/Z : ||wt|| <= 1/13}

be the closed danger comb.  For a coordinate r with w_r>13, put
w_r^-=w_r-13 and

    E_r(W)=(R/Z) \ union_(s != r) D_(w_s).

The exact essential-region identity says that the lowering is accepted iff

    E_r(W) subset D_(w_r^-).                               (1)

This script proves and audits a local sheet version of (1).  Write
t=(j+u)/13.  In a generic u-chamber the r-edge is

    Q_r(u)={A,B}={-k r^-1,-(k+1)r^-1},  k=floor(w_r u).

If W is accepted, the core-uncovered sheets are exactly the sheets uniquely
owned by r,

    U_r(u)={j in Q_r(u) : degree_W(j)=1}.                  (2)

Hence (1) is equivalent to U_r(u) subset Q_r^-(u) in every generic chamber.

There is a closed functional form for the lowered edge.  Put

    13u=m+theta,  0<theta<1,
    w_r u=k+alpha, 0<alpha<1.

Then

    a(u)=k-floor((w_r-13)u)=m+1_(alpha<theta),             (3)
    Q_r^-(u)=Q_r(u)+a(u)r^-1.                              (4)

Since Q_r is a two-point needle on the 13-cycle, its intersection with its
translate is

    a=0:  {A,B};       a=1:  {A};       a=-1: {B};
    otherwise: empty.                                      (5)

Equations (2)--(5) give the exact owner-stalk test:

    a=0  => any unique-owner subset is allowed;
    a=1  => U_r subset {A};
    a=12 => U_r subset {B};
    else => U_r is empty.                                  (6)

In particular, on the whole central cylinder

    2/13 < u < 11/13,

we have m in {2,...,10}, so a is in {2,...,11}; the old and lowered edges
are disjoint.  A lowerable coordinate therefore has *no unique ownership at
all* in that central band.  The only remaining obligations live on the four
fringe strips and retain a one-sided endpoint bit.

Endpoint completeness
---------------------
The essential region is open and D_(w_r^-) is closed.  If containment fails
at any wall or thirteenth boundary, the failure set E_r\D_(w_r^-) is open and
contains a generic point.  Thus the common refinement by all old and lowered
wall times loses no boundary case.  Equal wall times are not ordered.

The accepted hypothesis in (2) is essential.  For a loose W, a sheet can be
uncovered by the core and by r; it is then an essential obligation but not a
unique-owned sheet.  A literal two-lift counterexample is audited below.

Relation to the open bridge
---------------------------
PEHD13 would follow from: every accepted primitive full-residue packet with
some height >=13 has a high coordinate satisfying (6) in every chamber.
The central-band consequence is a sharper first subtarget, but not sufficient
without the oriented fringe conditions.  Dilation rays have central unique
ownership in every legal coordinate and are lowering-sinks, so primitivity is
again indispensable.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import gcd
from pathlib import Path
from random import Random


P = 13
DELTA = F(1, P)
RESIDUES = tuple(range(1, P))


def fmt(value: F | None) -> str:
    if value is None:
        return "none"
    return f"{value.numerator}/{value.denominator}"


def validate(packet: tuple[int, ...]) -> None:
    assert len(packet) == 12
    assert all(speed > 0 for speed in packet)
    assert all(speed % P == residue for residue, speed in zip(RESIDUES, packet))


def packet_from_heights(height_vector: tuple[int, ...]) -> tuple[int, ...]:
    packet = tuple(
        residue + P * height for residue, height in zip(RESIDUES, height_vector)
    )
    validate(packet)
    return packet


def dilation_packet(scale: int) -> tuple[int, ...]:
    assert scale > 0 and scale % P
    packet = [0] * 12
    for multiplier in RESIDUES:
        speed = scale * multiplier
        packet[speed % P - 1] = speed
    result = tuple(packet)
    validate(result)
    return result


def packet_gcd(packet: tuple[int, ...]) -> int:
    common = 0
    for speed in packet:
        common = gcd(common, speed)
    return common


def edge(speed: int, residue: int, u: F) -> tuple[int, int]:
    integer = (speed * u.numerator) // u.denominator
    inverse = pow(residue, -1, P)
    return ((-integer * inverse) % P, (-(integer + 1) * inverse) % P)


def distance_mod_one(value: F) -> F:
    phase = value % 1
    return min(phase, 1 - phase)


@dataclass(frozen=True)
class Failure:
    u: F
    strip: int
    shift: int
    core_uncovered: tuple[int, ...]
    unique_owned: tuple[int, ...]
    lowered_edge: tuple[int, int]


@dataclass(frozen=True)
class CoordinateAudit:
    original_accepted: bool
    lowering_accepted: bool
    unique_proxy_accepted: bool
    has_central_unique_owner: bool
    first_essential_failure: Failure | None
    first_unique_failure: Failure | None
    chambers: int


def coordinate_audit(packet: tuple[int, ...], coordinate: int) -> CoordinateAudit:
    """Audit essential containment, unique proxy, and formulas (3)--(6)."""
    validate(packet)
    assert 0 <= coordinate < 12 and packet[coordinate] > P
    residue = coordinate + 1
    speed = packet[coordinate]
    lowered_speed = speed - P
    candidate_speeds = (*packet, lowered_speed)
    events = sorted(
        {
            F(integer, candidate)
            for candidate in candidate_speeds
            for integer in range(1, candidate)
        }
    )
    boundaries = (F(0), *events, F(1))

    original_accepted = True
    lowering_accepted = True
    unique_proxy_accepted = True
    central_unique = False
    first_essential_failure = None
    first_unique_failure = None

    for left_boundary, right_boundary in zip(boundaries, boundaries[1:]):
        u = (left_boundary + right_boundary) / 2
        old_edges = tuple(
            edge(candidate, owner, u)
            for owner, candidate in zip(RESIDUES, packet)
        )
        degrees = [0] * P
        for old_edge in old_edges:
            for sheet in old_edge:
                degrees[sheet] += 1
        original_accepted &= min(degrees) >= 1

        core_uncovered = tuple(
            sheet
            for sheet in range(P)
            if all(
                sheet not in old_edges[other]
                for other in range(12)
                if other != coordinate
            )
        )
        unique_owned = tuple(
            sheet for sheet in old_edges[coordinate] if degrees[sheet] == 1
        )
        if min(degrees) >= 1:
            assert set(core_uncovered) == set(unique_owned)

        lowered_edge = edge(lowered_speed, residue, u)
        essential_ok = set(core_uncovered) <= set(lowered_edge)
        unique_ok = set(unique_owned) <= set(lowered_edge)
        lowering_accepted &= essential_ok
        unique_proxy_accepted &= unique_ok

        old_integer = (speed * u.numerator) // u.denominator
        new_integer = (lowered_speed * u.numerator) // u.denominator
        raw_shift = old_integer - new_integer
        strip = (P * u.numerator) // u.denominator
        theta = P * u - strip
        alpha = speed * u - old_integer
        assert 0 <= strip <= 12
        assert F(0) < theta < F(1)
        assert F(0) < alpha < F(1)
        assert raw_shift == strip + int(alpha < theta)
        shift = raw_shift % P

        inverse = pow(residue, -1, P)
        transported = tuple(
            (sheet + shift * inverse) % P for sheet in old_edges[coordinate]
        )
        assert set(transported) == set(lowered_edge)
        leading, trailing = old_edges[coordinate]
        if shift == 0:
            permitted = {leading, trailing}
        elif shift == 1:
            permitted = {leading}
        elif shift == P - 1:
            permitted = {trailing}
        else:
            permitted = set()
        assert (set(unique_owned) <= set(lowered_edge)) == (
            set(unique_owned) <= permitted
        )

        if 2 <= strip <= 10:
            assert shift not in (0, 1, P - 1)
            assert set(old_edges[coordinate]).isdisjoint(lowered_edge)
            central_unique |= bool(unique_owned)

        failure = Failure(
            u,
            strip,
            shift,
            core_uncovered,
            unique_owned,
            lowered_edge,
        )
        if not essential_ok and first_essential_failure is None:
            first_essential_failure = failure
        if not unique_ok and first_unique_failure is None:
            first_unique_failure = failure

    if original_accepted:
        assert lowering_accepted == unique_proxy_accepted
    return CoordinateAudit(
        original_accepted,
        lowering_accepted,
        unique_proxy_accepted,
        central_unique,
        first_essential_failure,
        first_unique_failure,
        len(boundaries) - 1,
    )


def audit_ray_sinks(scales: tuple[int, ...]) -> tuple[int, int, str, tuple]:
    legal = 0
    central_obstructions = 0
    digest = sha256()
    representatives = []
    for scale in scales:
        packet = dilation_packet(scale)
        assert packet_gcd(packet) == scale
        scale_representative = None
        for coordinate, speed in enumerate(packet):
            if speed <= P:
                continue
            legal += 1
            audit = coordinate_audit(packet, coordinate)
            assert audit.original_accepted
            assert not audit.lowering_accepted
            assert not audit.unique_proxy_accepted
            assert audit.has_central_unique_owner
            central_obstructions += 1
            failure = audit.first_unique_failure
            assert failure is not None
            digest.update(
                (
                    f"{scale}|{coordinate+1}|{speed}|{fmt(failure.u)}|"
                    f"{failure.strip}|{failure.shift}|{failure.unique_owned}|"
                    f"{failure.lowered_edge}\n"
                ).encode()
            )
            if coordinate == 0:
                scale_representative = (scale, coordinate + 1, speed, failure)
        if scale_representative is not None:
            representatives.append(scale_representative)
    return legal, central_obstructions, digest.hexdigest(), tuple(representatives)


def random_formula_audit(count: int, seed: int) -> int:
    rng = Random(seed)
    audited = 0
    for _ in range(count):
        height_vector = tuple(rng.randrange(0, 6) for _ in RESIDUES)
        legal = [index for index, height in enumerate(height_vector) if height]
        if not legal:
            continue
        coordinate = legal[rng.randrange(len(legal))]
        coordinate_audit(packet_from_heights(height_vector), coordinate)
        audited += 1
    return audited


def main() -> None:
    ray_scales = (2, 3, 4, 14, 17, 97)
    legal, central_obstructions, ray_digest, representatives = audit_ray_sinks(
        ray_scales
    )
    assert legal == central_obstructions == 59

    # Positive De Morgan/control case: lowering the only lift returns to [12].
    one_lift = packet_from_heights((1,) + (0,) * 11)
    positive_control = coordinate_audit(one_lift, 0)
    assert not positive_control.original_accepted
    assert positive_control.lowering_accepted
    assert positive_control.unique_proxy_accepted

    # Accepted is indispensable for replacing the full essential stalk by
    # "unique ownership".  Runner 12 never uniquely owns a violating sheet,
    # because the original loose packet already leaves that sheet uncovered.
    loose_proxy_packet = (1, 2, 3, 4, 5, 6, 7, 8, 9, 23, 11, 25)
    validate(loose_proxy_packet)
    loose_proxy = coordinate_audit(loose_proxy_packet, 11)
    assert not loose_proxy.original_accepted
    assert not loose_proxy.lowering_accepted
    assert loose_proxy.unique_proxy_accepted
    assert not loose_proxy.has_central_unique_owner
    proxy_failure = loose_proxy.first_essential_failure
    assert proxy_failure is not None
    assert proxy_failure.u == F(45, 506)
    assert proxy_failure.core_uncovered == (9,)
    assert proxy_failure.unique_owned == ()
    assert set(proxy_failure.lowered_edge) == {1, 2}

    random_count = 200
    random_seed = 20260718
    random_audited = random_formula_audit(random_count, random_seed)
    assert random_audited == random_count

    source_digest = sha256(Path(__file__).read_bytes()).hexdigest()
    print("A12 LOWERING OWNER-STALK / CENTRAL-CYLINDER LEMMA")
    print(f"source_sha256={source_digest}")
    print("domain=all shallow full-residue packets; unique-owner equivalence assumes original packet accepted")
    print("signs=D_w closed, E_r open; common generic chamber refinement is endpoint-complete")
    print("shift_law=raw_a=floor(wu)-floor((w-13)u)=m+1_{frac(wu)<frac(13u)}; bar_a=raw_a mod 13")
    print("local_rule=bar_a=0:any (raw_a=0 or 13); bar_a=1:leading endpoint only; bar_a=12:trailing endpoint only; otherwise no unique owner")
    print("central_rule=2/13<u<11/13 implies old/lowered edges disjoint, so a lowerable accepted colour has zero unique ownership")

    print("\nACCEPTED DILATION SINK AUDIT")
    print(
        f"scales={ray_scales} legal_lowerings={legal} lowerable=0 "
        f"central_unique_obstructions={central_obstructions} digest={ray_digest}"
    )
    for scale, residue, speed, failure in representatives:
        print(
            f"c={scale} r={residue} w={speed}: first owner failure "
            f"u={fmt(failure.u)} strip={failure.strip} shift={failure.shift} "
            f"unique={failure.unique_owned} lowered_edge={failure.lowered_edge}"
        )
    print("verdict=unqualified height descent is false at arbitrarily high accepted rays; primitive/dilation split is necessary")

    print("\nSTRICT HYPOTHESIS CONTROLS")
    print(
        f"positive one-lift W={one_lift}: original_accepted={positive_control.original_accepted} "
        f"lowering_accepted={positive_control.lowering_accepted} unique_proxy={positive_control.unique_proxy_accepted}"
    )
    print(
        f"loose proxy W={loose_proxy_packet}, lower r=12: original_accepted={loose_proxy.original_accepted} "
        f"lowering_accepted={loose_proxy.lowering_accepted} unique_proxy={loose_proxy.unique_proxy_accepted}"
    )
    print(
        f"proxy failure u={fmt(proxy_failure.u)} core_uncovered={proxy_failure.core_uncovered} "
        f"unique={proxy_failure.unique_owned} lowered_edge={proxy_failure.lowered_edge}"
    )
    print("verdict=without acceptedness, unique-owned sheets omit pre-existing holes and do not represent E_r")

    print("\nFORMULA REPLAY")
    print(
        f"deterministic random packet-coordinate audits={random_audited} seed={random_seed}; "
        "shift, transport, intersection, central-disjointness assertions all exact"
    )

    print("\nSHARPENED OPEN BRIDGE")
    print("PEHD13-owner form: accepted primitive max(h)>=13 => a high colour satisfying the local owner rule in every chamber.")
    print("central sublemma: such a packet has a high colour with no unique ownership on 2/13<u<11/13.")
    print("fringe sublemma: one central-inessential high colour obeys the leading/trailing endpoint rules on the four fringe strips.")
    print("warning=central sublemma alone is necessary, not sufficient; oriented fringe bits cannot be averaged away.")
    print("scope=would close only shallow winding with THM-770/795; THM-769 deep sheets and compact LRC14 remain open.")

    print("\nKAKEYA / TOURNAMENT CARRIER AUDIT")
    print("needles=the twelve labelled two-sheet edges moving on the 13-cut cylinder; unique trace is the essential stalk")
    print("pair_observable=which colour privately owns a cut in a chamber; switch=lowering translation bar_a*r^-1; tie path=chronological walls")
    print("preserves=exact E_r containment for accepted packets, including endpoint orientation")
    print("destroys=degree-only state loses owner label; central-owner count loses fringe side; runner tournament loses chamber chronology")


if __name__ == "__main__":
    main()
