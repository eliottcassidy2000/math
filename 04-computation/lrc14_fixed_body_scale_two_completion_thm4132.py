#!/usr/bin/env python3
"""Exact primary audit for THM-4132.

The proof joins THM-4129's closed U-safe interval to THM-3910's open
two-sheet quotient obstruction for the scale-two pair (1,9).
"""

from fractions import Fraction
from hashlib import sha256
from json import dumps


U = (1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22)
DELTA = Fraction(1, 14)
LEFT = Fraction(33, 70)
RIGHT = Fraction(27, 56)
J_LENGTH = RIGHT - LEFT
EXPECTED_C = (
    (Fraction(2, 21), Fraction(8, 63)),
    (Fraction(55, 63), Fraction(19, 21)),
)
EXPECTED_SEMANTIC = "aa6314ca8c25b5ad512bef0874dfb7262bd65cab8b200abba32d1cfe6b7710a6"


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def frac_pair(value):
    value = Fraction(value)
    return value.numerator, value.denominator


def semantic_digest(value):
    return sha256(
        dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def circle_distance(value):
    residue = value % 1
    return min(residue, 1 - residue)


def clearance(speeds, phase):
    return min(circle_distance(speed * phase) for speed in speeds)


def body_interval_audit():
    rows = []
    for speed in U:
        # If an integer occurred in speed*J, its distance would vanish.  With
        # no such integer, the concave tent has its minimum at an endpoint.
        first_integer = -(-(speed * LEFT).numerator // (speed * LEFT).denominator)
        last_integer = (speed * RIGHT).numerator // (speed * RIGHT).denominator
        require(first_integer > last_integer, f"integer crossing for body speed {speed}")
        left_gap = circle_distance(speed * LEFT)
        right_gap = circle_distance(speed * RIGHT)
        require(min(left_gap, right_gap) >= DELTA, f"unsafe body speed {speed}")
        rows.append((speed, frac_pair(left_gap), frac_pair(right_gap)))
    left_owners = tuple(speed for speed in U if circle_distance(speed * LEFT) == DELTA)
    right_owners = tuple(speed for speed in U if circle_distance(speed * RIGHT) == DELTA)
    require(left_owners == (15,), "left endpoint owner")
    require(right_owners == (4,), "right endpoint owner")
    require(J_LENGTH == Fraction(3, 280), "body interval length")
    return tuple(rows), left_owners, right_owners


def pair_bad(phase):
    return circle_distance(phase) < DELTA or circle_distance(9 * phase) < DELTA


def quotient_bad(w):
    """Both physical lifts z=w/2 and z+1/2 are pair-bad."""
    z = w / 2
    return pair_bad(z) and pair_bad(z + Fraction(1, 2))


def quotient_walls():
    walls = {Fraction(0), Fraction(1)}
    for shift in (Fraction(0), Fraction(1, 2)):
        for speed in (1, 9):
            for integer in range(speed):
                for sign in (-1, 1):
                    phase = (Fraction(integer) + sign * DELTA) / speed
                    walls.add((2 * (phase - shift)) % 1)
    return tuple(sorted(walls))


def quotient_components():
    walls = quotient_walls()
    active_cells = []
    for left, right in zip(walls, walls[1:]):
        midpoint = (left + right) / 2
        if quotient_bad(midpoint):
            active_cells.append((left, right))
    require(not quotient_bad(Fraction(0)), "quotient unexpectedly wraps at zero")
    components = []
    for cell in active_cells:
        if components and components[-1][1] == cell[0] and quotient_bad(cell[0]):
            components[-1] = (components[-1][0], cell[1])
        else:
            components.append(cell)
    for left, right in components:
        require(not quotient_bad(left) and not quotient_bad(right), "open endpoints")
        require(quotient_bad((left + right) / 2), "active component midpoint")
    return tuple(components)


def choose_literal_lift(t):
    """Construct a verified phase for one odd scale from the proof carrier."""
    require(t >= 3 and t % 2 == 1, "literal scale domain")
    candidates = {LEFT, RIGHT, (LEFT + RIGHT) / 2}
    for c_left, c_right in EXPECTED_C:
        for wall in (c_left, c_right):
            first = -(-(t * LEFT - wall).numerator // (t * LEFT - wall).denominator)
            last = (t * RIGHT - wall).numerator // (t * RIGHT - wall).denominator
            for integer in range(first - 1, last + 2):
                y = (wall + integer) / t
                if LEFT <= y <= RIGHT:
                    candidates.add(y)
    ordered = sorted(candidates)
    ordered += [
        (left + right) / 2
        for left, right in zip(ordered, ordered[1:])
    ]
    speeds = tuple(2 * speed for speed in U) + (t, 9 * t)
    for y in sorted(set(ordered)):
        if not (LEFT <= y <= RIGHT) or quotient_bad((t * y) % 1):
            continue
        for phase in (y / 2, (y + 1) / 2):
            if clearance(speeds, phase) >= DELTA:
                return y, phase, clearance(speeds, phase)
    raise RuntimeError(f"no literal lift at t={t}")


def main():
    body_rows, left_owners, right_owners = body_interval_audit()
    components = quotient_components()
    require(components == EXPECTED_C, "two-sheet quotient components")
    beta = max(right - left for left, right in components)
    require(beta == Fraction(2, 63), "quotient component maximum")

    first_scale_surplus = 3 * J_LENGTH - beta
    require(first_scale_surplus == Fraction(1, 2520), "first odd scale surplus")
    require(first_scale_surplus > 0, "compact-to-open gate")

    t1_phase = Fraction(1, 13)
    t1_speeds = tuple(2 * speed for speed in U) + (1, 9)
    t1_clearance = clearance(t1_speeds, t1_phase)
    require(t1_clearance == Fraction(1, 13), "t=1 clock")

    hostile_w = Fraction(1, 9)
    hostile_z = hostile_w / 2
    hostile = (
        circle_distance(hostile_z),
        circle_distance(9 * (hostile_z + Fraction(1, 2))),
    )
    require(hostile == (Fraction(1, 18), Fraction(0)), "fixed-sheet hostile")
    require(quotient_bad(hostile_w), "hostile not in quotient obstruction")

    literal_controls = []
    for t in range(3, 2002, 2):
        y, phase, gap = choose_literal_lift(t)
        require(gap >= DELTA, f"literal clearance at t={t}")
        require(len(set(tuple(2 * speed for speed in U) + (t, 9 * t))) == 13,
                f"distinct row at t={t}")
        literal_controls.append((t, frac_pair(y), frac_pair(phase), frac_pair(gap)))

    ledger = {
        "theorem": "THM-4132",
        "body": U,
        "body_interval": (frac_pair(LEFT), frac_pair(RIGHT), frac_pair(J_LENGTH)),
        "body_endpoint_owners": (left_owners, right_owners),
        "body_rows": body_rows,
        "quotient_components": tuple(
            (frac_pair(left), frac_pair(right)) for left, right in components
        ),
        "beta": frac_pair(beta),
        "first_odd_scale": 3,
        "first_scale_surplus": frac_pair(first_scale_surplus),
        "t1_clock": (frac_pair(t1_phase), frac_pair(t1_clearance)),
        "fixed_sheet_hostile": (frac_pair(hostile_w), tuple(map(frac_pair, hostile))),
        "literal_control_count": len(literal_controls),
        "literal_control_digest": semantic_digest(literal_controls),
        "scope": (
            "all positive odd t for 2U union {t,9t}; fixed U-body only; "
            "no arbitrary-body scale-two closure or LRC14"
        ),
    }
    digest = semantic_digest(ledger)
    require(digest == EXPECTED_SEMANTIC, "frozen semantic digest")
    print("status=PASS")
    print("theorem=THM-4132 fixed U-body exceptional scale-two completion")
    print(f"J=[{LEFT},{RIGHT}];length={J_LENGTH};owners={left_owners},{right_owners}")
    print(f"quotient_components={components};beta={beta}")
    print(f"first_odd_scale_surplus={first_scale_surplus}")
    print(f"t1_clock={t1_phase};clearance={t1_clearance}")
    print(f"fixed_sheet_hostile=w:{hostile_w};sheet_gaps:{hostile}")
    print(f"literal_controls={len(literal_controls)};range=odd 3..2001")
    print(f"literal_control_sha256={ledger['literal_control_digest']}")
    print("scope=fixed body 2U only; arbitrary bodies, physical entry, and LRC14 remain open")
    print(f"semantic_sha256={digest}")


if __name__ == "__main__":
    main()
