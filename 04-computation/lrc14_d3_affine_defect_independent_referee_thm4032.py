#!/usr/bin/env python3
"""Independent exact referee for the proposed d=3 affine-defect theorem.

This checker deliberately avoids the primary audit's wall/midpoint sweep and
its affine-centre box enumeration.  The semantic engine instead enumerates
the literal open y-intervals on which a named speed kills a named lift.  The
arithmetic engine enumerates only the proposed bounded defect triples and
then reconstructs centres by solving the two congruences modulo a.
"""

from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")
CHECKS = 0


def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def kill_components(speed, label):
    """Open y-intervals in [0,1] where speed kills lift label."""
    radius = Q(3, 14 * speed)
    answer = []
    # For 0<=y<=1 and 0<=label<=2, speed*(y+label)/3 lies in [0,speed].
    # The two spare integers are hostile guards on the range derivation.
    for integer in range(-1, speed + 2):
        centre = Q(3 * integer, speed) - label
        left, right = centre - radius, centre + radius
        if right > 0 and left < 1:
            answer.append((left, right, integer))
    return tuple(answer)


def assignment_witness(order):
    """Literal strict witness when order[label] is assigned to that lift."""
    banks = tuple(kill_components(speed, label) for label, speed in enumerate(order))
    for chosen in product(*banks):
        left = max(Q(0), *(item[0] for item in chosen))
        right = min(Q(1), *(item[1] for item in chosen))
        if left < right:
            return (left + right) / 2, tuple(item[2] for item in chosen)
    return None


def semantic_witness(speeds):
    """Full three-lift spoilage, independently expanded over all 3! owners."""
    for order in permutations(speeds):
        witness = assignment_witness(order)
        if witness is not None:
            return order, witness
    return None


def bound(x, y):
    return (3 * (x + y) - 1) // 14


def reconstruct(a, b, c, n_ab, n_ac, n_bc):
    """Solve the affine-centre equations without using the primary box scan."""
    k_ab = (n_ab - a * b) // 3
    k_ac = (n_ac - 2 * a * c) // 3
    k_bc = (n_bc - b * c) // 3
    for A in range(a):
        if (b * A - k_ab) % a or (c * A - k_ac) % a:
            continue
        B = (b * A - k_ab) // a
        C = (c * A - k_ac) // a
        require(c * B - b * C == k_bc, ("third equation", a, b, c))

        centres = (Q(A, a), Q(B, b) - Q(1, 3), Q(C, c) - Q(2, 3))
        radii = (Q(1, 14 * a), Q(1, 14 * b), Q(1, 14 * c))
        left = max(x - r for x, r in zip(centres, radii))
        right = min(x + r for x, r in zip(centres, radii))
        require(left < right, ("reconstructed intervals", a, b, c))

        # This proves the advertised finite box after fixing A modulo a.
        require(0 <= B <= 2 * b, ("B box", a, b, c, A, B, C))
        require(0 <= C <= 2 * c, ("C box", a, b, c, A, B, C))
        return A, B, C, (left + right) / 2
    raise AssertionError(("locally soluble defects failed to glue", a, b, c,
                          n_ab, n_ac, n_bc))


def oriented_certificate(a, b, c):
    for n_ab in range(-bound(a, b), bound(a, b) + 1):
        if n_ab % 3 != (a * b) % 3 or n_ab % gcd(a, b):
            continue
        for n_ac in range(-bound(a, c), bound(a, c) + 1):
            if n_ac % 3 != (2 * a * c) % 3 or n_ac % gcd(a, c):
                continue
            numerator = b * n_ac - c * n_ab
            if numerator % a:
                continue
            n_bc = numerator // a
            if abs(n_bc) > bound(b, c):
                continue
            if n_bc % 3 != (b * c) % 3 or n_bc % gcd(b, c):
                continue

            # Unit residues make nonzero defects and endpoint equality automatic.
            require(n_ab and n_ac and n_bc, ("automatic nonzero", a, b, c))
            require(14 * abs(n_ab) != 3 * (a + b), "ab endpoint tie")
            require(14 * abs(n_ac) != 3 * (a + c), "ac endpoint tie")
            require(14 * abs(n_bc) != 3 * (b + c), "bc endpoint tie")
            centre = reconstruct(a, b, c, n_ab, n_ac, n_bc)
            return (n_ab, n_ac, n_bc), centre
    return None


def certificate(speeds):
    a, b, c = speeds
    return oriented_certificate(a, b, c) or oriented_certificate(a, c, b)


def weakened_certificate(order, omitted):
    """First fake local triple when one proposed sidecar is omitted."""
    a, b, c = order

    def bank(x, y, residue):
        answer = []
        for n in range(-bound(x, y), bound(x, y) + 1):
            if not n or n % gcd(x, y):
                continue
            if omitted != "mod3" and n % 3 != residue % 3:
                continue
            answer.append(n)
        return answer

    ab = bank(a, b, a * b)
    ac = bank(a, c, 2 * a * c)
    bc = bank(b, c, b * c)
    for n_ab, n_ac, n_bc in product(ab, ac, bc):
        circuit = c * n_ab - b * n_ac + a * n_bc
        if omitted == "mod3" and circuit:
            continue
        # With the circuit omitted, every locally admissible triple is a
        # witness to the weakened predicate; record its nonzero defect.
        return n_ab, n_ac, n_bc, circuit
    return None


def first_hostile(omitted, ceiling=48):
    """Least (maximum speed, lexicographic profile, orientation, fake N)."""
    for maximum in range(1, ceiling + 1):
        values = tuple(x for x in range(1, maximum + 1) if x % 3)
        for speeds in combinations(values, 3):
            if speeds[-1] != maximum or semantic_witness(speeds) is not None:
                continue
            a, b, c = speeds
            for order in ((a, b, c), (a, c, b)):
                fake = weakened_certificate(order, omitted)
                if fake is not None:
                    return maximum, speeds, order, fake
    return None


def direct_bad_labels(y, speeds):
    labels = []
    for speed in speeds:
        bad = []
        for label in range(3):
            value = (Q(speed) * (y + label) / 3) % 1
            if min(value, 1 - value) < Q(1, 14):
                bad.append(label)
        require(len(bad) <= 1, ("one-label capacity", speed, y, bad))
        labels.append(tuple(bad))
    return tuple(labels)


def main():
    # Wider than the primary 1..23 audit and based on a different semantic engine.
    values = tuple(x for x in range(1, 33) if x % 3)
    profiles = 0
    positive = 0
    for speeds in combinations(values, 3):
        semantic = semantic_witness(speeds)
        arithmetic = certificate(speeds)
        require((semantic is not None) == (arithmetic is not None),
                ("semantic/arithmetic", speeds, semantic, arithmetic))
        if semantic is not None:
            positive += 1
            order, (y, _) = semantic
            killed = direct_bad_labels(y, order)
            require(tuple(item[0] for item in killed) == (0, 1, 2),
                    ("literal assigned labels", speeds, order, y, killed))
        profiles += 1

    require(semantic_witness((1, 4, 5)) is not None, "minimal positive control")
    require(certificate((1, 4, 5)) is not None, "minimal arithmetic control")
    require(semantic_witness((1, 4, 7)) is None, "gcd-gate false converse")
    require(certificate((1, 4, 7)) is None, "gcd-gate arithmetic hostile")

    no_mod3 = first_hostile("mod3")
    no_circuit = first_hostile("circuit")
    require(no_mod3 == (8, (2, 7, 8), (2, 7, 8), (-1, -2, -3, 0)),
            ("first no-mod3 hostile", no_mod3))
    require(no_circuit == (7, (1, 4, 7), (1, 4, 7), (1, -1, -2, 9)),
            ("first no-circuit hostile", no_circuit))

    print("LRC14_D3_AFFINE_DEFECT_INDEPENDENT_REFEREE")
    print("semantic_engine=all_6_literal_open_interval_assignments")
    print("arithmetic_engine=bounded_defects_plus_congruence_reconstruction")
    print(f"profiles={profiles};positive={positive};range=1..32_nonmultiples_of_3")
    print("box=A_mod_a,B_0_to_2b,C_0_to_2c:verified_on_every_reconstruction")
    print("mod3=>nonzero_and_no_endpoint_ties")
    print(f"first_omit_mod3_hostile={no_mod3}")
    print(f"first_omit_circuit_hostile={no_circuit}")
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
