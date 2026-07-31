#!/usr/bin/env python3
"""Exact semantic gate audit for the half-odometer and C221 cycles.

This dependency-free probe compares the locally packet-typed cycles in
THM-2698 and the transverse C221 scout with the literal THM-2305 endpoint
words on the canonical typed row.  It also checks the changed half-base
``B1(y)={13y+1/2}`` directly.

The conclusions are deliberately negative and typed:

* both displayed delayed cycles fail the literal endpoint word at the
  unshifted source-c1 safety factor;
* the physical event centres are not THM-2305 exclusive sources;
* all prescribed clocks 2, 4, 6 are even, and B1 agrees exactly with the
  ordinary base map at those clocks;
* every six-state literal Qa/Qb word under B1 is empty by two factor-sparse
  implications, although a strict four-state word exists;
* the C221 phase labels form the inverse-pair equality case of THM-2644,
  while the half-cycle projects to its sharp even-group hostile.

No scalar row, endpoint current, row exclusion, or LRC(14) conclusion is
asserted.
"""

from fractions import Fraction
from itertools import product
from math import gcd


P = 13
R = P**6
GUARD = 1
UNITS = (14, 27, 40, 53, 66)
C1 = P
CA = P**3
CB = 2 * P**5
CLOCKS = (2, 4, 6)
ROW = (GUARD,) + UNITS + (C1, CA, CB)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def centered(value):
    """Half-open representative in [-1/2,1/2)."""
    residue = value % 1
    return residue if residue < Fraction(1, 2) else residue - 1


def distance(value):
    return abs(centered(value))


def danger(speed, value, denominator=14):
    phase = centered(speed * value)
    radius = Fraction(1, denominator)
    return -radius <= phase < radius


def common_safe(value):
    return (
        not danger(GUARD, value, 7)
        and all(not danger(speed, value) for speed in UNITS)
    )


def exclusive_source(value, owner):
    require(owner in range(3), "unknown exclusive owner")
    blockers = (C1, CA, CB)
    return (
        common_safe(value)
        and danger(blockers[owner], value)
        and all(not danger(speed, value)
                for index, speed in enumerate(blockers) if index != owner)
    )


def terminal_word(value, source_owner):
    """Return the complete THM-2305 word, or None outside R_j."""
    require(source_owner in range(3), "unknown source owner")
    blockers = (C1, CA, CB)
    if not common_safe(value) or danger(blockers[source_owner], value):
        return None
    return tuple(
        index + 1 for index, speed in enumerate(blockers)
        if index != source_owner and danger(speed, value)
    ) or None


def literal_word(value, target):
    """Literal Qa=E_a or Qb=E_b from THM-2701."""
    require(target in ("a", "b"), "unknown literal target")
    if not common_safe(value) or danger(C1, value):
        return False
    a_hit = danger(CA, value)
    b_hit = danger(CB, value)
    return (a_hit and not b_hit) if target == "a" else (b_hit and not a_hit)


def literal_bits(value):
    return (
        not danger(GUARD, value, 7),
        tuple(not danger(speed, value) for speed in UNITS),
        not danger(C1, value),
        danger(CA, value),
        danger(CB, value),
    )


def B0(value):
    return (P * value) % 1


def B1(value):
    return (P * value + Fraction(1, 2)) % 1


def iterate(function, value, steps):
    for _ in range(steps):
        value = function(value)
    return value


def literal_margin(value, target):
    rows = [distance(value) - Fraction(1, 7)]
    rows.extend(distance(speed * value) - Fraction(1, 14)
                for speed in UNITS)
    rows.append(distance(C1 * value) - Fraction(1, 14))
    hit = CA if target == "a" else CB
    miss = CB if target == "a" else CA
    rows.append(Fraction(1, 14) - distance(hit * value))
    rows.append(distance(miss * value) - Fraction(1, 14))
    return min(rows)


def audit_displayed_cycles():
    half_points = (
        Fraction(55_232_507, 24 * R),
        Fraction(58_313_459, 24 * R),
    )
    half_lifts = (2_945_947, 9_005_121)
    c221_points = (
        Fraction(39_123_022, 17 * R),
        Fraction(41_305_508, 17 * R),
    )
    c221_lifts = (25_040_740, 76_541_689)

    for index in range(2):
        require(
            (P * half_points[index]
             + Fraction(half_lifts[index], 2 * R)) % 1
            == half_points[1 - index],
            "THM-2698 affine cycle changed",
        )
        require(
            (P * c221_points[index]
             + Fraction(c221_lifts[index], 17 * R)) % 1
            == c221_points[1 - index],
            "C221 affine cycle changed",
        )

    half_delayed = tuple((R * value) % 1 for value in half_points)
    c221_delayed = tuple((R * value) % 1 for value in c221_points)
    require(half_delayed == (Fraction(11, 24),) * 2,
            "THM-2698 delayed fixed phase changed")
    require(c221_delayed == (Fraction(4, 17), Fraction(13, 17)),
            "C221 delayed two-cycle changed")

    delayed_rows = []
    for family, values in (("half", half_delayed), ("C221", c221_delayed)):
        for value in values:
            bits = literal_bits(value)
            require(bits[0] and all(bits[1]) and not bits[2]
                    and bits[3] and not bits[4],
                    f"{family} delayed endpoint defect changed")
            require(not literal_word(value, "a")
                    and not literal_word(value, "b"),
                    f"{family} delayed endpoint became semantic")
            delayed_rows.append((family, value, distance(C1 * value), bits))

    physical_rows = []
    for family, values in (("half", half_points), ("C221", c221_points)):
        for event, value in enumerate(values):
            sources = tuple(owner + 1 for owner in range(3)
                            if exclusive_source(value, owner))
            require(not sources, f"{family} event entered an exclusive source")
            terminal = iterate(B0, value, 6)
            word = terminal_word(terminal, 2)
            require(word == (1, 2),
                    f"{family} deepest-clock fork endpoint changed")
            require(not exclusive_source(value, 2),
                    f"{family} deepest-clock source unexpectedly appeared")
            physical_rows.append((family, event, value, sources, terminal, word))

    return half_delayed, c221_delayed, tuple(delayed_rows), tuple(physical_rows)


def audit_half_semantic_language():
    # B1^n(y)=13^n y+(13^n-1)/24.  At every prescribed (even) clock,
    # the affine constant is integral, so B1^k=B0^k on the whole circle.
    clock_constants = tuple(
        (clock, Fraction(P**clock - 1, 24)) for clock in CLOCKS
    )
    require(all(constant.denominator == 1
                for _clock, constant in clock_constants),
            "half phase stopped cancelling at a prescribed clock")
    probes = (Fraction(0), Fraction(1, 17), Fraction(11, 24),
              Fraction(12894291, 80000000))
    for clock in CLOCKS:
        for value in probes:
            require(iterate(B1, value, clock) == iterate(B0, value, clock),
                    "even-clock endpoint identity failed")

    # Factor-sparse a obstruction.  B1^2=B0^2, hence c1 at state j+2
    # is exactly ca at state j.  A literal a therefore forbids every literal
    # endpoint two states later.
    require(C1 * P**2 == CA, "a/source two-clock identity changed")

    # Factor-sparse b obstruction.  For z=B1^5(y), the two target-b teeth at
    # y and B1(y) become D_2(z) and D_26(z).  The exact half-open contraction
    # from THM-2701 forces D_14(z), contradicting every literal endpoint at
    # state five.  The half shift disappears from all three even speeds.
    require(CB == 2 * P**5, "first b-to-z identity changed")
    require(CB * P == 26 * P**5, "second b-to-z identity changed")
    require(26 * Fraction(1, 28) == Fraction(13, 14),
            "b wrap threshold changed")
    require(Fraction(1, 26) < Fraction(1, 14),
            "b contraction stopped forcing speed-14 danger")

    certificates = {}
    for letters in product("ab", repeat=6):
        word = "".join(letters)
        early_a = next((j for j in range(4) if word[j] == "a"), None)
        if early_a is not None:
            certificates[word] = f"a@{early_a}->c1-danger@{early_a + 2}"
        else:
            require(word[:4] == "bbbb", "six-state dichotomy changed")
            certificates[word] = "b@0,1->speed14-danger@5"
    counts = {
        label: sum(value == label for value in certificates.values())
        for label in sorted(set(certificates.values()))
    }
    require(counts == {
        "a@0->c1-danger@2": 32,
        "a@1->c1-danger@3": 16,
        "a@2->c1-danger@4": 8,
        "a@3->c1-danger@5": 4,
        "b@0,1->speed14-danger@5": 4,
    }, "B1 length-six obstruction partition changed")

    # Strict positive control: the phase language is not empty merely at
    # one or two states.  This exact point realizes bbba with positive margin.
    witness = Fraction(12_894_291, 80_000_000)
    witness_states = []
    value = witness
    for step, target in enumerate("bbba"):
        require(literal_word(value, target),
                f"strict bbba witness failed at state {step}")
        margin = literal_margin(value, target)
        require(margin > 0, f"bbba witness hit a boundary at state {step}")
        witness_states.append((value, target, margin))
        value = B1(value)

    return clock_constants, counts, witness, tuple(witness_states)


def audit_thm2644_gate():
    # On the natural C221 phase-label projection, the two displayed lifts are
    # opposite.  Treating the labels as one nonnegative convolution weight
    # gives exactly THM-2644's inverse-pair equality boundary R=delta, not the
    # strict odd-torsor gate.  The actual affine stalk map has multiplier 13,
    # which is noninvertible modulo 221, so even this is only a sidecar audit.
    modulus = 221
    labels = (114, 107)
    require(labels[1] == (-labels[0]) % modulus,
            "C221 phase labels stopped being an inverse pair")
    mu = {label: 1 for label in labels}
    mass = sum(mu.values())
    energy = sum(value * value for value in mu.values())
    delta = mass * mass - energy
    oriented_return = sum(
        mu.get(g, 0) * mu.get((-g) % modulus, 0)
        for g in range(modulus)
    )
    identity_square_mass = mu.get(0, 0) ** 2
    require((mass, energy, delta, oriented_return, identity_square_mass)
            == (2, 2, 2, 2, 0),
            "C221 THM-2644 equality ledger changed")
    require(gcd(P, modulus) == P,
            "C221 affine multiplier unexpectedly became a group automorphism")

    single = labels[0]
    single_return = int((2 * single) % modulus == 0)
    require(single_return == 0 and (2 * single) % modulus == 7,
            "single C221 branch unexpectedly returned")

    # The half-cycle has the one nontrivial C2 parity label.  It is pure and
    # returns in two steps, but C2 is precisely THM-2644's sharp oddness
    # hostile: the branch is not the identity and A*=A.
    c2_label = 1
    c2_mass = c2_energy = c2_return = 1
    c2_identity_mass = 0
    require((-c2_label) % 2 == c2_label,
            "C2 orientation stopped coinciding with reversal")
    require((c2_mass, c2_energy, c2_return, c2_identity_mass)
            == (1, 1, 1, 0), "C2 hostile ledger changed")

    return (
        modulus, labels, mass, energy, delta, oriented_return,
        identity_square_mass, (2 * single) % modulus,
        (c2_mass, c2_energy, c2_return, c2_identity_mass),
    )


def main():
    half_delayed, c221_delayed, delayed_rows, physical_rows = (
        audit_displayed_cycles()
    )
    clock_constants, counts, witness, witness_states = (
        audit_half_semantic_language()
    )
    gate = audit_thm2644_gate()

    print("LRC14 half/C221 cycle semantic gate probe")
    print("status=VERIFIED-EXACT stopping certificate; canonical typed row")
    print(f"row={ROW}; clocks={CLOCKS}")
    print(f"half_delayed_cycle={half_delayed}; C221_delayed_cycle={c221_delayed}")
    print(f"delayed_endpoint_rows=(family,y,c1_distance,bits)={delayed_rows}")
    print("delayed_endpoint_verdict=all pass guard, five units, target-a, "
          "target-b-safe; all fail only the unshifted-c1-safe bit")
    print(f"physical_event_rows=(family,event,x,exclusive_sources,T6x,word)="
          f"{physical_rows}")
    print("physical_verdict=no event centre is an exclusive source; T6 reaches "
          "the deepest fork Q_(3,{1,2}) but not along an E3 source branch")
    print(f"B1_even_clock_constants=(k,(13^k-1)/24)={clock_constants}")
    print("prescribed_clock_verdict=B1^k=B0^k for k=2,4,6; the half phase "
          "cannot change any THM2305 endpoint span")
    print(f"B1_length6_certificate_partition={counts}")
    print("B1_factor_sparse=a@j forces c1-danger@j+2; "
          "b@j,b@j+1 force speed14-danger@j+5")
    print(f"B1_strict_word=bbba:witness={witness}:states={witness_states}")
    print("B1_scope=all length-six literal words empty and no recurrent SCC; "
          "this probe does not decide whether some length-five word exists")
    print("THM2644_C221=(modulus,labels,M,E,delta,R,J,2a)="
          f"{gate[:-1]}")
    print(f"THM2644_C2=(M,E,R,identity_mass)={gate[-1]}")
    print("THM2644_verdict=C221 is the inverse-pair equality/backtracking "
          "boundary; C2 has A*=A and is the sharp even-group hostile; neither "
          "supplies one odd-torsor same-oriented semantic transition")
    print("SCOPE: no semantic cospan, endpoint current, scalar row, row "
          "exclusion, or LRC14 conclusion")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
