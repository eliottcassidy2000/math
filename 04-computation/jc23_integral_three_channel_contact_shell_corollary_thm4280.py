#!/usr/bin/env python3
"""Finite exact controls for THM-4280's contact-shell corollary."""

from itertools import product


def norm(a: int, b: int) -> int:
    """Eisenstein norm of a+b*omega."""
    return a * a - a * b + b * b


def require(condition: bool, message: str) -> None:
    """Optimization-safe audit gate."""
    if not condition:
        raise RuntimeError(message)


residue_norms_mod2 = {
    (a, b): norm(a, b) % 2 for a, b in product(range(2), repeat=2)
}
require(
    residue_norms_mod2
    == {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 1},
    "mod-two norm table changed",
)

# Hostile bounded control of the proved valuation statement.
for a, b in product(range(-128, 129), repeat=2):
    n = norm(a, b)
    if n and n % 2 == 0:
        valuation = 0
        while n % 2 == 0:
            n //= 2
            valuation += 1
        require(valuation % 2 == 0, "odd two-adic valuation found")

examples = [
    (8, (1, 0), (1, 0)),
    (20, (1, 0), (2, 0)),
    (24, (1, -1), (1, -1)),
    (32, (2, 0), (2, 0)),
]
for degree, a, e in examples:
    require(degree == 4 * (norm(*a) + norm(*e)), "degree witness changed")
    require(a != (0, 0), "ramification-two witness lost")

for s in range(1, 129):
    target = 2 * s * s
    # This bounded search is a hostile control, not the proof of
    # nonrepresentation; the theorem proves it from the mod-two descent.
    bound = 2 * s + 2
    require(
        all(
            norm(a, b) != target
            for a, b in product(range(-bound, bound + 1), repeat=2)
        ),
        "2s^2 unexpectedly represented",
    )
    require(
        8 * s * s == 4 * (norm(s, 0) + norm(s, 0)),
        "visible degree-8s^2 witness changed",
    )

print("THM4280_PER_CONSUMER_CONTACT_EXACT_AUDIT_V1")
print("NORM_MOD2", residue_norms_mod2)
print("VALUATION_HOSTILE_BOX", -128, 128, "PASS")
print("SHARP_THREE_EXAMPLES", " ".join(str(d) for d, _, _ in examples))
print("DEGREE_8S2_HOSTILE_RANGE", "1..128", "PASS")
print("VERDICT PASS FINITE_CONTROLS; SYMBOLIC_PROOF_IN_THM4280_SECTION_6A")
