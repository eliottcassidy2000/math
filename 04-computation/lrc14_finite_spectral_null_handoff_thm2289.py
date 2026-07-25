#!/usr/bin/env python3
"""Exact arithmetic audit for THM-2289."""

from fractions import Fraction


THM = "THM-2289"
ALPHAS = {
    "strict": Fraction(15_041_431, 593_783_190),
    "repeated": Fraction(5_229_541, 593_783_190),
}
BETA = Fraction(2_593, 90_090)
SOURCE_COORDS = 9
TARGET_COORDS = 7


def jackson_ledger_margin(alpha: Fraction, n: int) -> Fraction:
    epsilon = Fraction(3, 2 * n)
    return (
        (alpha - SOURCE_COORDS * epsilon)
        * (BETA - TARGET_COORDS * epsilon)
        - (SOURCE_COORDS + TARGET_COORDS) * epsilon
    )


def jackson_ledger_closes(alpha: Fraction, n: int) -> bool:
    epsilon = Fraction(3, 2 * n)
    return (
        alpha - SOURCE_COORDS * epsilon > 0
        and BETA - TARGET_COORDS * epsilon > 0
        and jackson_ledger_margin(alpha, n) > 0
    )


def first_positive_n(alpha: Fraction) -> int:
    lo, hi = 2, 2
    while not jackson_ledger_closes(alpha, hi):
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if jackson_ledger_closes(alpha, mid):
            hi = mid
        else:
            lo = mid
    return hi


EXPECTED = {
    "strict": {
        "n": 33_810,
        "fail": Fraction(
            -91_689_987_875_828, 15_286_538_167_789_662_548_775
        ),
        "passed": Fraction(
            12_331_416_859, 792_352_055_420_125_200
        ),
        "source_factor": Fraction(4_766_997_229, 191_198_187_180),
        "target_factor": Fraction(117_991, 4_144_140),
        "h": 67_618,
        "h_over_7": 9_659,
        "q": 115_919,
        "m": int(
            "531427494109850809274382490322199270940210549004984817720024"
            "492918545763422455682"
        ),
    },
    "repeated": {
        "n": 96_570,
        "fail": Fraction(
            -21_997_372_328, 8_518_227_246_964_664_775
        ),
        "passed": Fraction(
            10_450_627_633, 246_356_440_619_713_023_600
        ),
        "source_factor": Fraction(11_044_460_029, 1_274_258_725_740),
        "target_factor": Fraction(5_543_557, 193_333_140),
        "h": 193_138,
        "h_over_7": 27_591,
        "q": 331_095,
        "m": int(
            "104282626719224343930980394141063003241014576694233744028445"
            "39686260225238346754941835938"
        ),
    },
}

RESULTS = {}
for branch, alpha in ALPHAS.items():
    n = first_positive_n(alpha)
    h = 2 * n - 2
    q = 2 * h + 1 - 2 * (h // 7)
    m = ((q**SOURCE_COORDS - 1) * (q**TARGET_COORDS - 1) // 2)
    epsilon = Fraction(3, 2 * n)
    result = {
        "n": n,
        "fail": jackson_ledger_margin(alpha, n - 1),
        "passed": jackson_ledger_margin(alpha, n),
        "source_factor": alpha - SOURCE_COORDS * epsilon,
        "target_factor": BETA - TARGET_COORDS * epsilon,
        "h": h,
        "h_over_7": h // 7,
        "q": q,
        "m": m,
    }
    assert result == EXPECTED[branch]
    RESULTS[branch] = result

assert RESULTS["strict"]["m"] < RESULTS["repeated"]["m"]

STRICT_PROFILES = [
    (1, b, c)
    for c in range(5, 20)
    for b in range(2, c)
]
REPEATED_PROFILES = [(1, 1, c) for c in range(5, 20)]
assert len(STRICT_PROFILES) == 150
assert len(REPEATED_PROFILES) == 15
assert set(STRICT_PROFILES).isdisjoint(REPEATED_PROFILES)

# Arithmetic hostile controls: bounded coefficients can certify an
# arbitrarily delayed scale, so the all-null-time injection is essential.
for k in (1, 2, 5, 19, 64, 127):
    h = 1
    q = 13**k + 1
    source_partial_sum = -h + q
    target_partial_sum = -h
    assert q % 13 == 1
    assert source_partial_sum != 0
    assert target_partial_sum != 0
    assert source_partial_sum + 13**k * target_partial_sum == 0

print(f"theorem={THM}")
print(f"beta={BETA}")
print(f"source_coordinates={SOURCE_COORDS}")
print(f"target_coordinates={TARGET_COORDS}")
print(f"strict_profiles={len(STRICT_PROFILES)}")
print(f"repeated_profiles={len(REPEATED_PROFILES)}")
for branch in ("strict", "repeated"):
    result = RESULTS[branch]
    print(f"{branch}_alpha={ALPHAS[branch]}")
    print(f"{branch}_N_fail={result['n'] - 1}")
    print(f"{branch}_margin_fail={result['fail']}")
    print(f"{branch}_N_first={result['n']}")
    print(f"{branch}_source_factor={result['source_factor']}")
    print(f"{branch}_target_factor={result['target_factor']}")
    print(f"{branch}_margin_pass={result['passed']}")
    print(f"{branch}_degree={result['h']}")
    print(f"{branch}_floor_degree_over_7={result['h_over_7']}")
    print(f"{branch}_coefficient_alphabet={result['q']}")
    print(f"{branch}_signed_certificate_classes={result['m']}")
    print(f"{branch}_window_cardinality={result['m'] + 1}")
print("late_certificate_controls=1,2,5,19,64,127")
print("status=PASS")
