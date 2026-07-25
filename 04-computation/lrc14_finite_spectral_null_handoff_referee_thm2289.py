#!/usr/bin/env python3
"""Independent cleared-integer referee for THM-2289."""

from math import gcd


B_NUM, B_DEN = 2_593, 90_090
R, S = 9, 7


def cleared_margin(a_num: int, a_den: int, n: int) -> tuple[int, int]:
    """Return numerator/denominator of the ledger after direct clearing."""
    # (A-3R/(2n))(B-3S/(2n))-3(R+S)/(2n)
    x_num = 2 * n * a_num - 3 * R * a_den
    y_num = 2 * n * B_NUM - 3 * S * B_DEN
    den = 4 * n * n * a_den * B_DEN
    num = x_num * y_num - 6 * n * (R + S) * a_den * B_DEN
    g = gcd(abs(num), den)
    return num // g, den // g


CASES = {
    "strict": {
        "alpha": (15_041_431, 593_783_190),
        "n": 33_810,
        "fail": (
            -91_689_987_875_828,
            15_286_538_167_789_662_548_775,
        ),
        "passed": (12_331_416_859, 792_352_055_420_125_200),
        "h": 67_618,
        "q": 115_919,
        "m": int(
            "531427494109850809274382490322199270940210549004984817720024"
            "492918545763422455682"
        ),
    },
    "repeated": {
        "alpha": (5_229_541, 593_783_190),
        "n": 96_570,
        "fail": (-21_997_372_328, 8_518_227_246_964_664_775),
        "passed": (10_450_627_633, 246_356_440_619_713_023_600),
        "h": 193_138,
        "q": 331_095,
        "m": int(
            "104282626719224343930980394141063003241014576694233744028445"
            "39686260225238346754941835938"
        ),
    },
}

AUDIT = {}
for branch, case in CASES.items():
    a_num, a_den = case["alpha"]
    fail = cleared_margin(a_num, a_den, case["n"] - 1)
    passed = cleared_margin(a_num, a_den, case["n"])
    assert fail == case["fail"]
    assert passed == case["passed"]
    assert fail[0] < 0 < passed[0]

    h = case["h"]
    xi = [
        index
        for index in range(-h, h + 1)
        if index == 0 or index % 7 != 0
    ]
    assert len(xi) == case["q"]
    assert xi[0] == -h
    assert xi[-1] == h
    assert all(index == 0 or index % 7 for index in xi)

    q = len(xi)
    ordered_nonzero_pairs = (q**9 - 1) * (q**7 - 1)
    assert ordered_nonzero_pairs % 2 == 0
    m = ordered_nonzero_pairs // 2
    assert m == case["m"]
    AUDIT[branch] = (fail, passed, len(xi), m)

# Directly referee the algebraic injection: a repeated certificate at two
# sample distinct times forces its target partial sum to vanish.
for k, ell in ((0, 1), (2, 9), (19, 64)):
    factor = 13**ell - 13**k
    assert factor != 0
    # Over Z, factor*B=0 implies B=0.
    for b_partial in (-17, -1, 1, 23):
        assert factor * b_partial != 0

print("theorem=THM-2289")
print("method=cleared-integer-margin-plus-literal-alphabet")
for branch in ("strict", "repeated"):
    fail, passed, q, m = AUDIT[branch]
    print(f"{branch}_fail_margin={fail[0]}/{fail[1]}")
    print(f"{branch}_pass_margin={passed[0]}/{passed[1]}")
    print(f"{branch}_literal_alphabet_size={q}")
    print(f"{branch}_signed_certificate_classes={m}")
print("certificate_reuse=impossible_at_distinct_times_when_B_nonzero")
print("status=PASS")
