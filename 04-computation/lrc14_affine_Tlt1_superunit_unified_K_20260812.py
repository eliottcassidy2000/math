#!/usr/bin/env python3
"""Exact arithmetic audit for the repaired unified T<1 tail envelope.

This repairs the omitted chamber T<1<T_infinity.  The companion note derives
the geometric inequalities; this verifier checks all constants, finite d
monotonicity, the two-turn ceiling, and the integer tail starts.
"""
from fractions import Fraction as F


DMAX = F(186_636_088_362, 11_773_143_757_375)
TARGET = DMAX / 5
JMIN = F(709, 48_048)
ZETA = F(3167, 3168)
OMEGA = F(3155, 3168)
KF = F(6864, 22085)
TINF_CAP = F(679, 672)
DRIFT_CAP = F(97, 96)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def tau(d):
    delta = F(d - 1) + F(1, 12)
    return 1 + delta / (264 * ZETA)


def kstar(d):
    """Safe unified nonzero-step error constant.

    Both phase-drift coefficients carry the superunit cap.  The Peano term
    is instead bounded directly from the finite hypothesis T<1; see the
    companion note.  Its displayed value is deliberately conservative.
    """
    delta = F(d - 1) + F(1, 12)
    t = tau(d)
    kp = d * KF + DRIFT_CAP * (F(d * (d - 1), 2 * ZETA) + F(d, ZETA))
    return (
        t * kp / d
        + 2 * t * delta / (7 * ZETA * OMEGA)
        + t * d / (42 * OMEGA)
        + F(d, 14 * OMEGA)
        + F(d * d, 264 * ZETA * OMEGA)
    )


def kzero(d):
    kp = d * KF + F(d * (d - 1), 2 * ZETA) + F(d, ZETA)
    return kp / d + F(d - 1) * F(2, 7 * OMEGA)


def main():
    require(ZETA == 1 - F(1, 3168), ZETA)
    require(OMEGA == 1 - F(13, 3168), OMEGA)

    # For p>=679 and 0<=r<=d-1, p/(p-r) is maximal at p=679,
    # r=d-1, then at d=8.
    caps = tuple(F(679, 680 - d) for d in range(1, 9))
    require(all(a < b for a, b in zip(caps, caps[1:])), caps)
    require(caps[-1] == TINF_CAP, caps[-1])

    taus = tuple(tau(d) for d in range(1, 9))
    require(all(a < b for a, b in zip(taus, taus[1:])), taus)
    require(taus[-1] == F(3252, 3167), taus[-1])
    comparison_caps = tuple(caps[d - 1] * taus[d - 1] for d in range(1, 9))
    require(all(x < 2 for x in comparison_caps), comparison_caps)
    require(comparison_caps[-1] == F(26287, 25336), comparison_caps[-1])

    constants = tuple(kstar(d) for d in range(1, 9))
    require(all(a < b for a, b in zip(constants, constants[1:])), constants)
    kmax = constants[-1]
    require(kmax == F(1_792_138_785_426, 221_510_098_565), kmax)
    threshold = kmax / (JMIN - TARGET)
    require(
        threshold
        == F(
            144_679_594_655_915_839_062_972_000,
            207_178_827_309_738_451_742_041,
        ),
        threshold,
    )
    require(698 < threshold < 699, threshold)
    gap = JMIN - kmax / 699 - TARGET
    require(gap > 0, gap)

    zero_constants = tuple(kzero(d) for d in range(1, 9))
    require(all(a < b for a, b in zip(zero_constants, zero_constants[1:])), zero_constants)
    zero_max = zero_constants[-1]
    require(zero_max == F(477_044_832, 69_943_195), zero_max)
    zero_threshold = zero_max / (JMIN - TARGET)
    require(
        zero_threshold
        == F(38_511_890_645_820_381_504_000, 65_418_006_728_682_807_623),
        zero_threshold,
    )
    require(588 < zero_threshold < 589, zero_threshold)
    zero_gap = JMIN - zero_max / 589 - TARGET
    require(zero_gap > 0, zero_gap)

    print("LRC14 UNIFIED T<1 CONTINUUM-ERROR EXACT AUDIT")
    print("Tinf_caps_by_d=", caps, sep="")
    print("Tinf_global_cap=", TINF_CAP, sep="")
    print("tau_by_d=", taus, sep="")
    print("comparison_path_caps_by_d=", comparison_caps, sep="")
    print("comparison_path_max=", comparison_caps[-1], ";less_than_2=True", sep="")
    print("Kstar_by_d=", constants, sep="")
    print("Kstar_max=", kmax, ";threshold=", threshold, ";integer_tail_start=699", sep="")
    print("p699_margin=", gap, sep="")
    print("A0_K_by_d=", zero_constants, sep="")
    print(
        "A0_Kmax=", zero_max, ";A0_threshold=", zero_threshold,
        ";A0_integer_tail_start=589;A0_p589_margin=", zero_gap, sep=""
    )
    print("status=exact arithmetic audit; geometric inequalities are in companion note")


if __name__ == "__main__":
    main()
