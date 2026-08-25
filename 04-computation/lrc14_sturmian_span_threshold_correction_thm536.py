#!/usr/bin/env python3
"""Two-engine exact correction audit for THM-536 (MISTAKE-489).

The original x-wall engine and the theta-space mechanical-word engine share
no breakpoint construction.  Exact adjacent values certify every finite cap
crossing; an all-m failure interval proves the extended threshold at cap one.
"""

from __future__ import annotations

from fractions import Fraction as F
import sys


sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def direct_x_measure(m: int) -> F:
    breakpoints = {F(0), F(1)}
    for speed in range(1, m):
        breakpoints.update(F(j, 7 * speed) for j in range(7 * speed + 1))
    ordered = sorted(breakpoints)
    total = F(0)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        residues = {int(7 * speed * midpoint) % 7 for speed in range(m)}
        if len(residues) == 7:
            total += right - left
    return total


def mechanical_theta_measure(m: int) -> F:
    breakpoints = {F(0), F(7)}
    for speed in range(1, m):
        breakpoints.update(F(j, speed) for j in range(7 * speed + 1))
    ordered = sorted(value for value in breakpoints if 0 <= value <= 7)
    total = F(0)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        residues = {int(speed * midpoint) % 7 for speed in range(m)}
        if len(residues) == 7:
            total += right - left
    return total / 7


def main() -> None:
    caps = {
        8: F(2243, 5880),
        9: F(1979, 4004),
        10: F(55, 91),
        11: F(66, 91),
        12: F(6, 7),
        13: F(1),
    }
    finite_thresholds = {8: 7, 9: 8, 10: 10, 11: 13, 12: 26}
    audit_m = sorted(
        {22, 23}
        | {threshold + 1 for threshold in finite_thresholds.values()}
        | {threshold + 2 for threshold in finite_thresholds.values()}
    )
    values = {}
    for m in audit_m:
        direct = direct_x_measure(m)
        mechanical = mechanical_theta_measure(m)
        require(direct == mechanical, f"two-engine mismatch m={m}")
        values[m] = direct

    # AP_m subset AP_(m+1) proves global monotonicity.  Adjacent exact values
    # therefore identify each finite supremum without any tail interpolation.
    for k, threshold in finite_thresholds.items():
        last = values[threshold + 1]
        first_failure = values[threshold + 2]
        require(last <= caps[k], f"last certified k={k}")
        require(first_failure > caps[k], f"first failure k={k}")

    require(values[22] == F(6155, 7497), "canonical minimal witness")
    require(values[22] < caps[12], "canonical N*=20 refuted")
    require(values[23] == F(4333193, 5222910), "transcript minimal witness")
    require(values[23] < caps[12], "transcript N*=21 refuted")

    # For m>=2, theta in [0,1/(m-1)) makes every floor(e theta) zero.
    # Its x-measure is 1/(7(m-1)), so a(m)<1 at every finite m.
    for m, value in values.items():
        if m >= 2:
            require(value <= 1 - F(1, 7 * (m - 1)), f"failure interval m={m}")

    print("THM536_STURMIAN_SPAN_THRESHOLD_CORRECTION_EXACT")
    print("definition=a(m)=meas{x:{floor(7ex)_mod_7:0<=e<m}=Z/7}")
    print("global_monotonicity=AP_m_subset_AP_(m+1)")
    for k, threshold in finite_thresholds.items():
        last = values[threshold + 1]
        first_failure = values[threshold + 2]
        print(
            f"k={k};Nstar={threshold};a(Nstar+1)={last};"
            f"cap_minus_last={caps[k]-last};a(Nstar+2)={first_failure};"
            f"first_minus_cap={first_failure-caps[k]}"
        )
    print("k=13;Nstar=infinity;a(m)<=1-1/(7(m-1))<1_for_every_finite_m>=2")
    print(f"canonical_minimal_witness=k12,N21,a22={values[22]}")
    print(f"transcript_minimal_witness=k12,N22,a23={values[23]}")
    print(f"two_engine_values={tuple((m, values[m]) for m in audit_m)}")
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
