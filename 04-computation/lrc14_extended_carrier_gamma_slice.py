#!/usr/bin/env python3
"""Exact THM-2763 one-fibre probe of the carrier-gauge character gamma."""

from hashlib import sha256
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
HELPER = HERE / "lrc14_extended_carrier_endpoint_lib.py"
HELPER_SHA256 = "4b3f9f195b1634e1e84a1bc8bccb878a1c8c44aec13f24d197f92547c9e36c57"
if sha256(HELPER.read_bytes().replace(b"\r\n", b"\n")).hexdigest() != HELPER_SHA256:
    raise RuntimeError("extended-carrier helper changed")
import lrc14_extended_carrier_endpoint_lib as base


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def dft(values, zeta, prime):
    return tuple(
        sum(
            value * pow(zeta, (-frequency * gamma) % base.P, prime)
            for gamma, value in enumerate(values)
        ) % prime
        for frequency in range(base.P)
    )


def main():
    source, target, terminal = base.build_carriers()["base"]
    prime, root = base.endpoint.MODS[0]
    zeta = pow(root, base.NRED // base.P, prime)
    ell = base.REPS[(0, 0)]
    present = tuple(base.endpoint.build_set(base.PAT_E3, ell))
    q_starts = [left for left, _right in terminal]
    tabs = base.endpoint.make_tabs(
        terminal, base.X0, ((prime, root),)
    )

    p_values = []
    q_values = []
    overlap_values = []
    source_present_masses = []
    target_present_masses = []
    for gamma in range(base.P):
        shift = -gamma * (base.T // base.P)
        source_gamma = base.marked.clutch.shift_union(source, shift)
        target_gamma = base.marked.clutch.shift_union(target, shift)
        left_carrier = base.intersect_sorted(present, source_gamma)
        right_carrier = base.intersect_sorted(present, target_gamma)
        source_present_masses.append(base.interval_mass(left_carrier))
        target_present_masses.append(base.interval_mass(right_carrier))
        p_value, overlap = base.endpoint.x_sweep(
            left_carrier, terminal, q_starts, base.X0,
            ((prime, root),), tabs,
        )
        q_value = base.endpoint.endpoint_sum(
            right_carrier, -base.Y0, ((prime, root),)
        )
        p_values.append(p_value[0])
        q_values.append(q_value[0])
        overlap_values.append(overlap)
    p_values = tuple(p_values)
    q_values = tuple(q_values)
    source_present_masses = tuple(source_present_masses)
    target_present_masses = tuple(target_present_masses)
    h_values = tuple(
        p_values[index] * q_values[index] % prime
        for index in range(base.P)
    )
    p_hat = dft(p_values, zeta, prime)
    q_hat = dft(q_values, zeta, prime)
    h_hat = dft(h_values, zeta, prime)

    expected_source_masses = (
        6_320_326_320, 0, 4_072_511_520, 0, 0, 0, 0,
        0, 0, 0, 6_320_326_320, 0, 0,
    )
    expected_target_masses = (
        6_320_326_320, 0, 4_046_066_640, 0, 0, 0, 0,
        0, 0, 0, 6_320_326_320, 0, 0,
    )
    require(
        source_present_masses == expected_source_masses
        and target_present_masses == expected_target_masses,
        "physical gamma-support masses changed",
    )
    expected_overlap_values = (
        374_977_060_157_700, 0, 241_617_017_842_200, 0, 0, 0, 0,
        0, 0, 0, 374_977_060_157_700, 0, 0,
    )
    require(tuple(overlap_values) == expected_overlap_values,
            "gamma overlap numerators changed")
    expected_raw_support = (0, 2, 10)
    require(
        tuple(index for index, value in enumerate(p_values) if value)
        == expected_raw_support
        and tuple(index for index, value in enumerate(q_values) if value)
        == expected_raw_support
        and tuple(index for index, value in enumerate(h_values) if value)
        == expected_raw_support,
        "raw dual-twist support changed",
    )
    require(q_values[-1] == 0 and q_values[0] != 0,
            "gauge hostile gamma=-1/0 mismatch")
    require(all(p_hat) and all(q_hat) and all(h_hat),
            "a primal carrier-imbalance residue vanished")

    print("THM-2763 extended carrier-bank diagonal gamma slice")
    print(f"helper_sha256={HELPER_SHA256}")
    print("slice=ell=(0,0), gamma=0,...,12; P uses source carrier/full terminal, Q uses target bare endpoint")
    print(f"P_gamma={p_values}")
    print(f"Q_gamma={q_values}")
    print(f"H_gamma=P_gamma*Q_gamma={h_values}")
    print(f"source_present_masses={source_present_masses}")
    print(f"target_present_masses={target_present_masses}")
    print(f"overlap_numerators={tuple(overlap_values)}")
    print(f"support(P,Q,H)={(sum(bool(v) for v in p_values), sum(bool(v) for v in q_values), sum(bool(v) for v in h_values))}")
    print(f"support_gamma_DFT(P,Q,H)={(sum(bool(v) for v in p_hat), sum(bool(v) for v in q_hat), sum(bool(v) for v in h_hat))}")
    print(f"P_gamma_DFT={p_hat}")
    print(f"Q_gamma_DFT={q_hat}")
    print(f"H_gamma_DFT={h_hat}")
    print("CONCLUSION: the extra carrier coordinate is active already on one quotient fibre.  The raw diagonal dual-twist bank is supported exactly at gamma={0,2,10}, while its gamma DFT is nonzero in every primal carrier-imbalance residue.  For two-sided insertion this is the (a,b)=(gamma,-gamma) slice of the full 28561-twist extension; it is not an exact-address or LRC conclusion.")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
