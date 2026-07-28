#!/usr/bin/env python3
"""Exact root-six address-selector / fixed-X spectral-cospan probe.

The THM-2712 semantic addresses are n=13*m.  Their reduced address atlas has
304 points in Z/(13^5), arranged in the seven-clock/five-tooth form recorded
by the based-macro companion.  This probe establishes three separate facts:

1. the atlas has no nontrivial translation stabilizer and its indicator has
   no zero at any nontrivial 13-power character;
2. after projection m mod 13, the unique maximal Fourier-energy pair is
   h=+/-6, matching the frozen following atom's deep root up to orientation;
3. at the oriented physical frequency X=12*13^4, for which the lift phase is
   h=6, the canonical THM-2625 marked current remains nonzero in every one of
   its 169 target marginals.  The two-step atom macro contributes a nonzero
   coherent multiplier 3042*S_6*S_{-6}.

Item 3 is a tensor-product cospan certificate.  It does not identify the
THM-2625 word/current with the THM-2680 semantic atom, make a packet address a
THM-2350 target residue, construct chronology, or prove LRC(14).
"""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_canonical_endpoint_current_thm2625 as current


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def polynomial_value(coefficients, point):
    value = Fraction(0)
    for coefficient in coefficients:
        value = value * point + coefficient
    return value


def interval_multiply(left, right):
    products = (
        left[0] * right[0],
        left[0] * right[1],
        left[1] * right[0],
        left[1] * right[1],
    )
    return min(products), max(products)


def polynomial_interval(coefficients, interval):
    value = (Fraction(coefficients[0]), Fraction(coefficients[0]))
    for coefficient in coefficients[1:]:
        value = interval_multiply(value, interval)
        value = (value[0] + coefficient, value[1] + coefficient)
    return value


def digest(bank, keys):
    payload = ",".join(str(bank[key]) for key in keys).encode()
    return sha256(payload).hexdigest()


def build_atlas():
    p = 13
    modulus = p**5
    centres = tuple(j * modulus // 7 for j in range(7))
    offsets = (
        ((0, 8), (13, 32), (37, 43), (319, 322), (326, 326)),
        ((0, 9), (13, 33), (38, 44), (319, 323), (327, 327)),
        ((0, 9), (13, 33), (37, 44), (319, 323), (327, 327)),
        ((0, 9), (13, 33), (37, 44), (319, 323), (327, 327)),
        ((0, 9), (13, 33), (37, 44), (319, 322), (327, 327)),
        ((0, 8), (13, 33), (37, 44), (319, 322), (327, 327)),
        ((0, 8), (13, 32), (37, 44), (319, 322), (327, 327)),
    )
    atlas = tuple(
        centre + value
        for centre, row in zip(centres, offsets)
        for left, right in row
        for value in range(left, right + 1)
    )
    require(atlas == tuple(sorted(set(atlas))), "atlas order/collision")
    require(len(atlas) == 304 and len(atlas) % p == 5, "atlas census")
    require(atlas[0] == 0 and atlas[-1] < modulus, "atlas range")
    return atlas


def main():
    p = 13
    address_modulus = p**5
    atlas = build_atlas()
    atlas_set = set(atlas)

    # If A+d=A, then d is the image of 0 and hence d is in A.  This makes an
    # exact 304-candidate hostile scan sufficient.  Independently, every
    # nonzero translation orbit in a cyclic 13-group has size divisible by
    # thirteen, whereas |A|=304 is not divisible by thirteen.
    stabilizer = tuple(
        shift
        for shift in atlas
        if {(value + shift) % address_modulus for value in atlas} == atlas_set
    )
    require(stabilizer == (0,), "nontrivial atlas translation stabilizer")

    counts = tuple(Counter(value % p for value in atlas)[residue] for residue in range(p))
    require(
        counts == (23, 27, 22, 27, 22, 23, 22, 23, 23, 23, 23, 23, 23),
        "inner residue histogram",
    )
    affine_stabilizer = tuple(
        (unit, shift)
        for unit in range(1, p)
        for shift in range(p)
        if all(counts[r] == counts[(unit * r + shift) % p] for r in range(p))
    )
    require(affine_stabilizer == ((1, 0),), "inner histogram affine symmetry")

    # For h != 0 the constant baseline 23 cancels, leaving
    # D(z)=4z-z^2+4z^3-z^4-z^6.  With x=z+z^-1,
    # |D(z)|^2=g(x).
    excess = tuple(value - 23 for value in counts)
    require(excess == (0, 4, -1, 4, -1, 0, -1) + (0,) * 6, "excess word")
    f_coefficients = (1, 1, -5, -4, 6, 3, -1)
    g_coefficients = (-4, 1, 12, 14, -8, 1)
    root_intervals = (
        (Fraction(-195, 100), Fraction(-194, 100)),
        (Fraction(-150, 100), Fraction(-149, 100)),
        (Fraction(-71, 100), Fraction(-70, 100)),
        (Fraction(24, 100), Fraction(25, 100)),
        (Fraction(113, 100), Fraction(114, 100)),
        (Fraction(177, 100), Fraction(178, 100)),
    )
    require(
        all(
            polynomial_value(f_coefficients, left)
            * polynomial_value(f_coefficients, right)
            < 0
            for left, right in root_intervals
        ),
        "real-cyclotomic root isolation",
    )
    # Six disjoint sign-changing intervals for a degree-six polynomial account
    # for all roots.  They are x_h=2*cos(2*pi*h/13), ordered h=6,...,1.
    energy_intervals = tuple(
        polynomial_interval(g_coefficients, interval) for interval in root_intervals
    )
    root_six_lower = energy_intervals[0][0]
    competing_upper = max(interval[1] for interval in energy_intervals[1:])
    require(root_six_lower > 100, "root-six energy lower bound")
    require(competing_upper < 40, "competing energy upper bound")
    require(root_six_lower > competing_upper, "root-six selector gap")

    # Cyclotomic nonvanishing is all-depth and needs no floating point.  If an
    # integer indicator polynomial with 304 terms vanished at a nontrivial
    # 13-power root, divisibility by the appropriate Phi_(13^r) would force its
    # value at 1 to be divisible by Phi_(13^r)(1)=13.  But 304=5 mod 13.
    require(len(atlas) % p, "13-divisibility nonvanishing certificate")

    # Choose the physical frequency whose phase on n=13*m is the oriented
    # inner character h=6:
    #   7*X*(13m)/13^6 = 6m/13 (mod 1).
    # The marked triangle Y=X+c3 has the identical address character because
    # c3=2*13^5.
    x_frequency = 12 * p**4
    y_frequency = x_frequency + current.W[current.TB]
    require(y_frequency == 38 * p**4, "selector triangle")
    require(current.nu13(x_frequency) == current.nu13(y_frequency) == 4, "selector valuation")
    require(7 * (x_frequency // p**4) % p == 6, "left selector phase")
    require(7 * (y_frequency // p**4) % p == 6, "right selector phase")

    # One certified finite-field embedding suffices: a nonzero image proves a
    # characteristic-zero cyclotomic integer nonzero.  Primality and exact
    # root order are checked at import by the THM-2625 companion.
    mods = (current.MODS[0],)
    prime, root = mods[0]
    zeta_13 = pow(root, current.NN // p, prime)
    require(pow(zeta_13, p, prime) == 1 and zeta_13 != 1, "13th root order")

    q_intervals = current.build_set(current.PAT_QA, current.ZERO_ELL)
    q_starts = [left for left, _ in q_intervals]
    tabs = current.make_tabs(q_intervals, x_frequency, mods)
    p_bank = {}
    q_bank = {}
    for address, ell in current.REPS.items():
        e_intervals = current.build_set(current.PAT_E, ell)
        left_values, _ = current.x_sweep(
            e_intervals,
            q_intervals,
            q_starts,
            x_frequency,
            mods,
            tabs,
        )
        right_values = current.endpoint_sum(e_intervals, -y_frequency, mods)
        deepest_phase = pow(zeta_13, current.M_DEEP * ell[current.TB] % p, prime)
        p_bank[address] = deepest_phase * left_values[0] % prime
        q_bank[address] = right_values[0] % prime

    # Separate quotient-gauge descent remains lawful at the new triangle.
    for address in ((0, 0), (1, 0), (0, 1), (3, 7)):
        ell = current.REPS[address]
        shifted = tuple((ell[i] + current.WMOD[i]) % p for i in range(9))
        e_intervals = current.build_set(current.PAT_E, shifted)
        left_values, _ = current.x_sweep(
            e_intervals,
            q_intervals,
            q_starts,
            x_frequency,
            mods,
            tabs,
        )
        right_values = current.endpoint_sum(e_intervals, -y_frequency, mods)
        deepest_phase = pow(
            zeta_13, current.M_DEEP * shifted[current.TB] % p, prime
        )
        require(
            deepest_phase * left_values[0] % prime == p_bank[address],
            "left representative descent",
        )
        require(right_values[0] % prime == q_bank[address], "right representative descent")

    endpoint_keys = tuple(current.REPS)
    powers = tuple(pow(zeta_13, exponent, prime) for exponent in range(p))
    left = {}
    right = {}
    target = {}
    for point in endpoint_keys:
        left[point] = sum(
            p_bank[address]
            * powers[
                -(
                    current.TAU[address][0] * point[0]
                    + current.TAU[address][1] * point[1]
                )
                % p
            ]
            for address in endpoint_keys
        ) % prime
        right[point] = sum(
            q_bank[address]
            * powers[
                (
                    current.TAU[address][0] * point[0]
                    + current.TAU[address][1] * point[1]
                )
                % p
            ]
            for address in endpoint_keys
        ) % prime
    for target_point in endpoint_keys:
        target[target_point] = sum(
            p_bank[address]
            * q_bank[address]
            * powers[
                -(
                    current.TAU[address][0] * target_point[0]
                    + current.TAU[address][1] * target_point[1]
                )
                % p
            ]
            for address in endpoint_keys
        ) % prime

    supports = tuple(
        sum(value != 0 for value in bank.values())
        for bank in (p_bank, q_bank, left, right, target)
    )
    require(supports == (169, 169, 169, 169, 169), "selector-current support hole")

    s_plus = sum(powers[6 * (value % p) % p] for value in atlas) % prime
    s_minus = sum(powers[-6 * (value % p) % p] for value in atlas) % prime
    require(s_plus and s_minus, "root-six atlas transform vanished")
    macro_multiplier = 3042 * s_plus * s_minus % prime
    require(macro_multiplier, "macro cospan multiplier vanished")
    require(
        all(value * macro_multiplier % prime for value in target.values()),
        "a target marginal died after address-cospan tensoring",
    )

    digests = tuple(
        digest(bank, endpoint_keys) for bank in (p_bank, q_bank, left, right, target)
    )
    expected_digests = (
        "6242def39eb7c33fb93520448792c695bade375e421c48cb5a647fca4ea547f0",
        "8cbcffdd371b4dbe559babf9a873e6b0be704dbe9e49975140810f3579f3bc53",
        "26d4666490a9eb078ecad81a1fbdcbd89b7641eeb09429dae9a1d981ab3b1d8b",
        "829afc83f04df7f1597d785790e56d504c6ddfc90ffbb090c1a26ff49602bb87",
        "ff2bbf7c4dc4a89ae7019aa82632b00f253bbfbbc334ce1b8da39c928a299c45",
    )
    require(digests == expected_digests, "selector-current digest drift")

    print("LRC14 root-six atom-atlas / fixed-X spectral-cospan probe")
    print(f"atlas: |A|={len(atlas)}=5 mod13; translation stabilizer={stabilizer}")
    print(f"inner residue counts={counts}; AGL(1,13) stabilizer={affine_stabilizer}")
    print("all 13-power address characters survive: Phi_(13^r)(1)=13 does not divide 304")
    print(
        "inner Fourier energy: unique maximal pair {+/-6}; "
        f"E_6>{root_six_lower}, every competitor <{competing_upper}"
    )
    print(f"selector triangle: (X,m,Y)=({x_frequency},1,{y_frequency}); address phase h=6 on both legs")
    print(f"exact field: p={prime}; support(P,Q,L*,R*,A)={supports}")
    print(f"root-six atlas images: S_6={s_plus}; S_-6={s_minus}")
    print(f"3042-path coherent multiplier={macro_multiplier} != 0")
    print(f"sha256(P,Q,L,R,A)={digests}")
    print("scope: nonzero tensor cospan; no semantic-current intertwiner, target/address identification, chronology, row exclusion, or LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
