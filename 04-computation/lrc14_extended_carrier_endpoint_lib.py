#!/usr/bin/env python3
"""Exact helper library for the THM-2749 -> THM-2334/2625 carrier gate.

This exact THM-2763 support library inserts the rail-eight, clock-one
THM-2749 common carrier into the deepest-owner R=13^6 marked
triangle

    (X,m,Y)=(12*13^4,1,38*13^4)

and computes the 169 separately allocated P_xi/Q_xi endpoint banks before
forming either the target difference or the determinant/Radon quotient.

Two target constructions are compared:

* carrier_only: translate only the THM-2749 carrier, leaving the canonical
  THM-2334 endpoint factors fixed;
* full_rechart: translate the carrier and every present endpoint factor.
  The delayed word is unchanged because 13^6*(7/13^6)=7.  The separate
  deepest Fourier leg is translated as well.

All decisions are made in one exact finite-field specialization of the
cyclotomic integer bank.  A nonzero image certifies characteristic-zero
nonvanishing; failure of a proposed covariance identity in this image
certifies failure of that identity in characteristic zero.
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

PINNED_DEPENDENCIES = {
    COMP / "lrc14_canonical_endpoint_current_thm2625.py":
        "eb89eb4753f95b00ba902e1ff7c691326da62ae52501817180b81c13ed6bc62f",
    COMP / "lrc14_fully_marked_root_zero_clutch_thm2749.py":
        "93b46b2701db8f72d00fa2ae131f9d9afd3200f32998959af3bb2e1fa2f56841",
}
for dependency, expected_hash in PINNED_DEPENDENCIES.items():
    actual_hash = sha256(
        dependency.read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    if actual_hash != expected_hash:
        raise RuntimeError(
            f"pinned dependency changed: {dependency.name}: {actual_hash}"
        )

import lrc14_canonical_endpoint_current_thm2625 as endpoint
import lrc14_fully_marked_root_zero_clutch_thm2749 as marked


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


P = 13
R = P**6
T = endpoint.T_DEN
SHIFT = marked.SHIFT
X = 12 * P**4
Y = 38 * P**4
X0 = X // P**4
Y0 = Y // P**4
DEEP = 8
OWNER = 8
TA = 6
TB = 7
CLOCK = 1
RAIL = 8
ZERO = (0,) * 9

# The available THM-2625 field has exact order 169*T.  Since X and Y are
# divisible by 13^4, it is exactly the effective cyclotomic order for the
# R*T endpoint phases at this triangle.
NRED = endpoint.NN
require(R * T // P**4 == NRED, "effective cyclotomic order mismatch")
require((X, Y, Y - X) == (342732, 1085318, endpoint.W[DEEP]),
        "fixed triangle mismatch")
require(R * SHIFT == 7 * T, "chart shift mismatch")

# Process-local specialization of the generic endpoint sweep to dilation R.
# Frequencies are divided by 13^4 at the same time, so the old exact roots of
# order NRED remain the correct roots.
endpoint.RDIL = R
endpoint.NN = NRED


def intersect_sorted(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            out.append((a, b))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def support_of(weighted):
    return marked.clutch.merge_intervals((a, b) for a, b, _weight in weighted)


def interval_mass(intervals):
    return sum(b - a for a, b in intervals)


PAT_E3 = {
    0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
    6: "out", 7: "out", 8: "in",
}
PAT_Q12 = {
    0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
    6: "in", 7: "in", 8: "out",
}


def generic_reps():
    """Deep-owner version of the THM-2309 two-target annihilator."""
    u0 = 5
    pivots = (OWNER, 0, 1, 2, 3, 4)
    rows = []
    for pivot in pivots:
        row = [0] * 9
        row[u0] = endpoint.W[pivot] % P
        row[pivot] = (-endpoint.W[u0]) % P
        rows.append(row)
    row_a = rows[pivots.index(1)]
    row_a[u0] = (row_a[u0] + endpoint.W[TA]) % P
    row_a[TA] = (row_a[TA] - endpoint.W[u0]) % P
    row_b = rows[pivots.index(2)]
    row_b[u0] = (row_b[u0] + endpoint.W[TB]) % P
    row_b[TB] = (row_b[TB] - endpoint.W[u0]) % P

    wmod = tuple(value % P for value in endpoint.W)
    v1 = tuple(-1 if i == 1 else 1 if i == TA else 0 for i in range(9))
    v2 = tuple(-1 if i == 2 else 1 if i == TB else 0 for i in range(9))
    v1 = tuple(value % P for value in v1)
    v2 = tuple(value % P for value in v2)
    require(endpoint.gf13_rank(rows) == 6, "deep-owner relation rank")
    for vector in (wmod, v1, v2):
        require(all(sum(a * b for a, b in zip(row, vector)) % P == 0
                    for row in rows), "annihilator vector mismatch")
    require(endpoint.gf13_rank((wmod, v1, v2)) == 3,
            "deep-owner quotient basis")
    reps = {
        (alpha, beta): tuple(
            (alpha * v1[i] + beta * v2[i]) % P for i in range(9)
        )
        for alpha in range(P) for beta in range(P)
    }
    require(all((ell[TA], ell[TB]) == key for key, ell in reps.items()),
            "deep-owner target coordinates")
    return reps, tuple(rows), wmod


REPS, RELATION_ROWS, WMOD = generic_reps()
KEYS = tuple(REPS)

_PRESENT_CACHE = None
_SHIFTED_PRESENT_CACHE = None


def present_cache(*, shifted=False):
    """Build the 169 expensive Boolean present sets only once per run."""
    global _PRESENT_CACHE, _SHIFTED_PRESENT_CACHE
    if _PRESENT_CACHE is None:
        _PRESENT_CACHE = {
            address: tuple(endpoint.build_set(PAT_E3, ell))
            for address, ell in REPS.items()
        }
    if not shifted:
        return _PRESENT_CACHE
    if _SHIFTED_PRESENT_CACHE is None:
        _SHIFTED_PRESENT_CACHE = {
            address: marked.clutch.shift_union(intervals, SHIFT)
            for address, intervals in _PRESENT_CACHE.items()
        }
    return _SHIFTED_PRESENT_CACHE


def build_carriers():
    module, _pair_prefixes, _, _, rails, present, _starts = (
        marked.clutch.relative.lift.m.core.build_carrier_data()
    )
    source_e3 = marked.semantic.exclusive_source(module, 3)
    terminal = marked.semantic.deepest_fork(module)
    prefixes = marked.build_semantic_prefixes(module, terminal)
    sections = marked.semantic_sections(module, source_e3, *marked.LABEL)
    source_vector, target_vector, details = marked.fully_marked_vectors(
        module, rails, present, prefixes, sections, RAIL
    )
    require(source_vector == target_vector and source_vector[CLOCK],
            "chosen THM-2749 fibre is not live")
    source = support_of(details[CLOCK][0])
    target = support_of(details[CLOCK][1])
    source_carry = support_of(details[CLOCK][2])
    target_carry = support_of(details[CLOCK][3])
    require(marked.clutch.shift_union(source, SHIFT) == target,
            "carrier translation")
    require(marked.clutch.shift_union(source_carry, SHIFT) == target_carry,
            "carry-carrier translation")
    terminal_full = tuple(endpoint.build_set(PAT_Q12, ZERO))
    require(tuple(terminal) == terminal_full, "terminal Q12 definitions differ")
    terminal_atom = marked.prefix_intervals(prefixes[0][CLOCK][marked.KAPPA])
    return {
        "base": (source, target, terminal_full),
        "coefficient_atom": (source_carry, target_carry, terminal_atom),
    }


def build_bank(carrier, terminal, *, shift_present=False,
               translate_deep_leg=False):
    """Return P_xi,Q_xi and their separately allocated endpoint DFTs."""
    prime, root = endpoint.MODS[0]
    z13 = pow(root, NRED // P, prime)
    q_starts = [a for a, _b in terminal]
    tabs = endpoint.make_tabs(terminal, X0, ((prime, root),))
    p_bank = {}
    q_bank = {}
    overlaps = {}
    present_sets = present_cache(shifted=shift_present)
    for address, ell in REPS.items():
        present = present_sets[address]
        restricted = intersect_sorted(present, carrier)
        left, overlap = endpoint.x_sweep(
            restricted, terminal, q_starts, X0,
            ((prime, root),), tabs,
        )
        right = endpoint.endpoint_sum(restricted, -Y0, ((prime, root),))
        deep_phase = pow(z13, ell[DEEP], prime)
        if translate_deep_leg:
            # Translation by tau sends the m=1 deepest coefficient through
            # exp(-2*pi*i*c3*tau)=zeta_13^(-1).
            deep_phase = deep_phase * pow(z13, -1, prime) % prime
        p_bank[address] = deep_phase * left[0] % prime
        q_bank[address] = right[0] % prime
        overlaps[address] = overlap

    powers = tuple(pow(z13, exponent, prime) for exponent in range(P))
    left_dft = {}
    right_dft = {}
    for point in KEYS:
        left_dft[point] = sum(
            p_bank[address]
            * powers[-(address[0] * point[0] + address[1] * point[1]) % P]
            for address in KEYS
        ) % prime
        right_dft[point] = sum(
            q_bank[address]
            * powers[(address[0] * point[0] + address[1] * point[1]) % P]
            for address in KEYS
        ) % prime

    # Inverse-transform controls keep the endpoint allocation before Radon.
    for address in KEYS:
        recover_p = sum(
            left_dft[point]
            * powers[(address[0] * point[0] + address[1] * point[1]) % P]
            for point in KEYS
        ) % prime
        recover_q = sum(
            right_dft[point]
            * powers[-(address[0] * point[0] + address[1] * point[1]) % P]
            for point in KEYS
        ) % prime
        require(recover_p == P**2 * p_bank[address] % prime,
                "left DFT inversion")
        require(recover_q == P**2 * q_bank[address] % prime,
                "right DFT inversion")
    return p_bank, q_bank, left_dft, right_dft, overlaps


def support(bank):
    return sum(value != 0 for value in bank.values())


def digest(bank):
    return sha256(",".join(str(bank[key]) for key in KEYS).encode()).hexdigest()


def product_bank(left, right):
    return {key: left[key] * right[key] % endpoint.MODS[0][0] for key in KEYS}


def affine_covariances(source, target, *, limit=8):
    """Find B(x)=c*zeta^(a.x)*A(x+d) on F_13^2."""
    prime, root = endpoint.MODS[0]
    z13 = pow(root, NRED // P, prime)
    powers = tuple(pow(z13, exponent, prime) for exponent in range(P))
    answers = []
    for d0 in range(P):
        for d1 in range(P):
            if any(
                (target[x] == 0) !=
                (source[((x[0] + d0) % P, (x[1] + d1) % P)] == 0)
                for x in KEYS
            ):
                continue
            first = next((x for x in KEYS if target[x]), None)
            if first is None:
                if all(not source[x] for x in KEYS):
                    answers.append(((d0, d1), (0, 0), 0))
                continue
            shifted_first = ((first[0] + d0) % P, (first[1] + d1) % P)
            for a0 in range(P):
                for a1 in range(P):
                    char_first = powers[(a0 * first[0] + a1 * first[1]) % P]
                    scalar = (
                        target[first]
                        * pow(source[shifted_first] * char_first % prime,
                              prime - 2, prime)
                    ) % prime
                    good = True
                    for x in KEYS:
                        shifted = ((x[0] + d0) % P, (x[1] + d1) % P)
                        expected = (
                            scalar * source[shifted]
                            * powers[(a0 * x[0] + a1 * x[1]) % P]
                        ) % prime
                        if target[x] != expected:
                            good = False
                            break
                    if good:
                        answers.append(((d0, d1), (a0, a1), scalar))
                        if len(answers) >= limit:
                            return tuple(answers)
    return tuple(answers)


def scalar_ratio(source, target):
    prime = endpoint.MODS[0][0]
    first = next((key for key in KEYS if source[key]), None)
    if first is None or (target[first] == 0):
        return None
    ratio = target[first] * pow(source[first], prime - 2, prime) % prime
    if all(target[key] == ratio * source[key] % prime for key in KEYS):
        return ratio
    return None


def main():
    carriers = build_carriers()
    prime, root = endpoint.MODS[0]
    z13 = pow(root, NRED // P, prime)
    phase_left = pow(root, (-X0 * R * SHIFT) % NRED, prime)
    phase_right = pow(root, (Y0 * R * SHIFT) % NRED, prime)
    expected_p_scalar = phase_left * pow(z13, -1, prime) % prime
    expected_q_scalar = phase_right
    require(expected_p_scalar * expected_q_scalar % prime == 1,
            "full marked-current translation phase")

    print("THM-2749 fixed-triangle endpoint-harmonic private probe")
    print(f"triangle=(X,m,Y)=({X},1,{Y}); R={R}; tau_grid={SHIFT}")
    print(f"deep-owner relation rows={RELATION_ROWS}")
    print("typing: 169 twists resolve only pi(r mod 13) in F13^2; "
          "they do not retain an exact r in the full relation lattice")

    for name, (source_carrier, target_carrier, terminal) in carriers.items():
        source = build_bank(source_carrier, terminal)
        carrier_only = build_bank(target_carrier, terminal)
        full_rechart = build_bank(
            target_carrier, terminal,
            shift_present=True, translate_deep_leg=True,
        )
        labels = ("P", "Q", "Lstar", "Rstar")
        source_banks = source[:4]
        carrier_banks = carrier_only[:4]
        rechart_banks = full_rechart[:4]
        source_h = product_bank(source[0], source[1])
        carrier_h = product_bank(carrier_only[0], carrier_only[1])
        rechart_h = product_bank(full_rechart[0], full_rechart[1])

        require(scalar_ratio(source[0], rechart_banks[0]) == expected_p_scalar,
                f"{name}: full-rechart P phase")
        require(scalar_ratio(source[1], rechart_banks[1]) == expected_q_scalar,
                f"{name}: full-rechart Q phase")
        require(scalar_ratio(source_h, rechart_h) == 1,
                f"{name}: translated full current")

        print(f"configuration={name}")
        print(f"  carrier_grid_mass(source,target)="
              f"({interval_mass(source_carrier)},{interval_mass(target_carrier)}); "
              f"terminal_grid_mass={interval_mass(terminal)}")
        print(f"  source_support(P,Q,L,R,H)="
              f"{tuple(support(bank) for bank in source_banks + (source_h,))}")
        print(f"  carrier_only_support(P,Q,L,R,H)="
              f"{tuple(support(bank) for bank in carrier_banks + (carrier_h,))}")
        print(f"  source_digests(P,Q,L,R,H)="
              f"{tuple(digest(bank) for bank in source_banks + (source_h,))}")
        print(f"  carrier_only_digests(P,Q,L,R,H)="
              f"{tuple(digest(bank) for bank in carrier_banks + (carrier_h,))}")
        for label, before, after in zip(labels, source_banks, carrier_banks):
            print(f"  carrier_only_affine_covariance[{label}]="
                  f"{affine_covariances(before, after)}")
        print(f"  carrier_only_affine_covariance[H]="
              f"{affine_covariances(source_h, carrier_h)}")
        print(f"  full_rechart_scalars(P,Q,H)="
              f"{(scalar_ratio(source[0], rechart_banks[0]), scalar_ratio(source[1], rechart_banks[1]), scalar_ratio(source_h, rechart_h))}")
        print(f"  expected_full_rechart_scalars(P,Q,H)="
              f"{(expected_p_scalar, expected_q_scalar, 1)}")
        print(f"  overlap_range(source)="
              f"({min(source[4].values())},{max(source[4].values())}); "
              f"overlap_range(carrier_only)="
              f"({min(carrier_only[4].values())},{max(carrier_only[4].values())})")

    print("CONCLUSION: translating only the THM-2749 carrier is not an "
          "affine shift/character covariance of either separate endpoint "
          "bank or their product.  Translating the entire present endpoint "
          "packet and the separate deepest leg restores exact P/Q covariance; "
          "their reciprocal 169th-root phases cancel before Radon.")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
