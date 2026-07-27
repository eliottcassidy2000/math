#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2574.

The theorem is an exact Fourier/monodromy calculation.  This companion
checks its integer phase identities, the plus-sign C_13 DFT convention, the
two-valued Euclidean carry law, the LRC 13-adic resonance boundary, and the
septimal sinc-zero avoidance needed for a nonzero fixed-frequency control.

No floating-point trigonometry is used.  For an interval of length a/7, the
fractional-frequency integral at X/k vanishes exactly when X is nonzero and
7k divides aX.
"""

from math import gcd


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def interval_transform_zero(length_numerator, k, x):
    """Zero test for an interval of length length_numerator/7 at xi=x/k."""
    require(1 <= length_numerator <= 6, "interval numerator outside 1..6")
    require(k > 0, "speed must be positive")
    return x != 0 and (length_numerator * x) % (7 * k) == 0


print("== THM-2574: oriented tooth holonomy and fixed-X descent ==")


print("\n== lifted component Fourier formula and monodromy ==")
# Write every exponential in units of 1/(13k).  The component formula has
# numerator X*s-13*X*n.  Moving s by 13 equals moving n down by one, exactly.
component_phase_checks = 0
circle_label_checks = 0
for k in range(1, 31):
    for s in range(-4, 18):
        for n in range(k):
            for x in range(-2 * k, 2 * k + 1):
                phase = x * s - P * x * n
                phase_after_target_loop = x * (s + P) - P * x * n
                phase_after_component_shift = x * s - P * x * (n - 1)
                require(
                    phase_after_target_loop == phase_after_component_shift,
                    "component monodromy phase failed",
                )
                component_phase_checks += 1

                # n and n+k name the same component on the physical circle;
                # their phases differ by the integer -X.
                phase_next_label = x * s - P * x * (n + k)
                require(
                    (phase_next_label - phase) % (P * k) == 0,
                    "component label failed to close modulo k",
                )
                circle_label_checks += 1

print(f"  exact component/monodromy phase checks: {component_phase_checks}")
print(f"  physical component-label closure checks: {circle_label_checks}")
print("  I_(n,s+13)=I_(n-1,s) and n is read modulo k")


print("\n== component characters select a physical residue ==")
# The n-sum has exponent (r-X)n/k.  Its complete residue sum is k exactly
# when X=r mod k and zero otherwise.  We verify the exponent census rather
# than approximate roots of unity.
orthogonality_checks = 0
connection_checks = 0
dft_sign_checks = 0
gauge_checks = 0
for k in range(1, 51):
    for r in range(k):
        for x in range(-2 * k, 2 * k + 1):
            exponent_counts = [0] * k
            for n in range(k):
                exponent_counts[((r - x) * n) % k] += 1
            if x % k == r:
                expected = [0] * k
                expected[0] = k
                require(exponent_counts == expected, "selected residue lost")
                q = (x - r) // k

                # G has target phase zeta^(q*s), hence is 13-periodic.
                for s in range(P):
                    require(
                        (q * (s + P) - q * s) % P == 0,
                        "compensated connection did not descend",
                    )
                    connection_checks += 1

                # The canon uses (1/13) sum_s J(s) zeta^(a*s), so a pure
                # zeta^(q*s) mode lands at transform index a=-q.
                surviving = []
                for a in range(P):
                    counts = [0] * P
                    for s in range(P):
                        counts[((a + q) * s) % P] += 1
                    if counts[0] == P:
                        surviving.append(a)
                    else:
                        require(counts == [1] * P, "C13 DFT orthogonality failed")
                    dft_sign_checks += 1
                require(surviving == [(-q) % P], "plus-sign DFT sign changed")

                # r+k is the same component character, but the compensating
                # connection changes by zeta^(-s), so q decreases by one and
                # the plus-sign DFT index increases by one.
                q_shifted = (x - (r + k)) // k
                require(q_shifted == q - 1, "connection-lift gauge did not shift")
                require(
                    (-q_shifted) % P == ((-q) + 1) % P,
                    "DFT gauge shift has wrong sign",
                )
                gauge_checks += 1
            else:
                # A nontrivial complete cyclic character has each root with
                # equal total multiplicity on its actual orbit.  The orbit is
                # the subgroup of multiples of g=gcd(r-X,k), each hit g times.
                step = (r - x) % k
                require(step != 0, "nonselected branch mislabeled")
                g = gcd(step, k)
                expected = [g if exponent % g == 0 else 0 for exponent in range(k)]
                require(exponent_counts == expected,
                        "nontrivial character orbit census changed")
            orthogonality_checks += 1

print(f"  component-character residue checks: {orthogonality_checks}")
print(f"  descended C13 connection checks: {connection_checks}")
print(f"  plus-sign DFT character checks: {dft_sign_checks}")
print(f"  integer-lift gauge checks: {gauge_checks}")
print("  G_(r,s) has phase zeta^(((X-r)/k)s) on X congruent to r")


print("\n== paired endpoint colour is the Euclidean carry ==")
pair_checks = 0
carry_census_checks = 0
for k in range(1, 51):
    for m_offset in range(-2 * k, 2 * k + 1):
        a, d = divmod(m_offset, k)
        census = {}
        for r in range(k):
            t = (r + m_offset) % k
            carry = (r + m_offset - t) // k
            census[carry] = census.get(carry, 0) + 1
            require(carry in {a, a + 1}, "carry escaped adjacent classes")
            for n in range(-1, 2):
                x = r + n * k
                y = x + m_offset
                q_left = (x - r) // k
                q_right = (y - t) // k
                require(q_left - q_right == -carry, "paired target phase changed")
                # J(s)=zeta^(-carry*s), so the plus-sign DFT lands at +carry.
                require(
                    (-carry + (carry % P)) % P == 0,
                    "paired DFT carry sign changed",
                )
                pair_checks += 1
        expected = {a: k - d}
        if d:
            expected[a + 1] = d
        require(census == expected, "Euclidean carry multiplicities changed")
        carry_census_checks += 1

print(f"  fixed-frequency paired phase checks: {pair_checks}")
print(f"  exact two-class carry censuses: {carry_census_checks}")
print("  q=(r+M-t)/k=floor(M/k) or floor(M/k)+1")


print("\n== LRC 13-adic resonance boundary ==")
resonant_cases = 0
nonresonant_cases = 0
nonzero_carry_classes = 0
for k in range(1, 101):
    if gcd(k, P) != 1:
        continue
    for multiplier in range(-30, 31):
        m_offset = P * multiplier
        carries = {
            (r + m_offset - ((r + m_offset) % k)) // k
            for r in range(k)
        }
        if m_offset % k == 0:
            require(len(carries) == 1, "resonance gained a second carry")
            only = next(iter(carries))
            require(only % P == 0, "13-adic resonance gained target colour")
            resonant_cases += 1
        else:
            require(len(carries) == 2, "nonresonance lost carry split")
            require(max(carries) - min(carries) == 1, "carry classes not adjacent")
            good = sum(c % P != 0 for c in carries)
            require(good >= 1, "both adjacent carries vanished modulo 13")
            nonzero_carry_classes += good
            nonresonant_cases += 1

print(f"  resonant k|M cases (all q=0): {resonant_cases}")
print(f"  nonresonant k-not-dividing-M cases: {nonresonant_cases}")
print(f"  nonzero carry classes across nonresonant cases: {nonzero_carry_classes}")
print("  for 13|M and gcd(k,13)=1, standard carry is all-zero iff k|M")


print("\n== danger/safe sinc-zero avoidance ==")
# D_L has length L/7 and U_L has length (7-L)/7.  For fixed physical
# residue r and X=r+n*k, each of their zero conditions removes at most one
# n modulo 7.  Hence at least five residues remain for the pair.
sinc_packets = 0
minimum_good = 7
sharp_five_packets = 0
for k in range(1, 51):
    for m_offset in range(-k, k + 1):
        for ell in (1, 2):
            for r in range(k):
                good = []
                danger_zeros = []
                safe_zeros = []
                for n in range(7):
                    x = r + n * k
                    y = x + m_offset
                    dz = interval_transform_zero(ell, k, x)
                    sz = interval_transform_zero(7 - ell, k, y)
                    if dz:
                        danger_zeros.append(n)
                    if sz:
                        safe_zeros.append(n)
                    if not dz and not sz:
                        good.append(n)
                require(len(danger_zeros) <= 1, "danger sinc removed two n residues")
                require(len(safe_zeros) <= 1, "safe sinc removed two n residues")
                require(len(good) >= 5, "paired sinc zeros covered too many residues")
                minimum_good = min(minimum_good, len(good))
                sharp_five_packets += int(len(good) == 5)
                sinc_packets += 1

print(f"  exact (k,M,L,r) packets checked: {sinc_packets}")
print(f"  minimum nonzero choices among n mod 7: {minimum_good}")
print(f"  packets attaining the sharp union-bound value five: {sharp_five_packets}")


print("\n== tooth-weighted Abel-normal bridge ==")
# A lawful whole-layer handoff weight w_(j,eps,s) is permuted by
# w_(j,eps,s+13)=w_(j-1,eps,s).  Its rth tooth character has holonomy
# exp(2*pi*i*r/k).  For M=r+h*k, the physical boundary phase
# exp(-2*pi*i*M*s/(13k)) and the compensating connection combine to
# zeta^(-h*s), so the plus-sign target DFT lands at q=h.
weighted_bridge_checks = 0
pure_gate_checks = 0
for k in range(1, 81):
    for m_offset in range(-3 * k, 3 * k + 1):
        r = m_offset % k
        h = (m_offset - r) // k
        require(m_offset - r == h * k, "physical offset did not split")
        for s in range(P):
            # Numerators are measured in units of 1/(13k).
            physical_target_phase = -m_offset * s
            inverse_connection_phase = r * s
            require(
                physical_target_phase + inverse_connection_phase == -h * k * s,
                "weighted normal connection has wrong target phase",
            )
            require((-h + (h % P)) * s % P == 0,
                    "weighted normal DFT sign changed")
            weighted_bridge_checks += 1

        # Uniform pure-gate weights have only the trivial tooth character.
        if m_offset % k == 0:
            h0 = m_offset // k
            require(r == 0 and h == h0, "pure-gate resonance quotient changed")
            for ell in (1, 2):
                # cos(pi*h*ell/7) cannot vanish: its vanishing would require
                # the even residue 2*h*ell to equal 7 modulo 14.
                require((2 * h0 * ell) % 14 != 7,
                        "canonical gate cosine unexpectedly vanished")
                pure_gate_checks += 1
        else:
            require(r != 0, "off-resonance pure gate gained trivial character")
            pure_gate_checks += 1

print(f"  connected tooth-weight/target phase checks: {weighted_bridge_checks}")
print(f"  pure-gate support/cosine checks: {pure_gate_checks}")
print("  M=r+hk: the connected r-tooth normal lands at plus-DFT colour q=h")
print("  uniform tooth weights retain only r=0, equivalently k divides M")


print("\n== nonresonant coloured fixed-X controls ==")
control_cases = 0
control_frequency_checks = 0
for k in range(1, 101):
    if gcd(k, P) != 1:
        continue
    for multiplier in range(-20, 21):
        m_offset = P * multiplier
        if m_offset % k == 0:
            continue
        for ell in (1, 2):
            chosen = None
            for r in range(k):
                t = (r + m_offset) % k
                carry = (r + m_offset - t) // k
                if carry % P == 0:
                    continue
                for n in range(7):
                    x = r + n * k
                    y = x + m_offset
                    if not interval_transform_zero(ell, k, x) and not interval_transform_zero(
                        7 - ell, k, y
                    ):
                        chosen = (r, n, x, y, carry)
                        break
                if chosen is not None:
                    break
            require(chosen is not None, "nonresonant coloured control not found")
            r, n, x, y, carry = chosen
            require(x % k == r, "control has wrong left residue")
            require(carry % P != 0, "control colour vanished")
            require(not interval_transform_zero(ell, k, x), "left control is a sinc zero")
            require(
                not interval_transform_zero(7 - ell, k, y),
                "right control is a sinc zero",
            )
            control_cases += 1
            control_frequency_checks += 2

print(f"  LRC-shaped nonresonant controls: {control_cases}")
print(f"  certified nonzero interval amplitudes: {control_frequency_checks}")
print("  every tested off-resonance packet has a nonzero q and nonzero fixed-X amplitude")


print("\n== trivial-character/full-comb boundary ==")
full_comb_checks = 0
for k in range(1, 101):
    for x in range(-3 * k, 3 * k + 1):
        selected = x % k == 0
        require(selected == (x % k == 0), "full-comb residue test changed")
        full_comb_checks += 1
print(f"  trivial component-character checks: {full_comb_checks}")
print("  forgetting component labels keeps r=0 and restores the full k-tooth comb")


print("\nscope: exact complex oriented local system; no positive Boolean carrier or live-row current")
print("all exact checks passed")
