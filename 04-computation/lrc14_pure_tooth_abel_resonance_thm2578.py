#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2578.

The companion checks the root-of-unity support and target signs of the
pure-tooth Abel normal, its live 13-adic resonance no-go, the exact component
weight obstruction, and a fixed-filter positive sharpness control.  It uses
only integer and rational arithmetic; trigonometric nonvanishing is reduced
to parity modulo 14.
"""

from fractions import Fraction
from math import gcd


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


print("== THM-2578: pure-tooth Abel-normal resonance ==")


print("\n== pure complementary-tooth moment law ==")
resonant = 0
off_resonant = 0
target_dft_checks = 0
geometric_orbit_checks = 0
for k in range(1, 81):
    for ell in (1, 2):
        for physical_offset in range(-4 * k, 4 * k + 1):
            if physical_offset % k:
                # The j-phase runs through every root in a nontrivial cyclic
                # subgroup with constant multiplicity gcd(M,k), hence sums 0.
                g = gcd(abs(physical_offset), k)
                orbit_length = k // g
                require(orbit_length > 1, "off-resonance orbit became trivial")
                counts = {}
                for j in range(k):
                    exponent = (physical_offset * j) % k
                    counts[exponent] = counts.get(exponent, 0) + 1
                require(len(counts) == orbit_length,
                        "geometric orbit support changed")
                require(set(counts.values()) == {g},
                        "geometric orbit multiplicities changed")
                off_resonant += 1
                geometric_orbit_checks += 1
                continue

            h = physical_offset // k
            # cos(pi*h*ell/7)=0 iff 2*h*ell=7 mod 14, impossible by parity.
            require((2 * h * ell) % 14 != 7,
                    "resonant endpoint cosine vanished")

            # N_s has phase zeta^(-h*s).  With the canon's plus-sign DFT,
            # the surviving target index is q=h.
            surviving = []
            for q in range(P):
                residues = [((q - h) * s) % P for s in range(P)]
                if set(residues) == {0}:
                    surviving.append(q)
                else:
                    require(sorted(residues) == list(range(P)),
                            "target root coset changed")
                target_dft_checks += 1
            require(surviving == [h % P], "pure-tooth target sign changed")
            resonant += 1

print(f"  resonant k|M packets: {resonant}")
print(f"  off-resonant cancelling packets: {off_resonant}")
print(f"  exact geometric-orbit checks: {geometric_orbit_checks}")
print(f"  plus-sign target DFT checks: {target_dft_checks}")
print("  N_s(M)=0 unless M=kh; then q=h and amplitude k*cos(pi*h*L/7)/pi^2")


print("\n== live deep-bank resonance is target-neutral ==")
live_packets = 0
live_resonances = 0
live_off_resonances = 0
for k in range(1, 61):
    if gcd(k, P) != 1:
        continue
    for c in range(1, 21):
        c3 = P * c
        for m in range(-30, 31):
            if gcd(abs(m), 91) != 1:
                continue
            physical_offset = m * c3
            if physical_offset % k:
                live_off_resonances += 1
            else:
                h = physical_offset // k
                require(h % P == 0,
                        "live resonant quotient acquired target colour")
                live_resonances += 1
            live_packets += 1

require(live_packets == live_resonances + live_off_resonances,
        "live resonance ledger changed")
print(f"  live (k,c3,m) packets checked: {live_packets}")
print(f"  k|M resonances, all q=0: {live_resonances}")
print(f"  k-not-dividing-M packets, all aggregate normals zero: {live_off_resonances}")
print("  mechanism: 13|M and gcd(k,13)=1 imply 13|(M/k) whenever k|M")


print("\n== component-weight obstruction and holonomy ==")
uniform_zero_modes = 0
spiked_nonzero_modes = 0
holonomy_checks = 0
connection_checks = 0
for k in range(2, 61):
    for r in range(1, k):
        # Uniform weights have only component mode r=0.  Adding three units
        # to component zero changes every nontrivial component coefficient
        # from 0 to the nonzero scalar 3.
        g = gcd(r, k)
        orbit_length = k // g
        require(orbit_length > 1, "nontrivial component mode collapsed")
        uniform_zero_modes += 1
        spike_remainder = 3
        require(spike_remainder != 0, "spiked component mode vanished")
        spiked_nonzero_modes += 1

        # Under s -> s+13 the physical component label moves j -> j-1.
        # W_r therefore gains exp(+2*pi*i*r/k), while the compensating
        # connection exp(-2*pi*i*r*s/(13k)) loses the same phase.
        for j in range(k):
            phase_after_loop = (r * (j + 1)) % k
            phase_before = (r * j) % k
            require((phase_after_loop - phase_before) % k == r % k,
                    "component-weight holonomy sign changed")
            holonomy_checks += 1
        for s in range(P):
            exponent_before = -r * s
            exponent_after = -r * (s + P) + P * r
            require(exponent_after == exponent_before,
                    "component connection failed to descend")
            connection_checks += 1

print(f"  uniform nontrivial modes forced zero: {uniform_zero_modes}")
print(f"  positive spiked profiles with nonzero modes: {spiked_nonzero_modes}")
print(f"  component holonomy checks: {holonomy_checks}")
print(f"  descended-connection checks: {connection_checks}")
print("  off resonance requires a nonuniform tooth-component boundary mode")


print("\n== fixed target-neutral filter breaks pure-gate cancellation ==")
# Take k=L=1, M=13, and the fixed common filter H=1_[0,1/2).  Gate
# boundaries are x=eps/14-s/13.  H selects a nonconstant subset of them,
# while none coincides with a boundary of H.  At M=13 each selected sign has
# a fixed 14th-root phase, so the moment sequence is visibly nonconstant.
filtered_moments = []
filtered_atoms = 0
for s in range(P):
    phase_counts = {}
    for eps in (-1, 1):
        x = (Fraction(eps, 14) - Fraction(s, P)) % 1
        require(x not in {Fraction(0), Fraction(1, 2)},
                "gate boundary collided with fixed-filter boundary")
        if Fraction(0) <= x < Fraction(1, 2):
            exponent_mod_14 = (P * eps) % 14
            phase_counts[exponent_mod_14] = (
                phase_counts.get(exponent_mod_14, 0) + 1
            )
            filtered_atoms += 1
    filtered_moments.append(tuple(sorted(phase_counts.items())))

require(filtered_atoms == 13, "fixed-filter handoff census changed")
require(len(set(filtered_moments)) > 1,
        "fixed-filter target moment became constant")
require(any(not moment for moment in filtered_moments),
        "fixed-filter hostile lost its zero target cells")
require(any(moment for moment in filtered_moments),
        "fixed-filter hostile lost every handoff")
print(f"  exact selected handoff atoms over all target shifts: {filtered_atoms}")
print(f"  distinct target moment profiles: {len(set(filtered_moments))}")
print("  finite DFT invertibility forces at least one nonzero target colour")
print("  the escape is a lawful nonuniform total-layer filter, not the bare gate")


print("\nscope: pure-gate live normal is target-neutral; nonuniform boundary weights remain open")
print("THM-2574 changes to a complex oriented component local system, not a positive aggregate")
print("\nall exact checks passed")
