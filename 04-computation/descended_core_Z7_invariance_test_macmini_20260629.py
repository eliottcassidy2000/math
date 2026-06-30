#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Working the DESCENDED CORE: execute klein-S7b's test -- are the 2-adic-descent odd cores O_j
Z_7-cyclic invariant (so the pure cyclotomic SOS closes rho_j), or is Gamma_0(14) averaging needed?
(mac-mini-2026-06-29-S28)

Setup: 2-adic descent (THM-580): S = O u E (odd/even); S' = E/2; recurse; cores [O_0,O_1,...].
The apex is 7 (14=2.7). klein-S7b/HYP-3566: rho_j>=c = the Z_7 cyclotomic Gram gap, IF the core is
Z_7-cyclic invariant; else Gamma_0(14) (HYP-3553) supplies the symmetry. klein-S8 (HYP-3571) confirmed
the floor is set-independent EMPIRICALLY (inf R'=0.344 >= 1/(2 zeta(2))); this tests the MECHANISM.

TEST per core O_j: (1) residues {v mod 7}; (2) is {v mod 7} a union of Z_7^* multiplier orbits
(=> Z_7-invariant => pure cyclotomic SOS)? (3) if not, AVERAGE the safe-set autocorrelation over
the Z_7^* multiplier (= the Gamma_0(14)/congruence symmetrization) -> the symmetrized Gram IS flat/
cyclotomic and SET-INDEPENDENT. The honest question: pure Z_7 or congruence-averaged?
"""
from __future__ import annotations
import functools, math, cmath
print = functools.partial(print, flush=True)


def descend(S):
    """2-adic descent: return the list of odd cores [O_0, O_1, ...]."""
    cores = []
    cur = sorted(set(S))
    seen = 0
    while cur and seen < 20:
        O = [v for v in cur if v % 2 == 1]
        E = [v for v in cur if v % 2 == 0]
        cores.append(O)
        if not E: break
        cur = sorted(set(v // 2 for v in E))
        seen += 1
    return cores


def Z7_orbit_union(residues):
    """Is the residue set (subset of Z_7) a union of Z_7^* multiplier orbits (=> Z_7-invariant)?
    Orbits of Z_7^* on Z_7: {0} and {1..6} (Z_7^* acts transitively on nonzero)."""
    R = set(residues)
    nz = R - {0}
    # invariant under mult by every unit u
    inv = all({(u * r) % 7 for r in R} == R for u in range(1, 7))
    return inv, sorted(R)


def safe_autocorr_spectrum(residues):
    """Treat the core's residues mod 7 as the danger-support; the SAFE set = Z_7 minus the
    danger residues' multiples. Compute the autocorrelation Gram spectrum of the safe set
    (the apex skeleton), and the Z_7^*-AVERAGED (congruence-symmetrized) spectrum."""
    w = cmath.exp(2j * math.pi / 7)
    # safe set: a in Z_7 with v*a != 0 mod 7 for all v in core with 7 nmid v; if some 7|v, that v
    # kills the apex-7 grid -> safe set on the 7-grid is empty, BUT the cyclotomic harmonics persist.
    # Use the residue INDICATOR's autocorrelation as the apex skeleton (the comb's Z_7 Fourier content).
    ind = [1 if (a % 7) in set(r % 7 for r in residues) else 0 for a in range(7)]
    spec = [abs(sum(ind[a] * w**(k*a) for a in range(7)))**2 for k in range(7)]
    # Z_7^*-averaged indicator (symmetrize over the multiplier orbit) -> the congruence average
    avg_ind = [0.0]*7
    for u in range(1, 7):
        for a in range(7):
            avg_ind[(u*a) % 7] += ind[a] / 6
    avg_spec = [abs(sum(avg_ind[a] * w**(k*a) for a in range(7)))**2 for k in range(7)]
    return spec, avg_spec


def main():
    print("=" * 80)
    print("Descended core Z_7-invariance test (klein-S7b) -- pure cyclotomic or Gamma_0(14)?")
    print("=" * 80)
    configs = {
        "tightest {1..12,182}": list(range(1, 13)) + [182],
        "consec {1..13}":       list(range(1, 14)),
        "skip-12 {1..11,13,84}":list(range(1, 12)) + [13, 84],
        "even-heavy {2,4,6,8,10,12,14,1,3,5,7,9,11}": [2,4,6,8,10,12,14,1,3,5,7,9,11],
    }
    n_inv = 0; n_tot = 0
    for name, S in configs.items():
        cores = descend(S)
        print(f"\n{name}:  cores (odd parts per descent level):")
        for j, O in enumerate(cores):
            res = sorted(set(v % 7 for v in O))
            inv, R = Z7_orbit_union(res)
            n_tot += 1; n_inv += int(inv)
            print(f"   O_{j} = {O}")
            print(f"        residues mod 7 = {R};  Z_7^*-invariant: {inv}"
                  + ("" if inv else "  (NOT a multiplier-orbit union)"))

    print("\n" + "=" * 80)
    print(f"Z_7-invariant cores: {n_inv}/{n_tot} -- the descended cores are GENERALLY NOT Z_7^*-invariant")
    print("(their residues mod 7 are not multiplier-orbit unions). So the PURE Z_7 cyclotomic SOS does")
    print("NOT directly apply -- the transitive symmetry must be MANUFACTURED by Z_7^*-AVERAGING = the")
    print("Gamma_0(14) congruence (HYP-3553). Demonstrate the averaging symmetrizes the apex Gram:")
    # demonstrate on a genuinely non-invariant, non-flat core (O_1 of consec = {1,3,5} mod 7)
    res = [1, 3, 5]
    spec, avg_spec = safe_autocorr_spectrum(res)
    print(f"\n   example NON-invariant core residues {res} (= O_1 of consec):")
    print(f"     raw apex Gram spectrum:        {[round(s,3) for s in spec]} (NOT flat => SET-DEPENDENT)")
    print(f"     Z_7^*-averaged (Gamma_0) spec: {[round(s,3) for s in avg_spec]} (FLAT off-0 => SET-INDEP)")
    print("   => the raw Gram of a non-invariant core is non-flat/set-dependent; averaging over the")
    print("   multiplier (the Gamma_0(14) congruence) FLATTENS it to the cyclotomic form, SET-INDEPENDENT")
    print("   -- 'manufacture the transitive symmetry it lacks'. (Caveat: complement-of-a-point cores")
    print("   like {0..5} are already flat -- it is set-INVARIANCE, not flatness per se, that fails.)")

    print("\n" + "=" * 80)
    print("CONCLUSION: the descended cores are NOT pure-Z_7-invariant (covering kills the apex-7 grid,")
    print("residues not multiplier-closed). So rho_j>=c is NOT closed by a bare cyclotomic SOS; it needs")
    print("the Gamma_0(14) congruence AVERAGING over Z_7^* to MANUFACTURE the transitive symmetry, after")
    print("which the averaged Gram is flat/cyclotomic and SET-INDEPENDENT (= klein-S8's inf R'=0.344 >=")
    print("1/(2 zeta(2))). The mechanism is congruence-averaging (Gamma_0(14)), not bare Z_7 -- this")
    print("SHARPENS the proof path: the next step is the Han-Lee Gamma_0(14) 2nd-moment constant.")
    print("=" * 80)


if __name__ == "__main__":
    main()
