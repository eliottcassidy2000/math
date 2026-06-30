#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THE LEAST-EIGENVALUE CERTIFICATE for the LRC apex floor. mac-mini-2026-06-30-S37.

The floor's binding value is 2 + lambda_min(A(C_p)) = 4cos^2(3pi/7) for p=7 (HYP-3590). This session makes
it a clean POSITIVE-DEFINITENESS (Bochner) certificate of an explicit matrix, identifies the eigenvector
(worst mode), proves it set-independent, and pins WHY it is positive: the apex cycle is ODD = NON-BIPARTITE.

THE CERTIFICATE. For an apex core O subset Z_p (p the apex prime, odd), the autocorrelation Gram is the
circulant whose eigenvalues are |sum_{x in O} w^{kx}|^2 (k=0..p-1). The k=0 (mean) mode is |O|^2; the
floor is the least NONZERO-mode eigenvalue:
    gap(O) = min_{k!=0} |sum_{x in O} w^{kx}|^2.
For the binding DOUBLET O={a,b} this Gram is exactly 2I + A(C_p) (C_p = Cay(Z_p,{+-(b-a)}) = the p-cycle),
and:
    lambda_min(2I + A(C_p)) = 2 + lambda_min(A(C_p)) = 2 - 2cos(pi/p) = 4 sin^2(pi/2p).
POSITIVE iff lambda_min(A(C_p)) > -2 iff C_p has NO eigenvalue -2 iff C_p is NON-BIPARTITE iff p ODD.
So: the least-eigenvalue certificate = the ODD-cycle (non-bipartite) certificate. It is set-INDEPENDENT
(depends only on p). klein-S18: this is the sigma-ODD / Bochner / DISCRETE side (the apex Gram gap, binds at
the doublet) -- the EXISTENCE-supporting content, NOT the measure-rho_j (which binds at the cusp).
"""
from __future__ import annotations
import functools, math, cmath, itertools
print = functools.partial(print, flush=True)


def cycle_adj_eigs(p):
    return [2 * math.cos(2 * math.pi * k / p) for k in range(p)]


def gram_min_nonzero(O, p):
    W = cmath.exp(2j * math.pi / p)
    if not O or set(x % p for x in O) == set(range(p)):
        return 0.0
    return min(abs(sum(W**(k * x) for x in O))**2 for k in range(1, p))


def main():
    print("=" * 84)
    print("THE LEAST-EIGENVALUE CERTIFICATE for the LRC apex floor (mac-mini-S37)")
    print("=" * 84)

    # ---- (1) the certificate: 2I + A(C_p) >0 iff p odd (non-bipartite) ----
    print("\n[1] CERTIFICATE -- the doublet Gram 2I + A(C_p) is POSITIVE DEFINITE iff p is ODD:")
    print(f"    {'p':>3} {'parity':>5} {'lambda_min A(C_p)':>18} {'lambda_min 2I+A':>16} {'C_p bipartite?':>15} {'cert':>9}")
    for p in (3, 5, 7, 11, 13, 4, 6):
        eigs = cycle_adj_eigs(p)
        la = min(eigs)
        lg = 2 + la
        bip = (p % 2 == 0)
        cert = "PSD>0" if lg > 1e-9 else "SINGULAR"
        print(f"    {p:>3} {('odd' if p%2 else 'even'):>5} {la:>18.5f} {lg:>16.5f} "
              f"{('YES' if bip else 'no'):>15} {cert:>9}")
    print("    A(C_p) eigenvalues are 2cos(2pi k/p) in (-2,2]; the bottom -2 is reached IFF p even")
    print("    (k=p/2 exists) = C_p bipartite. p ODD => lambda_min(A) = -2cos(pi/p) > -2 => 2I+A > 0.")
    print(f"    p=7: lambda_min(2I+A(C_7)) = 2-2cos(pi/7) = 4sin^2(pi/14) = {4*math.sin(math.pi/14)**2:.5f} "
          f"= 4cos^2(3pi/7) (HYP-3590, THM-590).")

    # ---- (2) the worst eigenvector (the Fejer-Bochner minorant direction) ----
    print("\n[2] THE WORST MODE (least eigenvector) -- the binding resonance direction:")
    for p in (5, 7, 11):
        eigs = cycle_adj_eigs(p)
        kmin = min(range(p), key=lambda k: eigs[k])
        print(f"    p={p}: A(C_p) min eigenvalue at Fourier mode k={kmin} (= (p-1)/2 or (p+1)/2, the near-pi/2 "
              f"frequency); eigenvector v_j = cos(2pi*{kmin}*j/p).")
    print("    => the floor binds in the MIDDLE Fourier mode (k ~ p/2) -- the most-oscillatory resonance,")
    print("    which for ODD p never hits the -1 phase exactly (that is the non-bipartiteness).")

    # ---- (3) the certificate is uniform over ALL cores (the Bochner minorant) ----
    print("\n[3] UNIFORM over all 127 non-full Z_7 cores -- the least-eigenvalue minorant (THM-590):")
    vals = sorted(set(round(gram_min_nonzero(O, 7), 5)
                      for sz in range(1, 7) for O in itertools.combinations(range(7), sz)))
    pos = [v for v in vals if v > 1e-9]
    print(f"    gap values over cores: {vals}; the GLOBAL least nonzero-mode eigenvalue = {min(pos):.5f}")
    print(f"    = 4cos^2(3pi/7) at the doublets. So inf over the finite family of the least eigenvalue is the")
    print(f"    DOUBLET's 4sin^2(pi/14) -- a single set-INDEPENDENT positive number (depends only on p=7).")

    # ---- (4) set-independence + the general apex-prime formula ----
    print("\n[4] SET-INDEPENDENCE -- the certificate depends ONLY on the apex prime p:")
    print(f"    {'apex p':>7} {'floor = 4sin^2(pi/2p)':>22} {'= 2-2cos(pi/p)':>15}")
    for p in (3, 5, 7, 11, 13):
        print(f"    {p:>7} {4*math.sin(math.pi/(2*p))**2:>22.5f} {2-2*math.cos(math.pi/p):>15.5f}")
    print("    The covering set never enters; only p does. This is the set-independent floor the gatekeeper")
    print("    needs -- a least-eigenvalue (Bochner/SOS) certificate, NOT a per-set estimate.")

    print("\n" + "=" * 84)
    print("CERTIFICATE (clean form): for the apex prime p (odd), the doublet Gram 2I+A(C_p) is positive")
    print("definite with lambda_min = 4sin^2(pi/2p) > 0, BECAUSE C_p is non-bipartite (p odd). Over all")
    print("non-full Z_p cores the least nonzero-mode eigenvalue is >= this (THM-590). It is set-independent")
    print("(only p), and it certifies the sigma-ODD / Bochner / DISCRETE content (the odd cycle is present &")
    print("non-degenerate = EXISTENCE, klein-S18), NOT the measure (which binds at the cusp). The floor's")
    print("positivity IS the apex cycle's odd length.")
    print("=" * 84)


if __name__ == "__main__":
    main()
