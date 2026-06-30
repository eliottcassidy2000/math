#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Three coordinated routes to the cyclic-Kershner covering radius (construction = LRC covering-min, n>=7).
mac-mini-2026-06-30-S45. (1) Torus lift; (2) three-distance gaps; (3) spectral/LP.

Construction {1,..,n-2, n(n-1)}, M = n/Phi_6(n). t* = n/Phi_6(n) = zeta_6 (the 60deg hexagonal rotation).
"""
from __future__ import annotations
import functools, math, cmath
print = functools.partial(print, flush=True)


def main():
    print("=" * 82)
    print("THREE ROUTES TO THE zeta_6-LINE COVERING RADIUS (cyclic-Kershner) -- mac-mini-S45")
    print("=" * 82)

    # ---- ROUTE 1: TORUS LIFT  Z[w]/(n-w) ~ Z/Phi6, w = x n ----
    print("\n[1] TORUS LIFT: Z[omega]/(n-omega) ~ Z/Phi_6(n), omega (60deg hex rotation) = mult-by-n:")
    print(f"    {'n':>3} {'Phi6':>5} {'n^2 mod Phi6':>12} {'n^3 mod Phi6':>12} {'ord(n)':>7} {'min|v n mod Phi6|':>18}")
    for n in (7, 11, 13, 14):
        P = n*n-n+1
        n2, n3 = n*n % P, n*n*n % P
        o = next(k for k in range(1, 2*P) if pow(n, k, P) == 1)
        S = list(range(1, n-1))+[n*(n-1)]
        mc = min(min((v*n) % P, P-(v*n) % P) for v in S)
        print(f"    {n:>3} {P:>5} {n2:>12} {n3 if n3<P else n3:>12} {o:>7} {mc:>18}  (= n; M=n/Phi6)")
    print("    n^2 ≡ n-1, n^3 ≡ -1 (mod Phi6) => ord(n)=6 = the hexagonal rotation. The construction's")
    print("    residues v*n mod Phi6 = a hexagonal AP; covering radius = n/Phi6. Reduces to KERSHNER 1939")
    print("    (hexagonal = thinnest LATTICE covering, PROVEN) for lattice/AP coverings.")

    # ---- ROUTE 2: THREE-DISTANCE ----
    print("\n[2] THREE-DISTANCE: the residue AP has gaps EXACTLY {1, n, 2n}; M = min|residue|/Phi6 = n/Phi6:")
    for n in (7, 14):
        P = n*n-n+1
        S = list(range(1, n-1))+[n*(n-1)]
        R = sorted((v*n) % P for v in S)
        gaps = sorted(set((R[(i+1) % len(R)]-R[i]) % P for i in range(len(R))))
        print(f"    n={n}: gaps={gaps} = {{1,{n},{2*n}}}; the missing origin (k=0) => nearest residue +-n => M=n/Phi6.")
    print("    The k=0 'lcm-killer' (the -1 mod Phi6 closing the AP) gives the gap-1; the origin-straddle gives 2n.")
    print("    For any {1}+AP covering, three-distance gives M in closed form; min max-gap among APs -> construction.")

    # ---- ROUTE 3: SPECTRAL / LP (Delsarte / Cohn-Elkies) ----
    print("\n[3] SPECTRAL/LP: covering radius = sup{lonely t} = inf{Fourier-positive certificates}:")
    n = 14; P = n*n-n+1; S = list(range(1, n-1))+[n*(n-1)]; t = n/P
    binding = [v for v in S if abs(min((v*t) % 1, 1-((v*t) % 1)) - n/P) < 1e-9]
    print(f"    n=14: at t*=n/Phi6={n}/{P}, min_v||v t*|| = {n/P:.5f}; BINDING PAIR = {binding} (smallest & outlier).")
    W = cmath.exp(2j*math.pi/P); R = [(v*n) % P for v in S]
    mags = [abs(sum(W**(j*r) for r in R))**2 for j in range(1, P)]
    print(f"    R-hat(j)=sum_R w^(jr): |.|^2 in [{min(mags):.2f},{max(mags):.2f}] -- AP => Dirichlet/Gauss-sum (peaked).")
    print("    The covering radius (max gap of R) is a Delsarte/Cohn-Elkies LP bound; the EISENSTEIN-SYMMETRIC")
    print("    (omega = x n - invariant) Fourier-positive optimizer is the open node (Viazovska-style); it handles")
    print("    ALL coverings, not just lattices/APs.")

    print("\n" + "=" * 82)
    print("SYNTHESIS: the 3 routes reduce the cyclic-Kershner (construction = covering-min) to ONE")
    print("Eisenstein-symmetric Fourier-positivity certificate. (1) lift = the geometry (hexagonal torus,")
    print("Kershner-lattice handles APs); (2) three-distance = the explicit M for APs; (3) spectral/LP = the")
    print("optimality machinery for ALL coverings, with the omega-symmetric extremal function the open node.")
    print("=" * 82)


if __name__ == "__main__":
    main()
