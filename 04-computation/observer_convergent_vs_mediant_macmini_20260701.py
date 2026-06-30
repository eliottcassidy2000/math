#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The observer's escape: convergent (n>=7) vs mediant (n<=6); and the spectral/Jacobi node.
mac-mini-2026-07-01-S46.

The covering-min witness t* (= M, speed-1 binding) is a continued fraction:
  n<=6: the MEDIANT  2/(2n-1) = [0; n-1, 2]  (drop-2 / projective regime, while PG(2,n-1) exists)
  n>=7: the CONVERGENT n/Phi_6(n) = [0; n-1, n]  (hexagonal/Eisenstein regime)
Both lie on the Stern-Brocot path from 1/(n-1): [0;n-1,2],[0;n-1,3],..,[0;n-1,n]=convergent.
The convergent's CF-recurrence matrix M=[[Phi6,n],[n-1,1]] is SL2 (det 1) -- the Jacobi/Stieltjes spectral
object; covering radius = its convergent = the A2-Catalan (semicircle) Jacobi spectral edge (klein).
Certificate = Lagrange (convergents are best approximations) = Markov-Stieltjes (Jacobi convergent optimal).
"""
from __future__ import annotations
import functools
from fractions import Fraction as F
import numpy as np
print = functools.partial(print, flush=True)


def cf(p, q):
    a = []
    while q:
        a.append(p // q); p, q = q, p - (p // q) * q
    return a


def main():
    print("=" * 80)
    print("OBSERVER'S ESCAPE: convergent (n>=7) vs mediant (n<=6) + the spectral/Jacobi node")
    print("=" * 80)

    print("\n[1] the escape t* (=M) along the Stern-Brocot path from 1/(n-1):")
    for n in (5, 6, 7, 8, 11, 14):
        med, con = F(2, 2*n-1), F(n, n*n-n+1)
        reg = "MEDIANT  [0;n-1,2]" if n <= 6 else "CONVERGENT [0;n-1,n]"
        print(f"    n={n:>2}: mediant {med}=[0;{cf(2,2*n-1)[1:]}]  convergent {con}=[0;{cf(n,n*n-n+1)[1:]}]  "
              f"covering-min = {reg}")
    print("    2/(2n-1) = mediant(1/n, 1/(n-1)) = the shallow escape (partial quotient 2).")
    print("    n/Phi_6 = [0;n-1,n] = the full convergent = the DEEPEST escape (partial quotient n).")
    print("    transition at n=7 = first PG(2,n-1) failure PG(2,6) (opus): the projective/mediant escape is")
    print("    blocked, the observer descends to the hexagonal/Eisenstein CONVERGENT.")

    print("\n[2] the spectral/Jacobi node: the convergent = the SL2 CF-recurrence matrix (det 1):")
    for n in (7, 14):
        A = np.array([[n, 1], [1, 0]]); B = np.array([[n-1, 1], [1, 0]]); M = A @ B
        P = n*n-n+1
        print(f"    n={n}: M=[[n,1],[1,0]][[n-1,1],[1,0]] = {M.tolist()} = [[Phi6={P},n],[n-1,1]], "
              f"det={int(round(np.linalg.det(M)))}, slope=n/Phi6={n}/{P}")
    print("    det M = Phi6 - n(n-1) = 1 (SL2/unimodular). The CF recurrence q_k=a_k q_{k-1}+q_{k-2} is")
    print("    tridiagonal (Jacobi); the convergent = the Pade/Gauss-quadrature spectral edge. klein: the A2")
    print("    'tridiagonalized' Jacobi operator has Catalan moments (semicircle); covering radius = its edge.")
    print("    CERTIFICATE = LAGRANGE (convergents = best approximations) = MARKOV-STIELTJES (the Jacobi")
    print("    convergent is the optimal spectral/quadrature bound) -- CLASSICAL.")

    print("\n" + "=" * 80)
    print("UNIFICATION: the observer's escape (the lonely runner's best-approximation lonely time) IS the")
    print("CF convergent n/Phi_6 = [0;n-1,n] for n>=7 (the mediant 2/(2n-1) for n<=6), and the convergent IS")
    print("the SL2/Jacobi spectral edge. The proof reduces to: (a) the LRC->CF reduction (the covering forces")
    print("the convergent escape past n=7 = the PG(2,6) failure = klein-S27 LRC->2D) -- the OPEN piece; plus")
    print("(b) convergent-is-best (Lagrange/Markov-Stieltjes) -- CLASSICAL. The one node = the LRC->CF reduction.")
    print("=" * 80)


if __name__ == "__main__":
    main()
