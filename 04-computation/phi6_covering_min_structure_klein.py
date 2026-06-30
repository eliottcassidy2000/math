#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ATTACK n/Phi_6(n) >= 1/n: the Q(sqrt-3) covering-min structure (klein-S23).

Phi_6(n) = n^2 - n + 1 (6th cyclotomic at n). Three identities:
 (A) Phi_6(n) = |n - zeta_6|^2 = the EISENSTEIN NORM (zeta_6 = e^{i pi/3}); so Q(sqrt-3).
 (B) Phi_6(n) = (n-1)^2 + (n-1) + 1 = #points of PG(2, n-1) (projective plane of order n-1).
 (C) thus a (Phi_6(n), n, 1) SINGER difference set exists when n-1 is a prime power.
The inequality n/Phi_6(n) >= 1/n  <=>  n^2 >= n^2-n+1  <=>  n >= 1 (TRIVIAL). So the content is NOT
the inequality but whether n/Phi_6(n) is the genuine COVERING-MIN (the construction is optimal).
"""
import math
from sympy import factorint, isprime, primefactors

def phi6(n): return n*n - n + 1
def is_prime_power(q):
    if q<2: return False
    p=primefactors(q); return len(p)==1

print("="*92)
print(" Phi_6(n) = n^2-n+1: Eisenstein norm | projective-plane order n-1 | prime factors | n/Phi6 vs 1/n")
print("="*92)
print(f" {'n':>3}{'Phi6(n)':>9}{'factors':>14}{'q=n-1 prime-pwr?':>17}{'n/Phi6':>10}{'1/n':>9}{'gap=(n-1)/nPhi6':>16}")
for n in list(range(3,21))+[14]:
    P=phi6(n); fac=factorint(P); facs="*".join(f"{p}^{e}" if e>1 else str(p) for p,e in fac.items())
    q=n-1; pp=is_prime_power(q)
    nP=n/P; inv=1/n; gap=(n-1)/(n*P)
    star=" <== LRC(14)" if n==14 else ""
    print(f" {n:>3}{P:>9}{facs:>14}{str(pp):>17}{nP:>10.5f}{inv:>9.5f}{gap:>16.6f}{star}")

print("\n"+"="*92)
print(" (A) Eisenstein-norm check: Phi_6(n) = |n - zeta_6|^2, zeta_6 = (1+i sqrt3)/2  [Q(sqrt-3)]")
print("="*92)
z=complex(0.5, math.sqrt(3)/2)  # primitive 6th root
for n in [3,7,14]:
    print(f"   n={n}: |n - zeta_6|^2 = {abs(n-z)**2:.4f}   Phi_6(n) = {phi6(n)}   match={abs(abs(n-z)**2-phi6(n))<1e-9}")

print("\n"+"="*92)
print(" prime factors of Phi_6(n) are =1 mod 6 (Eisenstein SPLIT primes), except 3 (ramified):")
print("="*92)
for n in [3,7,14,13,8]:
    P=phi6(n); pf=primefactors(P); cls=[(p, '3(ramified)' if p==3 else f'{p%6} mod6') for p in pf]
    print(f"   n={n}: Phi_6={P} -> primes {pf}; classes {cls}")
print("   => every prime p|Phi_6(n), p!=3, satisfies p=1 mod 6 (p splits in Z[zeta_6], the Eisenstein integers).")

print("\n"+"="*92)
print(" THE CONTENT (since the inequality is trivial): is n/Phi_6(n) the COVERING-MIN?")
print("="*92)
n=14; P=phi6(n)
print(f"   n=14: Phi_6=183=3*61 (61 = 1 mod 6, splits; 3 ramifies). q=13 is PRIME, so PG(2,13) EXISTS")
print(f"   => a (183, 14, 1) SINGER difference set of 14 speeds mod 183 (every nonzero diff once).")
print(f"   n/Phi_6 = 14/183 = {n/P:.6f};  1/n = 1/14 = {1/n:.6f};  margin = {(n-1)/(n*P):.6f}")
print(f"   The inequality 14/183 >= 1/14 is n>=1 (trivial). OPEN: is 14/183 the genuine covering-min,")
print(f"   i.e. is the (183,14,1) projective-plane covering OPTIMAL (no covering beats it)? -- a")
print(f"   design-theoretic optimality question, the Q(sqrt-3)/Eisenstein 'existence-column'.")
