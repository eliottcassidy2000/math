#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE ι-ODD CERTIFICATE IN Q(sqrt-3, sqrt-7): it LOCALIZES to sqrt21 -- exhibiting sqrt21 explicitly.

kind-pasteur-2026-07-01-S22. Next step from S22: compute the ι-odd (Gauss) certificate in the compositum
Q(sqrt-3, sqrt-7) = the biquadratic bridge of the Eisenstein COVERING field (sqrt-3, the p=3/LRC(6) apex) and
the certification field (sqrt-7, the p=7/LRC(14) apex).
 * g_3 = i*sqrt3, g_7 = i*sqrt7 (quadratic Gauss sums, both p=3mod4 => imaginary).
 * PRODUCT g_3 * g_7 = i^2 * sqrt21 = -sqrt21  -- a SINGLE REAL number in Q(sqrt21).
 * Hasse-Davenport composite: g_21 = chi_3(7)*chi_7(3)*g_3*g_7 = +sqrt21 (D=21=1 mod4 => real Gauss sum).
 * So the MIXED part of the compositum certificate LOCALIZES to the REAL QUADRATIC subfield Q(sqrt21) -- degree
   2, not degree 4: the bridge costs a field but NOT spectral complexity. And sqrt21 = g_{LRC6-apex} * g_{LRC14-apex}.
Verify all, and exhibit sqrt21 in the covering deep-well D_14 = 21*403.
"""
import sys, cmath, math
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def legendre(a,p):
    a%=p
    return 0 if a==0 else (1 if pow(a,(p-1)//2,p)==1 else -1)
def gauss(p):  # sum_a (a|p) e(a/p)
    return sum(legendre(a,p)*cmath.exp(2j*math.pi*a/p) for a in range(1,p))
def jacobi(a,n):  # a odd, n odd positive
    a%=n; r=1
    while a:
        while a%2==0:
            a//=2
            if n%8 in (3,5): r=-r
        a,n=n,a
        if a%4==3 and n%4==3: r=-r
        a%=n
    return r if n==1 else 0
def gauss_mod(N):  # sum_{a mod N} (a|N)_Jacobi e(a/N)  (quadratic Gauss sum for modulus N)
    return sum(jacobi(a,N)*cmath.exp(2j*math.pi*a/N) for a in range(1,N))

print("="*94); print(" THE COMPOSITUM CERTIFICATE g_3 * g_7 = -sqrt21 (localizes to Q(sqrt21))"); print("="*94)
g3=gauss(3); g7=gauss(7)
print(f"  g_3 = {g3.real:+.6f}{g3.imag:+.6f}i   (= i*sqrt3, |g3|={abs(g3):.6f}, sqrt3={math.sqrt(3):.6f})")
print(f"  g_7 = {g7.real:+.6f}{g7.imag:+.6f}i   (= i*sqrt7, |g7|={abs(g7):.6f}, sqrt7={math.sqrt(7):.6f})")
prod=g3*g7
print(f"  g_3 * g_7 = {prod.real:+.6f}{prod.imag:+.6f}i   ; -sqrt21 = {-math.sqrt(21):.6f}")
print(f"    => g_3*g_7 = -sqrt21 (REAL, single number in Q(sqrt21)): "
      f"{abs(prod.imag)<1e-9 and abs(prod.real+math.sqrt(21))<1e-9}")

print("\n"+"="*94); print(" HASSE-DAVENPORT composite: g_21 = chi_3(7) chi_7(3) g_3 g_7 = +sqrt21 (real Gauss sum, D=21=1mod4)"); print("="*94)
c37=legendre(7,3); c73=legendre(3,7)
pred=c37*c73*prod
g21=gauss_mod(21)
print(f"  chi_3(7)=(7|3)={c37}, chi_7(3)=(3|7)={c73}  => chi_3(7)chi_7(3) = {c37*c73}")
print(f"  predicted g_21 = {c37*c73}*(-sqrt21) = {pred.real:+.6f}{pred.imag:+.6f}i")
print(f"  direct g_21 = sum_a (a|21) e(a/21) = {g21.real:+.6f}{g21.imag:+.6f}i  ; +sqrt21 = {math.sqrt(21):+.6f}")
print(f"    => g_21 = +sqrt21 (D=21=1mod4 real), matches Hasse-Davenport: "
      f"{abs(g21.imag)<1e-9 and abs(g21.real-math.sqrt(21))<1e-9}")

print("\n"+"="*94); print(" THE BIQUADRATIC Q(sqrt-3, sqrt-7): subfields & localization"); print("="*94)
print("  Gal(Q(sqrt-3,sqrt-7)/Q) = (Z/2)^2; three quadratic subfields:")
print("    Q(sqrt-3)  [Eisenstein COVERING, p=3/LRC(6) apex]  -- imaginary")
print("    Q(sqrt-7)  [Gauss CERTIFICATION, p=7/LRC(14) apex] -- imaginary")
print("    Q(sqrt21)  [the REAL compositum, = g_3*g_7]        -- REAL, DEGREE 2")
print("  => the ι-odd certificate's MIXED part (neither pure covering nor pure cert) LOCALIZES to Q(sqrt21):")
print("     it is a SINGLE real quadratic number, NOT degree 4. The bridge costs a field, not spectral degree.")

print("\n"+"="*94); print(" sqrt21 EXHIBITED in the covering deep-well"); print("="*94)
n=14; D14=n*(n-1)*(n*n-n+4)
def squarefree(D):
    D=abs(D); out=1; d=2
    while d*d<=D:
        c=0
        while D%d==0: D//=d; c+=1
        if c%2: out*=d
        d+=1
    return out*(D if D>1 else 1)
sf=squarefree(D14)
print(f"  deep-well D_14 = n(n-1)(n^2-n+4) = {D14}, squarefree {sf} = 21*{sf//21}  => sqrt21 | sqrt(D_14): {sf%21==0}")
print(f"  covering-min = 14/183, Phi6(14)=183=3*61.  The apex-3 (Eisenstein) & apex-7 (Gauss) meet as sqrt21.")
print("  INTERPRETATION: sqrt21 = g_3 * g_7 = (the LRC(6) certificate) x (the LRC(14) certificate). The bridge")
print("   certificate for LRC(14) is the PRODUCT of the proved case's (p=3) and the open case's (p=7) Gauss sums,")
print("   living in the real quadratic Q(sqrt21) -- the certification remains QUADRATIC across the bridge.")
print("DONE.")
