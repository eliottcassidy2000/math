#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE NONABELIAN EXPANDER (PSL_2(F_p) / Fano), SUM-OF-TWO-SQUARES = p mod 4, and the FACILITY-LOCATION POTENTIAL/PoA.

kind-pasteur-2026-07-01-S24. Chasing the LTC lead's next step + the owner's jumping-off concepts.
 (1) EXPANDER: PSL_2(F_p) is the nonabelian apex group; |PSL_2(F_7)|=168=|Aut(Fano plane)|=|GL_3(F_2)|. The LPS
     Ramanujan (p+1)-regular expander has spectral gap |lambda|<=2*sqrt(p) -- the SAME sqrt(p) as the Gauss-sum
     certificate g_p=i*sqrt(p) (S20). So the expander's spectral gap IS the certificate's arithmetic.
 (2) SUM OF TWO SQUARES (owner's link) = the p mod 4 split = Brouwer/Borsuk-Ulam (S18): p=a^2+b^2 <=> p=1 mod4
     <=> -1 is a QR mod p <=> Gauss sum g_p REAL (sqrt p) <=> complement=automorphism <=> BROUWER/easy.
     p=3 mod4 (NOT sum of two squares) <=> g_p IMAGINARY (i sqrt p) <=> free Z_2 <=> BORSUK-ULAM/hard.
 (3) FACILITY-LOCATION GAME: runners=facilities, lonely observer=farthest point, covering-min=adversary's
     min-over-configs meas(L_C). OPTIMUM (equidistributed/independent) = (6/7)^k; adversarial min = extremizer.
     POTENTIAL Phi(S) = (6/7)^k - meas(L_C) = the resonance/additive-energy the adversary induces (a Koksma-
     Hlawka discrepancy term). PoA = optimum/adversarial. Bounding Phi (covering constraint) => inf meas > 0.
"""
import sys, math
from fractions import Fraction as Fr
from math import gcd
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def isprime(n):
    if n<2: return False
    for i in range(2,int(n**.5)+1):
        if n%i==0: return False
    return True

print("="*98); print(" (1) THE NONABELIAN EXPANDER: PSL_2(F_p), Ramanujan gap 2*sqrt(p) = Gauss-sum sqrt(p) arithmetic"); print("="*98)
print(f"  {'p':>3} {'|PSL_2(F_p)|=p(p^2-1)/2':>22} {'Ramanujan gap 2sqrt(p)':>22} {'Gauss |g_p|=sqrt(p)':>18}")
for p in [3,5,7,11,13,23,31]:
    if not isprime(p): continue
    order=p*(p*p-1)//2
    tag=' <== apex7 = Aut(Fano)=GL_3(F_2)=168' if p==7 else ''
    print(f"  {p:>3} {order:>22} {2*math.sqrt(p):>22.4f} {math.sqrt(p):>18.4f}{tag}")
print("  => the LPS Ramanujan expander on PSL_2(F_p) has spectral gap 2sqrt(p); the certificate is g_p=i sqrt(p).")
print("     The expander's arithmetic = the Gauss-sum arithmetic. Apex-7: PSL_2(F_7)=Aut(Fano)=QR_7/octonion group.")

print("\n"+"="*98); print(" (2) SUM OF TWO SQUARES = p mod 4 = the Brouwer/Borsuk-Ulam (easy/hard) split"); print("="*98)
def sum2sq(p):  # is p a sum of two squares?
    for a in range(int(p**.5)+1):
        b2=p-a*a
        if b2>=0 and int(b2**.5)**2==b2: return (a,int(b2**.5))
    return None
print(f"  {'p':>3} {'p mod4':>6} {'sum of two squares?':>20} {'-1 QR mod p?':>12} {'Gauss g_p':>12} {'regime':>16}")
for p in [3,5,7,11,13,17,19,23,31]:
    if not isprime(p): continue
    s=sum2sq(p); qr=(pow(p-1,(p-1)//2,p)==1)
    gp='real sqrt%d'%p if p%4==1 else 'i*sqrt%d'%p
    reg='BROUWER/easy' if p%4==1 else 'BORSUK-ULAM/hard'
    print(f"  {p:>3} {p%4:>6} {str(s) if s else 'NO':>20} {str(qr):>12} {gp:>12} {reg:>16}")
print("  => p=1mod4 <=> sum of two squares <=> -1 is QR <=> Gauss real <=> BROUWER (easy). p=3mod4 <=> not a sum")
print("     of two squares <=> Gauss imaginary i sqrt p <=> free Z_2 <=> BORSUK-ULAM (hard). apex 7 = the hard side.")

print("\n"+"="*98); print(" (3) FACILITY-LOCATION GAME: potential Phi = (6/7)^k - meas; PoA = optimum/adversarial"); print("="*98)
k=13  # 13 speeds at n=14
opt=(Fr(6,7))**k
pent=Fr(313,9702)     # pentagon (Z/10)* extremizer (kps S8/S9)
spor=Fr(389,12012)    # sporadic (Z/19) two-clash
print(f"  OPTIMUM (equidistributed/independent, each runner safe w.p. 6/7): (6/7)^{k} = {float(opt):.6f}")
print(f"  ADVERSARIAL min (extremizer): pentagon meas = {float(pent):.6f}; sporadic = {float(spor):.6f}")
print(f"  POTENTIAL Phi(S) = (6/7)^k - meas(L_C): pentagon Phi = {float(opt-pent):.6f} (the resonance penalty)")
print(f"  PRICE OF ANARCHY = optimum / adversarial-min = {float(opt/pent):.3f}  (pentagon binds)")
print(f"  => meas(L_C) = (6/7)^k - Phi(S); Phi = the pairwise-resonance / additive-energy (Koksma-Hlawka")
print(f"     discrepancy) the covering adversary induces. inf meas = (6/7)^k / PoA. Bounding PoA <= ~4.2 (i.e.")
print(f"     Phi <= (6/7)^k(1-1/PoA)) gives inf meas >= (6/7)^k/PoA_max > 0 -- the potential-function lower bound.")
print(f"  target 1/36={float(Fr(1,36)):.6f}: extremizer/target = {float(pent*36):.3f}x (TIGHT); PoA to reach 1/36 = {float(opt*36):.3f}")

print("\n"+"="*98); print(" TANGENTIAL (owner's jumping-off): Pochhammer fiber, integer complexity, Hlawka-Koksma"); print("="*98)
# Pochhammer fiber fraction f(n) = (1/2)_{n-2} / (n-2)!
def pochhammer(a,m):
    r=Fr(1)
    for i in range(m): r*=(a+i)
    return r
n=14; f=pochhammer(Fr(1,2),n-2)/math.factorial(n-2)
print(f"  POCHHAMMER fiber fraction f({n}) = (1/2)_{n-2}/({n-2})! = {float(f):.6f} ~ 1/sqrt(pi*n)={1/math.sqrt(math.pi*n):.6f}")
def intcomplexity(N):  # Mahler-Popken integer complexity (min 1's with +,*)
    C={1:1}
    for m in range(2,N+1):
        best=m  # all ones
        for a in range(1,m):
            if a in C:
                if m-a in C: best=min(best,C[a]+C[m-a])
                if m%a==0 and m//a in C: best=min(best,C[a]+C[m//a])
        C[m]=best
    return C[N]
print(f"  INTEGER COMPLEXITY (Mahler-Popken): ||7||={intcomplexity(7)}, ||14||={intcomplexity(14)}, ||21||={intcomplexity(21)} "
      f"(21=3*7 => ||21||<=||3||+||7||); tangential -- 'how hard to build the apex from 1s'.")
print("  HLAWKA-KOKSMA: |meas(L_C) - (6/7)^k| <= V * D*  (variation x star-discrepancy of the runner phases) =")
print("   the rigorous form of Phi(S) above: the adversary MAXIMIZES discrepancy D*(resonance), the covering")
print("   constraint CAPS it => the PoA/potential bound is a Koksma-Hlawka discrepancy bound.")
print("DONE.")
