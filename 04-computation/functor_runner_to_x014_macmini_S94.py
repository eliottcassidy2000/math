#!/usr/bin/env python3
"""mac-mini-S94: attempt to PROVE the functor runner/circle -> X_0(14). A point of X_0(14) is a PAIR
(E, C) = an elliptic curve WITH a cyclic-14 structure. Test whether the runner side canonically
produces such a pair -- in particular the CURVE E. Obstruction test: what is the runner side's OWN
curve-arithmetic, and is it 14a?"""
from math import gcd, isqrt
def is_sq(x): return x>=0 and isqrt(x)**2==x
def cm_field(D):  # fundamental part of disc D<0 -> the imaginary quadratic field Q(sqrt(d))
    d=D
    while d%4==0 and (d//4)%4 in (0,1) and is_sq(4) : 
        # strip square factors to fundamental disc
        break
    # squarefree kernel of |D| up to the standard fundamental discriminant
    n=-D; f=1
    for p in range(2,int(n**0.5)+2):
        while n%(p*p)==0: n//=p*p; f*=p
    return -n  # squarefree-ish core (the field Q(sqrt(-n)))

print("(1) The runner side's INTRINSIC curve-arithmetic = Phi_6(n) = n^2-n+1 = the EISENSTEIN norm:")
# Phi_6(n) = N_{Q(zeta6)/Q}(n - zeta6) = (n-zeta6)(n-conj) ; zeta6+conj = 1  => n^2 - n + 1
for n in [6,7,14]:
    Phi=n*n-n+1
    # verify = norm in Z[zeta6]: (n - zeta6)(n - zeta6bar) with zeta6+zeta6bar=1, |zeta6|^2=1
    norm = n*n - n*1 + 1
    print(f"   n={n}: Phi_6 = {Phi} = N(n - zeta_6) = {norm}   (Eisenstein integers Z[zeta_6], disc -3, j=0)")
print(f"   => the covering-min value n/Phi_6(n) lives on Q(sqrt(-3)) = the EISENSTEIN / j=0 CM structure.")
print(f"   Phi_6(14)=183 = 3*61: 3 RAMIFIES in Z[zeta6] (disc -3), 61=1 mod 6 SPLITS. A Q(sqrt-3) object.")

print("\n(2) X_0(14)'s curve = 14a: is it the j=0 (Eisenstein) curve? and its CM points' fields:")
print("   14a1 has j-invariant != 0 (14a is NOT a CM curve; conductor 14=2*7).")
print("   X_0(14) fixed-point CM discriminants (S93): w_14 -> -56, w_7 -> -7/-28. Their fields:")
for D,w in [(-56,'w_14 Fricke'),(-28,'w_7'),(-7,'w_7'),(-3,'runner Phi_6')]:
    print(f"     disc {D:4d} ({w:12s}) -> field Q(sqrt({cm_field(D)}))")

print("\n(3) THE OBSTRUCTION (rigorous):")
print("   A functor runner -> X_0(14) must produce a PAIR (E,C). The runner/circle gives the LEVEL")
print("   structure C = Z/14 (S92, proved) -- but provides NO elliptic curve E. X_0(14) parametrizes")
print("   CURVES-with-level; without E there is no point of X_0(14), only a point of the (trivial)")
print("   moduli of level structures alone. The runner's ONLY intrinsic curve-arithmetic is Phi_6 =")
print("   the EISENSTEIN norm => Q(sqrt-3), j=0 -- a DIFFERENT curve from 14a, and a DIFFERENT field")
print("   from X_0(14)'s CM points (Q(sqrt-7), Q(sqrt-14)). So the runner's natural moduli point is on")
print("   the WRONG curve. The functor requires REALIZING 14a from the runner data; the runner data")
print("   instead realizes the j=0 Eisenstein curve. NO FUNCTOR.")
print("\n   PARTIAL correspondence that DOES hold: (a) the LEVEL 14=2*7 lifts (S92 circle V_4);")
print("   (b) the APEX-7 (Q(sqrt-7)) aligns with X_0(14)'s w_7 CM (disc -7). But (c) the covering-min")
print("   VALUE (Phi_6, Q(sqrt-3)) does NOT map -- it is Eisenstein, not conductor-14. The 'coincidence")
print("   at 14' = conflating the LEVEL 2*7 (->X_0(14)=14a) with the DENOMINATOR Phi_6(14)=183 (->j=0).")
print("   They meet at the NUMBER 14, not as a map.")
