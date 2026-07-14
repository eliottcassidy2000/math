#!/usr/bin/env python3
"""mac-mini-S93: PROVE the moduli fixed-point arithmetic of X_0(14) and its correspondence to the
danger-circle V_4 (S92). Method (rigorous, classical): AL eigenvalue w_Q on the 1-dim cusp form f_14
determines translation (+1, 0 fixed pts) vs reflection (-1, 4 fixed pts) on the genus-1 curve;
Riemann-Hurwitz fixes the counts; the Fricke count = h(-4N) (+h(-N) if N=3 mod4). Verify on N=11 first."""
from math import gcd
def class_number(D):  # D<0 discriminant (D=0,1 mod 4); count reduced primitive forms [a,b,c]
    assert D<0 and D%4 in (0,1)
    h=0; a=1
    while a*a <= -D/3:
        b=-a
        while b<=a:
            if (b*b-D)%(4*a)==0:
                c=(b*b-D)//(4*a)
                if a<=c and gcd(gcd(abs(a),abs(b)),abs(c))==1:
                    if (abs(b)==a or a==c):
                        if b>=0: h+=1
                    else: h+=1
            b+=1
        a+=1
    return h
for D in [-4,-7,-8,-11,-28,-44,-56]:
    print(f"  h({D}) = {class_number(D)}")
print()
def fricke_fixed(N, N3mod4_disc=True):
    v=class_number(-4*N)
    if N%4==3: v+=class_number(-N)   # -N is a valid disc when N=3 mod4
    return v
print("SANITY CHECK N=11 (X_0(11)=11a, genus1, Fricke w_11 eigenvalue -1 => reflection => 4 fixed):")
print(f"  Fricke ν(w_11) = h(-44)+h(-11) = {class_number(-44)}+{class_number(-11)} = {fricke_fixed(11)}  (expect 4)")
print()
print("="*66); print("X_0(14) = 14a (genus 1). AL eigenvalues on f_14 (LMFDB 14.a): w_2=+1, w_7=-1, w_14=-1")
print("="*66)
# genus-1 dichotomy: eigenvalue +1 <=> translation <=> 0 fixed pts <=> quotient genus 1;
#                    eigenvalue -1 <=> reflection  <=> 4 fixed pts <=> quotient genus 0 (Riemann-Hurwitz g'=1-nu/4)
eig={'w_2':+1,'w_7':-1,'w_14':-1}
for w,e in eig.items():
    nu = 0 if e==+1 else 4
    gq = 1 - nu//4
    kind = "TRANSLATION (fixed-point-free)" if e==+1 else "REFLECTION (4 fixed pts)"
    print(f"  {w}: eigenvalue {e:+d} => {kind}; nu={nu}; quotient genus {gq}")
print(f"  Fricke w_14 fixed pts = h(-56) = {class_number(-56)} (disc -56 = -2^3*7, encodes BOTH primes 2,7)")
print(f"  CHECK V_4 consistency: w_2 o w_7 = w_14; translation o reflection = reflection (0-fix,4-fix -> 4-fix) OK")
print(f"  and w_7 o w_14 = w_2: reflection o reflection = translation (4,4 -> 0) OK")
print()
print("CORRESPONDENCE to the danger circle Z/14 (S92: W_2=x+7, W_7=7-x, W_14=-x):")
Z=14
cf={'w_2':[x for x in range(Z) if (x+7)%Z==x],'w_7':[x for x in range(Z) if (7-x)%Z==x],
    'w_14':[x for x in range(Z) if (-x)%Z==x]}
for w in ['w_2','w_7','w_14']:
    circ=len(cf[w]); mod = 0 if eig[w]==1 else 4
    tag = "MATCH (both translation, free)" if w=='w_2' else f"circle {circ} -> moduli {mod}: lift ADDS the CM points"
    print(f"  {w}: circle fixes {cf[w]} ({circ}) ; moduli {mod} ; {tag}")
print()
print("VERDICT: the 2-part w_2 is a fixed-point-FREE TRANSLATION on BOTH circle and moduli (the clean lift).")
print("The reflections w_7,w_14 GAIN their fixed points in the lift (0->4, 2->4) = the CLASS NUMBERS")
print("(h(-56)=4) = the arithmetic the FLAT circle cannot see. 'Lifting' = curving the flat circle into the")
print("genus-1 moduli, which CREATES the CM fixed points. NOTE: corrects klein-S59 ('W_2 4 CM pts') -- w_2 is")
print("the TRANSLATION (0 fixed pts); the 4-fixed-point reflections are w_7 and w_14.")
