#!/usr/bin/env python3
"""mac-mini-S96: does the mod-6 Eisenstein congruence (S95) FORCE the covering-min value 14/183 (=n/Phi_6(n))?
HONEST TEST. The value is PROVED by THM-724 (metric/Farey, NO modular input). Test what the congruence
actually shares with the value: decompose the covering-min Dedekind sum into the E_2 anomaly (Eisenstein
backbone, forced by the mod-6 structure) + the finite Farey correction (NOT forced by the congruence)."""
from fractions import Fraction as F
n=14; Phi=n*n-n+1  # 183
print(f"covering-min value M = n/Phi_6(n) = {n}/{Phi} = {F(n,Phi)}  (=14/183, PROVED by THM-724, metric)")
print(f"margin = M - 1/n = {F(n,Phi)-F(1,n)} = (n-1)/(n*Phi) = {F(n-1,n*Phi)}")
# HYP-3768: margin = -12 s(n,Phi)/n^2 ; s(n,Phi) = -(Phi-1)/(12 Phi)
s = F(-(Phi-1),12*Phi)
print(f"\nDedekind sum s(n,Phi_6) = -(Phi-1)/(12 Phi) = {s} = {float(s):.6f}   (HYP-3768, proved by reciprocity)")
# DECOMPOSITION: s = -1/12 + 1/(12 Phi)
backbone = F(-1,12); correction = F(1,12*Phi)
print(f"DECOMPOSE:  s = -1/12 + 1/(12*Phi) = {backbone} + {correction} = {backbone+correction}  (equal: {backbone+correction==s})")
print(f"   -1/12 = -B_2/2 = the E_2 EISENSTEIN ANOMALY (n->infinity limit)  <-- the Eisenstein 'backbone'")
print(f"   +1/(12*Phi) = +1/{12*Phi} = the FINITE FAREY correction        <-- NOT the Eisenstein limit")
print(f"\nSo the covering-min value = 1/n + margin, margin = -12s/n^2 = -12(-1/12 + 1/(12Phi))/n^2:")
m_from_s = -12*s/F(n*n)
print(f"   margin = {m_from_s} = (Phi-1)/(n^2 Phi); matches (n-1)/(n Phi): {m_from_s==F(n-1,n*Phi)}")
print()
print("WHAT THE mod-6 CONGRUENCE FORCES vs NOT:")
print("  FORCES (shared): the E_2 EISENSTEIN backbone. f_14 = E_2 mod 6 (S95); s(n,Phi_6) -> -1/12 = the E_2")
print("    anomaly (the n->inf limit). Both are E_2 shadows. The '-1/12' is the Eisenstein leading term.")
print("  DOES NOT FORCE: the FINITE VALUE 14/183. The correction +1/(12*Phi_6) = the Farey/reciprocity data")
print("    (Phi_6 = 1 mod n, CF [0;n-1,n]), PROVED by THM-724 with ZERO modular input. The congruence is about")
print("    f_14's coefficients mod 6; the value is metric. Consistent with S94 (no functor).")
print()
print("THE TWO 'SIXES' ARE INDEPENDENT (no common cause):")
# 6a: order of n mod Phi_6(n) (why Phi_6 = 6th cyclotomic appears in the value)
def mult_order(a,m):
    a%=m; k=1; x=a%m
    while x!=1:
        x=(x*a)%m; k+=1
        if k>m: return None
    return k
print(f"  6a = ord({n} mod {Phi}) = {mult_order(n,Phi)}  (n^3={n**3}=-1 mod {Phi}; by CONSTRUCTION Phi_6 makes n order 6)")
print(f"  6b = torsion of 14a = Z/6  (a fact about the CURVE, giving the mod-6 congruence)")
print(f"  Same integer 6, DIFFERENT origins: 6a is the covering CF/Farey construction; 6b is the curve's arithmetic.")
print(f"  No rigorous common cause => the congruence (6b) does not force the value (6a/Farey).")
print()
print("VERDICT: the mod-6 congruence does NOT force the covering-min value. It forces the EISENSTEIN BACKBONE")
print("(the -1/12 E_2 anomaly = the covering-min's n->inf limit); the finite value 14/183 is FAREY-forced")
print("(THM-724), congruence-independent. 'Linked at the backbone' (S95), NOT 'forcing the value'.")
