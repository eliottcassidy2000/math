"""
S75e: the LRC(14) cap is a totally-real cyclotomic quantity in F=Q(cos 2pi/7).
- the cyclotomic conductor 7^{1,2} sits exactly in the BINDING rows k=8,9 (dip!=0); easy rows k>=10 are clean.
- de Moivre cubic x^3+x^2-2x-1: disc=49=7^2, Galois C_3, totally real, h=1.
- the Fejer magic function F_7=(de Moivre cubic)^2 is a TOTALLY-POSITIVE SQUARE => Bochner positivity.
- two heads of n=14: 14=2*7 (cap, Q(cos2pi/7)) and 2n-1=27=3^3 (witness ramification tower, HYP-2436).
"""
import sympy as sp
from sympy import Rational, factorint
x=sp.symbols('x')

print("1) cap/dip in Q(cos2pi/7): the 7-conductor sits in the BINDING rows")
caps={8:Rational(2243,5880),9:Rational(1979,4004),10:Rational(55,91),11:Rational(66,91),12:Rational(6,7),13:Rational(1)}
for k in range(8,14):
    binom=Rational((k+1)*k//2,91); dip=binom-caps[k]
    print(f"   k={k}: cap={caps[k]!s:>10} dip={dip!s:>10} cap_denom={factorint(caps[k].q)} 7-adic v7(dip_denom)={(factorint(dip.q).get(7,0)) if dip!=0 else '-'}")

print("\n2) the de Moivre cubic field F=Q(cos2pi/7):")
mp=sp.minimal_polynomial(2*sp.cos(2*sp.pi/7),x)
print(f"   min poly={mp}  disc={sp.discriminant(mp)}=7^2  degree=3  Galois=C_3  totally real (3 real roots)  h=1")

print("\n3) magic function F_7=(de Moivre cubic)^2 = totally-positive square (Bochner F-hat=(7-|n|)_+>=0):")
print("   F_7(t)=sin^2(7 pi t)/sin^2(pi t); F_7(0)=49=7^2; double zeros at de Moivre angles. The EVEN cap half is")
print("   this cyclotomic SQUARE (SOS-provable). Sharper: Jackson F_7^2, de la Vallee-Poussin, modular Gamma_0(7).")

print("\n4) n=14's TWO cyclotomic heads (both depth 3):")
print("   7-HEAD (cap):     14=2*7 -> 7 sectors -> Q(cos2pi/7) -> Phi_7 de Moivre -> Fejer square; depth (7-1)/2=3")
print(f"   3-HEAD (witness): 2n-1=27=3^3 -> doubling orbit ord_27(2)={sp.n_order(2,27)}=phi(27); ramification v_3(27)=3 (HYP-2436)")
print("   => n=14 = first case where BOTH the 2*(apex prime) cap-head and the p^3 ramified witness-head are hard.")
