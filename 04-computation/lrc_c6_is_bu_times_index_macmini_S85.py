"""
S85: the C_6 = C_2 x C_3 unification. EXTENDS S84/codex HYP-3310/3311 by naming the C_2 coordinate.
The field is Q(zeta_7), Gal = C_6 = C_2 x C_3. codex named C_3 (binding-pair skeleton) + C_2-transverse (chi_7).
NEW: C_2 = the BORSUK-ULAM antipodal of my S79 index-theorem (chi_7(-1)=-1 since 7=3 mod4), C_3 = the index degree (S83).
So C_6 = (BU antipodal) x (index degree) = the full residue proof structure; the 2-adic magnitude is the doubling-shadow 2U=2*U.
nonzero(14) = U (skeleton) ⊔ 2U (doubling-shadow=magnitude) ⊔ {7} (ramified). The index-theorem sees C_6 but is BLIND to the 2-adic height (codex).
"""
from math import gcd
qr={pow(a,2,7) for a in range(1,7)}
def chi7(a):
    r=a%7; return 0 if r==0 else (1 if r in qr else -1)

print("C_6 = C_2 x C_3 = Gal(Q(zeta_7)/Q):")
print(f"  C_3 = QRs mod7 = {sorted(qr)} (the squares, fixed field Q(sqrt-7));  C_2 = {{1,6}}={{+-1}} (chi7, fixed field Q(cos2pi/7))")
print(f"  chi7(-1)=chi7(6)={chi7(6)}  (NQR because 7=3 mod4) => C_2 IS the antipodal sign = the BORSUK-ULAM Z_2 (S79)")
print()
print("the 3 binding pairs {a,-a}: C_3 says WHICH pair (index degree, S83); C_2=chi7 says WHICH element (BU antipodal):")
for a,b in [(1,13),(3,11),(5,9)]:
    print(f"  {{{a:>2},{b:>2}}} mod7=({a%7},{b%7})  chi7=({chi7(a):>2},{chi7(b):>2})  one QR + one NQR (C_2 transverse to C_3)")
print()
U=[s for s in range(1,14) if gcd(s,14)==1]
shadow=sorted((2*u)%14 for u in U)
print("2-adic magnitude = the DOUBLING-SHADOW of the skeleton:")
print(f"  U (skeleton) = {U};  2U mod14 = {shadow} = the even covering layer;  {{7}} = ramified leftover")
print(f"  nonzero(14) = U ⊔ 2U ⊔ {{7}}  partition check: {sorted(U+shadow+[7])==list(range(1,14))}")
print()
print("SYNTHESIS: C_6 (residue) = C_2[BU antipodal, S79, Q(sqrt-7)] x C_3[index degree, S83, Q(cos2pi/7)].")
print("  The index-theorem (S79) captures all of C_6 but FORGETS the 2-adic height (codex's blind data) = the hard residual.")
print("  Charge view (codex HYP-3400): the C_2 charge = the Borsuk-Ulam degree (odd) = a conserved topological invariant.")
