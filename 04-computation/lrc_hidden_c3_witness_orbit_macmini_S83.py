"""
S83: the hidden C_3. The LRC(14) witness space is a SINGLE C_3-Galois orbit, and every '3' is this one C_3.
The 3 binding-pair witnesses {1,13},{3,11},{5,9} are cycled by x3 mod 14 (one C_3 orbit).
C_3 = (Z/14)*/{+-1} = (Z/7)*/{+-1} = Gal(Q(cos2pi/7)/Q) = #de Moivre angles = index (p-1)/2 = #QR mod 7.
LEAP: run the proof C_3-equivariantly -- HYP-2909 (one binding pair, PROVED) + C_3 generates the equioscillation;
cap = the C_3-trace (rational); rigidity & equidistribution = C_3-equivariant residuals on Q(cos2pi/7).
"""
from math import gcd
import sympy as sp

units=[a for a in range(1,14) if gcd(a,14)==1]
print(f"(Z/14)* = {units} = C_6 (gen 3: {[pow(3,k,14) for k in range(6)]})")
# the 3 binding pairs as a single C_3 orbit under x3
cur=frozenset([1,13]); orbit=[]
for _ in range(3):
    orbit.append(tuple(sorted(cur))); cur=frozenset([(x*3)%14 for x in cur])
target={frozenset([1,13]),frozenset([3,11]),frozenset([5,9])}
print(f"x3 orbit of {{1,13}}: {orbit}  single C_3 orbit = the 3 binding pairs? {set(map(frozenset,orbit))==target}")
print()
print("the SAME C_3, six ways (all = 3):")
print(f"  |(Z/14)*/{{+-1}}| = {len(units)//2}")
qr7=sorted({pow(a,2,7) for a in range(1,7)}); print(f"  #QR mod 7 = {len(qr7)} (QRs={qr7})")
x=sp.symbols('x'); mp=sp.minimal_polynomial(2*sp.cos(2*sp.pi/7),x)
print(f"  de Moivre cubic = {mp}, degree {sp.degree(mp)} = |Gal(Q(cos2pi/7))| = C_3")
print(f"  index (p-1)/2 for p=7 = {(7-1)//2}")
print()
print("LEAP: HYP-2909 (one binding pair, PROVED+Lean) = ONE C_3-orbit point; C_3 (x3) generates the other 2")
print("  => the full equioscillation at the 6 units (kps S255) from the proved single pair + the C_3 symmetry.")
print("  cap = C_3-trace in Q(cos2pi/7) (rational, disc 7^2, S75e). Rigidity & equidist = C_3-equivariant residuals.")
print("  DIRECTION: run the proof on the C_3 witness space / Q(cos2pi/7), not the 13-dim config space.")
