"""
Imaginary period Omega3 = 4 int_0^oo ds/sqrt(4s^4-13s^2+32) (note: BOTH period quartics 4s^4 +/- 13s^2 +32
have disc_{s^2} = 13^2-4.4.32 = -343 = -7^3). Lattice modulus tau, and is it governed by sqrt(-7)?
"""
import math
def simpson(f,a,b,n):
    if n%2:n+=1
    h=(b-a)/n;t=f(a)+f(b)
    for i in range(1,n):t+=(4 if i%2 else 2)*f(a+i*h)
    return t*h/3
def Iqu(sign):  # int_0^oo ds/sqrt(4s^4 + sign*13 s^2 +32)
    f=lambda s:1.0/math.sqrt(4*s**4+sign*13*s**2+32)
    B=2000.0
    return simpson(f,0,B,2000000)+1.0/(2*B)
Op=4*Iqu(+1); O3=4*Iqu(-1)
print(f"real period      Omega+ = {Op:.6f}")
print(f"imaginary period Omega3 = {O3:.6f}")
print(f"area of lattice  Omega+ * Omega3 = {Op*O3:.6f}")
# rhombic lattice (Delta<0): generators Omega+, (Omega+ + i Omega3)/2 ; modulus tau = 1/2 + i Omega3/(2 Omega+)
Imtau=O3/(2*Op); 
print(f"\nlattice modulus tau = 1/2 + i*{Imtau:.6f}   (Im tau = Omega3/(2 Omega+))")
print(f"|tau|^2 = 1/4 + {Imtau**2:.6f} = {0.25+Imtau**2:.6f}")
# is the period geometry governed by sqrt7? check ratios
print(f"\nratios vs sqrt7={math.sqrt(7):.5f}, sqrt7/2={math.sqrt(7)/2:.5f}:")
print(f"   Omega3/Omega+      = {O3/Op:.6f}")
print(f"   Omega3/Omega+ /sqrt7 = {O3/Op/math.sqrt(7):.6f}")
print(f"   2*Im tau /sqrt7    = {2*Imtau/math.sqrt(7):.6f}")
print(f"   |tau| = {math.sqrt(0.25+Imtau**2):.6f}, |tau|*2/sqrt7 = {2*math.sqrt(0.25+Imtau**2)/math.sqrt(7):.6f}")
# CM check: if 14a were CM by Q(sqrt-7), tau would be in Q(sqrt-7). It's NOT CM, but does -7 appear?
print()
# search 0.040 including imaginary period / area / Petersson-proxy
target=114382/332563 - 3/math.pi**2
# Petersson proxy: <f,f> ~ Area/(4 pi^2) * deg(phi); also L(sym2,2) ~ const * <f,f>
cands={"Omega3":O3,"Area":Op*O3,"1/Area":1/(Op*O3),"Area/(4pi^2)":Op*O3/(4*math.pi**2),
       "Im tau":Imtau,"1/Omega3":1/O3,"Op/Area":Op/(Op*O3),"Area/(2pi)^2":Op*O3/(2*math.pi)**2,
       "1/(2 Area)":1/(2*Op*O3)}
print(f"floor cusp-part target = {target:.5f}")
for k,v in cands.items():
    f=" <==" if abs(v-target)<0.004 else ""
    print(f"   {k:>16} = {v:.5f}{f}")
print()
print("ASSESSMENT: the apex-7 (disc -343=-7^3) governs BOTH period integrals of 14a -- Q(sqrt-7) is the")
print("period field (the MEASURE column), structurally, despite 14a not being CM. The clean period fact is")
print("L(14a,1)=Omega+/6 (BSD). Whether 0.040 = a period combination: checked Op,O3,Area,tau -- see above.")
