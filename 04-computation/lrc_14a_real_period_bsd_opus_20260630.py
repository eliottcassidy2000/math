"""
Periods of 14a (minimal model y^2+xy+y=x^3+4x-6). Complete the square: (2y+x+1)^2 = 4x^3+x^2+18x-23
= (x-1)(4x^2+5x+23); disc(4x^2+5x+23) = -343 = -7^3 (the apex-7! Q(sqrt-7)). Real period Omega+ =
2 int_1^oo dx/sqrt((x-1)(4x^2+5x+23)) = 4 int_0^oo ds/sqrt(4s^4+13s^2+32). BSD-check vs L(1)=0.330224.
"""
import math
def integrand(s): return 1.0/math.sqrt(4*s**4+13*s**2+32)
# integrate 0..inf via substitution s=tan(theta) or just adaptive to large bound
def simpson(f,a,b,n):
    if n%2: n+=1
    h=(b-a)/n; tot=f(a)+f(b)
    for i in range(1,n):
        tot+= (4 if i%2 else 2)*f(a+i*h)
    return tot*h/3
# tail beyond B ~ int 1/(2 s^2) = 1/(2B)
B=2000.0
I=simpson(integrand,0,B,2000000)+1.0/(2*B)
Omega_plus=4*I
print(f"real period Omega+ (14a) = {Omega_plus:.6f}")
L1=0.330224
print(f"L(14a,1) = {L1}")
print(f"L(1)/Omega+ = {L1/Omega_plus:.6f}   (BSD: = #Sha . prod c_p / #tors^2; #tors=6 => /36)")
for cp in [1,2,3,4,6,12]:
    val=cp/36
    if abs(L1/Omega_plus - val)<0.01:
        print(f"   matches prod c_p . #Sha = {cp}: L(1)/Omega+ = {cp}/36 = {val:.6f}  => Omega+ = 36 L(1)/{cp} = {36*L1/cp:.5f}")
# imaginary period: lattice covolume; for Delta<0 the second period is complex. Area = Omega+ * Omega_imag.
# Omega_imag via the complementary integral 2 int over the complex/oval. Use tau = Omega2/Omega1 from AGM not needed;
# estimate area via L-function / modular: covolume = 2*Im. Use Omega_minus from c_oo: real curve Delta<0 has c_oo=1.
print()
# search 0.040 among period quantities
floor_cusppart=114382/332563 - 3/math.pi**2
cands={"Omega+":Omega_plus,"1/Omega+":1/Omega_plus,"Omega+/(2pi)":Omega_plus/(2*math.pi),
       "1/(2 Omega+)":1/(2*Omega_plus),"L(1)":L1,"L(1)/Omega+":L1/Omega_plus,
       "1/Omega+^2":1/Omega_plus**2,"Omega+/N":Omega_plus/14,"1/(Omega+ pi)":1/(Omega_plus*math.pi)}
print(f"floor cusp-part target = {floor_cusppart:.5f}")
for k,v in cands.items():
    flag=" <== close" if abs(v-floor_cusppart)<0.005 else ""
    print(f"   {k:>14} = {v:.5f}{flag}")
print()
print(f"note: disc of the quadratic factor = -343 = -7^3 -> Q(sqrt-7) (the MEASURE column field) sits in")
print(f"      the period geometry of 14a, even though 14a is NOT CM. The apex-7 is in the period lattice.")
