from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps=40
def ek_from_roots(roots):
    # elementary symmetric e_0..e_d exactly from a list of rational roots
    e=[Fr(1)]
    for r in roots:
        ne=e+[Fr(0)]
        for i in range(len(e),0,-1): ne[i]=e[i-1]*r + (e[i] if i<len(e) else Fr(0))
        ne[0]=e[0]
        e=ne
    return e
def binom(n,k):
    from math import comb; return comb(n,k)
def analyze(roots, ks):
    d=len(roots); e=ek_from_roots(roots)
    h=[Fr(e[k],binom(d,k)) for k in range(d+1)]
    # moments of root measure
    m1=sum(roots)/d; m2=sum(r*r for r in roots)/d; m3=sum(r**3 for r in roots)/d
    x=float(m2/m1**2); z=float(m3/m1**3); curv=3*x*x-2*z-1
    out={}
    for k in ks:
        if k+1>d or k-2<0: continue
        Rk=Fr(h[k]**2,h[k-1]*h[k+1]); Rk1=Fr(h[k-1]**2,h[k-2]*h[k])
        val=float(mp.log(mp.mpf(Rk.numerator)/Rk.denominator)-mp.log(mp.mpf(Rk1.numerator)/Rk1.denominator))
        out[k]=val*d*d
    return curv,out

print("THM-3000 INDEPENDENT VERIFY + MOVING-EDGE PROBE (real-rooted, exact e_k, curvature = 3x^2-2z-1)")
print("="*92)
# measure: roots = arithmetic 1..d
for d in [200,400,800]:
    roots=[Fr(i) for i in range(1,d+1)]
    curv,out=analyze(roots,[2,3,4, d//8, d//4, d//2, (3*d)//4, d-2])
    print(f"\narithmetic roots 1..{d}:  curvature 3x^2-2z-1 = {curv:.5f}")
    print("   k       d^2 log(R_k/R_k-1)   k/d")
    for k in sorted(out):
        print(f"   {k:5d}   {out[k]:18.5f}   {k/d:.3f}   {'<- fixed-k universal' if k<=4 else '<- MOVING EDGE'}")
print("\n"+"="*92)
print("BOLD REFRAME PROBE: AMM-12592 certificate biases as a 2-point root measure")
# p_A=1285/2181, p_B=8847357/11821757 ; try root measures reflecting the biases
for name,(a,b,wa) in [("bias p_A=1285/2181 -> roots {q_A,p_A}",(Fr(896,2181),Fr(1285,2181),Fr(1,2))),
                      ("bias p_B two-point",(Fr(2974400,11821757),Fr(8847357,11821757),Fr(1,2)))]:
    m1=wa*a+(1-wa)*b; m2=wa*a*a+(1-wa)*b*b; m3=wa*a**3+(1-wa)*b**3
    x=float(m2/m1**2); z=float(m3/m1**3); curv=3*x*x-2*z-1
    print(f"  {name}: curvature = {curv:.6f}")
print(f"  targets to spot: 2457/6592={2457/6592:.6f}, 1/25=0.04, gamma-ish; log_B/log_A~3.02")
