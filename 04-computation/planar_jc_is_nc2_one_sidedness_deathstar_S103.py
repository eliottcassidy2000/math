# HISTORICAL / PARTIAL: read MISTAKE-237 and HYP-8905.  This script illustrates
# only the binary homogeneous family P=A z^d+B zbar^d.  It proves no
# NC2/GMC/JC implication.  THM-1300 refutes JC(n) for n>=3; JC(2) remains open.
import numpy as np
I=1j
def P(x,y,A,B,d): return A*(x+I*y)**d + B*(x-I*y)**d
def hess(f,x0,y0,h=1e-4):
    fxx=(f(x0+h,y0)-2*f(x0,y0)+f(x0-h,y0))/h**2
    fyy=(f(x0,y0+h)-2*f(x0,y0)+f(x0,y0-h))/h**2
    fxy=(f(x0+h,y0+h)-f(x0+h,y0-h)-f(x0-h,y0+h)+f(x0-h,y0-h))/(4*h**2)
    return np.array([[fxx,fxy],[fxy,fyy]])

print("HISTORICAL / PARTIAL: read MISTAKE-237 and HYP-8905; only the binary homogeneous formula survives.")
print("="*70)
print("[1] Binary family P=A z^d+B zbar^d: nilpotent Hessian iff A*B=0")
print("="*70)
x0,y0=0.37+0.11j, -0.19+0.23j
for d in [3,4,5]:
    for (A,B,tag) in [(1,0,"holomorphic (one-sided)"),(0,1,"antiholo (one-sided)"),(1,1,"TWO-sided")]:
        H=hess(lambda x,y:P(x,y,A,B,d),x0,y0)
        tr=H[0,0]+H[1,1]; det=H[0,0]*H[1,1]-H[0,1]**2
        nil = abs(tr)<1e-2 and abs(det)<1e-2
        print(f"  d={d} A={A} B={B} [{tag:24s}]: trace={abs(tr):.2e} det={abs(det):.2e} -> nilpotent={nil}")
print("  => inside this binary family only, nilpotence occurs exactly when A*B=0.")

print("\n"+"="*70)
print("[2] one-sided map F=(x,y)+grad P is ONE-FIBER-LINEAR (F2 - i F1 = -i z), hence TAME/invertible")
print("="*70)
d=4;A=1.0
def gradP(x,y): 
    h=1e-6
    return ((P(x+h,y,A,0,d)-P(x-h,y,A,0,d))/(2*h), (P(x,y+h,A,0,d)-P(x,y-h,A,0,d))/(2*h))
for (x0,y0) in [(0.3+0.2j,-0.1+0.4j),(1.0,0.5j)]:
    gx,gy=gradP(x0,y0); F1,F2=x0+gx,y0+gy
    combo=F2-I*F1; z=x0+I*y0
    print(f"  point z={z:.3f}: F2 - i F1 = {combo:.4f}  vs  -i z = {(-I*z):.4f}  (match => linear pencil member)")

print("\n"+"="*70)
print("[3] Zhao VC / moment: one-sided P is HARMONIC (holomorphic), Delta=0, so Delta^m(P^m)=0 trivially")
print("="*70)
def lap(f,x0,y0,h=1e-4): return (f(x0+h,y0)-2*f(x0,y0)+f(x0-h,y0))/h**2 + (f(x0,y0+h)-2*f(x0,y0)+f(x0,y0-h))/h**2
for d in [3,4]:
    val=lap(lambda x,y:(x+I*y)**d, 0.3+0.1j,-0.2+0.15j)
    print(f"  Delta((x+iy)^{d}) = {abs(val):.2e}  (=0: holomorphic is harmonic => Delta^m(P^m)=0, Zhao VC trivial one-sided)")
print("\n"+"="*70)
print("AUDIT BOUNDARY")
print("="*70)
print("  The exact symbolic proof is in HYP-8905; THM-2063 gives tameness for the resulting one-fiber-linear map.")
print("  The general symmetric reduction for JC(2) changes dimension; it is not this binary calculation.")
print("  Rank analogies and moment syntax supply no NC2/GMC-to-JC map. JC(2) remains open.")
