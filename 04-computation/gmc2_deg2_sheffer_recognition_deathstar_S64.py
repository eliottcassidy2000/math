from fractions import Fraction as Fr
from math import factorial
import cmath
# E[P^m]=E_r[L_m(alpha,beta)], alpha=a*c*r, beta=b, L_m=sum_k m!/(k!^2(m-2k)!)alpha^k beta^{m-2k}, E_r[r^j]=j!
def pmul(p,q):
    d={}
    for i,u in enumerate(p):
        for j,v in enumerate(q):
            if u and v: d[i+j]=d.get(i+j,0)+u*v
    n=max(d)+1 if d else 1; o=[0]*n
    for k,v in d.items(): o[k]=v
    return o
def padd(*ps):
    n=max(len(p) for p in ps); o=[0]*n
    for p in ps:
        for i,u in enumerate(p): o[i]+=u
    return o
def Er(p,fac): return sum(c*fac(j) for j,c in enumerate(p))
def EPm(a,b,c,M,fac=factorial):
    ac=pmul(a,c); alpha=[0]+ac; res=[]
    for m in range(1,M+1):
        tot=[0]
        for k in range(m//2+1):
            coef=Fr(factorial(m),factorial(k)**2*factorial(m-2*k))
            term=[1]
            for _ in range(k): term=pmul(term,alpha)
            for _ in range(m-2*k): term=pmul(term,b)
            tot=padd(tot,[coef*u for u in term])
        res.append(Er(tot,fac))
    return res

print("=== (1) REAL case: E[P^2]=E_r[b^2]+2ac > 0 for real a=c, ANY b (positivity, all degrees) ===")
for bb in ([Fr(1),Fr(0),Fr(1)],[Fr(0),Fr(0),Fr(1)],[Fr(-2),Fr(3),Fr(1)]):
    E=EPm([Fr(1)],bb,[Fr(1)],2)
    print(f"  b={[str(x) for x in bb]}, a=c=1: E[P^2]={E[1]} (>0: {E[1]>0})  -> killed by positivity regardless of deg")

print("\n=== (2) Does the NAIVE curve E[P^m]=E_r[s^m He_m(b/s)] hold? (expect NO for non-const b) ===")
def He(m,x):
    if m==0: return 1
    h0,h1=1,x
    for k in range(1,m): h0,h1=h1,x*h1-k*h0
    return h1
# compare E[P^m] vs E_r[s(r)^m He_m(b(r)/s(r))], s^2=-2*acr, via series in r then E_r (numerically hard w/ sqrt; do float sampling)
import numpy as np
a=c=1.0; b=[0.0,0.0,1.0]  # b=r^2
def bval(r): return b[0]+b[1]*r+b[2]*r*r
E_true=[float(x) for x in EPm([Fr(1)],[Fr(0),Fr(0),Fr(1)],[Fr(1)],6)]
# E_r[s^m He_m(b/s)] by quadrature, s^2=-2*a*c*r => s=sqrt(-2r) imaginary
import numpy as np
rs=np.linspace(1e-6,60,400000); w=np.exp(-rs)
naive=[]
for m in range(1,7):
    s=np.sqrt(-2*a*c*rs+0j)
    val=np.trapz((s**m)*np.array([He(m,complex(b[0]+b[1]*r+b[2]*r*r)/complex(sc)) for r,sc in zip(rs,s)])*w, rs)
    naive.append(val)
print("  E[P^m] (true) m=1..6:", [round(x,2) for x in E_true])
print("  E_r[s^m He_m(b/s)]   :", [complex(round(x.real,2),round(x.imag,2)) for x in naive])
print("  => naive-curve equals true?", all(abs(E_true[i]-naive[i])<1e-3 for i in range(6)))

print("\n=== (3) The EXACT Sheffer form: E[P^m]=sum_k m!/(k!(m-2k)!) gamma^k E_{Gamma(k+1)}[b^{m-2k}] ===")
def EGamma(kp1,poly):  # E_{Gamma(kp1)}[poly], moments r^i -> (kp1-1+i)!/(kp1-1)! = rising factorial
    k=kp1-1
    return sum(c*Fr(factorial(k+i),factorial(k)) for i,c in enumerate(poly))
def bpow(b,j):
    r=[1]
    for _ in range(j): r=pmul(r,b)
    return r
def EPm_sheffer(gamma,b,M):
    res=[]
    for m in range(1,M+1):
        tot=Fr(0)
        for k in range(m//2+1):
            tot+=Fr(factorial(m),factorial(k)*factorial(m-2*k))*gamma**k*EGamma(k+1,bpow(b,m-2*k))
        res.append(tot)
    return res
E1=EPm([Fr(1)],[Fr(0),Fr(1),Fr(1)],[Fr(1)],8)          # a=c=1(gamma=1), b=r+r^2
E2=EPm_sheffer(Fr(1),[Fr(0),Fr(1),Fr(1)],8)
print("  E[P^m] direct  :", [str(x) for x in E1])
print("  E[P^m] Sheffer :", [str(x) for x in E2])
print("  Sheffer form matches:", E1==E2)
