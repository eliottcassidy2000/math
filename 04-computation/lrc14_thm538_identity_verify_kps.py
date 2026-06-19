import sys, itertools, cmath
from cmath import exp, pi, sin
from fractions import Fraction
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def chat(T, n):
    if n == 0: return complex(1 - len(T)/7.0, 0)
    s = sin(pi*n/7)/(pi*n)
    return -sum(cmath.exp(-2j*pi*n*(j+0.5)/7)*s for j in T)
def K(nvec):
    tot=0+0j
    for r in range(0,7):
        for T in itertools.combinations(range(1,7), r):
            p=1+0j
            for nj in nvec: p*=chat(set(T), nj)
            tot+=((-1)**r)*p
    return tot
def meas_S7(E):  # exact engine
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        if len(hit)==7: tot+=hi-lo
    return float(tot)
def M7(k): return sum(((-1)**t)*__import__('math').comb(6,t)*((7-t)/7)**(k-1) for t in range(7))
# E = {0,1,2,3} : offsets (1,2,3), Lambda = {m in Z^3: m1+2m2+3m3=0}. Enumerate |m|<=L.
E=[0,1,2,3]; offs=[1,2,3]; k=len(E)
exact=meas_S7(E); m7=M7(k)
print(f"E={E}: exact meas_S7={exact:.6f}  M7({k})={m7:.6f}  target corr = {exact-m7:.6f}")
for L in [3,6,10,15,20,30]:
    s=0+0j; nlow=0; clow=0+0j
    for m1 in range(-L,L+1):
     for m2 in range(-L,L+1):
      for m3 in range(-L,L+1):
        if m1+2*m2+3*m3==0 and (m1,m2,m3)!=(0,0,0):
            # pad to 7? No -- k-1=3 offsets, but K expects a vector; use the 3 coords (chat over them)
            kk=K((m1,m2,m3))
            s+=kk
            supp=sum(1 for x in (m1,m2,m3) if x!=0)
            if supp<=5: clow+=kk; nlow+=1
    print(f"  |m|<=L={L:2d}: sum K = {s.real:+.6f} (target {exact-m7:+.6f})  | low-supp(<=5) contribution = {clow.real:+.6f} from {nlow} relations")
