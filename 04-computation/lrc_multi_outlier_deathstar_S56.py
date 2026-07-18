from fractions import Fraction as F
from math import gcd
def M_lower(fam,Qcap):
    best=F(0)
    for q in range(2,Qcap+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q)
    return best
def test(core,j,scale,Qcap=60):
    primes=[101,103,107,109,113,127,131]
    V=sorted(core+[scale*p for p in primes[:j]])
    MV=M_lower(V,Qcap); MC=M_lower(core,Qcap); pred=min(MC,F(1,2*j))
    print(f"  core|{len(core)}| j={j} scale={scale}: M(C)~{float(MC):.4f} 1/(2j)={float(F(1,2*j)):.4f} "
          f"pred={float(pred):.4f}  M(V)>={float(MV):.4f}  [>=pred:{MV>=pred}] M(V)>1/14:{MV>F(1,14)}")
print("core={1..12}(M=1/13); j far outliers -> M(V) >= min(1/13,1/(2j)):")
C=list(range(1,13))
for j in range(1,7): test(C,j,60)
print("\napex-7: core={1..6}, j=7 far outliers, 1/(2j)=1/14:")
for scale in (20,60,200): test(list(range(1,7)),7,scale)
