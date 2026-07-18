from fractions import Fraction as F
from math import gcd
def M_exact(fam):
    Q=2*max(fam)+2; best=F(0)
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q)
    return best
def check(name,fam):
    fam=sorted(fam); w=fam[-1]; v2=fam[-2]; Vp=fam[:-1]
    MV=M_exact(fam); MVp=M_exact(Vp)
    bound=F(w)*MVp/(w+v2)                 # claimed: M(V) >= w*M(V')/(w+v2)
    bound_lrc=F(w,(len(fam))*(w+v2))      # >= w/((n)*(w+v2))? no: M(V')>=1/(n-1); n-1=len-1... use 1/13
    b13=F(w,13*(w+v2))
    ok = MV>=bound
    print(f"{name}: M(V)={MV}={float(MV):.4f}  bound w*M(V')/(w+v2)={float(bound):.4f} [OK:{ok}]  "
          f"1/13-form w/(13(w+v2))={float(b13):.4f}  (M(V')={float(MVp):.4f})")

check("AP {1..13}", list(range(1,14)))
check("GW", [1,2,3,4,5,6,7,8,9,10,11,13,24])
check("deep well {1..12,182}", list(range(1,13))+[182])
check("far element {1..12,500}", list(range(1,13))+[500])
check("covering min [2,5,6,8,11,13,14,16,19,20,24,25,27]", [2,5,6,8,11,13,14,16,19,20,24,25,27])
check("2*AP", [2*k for k in range(1,14)])
check("comparable random", [3,5,7,8,11,13,17,19,23,25,29,31,37])
