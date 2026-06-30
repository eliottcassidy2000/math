"""Verify construction M = n/Phi_6(n) for n=4..20 (ductile picture), margin, asymptotics, n-1 prime-power flag."""
import math
from fractions import Fraction
def M_exact(S,Qmax):
    best=Fraction(0)
    for q in range(2,Qmax+1):
        bb=0
        for a in range(1,q):
            if math.gcd(a,q)!=1: continue
            m=q
            for s in S:
                r=(s*a)%q; d=r if r<=q-r else q-r
                if d<m:m=d
                if m<=bb:break
            if m>bb:bb=m
        v=Fraction(bb,q)
        if v>best:best=v
    return best
def is_pp(m):
    if m<2: return False
    for p in range(2,int(m**.5)+1):
        if m%p==0:
            while m%p==0: m//=p
            return m==1
    return True
def factor(m):
    f={};d=2
    while d*d<=m:
        while m%d==0: f[d]=f.get(d,0)+1; m//=d
        d+=1
    if m>1: f[m]=f.get(m,0)+1
    return f
print(f"{'n':>3} {'M_constr':>9} {'=n/Phi6?':>9} {'Phi6=n^2-n+1':>14} {'Phi6 factored':>16} {'n-1 pp(PG)?':>11} {'margin~1/n^2':>13}")
for n in range(4,21):
    Phi6=n*n-n+1
    S=list(range(1,n-1))+[(n-1)*n]
    M=M_exact(S,Phi6+5)
    ok = (M==Fraction(n,Phi6))
    fc="*".join(f"{p}^{e}" if e>1 else f"{p}" for p,e in factor(Phi6).items())
    margin=float(Fraction(n,Phi6)-Fraction(1,n))
    print(f"{n:>3} {str(M):>9} {str(ok):>9} {Phi6:>14} {fc:>16} {str(is_pp(n-1)):>11} {margin:>13.6f}")
print("\nBIG PICTURE: C(n) = 2/(2n-1) [drop-2] for n<=6; = n/Phi_6(n) [construction] for n>=7. ONE transition at n=7=apex.")
print("  margin C(n)-1/n = (n-1)/(n.Phi_6) ~ 1/n^2 -> 0 (the razor-thin covering margin, vanishing).")
print("  n.C(n) = n^2/Phi_6 = n^2/(n^2-n+1) -> 1^- (covering-min approaches 1/n from ABOVE).")
print("  PG(2,n-1) exists iff n-1 prime power; n=7 (n-1=6, NOT pp, PG(2,6) fails by Bruck-Ryser) IS the transition.")
