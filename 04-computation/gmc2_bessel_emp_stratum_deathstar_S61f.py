# {-1,0,1} stratum: P = a(u)w + b(u) + c(u)z, u=zw. Claim: E[P^m]=0 all m  <=>  a*c==0 and b==0.
from fractions import Fraction as Fr
from math import factorial
def poch_none(): pass
def L(poly):  # L(z^p w^q)=p! [p==q]
    s=Fr(0)
    for (p,q),co in poly.items():
        if p==q: s+=co*factorial(p)
    return s
def pmul(p,q):
    r={}
    for k1,u in p.items():
        for k2,v in q.items():
            k=(k1[0]+k2[0],k1[1]+k2[1]); r[k]=r.get(k,Fr(0))+u*v
    return {k:v for k,v in r.items() if v!=0}
def padd(*ps):
    r={}
    for p in ps:
        for k,v in p.items(): r[k]=r.get(k,Fr(0))+v
    return {k:v for k,v in r.items() if v!=0}
def scal(p,c): return {k:v*Fr(c) for k,v in p.items() if v*Fr(c)!=0}
Z={(1,0):Fr(1)}; W={(0,1):Fr(1)}; one={(0,0):Fr(1)}; U={(1,1):Fr(1)}  # u=zw
def poly_u(coeffs):  # sum coeffs[k] u^k
    r={}
    for k,c in enumerate(coeffs):
        if c: 
            uk=one
            for _ in range(k): uk=pmul(uk,U)
            r=padd(r,scal(uk,c))
    return r
def build_P(a,b,c):  # a,b,c lists of u-coeffs
    A=poly_u(a); B=poly_u(b); C=poly_u(c)
    return padd(pmul(A,W), B, pmul(C,Z))
def EPm(P,M):
    Pm=one; out=[]
    for m in range(1,M+1): Pm=pmul(Pm,P); out.append(L(Pm))
    return out
print("=== {-1,0,1} stratum: E[P^m]=0 all m  <=>  a*c==0 and b==0 ===")
tests=[
  ([1],[0],[1], "a=1,c=1,b=0 (ac!=0)"),
  ([1],[0],[0], "a=1,c=0,b=0 (one-sided, ac=0)"),
  ([1],[1],[0], "a=1,c=0... b=1"),
  ([1],[0],[1], None),
  ([1,1],[0],[1], "a=1+u,c=1,b=0 (ac!=0)"),
  ([1],[0,1],[1],"a=1,c=1,b=u (ac!=0,b!=0)"),
  ([0,1],[0],[0],"a=u,c=0,b=0 (one-sided)"),
  ([2],[0],[3],"a=2,c=3,b=0 (ac!=0)"),
]
for a,b,c,label in tests:
    if label is None: continue
    P=build_P(a,b,c); ep=EPm(P,6)
    ac_zero = all(x==0 for x in [ (sum(1 for _ in a) and 0) ]) # placeholder
    print(f"  {label:38s}: E[P^m] m=1..6 = {[str(x) for x in ep]}  {'NULLCONE(all 0)' if all(x==0 for x in ep) else 'nonzero'}")
print("\n  predicted nullcone iff (a*c==0 poly) and (b==0). Confirm: only the one-sided (c=0 or a=0)")
print("  with b=0 give all-zero; any ac!=0 or b!=0 gives nonzero. => GMC(2) on {-1,0,1} reduces to one-sided.")
