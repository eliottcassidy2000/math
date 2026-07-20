# HONEST RERUN of my S61h domination statistic to m=24 (kp-S128c120's challenge: my window m<=8
# is too small; does |top|/|mass| decay below 1/2 like klein's did, 0.67 -> 0.04%?).
from fractions import Fraction as Fr
from math import factorial
def pmul(p,q):
    r={}
    for (a,b),u in p.items():
        for (c,d),v in q.items():
            k=(a+c,b+d); r[k]=r.get(k,Fr(0))+u*v
    return r
def stats(P,M):
    print(f"  m : |top|/totalmass   top-dominates(|top|>|rest|)   E[P^m]?=0")
    Pm={(0,0):Fr(1)}
    for m in range(1,M+1):
        Pm=pmul(Pm,P)
        g={}
        for (a,b),c in Pm.items():
            if a==b: g[a]=g.get(a,Fr(0))+c
        if not g:
            print(f"  {m:2d}: (no charge-0)"); continue
        amax=max(g)
        terms={a:g[a]*factorial(a) for a in g}          # gamma_a * a!
        E=sum(terms.values())
        mass=sum(abs(v) for v in terms.values())
        top=abs(terms[amax]); rest=mass-top
        share=float(top/mass) if mass!=0 else 0.0
        dom = top>rest
        print(f"  {m:2d}: {share:.4f}            {str(dom):5s}                    {'0' if E==0 else '!=0'}")
# my S61h P (two-sided, non-constant coeff): Z^2 + (1+ZW)W = Z^2 + W + Z W^2  (charges +2,-1,-1)
print("P = Z^2 + W + Z*W^2  (my S61h example):")
stats({(2,0):Fr(1),(0,1):Fr(1),(1,2):Fr(1)},24)
# b-sweep on the {-1,0,1} stratum P = z + w + b   (charges +1,-1, and 0 via constant b)
for b in (Fr(1),Fr(2),Fr(3)):
    print(f"\nP = Z + W + {b}  (charge-0 weight b={b}, kp's b-sweep analog):")
    stats({(1,0):Fr(1),(0,1):Fr(1),(0,0):b},24)
