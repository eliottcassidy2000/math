#!/usr/bin/env python3
"""cont.51 part 2: find closed forms for the consec sector-distribution tail (p6,p5,p4).
p6 = 1/(7(k-1)) conjectured. Test p6 exactly, seek p5, p4 forms; simplify BUNCH = p5+3p6."""
from fractions import Fraction as F
def pdist(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1); p=[F(0)]*8
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        p[7-len(hit)]+=b-a
    return p
print("p6 = 1/(7(k-1))?  and  BUNCH(consec) = p5+3p6:")
for k in range(7,14):
    p=pdist(list(range(k)))
    p6c=F(1,7*(k-1))
    BUNCH=p[5]+3*p[6]
    print(f"  k={k}: p6={p[6]} vs 1/(7(k-1))={p6c} [{'MATCH' if p[6]==p6c else 'NO'}]  p5={p[5]}  BUNCH={BUNCH}={float(BUNCH):.4f}")
print()
print("klein BUNCH-max = 18/(7k-6); consec BUNCH vs that + seek p5 closed form:")
for k in range(7,14):
    p=pdist(list(range(k)))
    # try p5 = a/(7(k-1)(k-2))? or similar
    p5=p[5]
    guess1=F(3*k-8, 7*(k-1)*(k-2)) if k>2 else F(0)
    print(f"  k={k}: p5={p5}={float(p5):.5f}  BUNCH-max klein 18/(7k-6)={F(18,7*k-6)}={float(F(18,7*k-6)):.4f}")
print()
# the FULL distribution closed form attempt: p_j numerators over 7*(k-1)*...?
print("cumulative T_j = P(N>=j) for consec -- cleaner than p_j?")
for k in [8,9]:
    p=pdist(list(range(k)))
    T=[sum(p[n] for n in range(j,8)) for j in range(8)]
    print(f"  k={k}: T1..T6 = " + ", ".join(f"{float(T[j]):.4f}" for j in range(1,7)))
    print(f"        T6=p6={T[6]}=1/(7(k-1)); T5={T[5]}; T4={T[4]}")
