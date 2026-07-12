#!/usr/bin/env python3
"""cont.52: hunt closed forms for the full consec sector-distribution p_j(k). p6=1/(7(k-1))
found. Seek p5,p4,p3,... via structure. Test candidate forms; if J(consec_k) has a closed
form, THM-718's compact check collapses to a formula."""
from fractions import Fraction as F
def pdist(k):
    E=list(range(k))
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
print("p_j(consec_k) exact, seeking closed forms (p6=1/(7(k-1)) known):")
data={j:{} for j in range(8)}
for k in range(7,16):
    p=pdist(k)
    for j in range(8): data[j][k]=p[j]
# p6 = 1/(7(k-1)); test p5 candidates
print("\np5: exact values and ratio p5/p6 (=7(k-1)*p5):")
for k in range(7,16):
    p5=data[5][k]; p6=data[6][k]
    print(f"  k={k}: p5={p5}={float(p5):.5f}, p5/p6 = {p5/p6}")
print("\np5 * 7(k-1)(k-2) (seek integer/linear numerator):")
for k in range(7,16):
    val=data[5][k]*7*(k-1)*(k-2)
    print(f"  k={k}: p5*7(k-1)(k-2) = {val} = {float(val):.4f}")
print("\nJ(consec_k) exact + candidate closed form:")
for k in range(7,16):
    p=pdist(k); J=sum(p[j]*j*(7-j) for j in range(8))
    print(f"  k={k}: J = {J} = {float(J):.5f}")
# the KEY: does J have a rational closed form in k? fit J*7*(k-1) etc.
print("\nJ * common structure (seek polynomial/k pattern):")
for k in range(8,15):
    p=pdist(k); J=sum(p[j]*j*(7-j) for j in range(8))
    print(f"  k={k}: J={float(J):.4f}, J*(k-1)={float(J*(k-1)):.4f}, 7J={float(7*J):.4f}")
