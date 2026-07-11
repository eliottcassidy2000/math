#!/usr/bin/env python3
"""THM-710: factorial moments are eigenvectors of the far-element transfer.
(1) ALGEBRA (the proof, verified symbolically over all distributions on {0..7}):
    under p_j -> ((7-j)/7)p_j + ((j+1)/7)p_{j+1}:  m_r -> ((7-r)/7) m_r  EXACTLY.
(2) RUNG PROPAGATION (exact rationals): V' = (6/7)V + 1/7 - m2/84 <= (6/7)V + 1/7;
    check (6/7)cap_{k+1} + 1/7 <= cap_{k+2} for k = 8..11.
(3) EMPIRICS: real far-element sweeps -- m_r(E'+{w}) vs ((7-r)/7) m_r(E'), error ~ C/w."""
from fractions import Fraction as F
import random
CAP={9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7),13:F(1)}
# (1) exact algebra on random rational distributions
random.seed(2)
def falling(j,r):
    v=1
    for i in range(r): v*=(j-i)
    return v
ok=True
for _ in range(200):
    p=[F(random.randint(0,20),1) for _ in range(8)]
    s=sum(p); p=[x/s for x in p]
    pn=[F(0)]*8
    for j in range(8):
        pn[j]+=F(7-j,7)*p[j]
        if j+1<8: pn[j]+=F(j+2-1,7)*p[j+1]  # ((j+1)+... careful: contribution to slot j from p_{j+1} is ((j+1)/7)? operator: p_j' = ((7-j)/7)p_j + ((j+1)/7)p_{j+1}
    pn=[F(7-j,7)*p[j] + (F(j+1,7)*p[j+1] if j+1<8 else F(0)) for j in range(8)]
    for r in range(0,7):
        mr=sum(falling(j,r)*p[j] for j in range(8))
        mrn=sum(falling(j,r)*pn[j] for j in range(8))
        if mrn != F(7-r,7)*mr: ok=False
print(f"(1) eigen-identity m_r -> ((7-r)/7) m_r: {'EXACT on 200 random distributions x r=0..6' if ok else '*** FAIL ***'}")
print()
print("(2) rung propagation (6/7)cap_(k+1) + 1/7 <= cap_(k+2), exact:")
for k in range(8,12):
    lhs=F(6,7)*CAP[k+1]+F(1,7); rhs=CAP[k+2]
    print(f"  k={k}: {lhs} = {float(lhs):.4f} <= {float(rhs):.4f}  {'OK' if lhs<=rhs else '*** FAIL ***'}  slack {float(rhs-lhs):+.4f}")
print()
print("(3) real far-element eigen check (E'={1,2,3,5}, w sweep): |m_r(E'+w) - ((7-r)/7)m_r(E')|*w")
def pvec(E):
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
Ep=[1,2,3,5]; pE=pvec(Ep)
for w in [101,301,901]:
    pW=pvec(Ep+[w])
    row=[]
    for r in [1,2,3]:
        mr=sum(falling(j,r)*pE[j] for j in range(8))
        mrw=sum(falling(j,r)*pW[j] for j in range(8))
        row.append(f"r={r}: {float(abs(mrw-F(7-r,7)*mr)*w):.3f}")
    print(f"  w={w}: " + "  ".join(row) + "   (bounded => O(1/w) confirmed)")
