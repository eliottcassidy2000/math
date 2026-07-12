#!/usr/bin/env python3
"""Slice 3: the SHARP adversarial question. J = mu(7-mu) - Var. Consec {1..9} sits at
(mu=1.446, Var=2.967, J=5.062). Can any family beat it? Two attack vectors:
(A) LOW-mu families (best coverers) -- minimize mu, does Var stay high enough?
(B) simulated-annealing on J directly over 9-sets in a wide box."""
from fractions import Fraction as F
import random
def stats(E):
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
    mu=sum(j*p[j] for j in range(8)); m2=sum(j*(j-1)*p[j] for j in range(8))
    return float(mu), float(6*mu-m2)
Jmin=float(F(4465,882))
print(f"target: beat consec {{1..9}} J={Jmin:.4f} (mu=1.446)")
# (A) low-mu hill-climb: minimize mu, report J
random.seed(1)
def climb(objective, sign):
    E=sorted(random.sample(range(0,30),9)); cur=objective(E)
    for _ in range(600):
        i=random.randrange(9); nd=random.choice([-2,-1,1,2])
        E2=E[:]; E2[i]=E2[i]+nd
        if len(set(E2))!=9 or min(E2)<0 or max(E2)>60: continue
        v=objective(sorted(E2))
        if sign*v < sign*cur: cur=v; E=sorted(E2)
    return E,cur
print("\n(A) low-mu hill-climbs (minimize mu; report resulting J):")
worst_mu=None
for t in range(6):
    E,mu=climb(lambda E: stats(E)[0], +1)
    _,J=stats(E)
    print(f"  min-mu run {t}: mu={mu:.4f} J={J:.4f} at {E}  {'*** BEATS CONSEC' if J<Jmin else ''}")
    if worst_mu is None or J<worst_mu[0]: worst_mu=(J,E)
print("\n(B) direct J-minimization (simulated-annealing style):")
best=(Jmin,[1,2,3,4,5,6,7,8,9])
for t in range(8):
    E,J=climb(lambda E: stats(E)[1], +1)
    if J<best[0]: best=(J,E)
    print(f"  J-min run {t}: J={J:.4f} at {E}  {'*** BEATS CONSEC' if J<Jmin-1e-6 else ''}")
print(f"\nGLOBAL best found: J={best[0]:.4f} at {best[1]}")
print(f"CONSEC {{1..9}} SURVIVES as J-min: {best[0]>=Jmin-1e-6}  (margin over threshold 4.747: {best[0]-4.7473:+.4f})")
