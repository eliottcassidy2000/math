#!/usr/bin/env python3
"""cont.44 FREE EXPLORATION: the cross-functional extremal atlas. Evaluate the LRC(14)
functionals on the ALGEBRAIC families (klein's atlas: AP/consec, deep-well, mod-7 pole,
detuned-AP) across k. Look for a universal extremal + continued-fraction structure in the
extremal values (Ostrowski-ladder inspiration, mac-mini-S38)."""
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
def funcs(E):
    p=pdist(E)
    mu=sum(j*p[j] for j in range(8)); m2=sum(j*(j-1)*p[j] for j in range(8))
    J=6*mu-m2; T=[sum(p[n] for n in range(j,8)) for j in range(8)]
    POS=6*T[1]+4*T[2]+2*T[3]; BUNCH=2*(p[5]+3*p[6]); var=sum(j*j*p[j] for j in range(8))-mu*mu
    return dict(mu=mu,J=J,POS=POS,BUNCH=BUNCH,var=var,p0=p[0])
def cf(fr):  # continued fraction of a Fraction in [0,1]
    a=[]; x=fr
    for _ in range(8):
        if x==0: break
        inv=1/x; ai=inv.numerator//inv.denominator; a.append(ai); x=inv-ai
    return a
print("EXTREMAL ATLAS at k=9 (the binding row), algebraic families:")
fams={"consec/AP {0..8}":list(range(9)),
      "AP {1..9}":list(range(1,10)),
      "mod-7 pole {1,8..57}":[1,8,15,22,29,36,43,50,57],
      "detuned-AP {0..7,9}":[0,1,2,3,4,5,6,7,9],
      "deep-well-ish {0..6,8,20}":[0,1,2,3,4,5,6,8,20],
      "double {0..7,16}":[0,1,2,3,4,5,6,7,16]}
print(f"{'family':24s} {'J':>9s} {'POS':>9s} {'BUNCH':>8s} {'nu-var':>8s} {'p0':>7s}")
for nm,E in fams.items():
    f=funcs(E)
    print(f"  {nm:22s} {float(f['J']):9.4f} {float(f['POS']):9.4f} {float(f['BUNCH']):8.4f} {float(f['var']):8.4f} {float(f['p0']):7.4f}")
print()
print("J(consec {0..k-1}) across k -- continued-fraction / Ostrowski structure in the VALUES?")
for k in range(6,13):
    E=list(range(k)); f=funcs(E)
    Jv=f['J']
    print(f"  k={k}: J = {Jv} = {float(Jv):.5f}   cf(J-5)={cf(Jv-5) if Jv>5 else cf(Jv)}   p0={float(f['p0']):.5f}")
print()
print("p0(consec {0..k-1}) across k -- the best-coverer value; three-gap closed form?")
for k in range(6,13):
    E=list(range(k)); f=funcs(E)
    print(f"  k={k}: p0 = {f['p0']} = {float(f['p0']):.5f}   1-p0 (=T1=meas S7) = {1-f['p0']} = {float(1-f['p0']):.5f}")
