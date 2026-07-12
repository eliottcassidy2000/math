#!/usr/bin/env python3
"""cont.52 part 2: verify p5 parity closed forms, extend to p4, p3; assemble J(consec_k)."""
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
print("VERIFY closed forms:")
print("  p6 = 1/(7(k-1));  p5(even) = 5/(14(k-1));  p5(odd) = (5k-9)/(14(k-1)(k-2))")
ok6=ok5=True
for k in range(7,18):
    p=pdist(k)
    p6c=F(1,7*(k-1))
    p5c = F(5,14*(k-1)) if k%2==0 else F(5*k-9,14*(k-1)*(k-2))
    m6 = p[6]==p6c; m5 = p[5]==p5c
    ok6 &= m6; ok5 &= m5
    print(f"  k={k}: p6 {'OK' if m6 else 'NO '} | p5 {'OK' if m5 else 'NO %s vs %s'%(p[5],p5c)}")
print(f"\n  p6 closed form: {'CONFIRMED' if ok6 else 'FAILS'}; p5 parity closed form: {'CONFIRMED' if ok5 else 'FAILS'}")
print("\np4 hunt: p4 * 14(k-1)(k-2)(k-3) (seek pattern):")
for k in range(7,18):
    p=pdist(k)
    val=p[4]*14*(k-1)*(k-2)*(k-3)
    print(f"  k={k} ({'even' if k%2==0 else 'odd '}): p4*14(k-1)(k-2)(k-3) = {val} = {float(val):.3f}")
