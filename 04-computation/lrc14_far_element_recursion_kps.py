import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def dist_p(E):
    """p_t = meas{N=t}, N=#missed of 6 inner sectors."""
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p=[Fraction(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        p[sum(1 for j in range(1,7) if j not in hit)]+=hi-lo
    return p
def primitive(E): return reduce(gcd,E)==1
# Q(m) = max over m-element sets (0 in E) of [p_0 + (1/7) p_1] = the far-element plateau bound.
# (p_1 = P(miss exactly 1 of the 6 inner sectors); a far point fills it w.p. 1/7... but careful:
#  adding an independent uniform point to orbit O: P(O U {y} covers all 7) = p_0(O) + (1/7)*P(O misses exactly 1 of 7 total).
#  sector 0 always hit by e=0, so "misses exactly 1 of 7" = misses exactly 1 of the 6 inner = p_1.)
capf={7:0.2649,8:0.38153,9:0.49426,10:0.6044}
def Q(m, box):
    best=(Fraction(0),None)
    for tail in itertools.combinations(range(1,box+1),m-1):
        E=(0,)+tail
        if not primitive(E): continue
        p=dist_p(E); val=p[0]+Fraction(1,7)*p[1]
        if val>best[0]: best=(val,E)
    return best
print("=== RECURSION: Q(k-1)=max_{(k-1)-sets}[p0+(1/7)p1] vs cap_k (the far-element plateau bound) ===\n")
for k in [8,9,10]:
    m=k-1; box={7:14,8:14,9:13}[m]
    qval,qset=Q(m,box)
    cap=capf[k]
    consec=list(range(m)); pc=dist_p(consec); qc=pc[0]+Fraction(1,7)*pc[1]
    print(f"k={k}: Q(k-1={m})={float(qval):.5f} at {qset}  (consec_{m} gives {float(qc):.5f})  cap_{k}={cap}  "
          f"Q<cap? {float(qval)<cap} margin={cap-float(qval):.4f}")
print()
print("=== so: a k-set with a FAR element has p0 -> ~Q(k-1) < cap_k. Recursion DOWN in k. ===")
print("=== also verify the plateau formula matches the far-element data (k=9 core consec_8) ===")
core=list(range(8)); pc=dist_p(core)
plat=pc[0]+Fraction(1,7)*pc[1]
print(f"   core=consec_8: p0={float(pc[0]):.5f} p1={float(pc[1]):.5f} => plateau={float(plat):.5f} (data far-w plateau ~0.36) ")
