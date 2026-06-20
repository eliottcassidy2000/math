from fractions import Fraction as F
from math import gcd
def bps(E):
    E=sorted(set(E)); b={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*e+1): b.add(F(a,7*e))
    return sorted(x for x in b if 0<=x<=1)
def p0(E):
    E=sorted(set(E)); B=bps(E); s=F(0)
    for lo,hi in zip(B,B[1:]):
        if hi<=lo: continue
        mid=(lo+hi)/2
        if not (set(range(1,7))-set(int((e*mid)%1*7) for e in E)): s+=hi-lo
    return s
CAPS={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7)}
# the over-QB witness, exact
E=[0,19,20,21,22,23,24,25]
v=p0(E)
print("witness k=8",E,"p0=",v,"=",float(v),"cap=",float(CAPS[8]),"QB=0.19660")
print("  margin to cap:",float(CAPS[8]-v))
# scan singleton + consec_7 over all far scales g to find the TRUE finite-scale max
print("\n[singleton 0] + (g + consec_7), scan g: find max p0 over finite g")
worst=F(0); wg=None
for g in range(8,300):
    E=[0]+[g+i for i in range(7)]
    # primitive? gcd of (0,g,...,g+6) -> gcd includes consecutive => 1
    vv=p0(E)
    if vv>worst: worst=vv; wg=g
print("  max p0=",float(worst),"at g=",wg,"  (QB=0.19660, cap=%.5f)"%float(CAPS[8]))
print("  exceeds QB by",float(worst)-F(19660,100000)," ; under cap by",float(CAPS[8]-worst))
# also consec_7 + far singleton (other orientation): consec_7 U {g}
print("\nconsec_7 U {g} (singleton FAR), scan g:")
worst2=F(0); wg2=None
for g in range(7,300):
    E=list(range(7))+[g]
    vv=p0(E)
    if vv>worst2: worst2=vv; wg2=g
print("  max p0=",float(worst2),"at g=",wg2)
