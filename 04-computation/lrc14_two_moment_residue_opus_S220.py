"""
opus-2026-07-11-S220: the two-moment (pair-correlation) residue -- the SINGLE remaining lemma of LRC(14)-S3.

After kps THM-701 (wide recursion closes) + mac-mini THM-702/703, the whole program reduces to:
  Phi(F) = p0(F) + (1/3)p1(F) = E[g(N)] <= E[q(N)] = 1 - (2/3)m1 + (1/6)m2  <=  cap_{|F|+1},  for all bounded cores F,
where N = #missed inner sectors {1..6}, g(0)=1,g(1)=1/3,g(>=2)=0, q(t)=1-(2/3)t+(1/6)t(t-1) (quadratic majorant q>=g on {0..6}),
  m1 = E[N] = sum_{j=1}^6 meas{avoid sector j}       (empty-MEAN, a 1-arc avoidance sum),
  m2 = E[N(N-1)] = sum_{j!=j'} meas{avoid j,j'}       (empty-PAIR = PAIR-CORRELATION of arc-avoidance).
=> the residue is the two-moment inequality  4 m1 - m2 >= 6(1 - cap_{k+1}),  consec the extremizer.
This script: (1) verifies the majorant q>=g and Phi<=majorant; (2) computes m1,m2,majorant vs cap_{k+1};
(3) finds the extremal core (min of 4m1-m2); (4) tests whether consec is extremal + the margin.
"""
from fractions import Fraction as F
from itertools import combinations

def sector(e,x): return int(((e*x)%1)*7)
def brk(E):
    Es=[abs(e) for e in E if e!=0]; bps={F(0),F(1)}
    for e in Es:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    return sorted(b for b in bps if 0<=b<1)
def avoid(E,A):
    Aset=set(A); pts=brk(E)+[F(1)]; tot=F(0)
    for i in range(len(pts)-1):
        a,b=pts[i],pts[i+1]
        if b<=a: continue
        mid=(a+b)/2
        if all(sector(e,mid) not in Aset for e in E): tot+=(b-a)
    return tot
def missdist(E):
    """return (p0,p1,m1,m2): p0=P(N=0), p1=P(N=1), m1=E[N], m2=E[N(N-1)], N=#missed of sectors 1..6."""
    pts=brk(E)+[F(1)]; p=[F(0)]*7
    for i in range(len(pts)-1):
        a,b=pts[i],pts[i+1]
        if b<=a: continue
        mid=(a+b)/2
        occ=set(sector(e,mid) for e in E)
        N=len(set(range(1,7))-occ)           # missed inner sectors
        p[N]+=(b-a)
    m1=sum(F(n)*p[n] for n in range(7)); m2=sum(F(n*(n-1))*p[n] for n in range(7))
    return p[0],p[1],m1,m2

# caps (inner-sector convention, from THM-532/534)
cap={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7),13:F(1)}

def q(t): return 1 - F(2,3)*t + F(1,6)*t*(t-1)
def g(t): return {0:F(1),1:F(1,3)}.get(t,F(0))
print("majorant check q>=g on {0..6}:", all(q(t)>=g(t) for t in range(7)), " q =", [str(q(t)) for t in range(7)])
print()

# cores: consec_k vs perturbations (all include 0; k = |E|); the residue is on BOUNDED cores
def consec(k): return list(range(k))
fams = {
  "consec_8   ": consec(8),
  "consec_9   ": consec(9),
  "consec_10  ": consec(10),
  "AP2 k9     ": [0,2,4,6,8,10,12,14,16],
  "near-AP k9 ": [0,1,2,3,4,5,6,7,9],
  "gap k9     ": [0,1,2,3,4,5,6,7,13],
  "2blk k10   ": [0,1,2,3,4,10,11,12,13,14],
}
print(f"{'core':>12} {'k':>3} {'Phi':>8} {'majorant':>9} {'cap_{k+1}':>9} {'maj<=cap?':>9} {'m1':>7} {'m2':>7} {'4m1-m2':>8} {'6(1-cap)':>9}")
for name,E in fams.items():
    k=len(E); p0,p1,m1,m2=missdist(E)
    Phi=p0+F(1,3)*p1; maj=1-F(2,3)*m1+F(1,6)*m2; c=cap[k+1] if k+1 in cap else F(1)
    lhs=4*m1-m2; rhs=6*(1-c)
    ok="OK" if maj<=c else "FAIL"
    print(f"{name:>12} {k:>3} {float(Phi):>8.5f} {float(maj):>9.5f} {float(c):>9.5f} {ok:>9} {float(m1):>7.4f} {float(m2):>7.4f} {float(lhs):>8.4f} {float(rhs):>9.4f}")

print("\n=== is consec the extremizer (min 4m1-m2 <=> max majorant) among k=9 cores? scan diam<=11 primitive ===")
best=None
cnt=0
for extra in range(8,12):   # consec_8 + one element 'extra' in [8,11]
    E=list(range(8))+[extra]; k=len(E)
    p0,p1,m1,m2=missdist(E); maj=1-F(2,3)*m1+F(1,6)*m2
    cnt+=1
    if best is None or maj>best[0]: best=(maj,E,m1,m2)
print(f"among consec_8+[8..11]: max majorant = {float(best[0]):.5f} at {best[1]}  (cap_9={float(cap[9]):.5f}, consec_9 maj={float(1-F(2,3)*missdist(consec(9))[2]+F(1,6)*missdist(consec(9))[3]):.5f})")
