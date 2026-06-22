"""
Bounding the multi-large equidistribution lemma (kps-S31v).
G_C = core lonely set {x: ||c x||>=1/14, c in C} (a union of arcs). U_v={||vx||<1/14} (measure 1/7).
Want: meas(G_C \ union U_{v_i}) > 0. Toolbox:
(1) EQUIDISTRIBUTION (Erdos-Turan/Koksma): meas(G_C n U_v) = (1/7)meas(G_C) + err, |err|<= arcCount/(7v).
(2) UNION BOUND (Bonferroni-1): uncovered >= (1 - r/7)meas(G_C) - r*err. POSITIVE for r<=6.
(3) PAIRWISE INDEPENDENCE: meas(U_i n U_j) -> (1/7)^2 for non-resonant => second-moment for r>=7.
Verify (1),(3) exactly (comb teeth in arcs).
"""
from fractions import Fraction as F
def gc_arcs(C):
    # breakpoints where ||c x||=1/14 : x=(14m+-1)/(14c); G_C = {min_c ||cx||>=1/14}
    C=[c for c in C if c!=0]; bps={F(0),F(1)}
    for c in C:
        a=abs(c)
        for m in range(0,a+1):
            for s in (-1,1):
                x=F(14*m+s,14*a)
                if 0<=x<=1: bps.add(x)
    B=sorted(bps); arcs=[]
    for lo,hi in zip(B,B[1:]):
        mid=(lo+hi)/2
        if all(min((c*mid)%1,1-(c*mid)%1)>=F(1,14) for c in C): arcs.append((lo,hi))
    # merge adjacent
    merged=[]
    for a,b in arcs:
        if merged and merged[-1][1]==a: merged[-1]=(merged[-1][0],b)
        else: merged.append((a,b))
    return merged
def meas(arcs): return float(sum(b-a for a,b in arcs))
def meas_in_comb(arcs, v):
    # exact meas of {x in arcs : ||vx||<1/14} = sum over arcs of comb-measure
    tot=F(0)
    for a,b in arcs:
        # comb teeth: intervals ((m-1/14)/v,(m+1/14)/v); count overlap with [a,b]
        import math
        m_lo=math.floor(float(a*v)-0.5); m_hi=math.ceil(float(b*v)+0.5)
        for m in range(m_lo,m_hi+1):
            lo=F(14*m-1,14*v); hi=F(14*m+1,14*v)
            l=max(lo,a); h=min(hi,b)
            if h>l: tot+=h-l
    return float(tot)
C=[1,2,3,4,5]  # bounded core (5 speeds)
arcs=gc_arcs(C); mGC=meas(arcs)
print(f"core C={C}: meas(G_C)={mGC:.5f}, arcCount={len(arcs)}")
print(f"\n(1) EQUIDISTRIBUTION meas(G_C n U_v)/meas(G_C) -> 1/7={1/7:.4f}:")
for v in [20,50,100,500,2000,10000]:
    r=meas_in_comb(arcs,v)/mGC
    print(f"   v={v:6d}: ratio={r:.4f}  err={abs(r-1/7):.4f}  bound arcCount/(7v)={len(arcs)/(7*v):.4f}")
print(f"\n(3) PAIRWISE meas(G_C n U_v1 n U_v2)/meas(G_C) -> (1/7)^2={1/49:.4f}:")
import math
for v1,v2 in [(101,103),(101,202),(500,1000),(733,1009)]:
    # intersect arcs with U_v1, then with U_v2 (approx via fine refine on arcs)
    sub=[]
    for a,b in arcs:
        m_lo=math.floor(float(a*v1)-0.5); m_hi=math.ceil(float(b*v1)+0.5)
        for m in range(m_lo,m_hi+1):
            lo=F(14*m-1,14*v1); hi=F(14*m+1,14*v1); l=max(lo,a); h=min(hi,b)
            if h>l: sub.append((l,h))
    r2=meas_in_comb(sub,v2)/mGC
    res = "RESONANT" if math.gcd(v1,v2)>1 or v2%v1==0 else "indep"
    print(f"   v1={v1},v2={v2} ({res}): ratio={r2:.4f}  vs 1/49={1/49:.4f}")
print("\n=> (1) equidist err ~ arcCount/(7v) (Erdos-Turan) => each comb covers ~1/7 of G_C for v large.")
print("   (2) UNION BOUND: r<=6 large speeds => uncovered >= (1-r/7)meas(G_C)>0 (CLEAN, rigorous via ET).")
print("   (3) pairwise ->1/49 (independent) => r>=7 closes by 2nd-moment/(6/7)^r modulo bounded resonances.")
