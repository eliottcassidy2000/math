#!/usr/bin/env python3
"""
lrc14_k5_exact_closure_macmini.py   (mac-mini, ANGLE E)

GOAL: close LRC(14) case S3 for cluster size k=|L|=5 (the first open size; pigeonhole
margin c=|L|-1=4 gives 7/4-1=3/4<1, so pigeonhole fails -> need rho*).

OBJECT (THM-527):  offset set E, 0 in E, |E|=k.  For x in [0,1) put the k points
{frac(e*x): e in E} on R/Z.  good(x) <=> circular max-gap > 2/7 (fit in open arc < 5/7).
  mu(E)      = meas{x in [0,1):       good(x)}
  rho*(P,E)  = meas{x in G_P:         good(x)},  G_P={tau: ||p tau||>=1/14 for all p in P}.
THM-527: rho*(P,E)>0  ==>  M(P u L) >= 1/14  (case S3).  A uniform floor c0>0 proves S3.
For k=5: |L|=5 large runners, |P|=13-5=8 small runners P subset {1..13}.

EXACT (Fractions, no float decisions): maxgap{frac(e_i x)} is piecewise-linear in x with
rational breakpoints n/e_i (point crosses 0) and n/(e_i-e_j) (collision/order-swap). On each
open subinterval the cyclic order is fixed and every gap g(x)=alpha x+beta is linear; solve
g>2/7 exactly. Intersect with G_P (exact rational arcs) -> rho* exact Fraction.

USAGE:  python3 lrc14_k5_exact_closure_macmini.py [section]
  sections: sanity | enum SMAX | tail | rhofull SMAX | rhok K SMAX | all   (default: sanity)

RESULTS (mac-mini-2026-06-18-S2e, ANGLE E):
  - PURE mu floor k=5: min over ALL 5-shapes (spread<=20) = 9/14, attained at consecutive
    {0,1,2,3,4} (and scalings). For k=5 consecutive IS the global mu-extremizer (unlike k>=7).
    k=6 pure mu floor = 4/7 at {0,1,2,3,4,5}.
  - rho* floor k=5: min over ALL P (|P|=8) x ALL shapes (spread<=16 exhaustive; 1820 shapes
    x 1287 P; gcd=1 tail to spread>=17; random to spread 60) = 95/2548 ~= 0.03728, attained
    at E={0,2,4,6,8} (=2*consecutive), P={1,2,3,4,7,9,12,13}. NO rho*=0 anywhere.
  - rho* floor k=6: min over ALL P (|P|=7) x shapes (spread<=9 exhaustive) = 3488/63063
    ~= 0.05531 at consecutive {0,1,2,3,4,5}. NO zeros.
  - SCALING (PROVED exact): good(g*E)=preimage of good(E) under x->g*x. mu scale-invariant;
    rho* NOT (G_P fixed). Min over scalings is at SMALL g; rho*(P,g*E)->mu(E)*meas(G_P) as
    g->inf (Weyl), a POSITIVE limit (k=5: >=(9/14)(2479/17640)=0.0903). Runaway shapes give
    min_P rho* ~0.09-0.12, all >> floor.
  - END-TO-END: real covering 13-sets S={1,2,3,4,7,9,12,13}u{Vmax-e:e in {0,2,4,6,8}} have
    meas(G_S at 1/14)>0 (explicit witness tau) => M(S)>=1/14 directly.
  - HONEST GAP: rho*>=meas(G_P)-5/14 (IE bound) positive for only 6/1287 P -> positivity
    needs EXACT alignment, not measure-counting. Bounded spread/scaling = verified finite
    check; spread->inf & scaling g->inf tails controlled by the Weyl limit but need a
    quantitative rate (= OPEN-Q-108). k=5,6 CLOSED MODULO the same tail rate as THM-527.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

H   = F(1,14)
THR = F(2,7)          # gap threshold: good <=> maxgap > 2/7

# ---------- exact G_P machinery ----------
def danger_arcs(u, h=H):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else:   iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def safe_set(A, h=H):
    if not A: return [(F(0),F(1))]
    dz=merge([iv for u in A for iv in danger_arcs(u,h)])
    safe=[]; prev=F(0)
    for a,b in dz:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<F(1): safe.append((prev,F(1)))
    return safe
def meas(arcs): return sum((b-a) for a,b in arcs)
def intersect(A,B):
    out=[]; i=j=0; A=sorted(A); B=sorted(B)
    while i<len(A) and j<len(B):
        a0,a1=A[i]; b0,b1=B[j]
        lo=max(a0,b0); hi=min(a1,b1)
        if lo<hi: out.append((lo,hi))
        if a1<b1: i+=1
        else: j+=1
    return out

# ---------- EXACT good-set of an offset set E ----------
def good_set_exact(E):
    es=sorted(set(E))
    bps=set([F(0),F(1)])
    for e in es:
        for n in range(1,e): bps.add(F(n,e))
    for i in range(len(es)):
        ei=es[i]
        for j in range(i+1,len(es)):
            d=es[j]-ei
            for n in range(1,d): bps.add(F(n,d))
    bp=sorted(bps)
    good=[]
    for t in range(len(bp)-1):
        lo,hi=bp[t],bp[t+1]; mid=(lo+hi)/2
        valmap={}
        for e in es:
            v=(e*mid)%1
            if v not in valmap: valmap[v]=(e,(e*mid)-v)  # (rep e, floor)
        keys=sorted(valmap)
        slot_e=[valmap[v][0] for v in keys]
        slot_c=[valmap[v][1] for v in keys]
        m=len(keys)
        if m==1:
            good.append((lo,hi)); continue
        sub=[]
        for k in range(m):
            if k<m-1:
                alpha=slot_e[k+1]-slot_e[k]; beta=-(slot_c[k+1]-slot_c[k])
            else:
                alpha=slot_e[0]-slot_e[m-1]; beta=F(1)-(slot_c[0]-slot_c[m-1])
            if alpha==0:
                if beta>THR: sub.append((lo,hi))
            else:
                xb=(THR-beta)/alpha
                if alpha>0:
                    a=max(lo,xb);
                    if a<hi: sub.append((a,hi))
                else:
                    b=min(hi,xb)
                    if lo<b: sub.append((lo,b))
        for a,b in merge(sub):
            if a<b: good.append((a,b))
    return merge(good)

def mu_exact(E):           return meas(good_set_exact(E))
def rho_star_exact(P,E):
    g=good_set_exact(E); gp=safe_set(P)
    return meas(intersect(g,gp)), meas(gp)

# ===================================================================== SANITY
def sec_sanity():
    print("="*90)
    print("SANITY: pure consecutive mu_k (THM-527: mu_3=1,mu_4=19/21,mu_5=9/14,mu_13=829/4620)")
    print("="*90)
    exp={3:F(1),4:F(19,21),5:F(9,14),13:F(829,4620)}
    for k in range(3,14):
        m=mu_exact(list(range(k)))
        tag=""
        if k in exp: tag="  EXPECT "+str(exp[k])+("  OK" if m==exp[k] else "  *** MISMATCH")
        print(f"  k={k:2d}: mu_k = {str(m):>14s} = {float(m):.6f}{tag}")

# ===================================================================== ENUM (min mu over shapes)
def sec_enum(smax):
    print("="*90)
    print(f"k=5 PART 1: min pure mu(E) over ALL shapes (0 in E,|E|=5) by spread s=max(E), s<= {smax}")
    print("="*90)
    glob=(F(2),None); per={}
    for s in range(4, smax+1):
        loc=(F(2),None)
        for T in itertools.combinations(range(1,s),3):     # choose 3 from {1..s-1}; s itself forced in
            E=(0,)+T+(s,)
            m=mu_exact(E)
            if m<loc[0]: loc=(m,E)
        per[s]=loc
        if loc[0]<glob[0]: glob=loc
        print(f"  s={s:2d}: min mu = {str(loc[0]):>16s} = {float(loc[0]):.6f}  at E={loc[1]}")
    print(f"\n  *** GLOBAL min mu (spread<= {smax}) = {str(glob[0])} = {float(glob[0]):.6f}  at E={glob[1]} ***")
    return glob,per

# ===================================================================== RHO* over P
def sec_rho(worst_shapes):
    print("="*90)
    print("k=5 PART 2: rho*(P,E) EXACT over ALL P subset {1..13}, |P|=8, for worst shapes")
    print("  (THM-527: rho*>0 ==> M>=1/14). Find min rho* and whether any rho*=0.")
    print("="*90)
    overall=(F(2),None,None); zeros=0; total=0
    Ps=list(itertools.combinations(range(1,14),8))
    print(f"  #P (|P|=8) = {len(Ps)}; #worst shapes = {len(worst_shapes)}")
    for E in worst_shapes:
        mn=(F(2),None);
        gp_min=(F(2),None)
        for P in Ps:
            r,gp=rho_star_exact(list(P),list(E))
            total+=1
            if r==0: zeros+=1
            if r<mn[0]: mn=(r,P)
            if gp<gp_min[0]: gp_min=(gp,P)
        if mn[0]<overall[0]: overall=(mn[0],E,mn[1])
        print(f"  E={str(E):28s} mu={float(mu_exact(list(E))):.4f}  min_P rho* = {str(mn[0]):>16s}"
              f" = {float(mn[0]):.6f}  at P={mn[1]}  (min meas G_P={float(gp_min[0]):.5f})")
    print(f"\n  *** OVERALL min rho* over (worst shapes x all P) = {str(overall[0])} "
          f"= {float(overall[0]):.7f} ***")
    print(f"      at E={overall[1]} P={overall[2]}")
    print(f"      #(rho*=0) = {zeros} / {total}   ==>  {'POSITIVE FLOOR' if zeros==0 else 'HAS ZEROS!'}")
    return overall

# ===================================================================== TAIL (large spread)
def sec_tail():
    print("="*90)
    print("k=5 PART 3: large-spread tail.  Angle C claim: huge spread RAISES mu.")
    print("  Test families that GROW one offset to infinity; show mu stays bounded below.")
    print("="*90)
    # family A: {0,1,2,3,M}  -- 4 consecutive + 1 runaway
    print("  Family {0,1,2,3,M} (4 consecutive + runaway M):")
    for M in [4,5,6,8,12,20,50,200,1000,5000]:
        m=mu_exact([0,1,2,3,M]); print(f"     M={M:5d}: mu={str(m):>14s} = {float(m):.6f}")
    # family B: {0,1,M-1,M, 2M} growing AP-like
    print("  Family {0,2,4,6,M} (spread-out base + runaway):")
    for M in [8,10,14,20,50,200,1000]:
        m=mu_exact([0,2,4,6,M]); print(f"     M={M:5d}: mu={str(m):>14s} = {float(m):.6f}")
    # family C: scale the whole consecutive cluster {0,g,2g,3g,4g}
    print("  Family {0,g,2g,3g,4g} (consecutive scaled by g):")
    for g in [1,2,3,5,7,11,13,100]:
        m=mu_exact([0,g,2*g,3*g,4*g]); print(f"     g={g:4d}: mu={str(m):>14s} = {float(m):.6f}")

if __name__=="__main__":
    sec = sys.argv[1] if len(sys.argv)>1 else "sanity"
    if sec=="sanity": sec_sanity()
    elif sec=="enum":
        smax=int(sys.argv[2]) if len(sys.argv)>2 else 16
        sec_enum(smax)
    elif sec=="tail": sec_tail()
    elif sec=="all":
        sec_sanity()
        glob,per=sec_enum(int(sys.argv[2]) if len(sys.argv)>2 else 16)
        sec_tail()

# ===================================================================== RHO* full scan (called explicitly)
def sec_rhofull(smax_shapes):
    """For EVERY 5-shape with spread<=smax_shapes, and EVERY P subset {1..13} |P|=8,
       compute rho* exactly. Find global min and any zeros. This is the finite check."""
    print("="*90)
    print(f"k=5 PART 2 (FULL): min_P rho*(P,E) EXACT over ALL P (|P|=8) and ALL shapes spread<= {smax_shapes}")
    print("="*90)
    Ps=[list(p) for p in itertools.combinations(range(1,14),8)]
    print(f"  #P={len(Ps)}")
    overall=(F(2),None,None); zeros=[]; shape_cnt=0
    worst_per_shape=[]
    for s in range(4, smax_shapes+1):
        for T in itertools.combinations(range(1,s),3):
            E=(0,)+T+(s,)
            g=good_set_exact(list(E))
            mn=(F(2),None)
            for P in Ps:
                gp=safe_set(P)
                r=meas(intersect(g,gp))
                if r<mn[0]: mn=(r,tuple(P))
                if r==0: zeros.append((E,tuple(P)))
            shape_cnt+=1
            worst_per_shape.append((mn[0],E,mn[1]))
            if mn[0]<overall[0]: overall=(mn[0],E,mn[1])
    worst_per_shape.sort()
    print(f"  shapes checked: {shape_cnt}")
    print("  10 smallest min_P rho* (shape, P):")
    for r,E,P in worst_per_shape[:10]:
        print(f"    rho*={str(r):>16s} = {float(r):.7f}  E={E}  P={P}")
    print(f"\n  *** GLOBAL min rho* = {str(overall[0])} = {float(overall[0]):.7f} ***")
    print(f"      E={overall[1]}  P={overall[2]}")
    print(f"  #(rho*=0) = {len(zeros)}  ==> {'POSITIVE FLOOR (S3 k=5 closed on this range)' if not zeros else 'HAS ZEROS'}")
    for E,P in zeros[:8]: print(f"     ZERO at E={E} P={P}")
    return overall, zeros

if __name__=="__main__" and len(sys.argv)>1 and sys.argv[1]=="rhofull":
    sec_rhofull(int(sys.argv[2]) if len(sys.argv)>2 else 8)

# ===================================================================== generic-k rho* full scan
def sec_rhofull_k(kk, smax_shapes):
    """min_P rho* over ALL P subset{1..13} with |P|=13-kk, and ALL kk-shapes spread<=smax."""
    psz=13-kk
    print("="*90)
    print(f"k={kk} (FULL): min_P rho* EXACT over ALL P (|P|={psz}) and ALL shapes spread<= {smax_shapes}")
    print("="*90)
    Ps=[list(p) for p in itertools.combinations(range(1,14),psz)]
    print(f"  #P={len(Ps)}")
    overall=(F(2),None,None); zeros=[]; cnt=0; worst=[]
    for s in range(kk-1, smax_shapes+1):
        for T in itertools.combinations(range(1,s),kk-2):
            E=(0,)+T+(s,)
            g=good_set_exact(list(E))
            mn=(F(2),None)
            for P in Ps:
                r=meas(intersect(g,safe_set(P)))
                if r<mn[0]: mn=(r,tuple(P))
                if r==0: zeros.append((E,tuple(P)))
            cnt+=1; worst.append((mn[0],E,mn[1]))
            if mn[0]<overall[0]: overall=(mn[0],E,mn[1])
    worst.sort()
    print(f"  shapes checked: {cnt}")
    print("  8 smallest min_P rho*:")
    for r,E,P in worst[:8]:
        print(f"    rho*={str(r):>18s} = {float(r):.7f}  E={E}  P={P}")
    print(f"\n  *** GLOBAL min rho* (k={kk}) = {str(overall[0])} = {float(overall[0]):.7f} ***")
    print(f"      E={overall[1]} P={overall[2]}")
    print(f"  #(rho*=0)={len(zeros)} ==> {'POSITIVE FLOOR' if not zeros else 'HAS ZEROS'}")
    return overall, zeros

if __name__=="__main__" and len(sys.argv)>1 and sys.argv[1]=="rhok":
    sec_rhofull_k(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv)>3 else 9)
