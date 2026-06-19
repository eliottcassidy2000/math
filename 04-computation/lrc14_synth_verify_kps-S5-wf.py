#!/usr/bin/env python3
"""
lrc14_synth_verify_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5, SYNTHESIS VERIFIER)

Independent re-implementation of mu_theta and the decisive contested numbers.
Two INDEPENDENT exact methods for mu_theta(E):
  METHOD A (order-cell, affine-gap): like the canon engine.
  METHOD B (direct breakpoint union): for each ORDERED pair (i,j) and each gap-threshold
            crossing, build the full breakpoint set and integrate the indicator 1[maxgap>theta]
            by exact midpoint evaluation per atomic cell (maxgap is NOT affine globally but on
            each atomic cell the *winning* gap is affine, and the cell is small enough that
            1[maxgap>theta] is constant on the open cell IF we include gap=theta crossings).
We cross-check A==B on every E.

Then we test:
  (1) The SECOND-VERDICT counterexamples claiming mu(k=13) << 1/14.
  (2) Canon anchors mu_4=19/21, mu_5=9/14, mu_6=4/7, mu_13(consec)=829/4620.
  (3) A float sandwich with N large to confirm.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

TWO7 = F(2,7)
ONE7 = F(1,7)

def merge(iv):
    iv = sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]:
            out[-1]=(out[-1][0], max(out[-1][1], b))
        else:
            out.append((a,b))
    return out

def meas(arcs):
    return sum((b-a for a,b in arcs), F(0))

# ---------- METHOD A: order-cell affine-gap (independent reimpl) ----------
def good_A(E, theta):
    E=sorted(set(E)); k=len(E)
    if k==1: return [(F(0),F(1))]
    diffs=set()
    for a in range(k):
        for b in range(a+1,k):
            diffs.add(E[b]-E[a])
    bps={F(0),F(1)}
    for d in diffs:
        for m in range(0,d+1):
            bps.add(F(m,d))
    bps=sorted(x for x in bps if 0<=x<=1)
    good=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        pts=sorted(((E[t]*xm)%1, E[t]) for t in range(k))
        order=[e for _,e in pts]
        floors=[int((e*xm)//1) for e in order]
        for idx in range(k):
            e_cur=order[idx]; f_cur=floors[idx]
            if idx<k-1:
                e_nx=order[idx+1]; f_nx=floors[idx+1]; wrap=F(0)
            else:
                e_nx=order[0]; f_nx=floors[0]; wrap=F(1)
            A=F(e_nx-e_cur); Cc=F(f_cur-f_nx)+wrap
            if A==0:
                if Cc>theta: good.append((x0,x1))
                continue
            xb=(theta-Cc)/A
            if A>0: lo=max(x0,xb); hi=x1
            else:   lo=x0; hi=min(x1,xb)
            if lo<hi: good.append((lo,hi))
    return merge(good)

def mu_A(E, theta=TWO7):
    return meas(good_A(E, theta))

# ---------- METHOD B: full breakpoint set incl gap=theta, constant per atomic cell ----------
def good_B(E, theta):
    """Build ALL breakpoints: order-changes m/d AND every gap=theta crossing for every ordered
    pair as a potential consecutive-in-order pair. Then on each atomic cell the maxgap>theta
    indicator is CONSTANT (we evaluate at the midpoint exactly via rationals)."""
    E=sorted(set(E)); k=len(E)
    if k==1: return [(F(0),F(1))]
    diffs=set()
    for a in range(k):
        for b in range(a+1,k):
            diffs.add(E[b]-E[a])
    bps={F(0),F(1)}
    # order-change points
    for d in diffs:
        for m in range(0,d+1):
            bps.add(F(m,d))
    # gap=theta crossings: for every ordered pair (i,j) i!=j, the gap from point i to point j
    # (when they are cyclically consecutive) is (e_j-e_i)x + integer +maybe wrap = theta.
    # We add ALL candidate crossings (e_j-e_i)x + c = theta for the relevant integer range of c.
    for i in range(k):
        for j in range(k):
            if i==j: continue
            slope=E[j]-E[i]   # can be negative
            if slope==0: continue
            # gap value = frac(e_j x) - frac(e_i x) mod 1; as x ranges, candidate linear pieces
            # are slope*x + c for integer c. Solve slope*x + c = theta and slope*x + c = theta with wrap.
            # c ranges so that x in [0,1]: x=(theta-c)/slope in [0,1].
            # cover a generous integer range for c.
            lo_c = -abs(slope)-2
            hi_c = abs(slope)+2
            for c in range(lo_c, hi_c+1):
                denom=slope
                xb=(theta-F(c))/denom
                if 0<=xb<=1: bps.add(xb)
    bps=sorted(set(x for x in bps if 0<=x<=1))
    good=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        fr=sorted(set((e*xm)%1 for e in E))
        if len(fr)==1:
            mg=F(1)
        else:
            gaps=[fr[t+1]-fr[t] for t in range(len(fr)-1)]+[fr[0]+1-fr[-1]]
            mg=max(gaps)
        if mg>theta:
            good.append((x0,x1))
    return merge(good)

def mu_B(E, theta=TWO7):
    return meas(good_B(E, theta))

# ---------- float sandwich ----------
def mu_float(E, theta=2/7, N=4_000_000):
    import random
    theta=float(theta)
    cnt=0
    Ef=[float(e) for e in E]
    rnd=random.Random(12345)
    for _ in range(N):
        x=rnd.random()
        pts=sorted((e*x)%1.0 for e in Ef)
        if len(pts)==1:
            mg=1.0
        else:
            mg=0.0
            for t in range(len(pts)-1):
                d=pts[t+1]-pts[t]
                if d>mg: mg=d
            d=pts[0]+1.0-pts[-1]
            if d>mg: mg=d
        if mg>theta: cnt+=1
    return cnt/N

def main():
    print("=== ANCHOR CHECK (theta=2/7) ===")
    anchors={
        "mu4_consec":([0,1,2,3], F(19,21)),
        "mu5_consec":([0,1,2,3,4], F(9,14)),
        "mu6_consec":([0,1,2,3,4,5], F(4,7)),
        "mu7_perf":([0,2,3,4,5,6,8], F(13,35)),
        "mu13_consec":(list(range(13)), F(829,4620)),
        "L1_scale":([0,7,14,21,28,35], F(4,7)),
    }
    for name,(E,exp) in anchors.items():
        a=mu_A(E); b=mu_B(E)
        ok = (a==exp) and (b==exp)
        print(f"  {name:14s} A={a} B={b} expected={exp}  {'OK' if ok else '*** MISMATCH ***'}")

    print("\n=== CONTESTED SECOND-VERDICT COUNTEREXAMPLES (k=13, claim mu << 1/14) ===")
    print(f"  1/14 = {float(F(1,14)):.6f}")
    ces = {
        "ce1":[0,12,17,20,24,26,27,28,32,36,37,47,60],
        "ce2":[0,11,21,28,33,35,36,37,39,42,45,49,62],
        "ce3":[0,7,11,20,21,23,28,33,35,36,39,42,45],
        "ce4":[0,7,13,14,15,16,17,19,21,24,27,29,40],
    }
    claimed={"ce1":F(3736,85785),"ce2":F(1314101,28198716),"ce3":F(23059,412335),"ce4":F(3303,52780)}
    for name,E in ces.items():
        a=mu_A(E); b=mu_B(E)
        fl=mu_float(E, N=2_000_000)
        cl=claimed[name]
        agree_AB = (a==b)
        agree_claim = (a==cl)
        below = a < F(1,14)
        print(f"  {name}: A={a} (~{float(a):.6f}) B={b} float~{fl:.6f}")
        print(f"        claimed={cl} (~{float(cl):.6f})  A==B:{agree_AB}  A==claimed:{agree_claim}  A<1/14:{below}")

    print("\n=== mu_{1/7} of these contested sets (the CORRECT global-witness object) ===")
    for name,E in ces.items():
        a17=mu_A(E, ONE7); b17=mu_B(E, ONE7)
        print(f"  {name}: mu_1/7 A={a17} (~{float(a17):.6f}) B={b17}  A==B:{a17==b17}  >=consec13(477/1078={float(F(477,1078)):.4f}):{a17>=F(477,1078)}")

if __name__=="__main__":
    main()
