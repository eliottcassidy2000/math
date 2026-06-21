#!/usr/bin/env python3
r"""
lrc14_threadB_fastcheck_macmini_0620sB.py   (mac-mini, Thread B stage 5b -- FAST deliverable)

FAST exact exhaustive measS7<=cap over ALL primitive E (0 in E, |E|=k), span<=W.
Uses INTEGER arithmetic over a common denominator (no Fraction objects in the hot loop):
  - breakpoints of speed e are m/(7e); put everything over D = 7*lcm(speeds).
  - on each cell, the sector of e at x=p/D is floor(7*e*p/D) mod 7, computed in ints.
  - measS7 = (sum of cell widths where all 7 sectors hit) / D, kept as an integer
    numerator over D; compared to cap exactly via cross-multiplication.

This is the SAME exact value as the Fraction breakpoint engine (cross-checked here),
but ~50-100x faster, making W=20 (k=8), W=17 (k=9), W=16 (k=10) feasible.

DELIVERABLE: worst measS7 over the whole span-window, #over-cap (must be 0), and the
worst (tightest) shapes with their relation-height lambda_max.
"""
import sys, itertools, math
from math import gcd, comb
from fractions import Fraction as F
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def gcd_all(xs): return reduce(gcd, [abs(x) for x in xs if x != 0], 0)
def lcm(a,b): return a*b//gcd(a,b)
def banner(t): print("\n"+"="*82+f"\n{t}\n"+"="*82)

CAPS = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91)}

def measS7_int(speeds):
    """Exact measS7 of E={0}+speeds, returned as Fraction.  Integer-grid common-refinement.
       speeds: tuple of distinct positive ints."""
    # common denominator D = 7 * lcm(speeds)
    L = 1
    for e in speeds: L = lcm(L, e)
    D = 7*L
    # breakpoints (as integer numerators over D): for speed e, x=m/(7e)=m*L/(e)/D ... wait
    # m/(7e) over D=7L: numerator = m*(L/e). So breakpoints = { m*(L//e) : m=0..7e } for each e.
    bset = {0, D}
    for e in speeds:
        step = L//e   # since 7e * (m/(7e)) ; numerator over D is m*(7L/(7e))=m*(L/e)
        for m in range(0, 7*e+1):
            bset.add(m*step)
    bps = sorted(bset)
    num = 0
    sevenfull = (1<<7) - 1
    for i in range(len(bps)-1):
        a = bps[i]; b = bps[i+1]
        if b<=a: continue
        # midpoint numerator over 2D: pm/(2D) where pm = a+b
        pm = a+b   # over 2D
        mask = 1   # bit 0: the observer e=0 in E always sits in sector 0 (floor(7*frac(0))=0)
        # sector of speed e at x=pm/(2D): floor(7*e*pm/(2D)) mod 7
        for e in speeds:
            # 7*e*x = 7*e*pm/(2D) = 7*e*pm/(2*7*L) = e*pm/(2L)
            sec = (e*pm)//(2*L) % 7
            mask |= (1<<sec)
            if mask==sevenfull: break
        if mask==sevenfull:
            num += (b-a)
    return F(num, D)

# cross-check against Fraction engine
def measS7_frac(E):
    E=sorted(set(E)); bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for m in range(0,7*abs(e)+1): bps.add(F(m,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        secs={int(((e*xm)%1)*7) for e in E}
        if len(secs)==7: tot+=x1-x0
    return tot

# lambda_max for the tight set only
def mat_rank_Q(rows):
    M=[[F(x) for x in row] for row in rows]; r0=0; rank=0; ncol=len(M[0])
    for c in range(ncol):
        piv=None
        for r in range(r0,len(M)):
            if M[r][c]!=0: piv=r; break
        if piv is None: continue
        M[r0],M[piv]=M[piv],M[r0]; pv=M[r0][c]
        for r in range(len(M)):
            if r!=r0 and M[r][c]!=0:
                f=M[r][c]/pv; M[r]=[M[r][j]-f*M[r0][j] for j in range(ncol)]
        r0+=1; rank+=1
        if r0==len(M): break
    return rank
def lambda_max(speeds, Hcap=20):
    d=len(speeds); need=d-1; found=[]; mins=[]
    for H in range(1,Hcap+1):
        for n in itertools.product(range(-H,H+1),repeat=d):
            if max(abs(x) for x in n)!=H: continue
            if sum(n[i]*speeds[i] for i in range(d))!=0: continue
            if not found:
                if any(x!=0 for x in n): found.append([F(x) for x in n]); mins.append(H)
            else:
                if mat_rank_Q(found+[[F(x) for x in n]])>len(found):
                    found.append([F(x) for x in n]); mins.append(H)
            if len(found)>=need: return max(mins)
    return None

def crosscheck():
    banner("CROSS-CHECK: measS7_int == measS7_frac")
    ok=True
    for E in [[0,1,2,3,4,5,6,7],[0,2,3,5,8,11,13,17],[0,1,3,7,12,20],[0,1,2,3,4,5,6,8,9,10]]:
        sp=tuple(e for e in E if e!=0)
        a=measS7_int(sp); b=measS7_frac(E)
        if a!=b: ok=False; print(f"  MISMATCH {E}: int={a} frac={b}")
    print(f"  int engine matches Fraction engine on 4 shapes: {ok}")
    return ok

def exhaustive(k, W):
    banner(f"EXHAUSTIVE measS7<=cap : all primitive E, |E|={k}, span<=W={W}")
    cap=CAPS[k]; capnum,capden=cap.numerator,cap.denominator
    consec=measS7_int(tuple(range(1,k)))
    print(f"  cap_{k}={cap}={float(cap):.5f}; consec_{k}={float(consec):.5f}")
    import time; t0=time.time()
    cnt=0; over=0; worst=F(0); argworst=None; near=[]
    for combo in itertools.combinations(range(1,W+1), k-1):
        if gcd_all(combo)!=1: continue
        cnt+=1
        m=measS7_int(combo)
        if m>worst: worst=m; argworst=(0,)+combo
        if m>cap: over+=1
        if m*F(50) > cap*F(50)-F(1):  # within ~0.02 of cap
            if m > cap - F(1,50): near.append(((0,)+combo, m))
    dt=time.time()-t0
    print(f"  #primitive E (span<=W): {cnt}   ({dt:.0f}s)")
    print(f"  WORST measS7 = {worst} = {float(worst):.6f}  at {argworst}")
    print(f"  consec is the worst? {argworst==tuple(range(k))}")
    print(f"  # over cap_{k}: {over}   =>  measS7 <= cap_{k} for ALL: {over==0}")
    print(f"  margin cap - worst = {float(cap-worst):.6f}")
    print(f"  # shapes within 0.02 of cap: {len(near)}")
    for E,m in sorted(near,key=lambda t:-t[1])[:10]:
        lm=lambda_max(tuple(e for e in E if e!=0),Hcap=10)
        print(f"      measS7={float(m):.5f} span={max(E)} lambda_max={lm}  E={list(E)}")
    return cnt,over,worst

if __name__=="__main__":
    print("#"*82)
    print("# THREAD B stage 5b -- FAST exhaustive span-window finite check")
    print("#"*82)
    crosscheck()
    exhaustive(8, 20)
    exhaustive(9, 17)
    exhaustive(10, 16)
    print("\nDONE (Thread B stage 5b).")
