#!/usr/bin/env python3
"""
lrc14_sector_lp_moments_macmini_0618s7.py   (mac-mini-2026-06-18-S7)

ANGLE D continued.  The (S_1,S_2) LP and even-Bonferroni B_2/B_4 do NOT close the
dangerous rows k=8,9,10.  Question: with how many EXACT factorial moments S_1..S_R
does the moment-LP upper bound on p_0 fall <= cap_k?  And: which moments carry the
binding dual constraint (-> the certificate structure)?

p_t = meas{exactly t of sectors {1..6} are missed},  t=0..6.
S_r = E[C(N,r)] = sum_t C(t,r) p_t,  exactly computable.
moment-LP_R:  max p_0  s.t.  sum_t C(t,r) p_t = S_r (r=0..R),  p_t>=0.
As R->6 the constraints pin p exactly so LP-> meas(S7).  We solve EXACTLY via
vertex enumeration (support size <= R+1).

We tabulate moment-LP_R(consec_k) for R=2..6 vs cap_k, and identify the SMALLEST R
that closes each k.  Then we read the optimal DUAL y=(y_0..y_R):
   p_0 <= y_0 + sum_{r=1}^R y_r S_r     with    y_0 + sum_r y_r C(t,r) >= [t=0]  for all t.
That linear functional of the exact moments IS the certificate for that k (modulo a
monotonicity statement S_r(E)<=S_r(consec) which we also test).
"""
import sys, itertools
from math import comb
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
H=F(1,14)

def measJ(A, E):
    Aset=set(A); E=sorted(set(E))
    if 0 in Aset: return F(0)
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        if all(int(((e*xm)%1)*7) not in Aset for e in E): total+=x1-x0
    return total
def Svec(E):
    S=[F(0)]*7; S[0]=F(1)
    for r in range(1,7):
        S[r]=sum(measJ(A,E) for A in itertools.combinations(range(1,7), r))
    return S
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set(int(((e*xm)%1)*7) for e in E)
        if len(secs)==7: total+=x1-x0
    return total

def danger(u):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-H/u)%1; b=(c+H/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def mgmerge(iv):
    iv=sorted(iv); o=[]
    for a,b in iv:
        if o and a<=o[-1][1]: o[-1]=(o[-1][0],max(o[-1][1],b))
        else: o.append((a,b))
    return o
def measGP(P):
    if not P: return F(1)
    dz=mgmerge([iv for u in P for iv in danger(u)]); s=F(0); prev=F(0)
    for a,b in dz:
        if a>prev: s+=a-prev
        prev=max(prev,b)
    if prev<1: s+=1-prev
    return s
import functools
@functools.lru_cache(None)
def cap(k):
    psz=13-k
    if psz==0: return F(1)
    return min(measGP(P) for P in itertools.combinations(range(1,14),psz))

# ---- exact moment-LP: max p_0, support t in 0..6, moments S_r r=0..R fixed ----
def solve3(rows, b):
    n=len(rows)
    M=[rows[i][:]+[b[i]] for i in range(n)]
    for c in range(n):
        piv=None
        for r in range(c,n):
            if M[r][c]!=0: piv=r;break
        if piv is None: return None
        M[c],M[piv]=M[piv],M[c]
        pv=M[c][c]; M[c]=[x/pv for x in M[c]]
        for r in range(n):
            if r!=c and M[r][c]!=0:
                f=M[r][c]; M[r]=[M[r][j]-f*M[c][j] for j in range(n+1)]
    return [M[i][n] for i in range(n)]

def moment_lp_max_p0(S, R):
    """max p_0 over prob measures on t=0..6 with C(t,r) moments = S_r, r=0..R.
       Optimum at vertex with support size <= R+1.  Enumerate supports."""
    ts=list(range(7))
    best=None
    for supp in itertools.combinations(ts, R+1):
        rows=[[F(comb(t,r)) for t in supp] for r in range(R+1)]  # r=0 -> all 1
        b=[S[r] for r in range(R+1)]
        sol=solve3(rows,b)
        if sol is None: continue
        if any(x<F(0) for x in sol): continue
        p0=sum(sol[i] for i in range(len(supp)) if supp[i]==0)
        if best is None or p0>best: best=p0
    # also smaller supports (size<=R): pad by trying all sizes
    for sz in range(1,R+1):
        for supp in itertools.combinations(ts, sz):
            # need moments r=0..R consistent; overdetermined -> solve first sz, check rest
            rows=[[F(comb(t,r)) for t in supp] for r in range(sz)]
            b=[S[r] for r in range(sz)]
            sol=solve3(rows,b)
            if sol is None: continue
            if any(x<F(0) for x in sol): continue
            ok=all(sum(sol[i]*F(comb(supp[i],r)) for i in range(sz))==S[r] for r in range(sz,R+1))
            if not ok: continue
            p0=sum(sol[i] for i in range(sz) if supp[i]==0)
            if best is None or p0>best: best=p0
    return best

print("="*100)
print("moment-LP_R(consec_k) = exact max p_0 given S_1..S_R  vs cap_k.  Smallest R that closes k.")
print("="*100)
hdr=f"{'k':>3}" + "".join(f"{'LP_R='+str(R):>11}" for R in range(2,7)) + f"{'meas S7':>11}{'cap_k':>11}{'minR<=cap':>11}"
print(hdr)
results={}
for k in range(8,14):
    E=list(range(k)); S=Svec(E); ck=cap(k); s7=measS7(E)
    row=f"{k:>3}"
    minR=None
    for R in range(2,7):
        lp=moment_lp_max_p0(S,R)
        results[(k,R)]=lp
        row+=f"{float(lp):>11.5f}" if lp is not None else f"{'--':>11}"
        if minR is None and lp is not None and lp<=ck: minR=R
    row+=f"{float(s7):>11.5f}{float(ck):>11.5f}{str(minR):>11}"
    print(row)

print("\nInterpretation: minR = fewest exact factorial moments whose LP upper bound on")
print("meas(S7) already lies <= cap_k.  That is the degree of the Bonferroni/dual certificate.")

# ---- read the dual at the closing R for k=8 (the hardest) ----
print("\n"+"="*100)
print("DUAL certificate at the closing degree for each k (the binding Bonferroni-type bound)")
print("="*100)
def dual_from_primal(S, R):
    """Find dual y_0..y_R s.t. p_0 <= y_0+sum y_r S_r and it's tight at LP opt.
       Dual feasibility: g(t):=sum_r y_r C(t,r) >= [t==0] for t=0..6, minimize sum y_r S_r.
       Solve by enumerating which constraints are tight (support of optimal primal)."""
    # The optimal primal support gives tight dual constraints. Reconstruct from support.
    ts=list(range(7))
    best=None;bestval=None
    for supp in itertools.combinations(ts, R+1):
        rows=[[F(comb(t,r)) for t in supp] for r in range(R+1)]
        b=[S[r] for r in range(R+1)]
        sol=solve3(rows,b)
        if sol is None or any(x<F(0) for x in sol): continue
        p0=sum(sol[i] for i in range(len(supp)) if supp[i]==0)
        # dual: g(t)=[t==0] on support (complementary slackness), solve for y
        drows=[[F(comb(t,r)) for r in range(R+1)] for t in supp]
        db=[F(1) if t==0 else F(0) for t in supp]
        y=solve3(drows,db)
        if y is None: continue
        # check dual feasibility g(t)>=[t==0] all t
        feas=all(sum(y[r]*F(comb(t,r)) for r in range(R+1)) >= (F(1) if t==0 else F(0)) for t in ts)
        val=sum(y[r]*S[r] for r in range(R+1))
        if feas and (bestval is None or val<bestval):
            if val==p0:
                bestval=val;best=(y,p0)
    return best

for k in range(8,14):
    E=list(range(k)); S=Svec(E); ck=cap(k)
    # find minR closing
    minR=None
    for R in range(2,7):
        if results[(k,R)] is not None and results[(k,R)]<=ck: minR=R;break
    if minR is None:
        print(f"  k={k}: no R<=6 closes (should not happen at R=6).");continue
    d=dual_from_primal(S,minR)
    if d is None:
        print(f"  k={k}: minR={minR}, dual reconstruction failed (degenerate).");continue
    y,p0=d
    terms=" + ".join(f"({str(y[r])})*S_{r}" for r in range(minR+1))
    print(f"  k={k}: degree R={minR}.  meas(S7) <= {terms}")
    print(f"         = {float(p0):.5f} <= cap_k={float(ck):.5f}   [dual y = {[str(x) for x in y]}]")

print("\nDONE.")
