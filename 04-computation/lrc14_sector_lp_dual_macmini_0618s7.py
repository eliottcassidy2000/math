#!/usr/bin/env python3
"""
lrc14_sector_lp_dual_macmini_0618s7.py   (mac-mini-2026-06-18-S7)

ANGLE D, the actual LP + DUAL CERTIFICATE for  meas(S7(E)) <= cap_k.

IDEA (Bonferroni / Galambos LP).  meas(S7) = int prod_j [sector j hit].  Write
N(x) = # MISSED sectors = sum_{j=0}^6 1[sector j not hit at x].  Then
   S7 hit  <=>  N(x)=0.
meas(S7) = meas{N=0}.  Let p_t = meas{N(x)=t}, t=0..6 (a probability vector, since
total measure 1).  Define the binomial-moment (factorial-moment) data
   S_r = sum_{|A|=r, A subset {0..6}} meas{ all sectors in A missed }
       = E[ C(N, r) ]  = sum_t C(t,r) p_t.
These S_r are EXACTLY computable: S_r = sum_{|A|=r} J(A,E) and (since 0 in A kills it)
only A subset {1..6} survive, so S_r = sum_{|A|=r, A subset {1..6}} J(A,E).
We have S_0 = 1.  Sector 0 is ALWAYS hit (e=0), so N(x) <= 6 and p_t=0 for t s.t.
sector-0 missing -> effectively N counts missed sectors among {1..6}, N in {0..6}.

meas(S7) = p_0 = 1 - (p_1+...+p_6).  The S_r give linear constraints
   sum_t C(t,r) p_t = S_r,  p_t>=0,  sum p_t = 1.
LP:  maximize p_0  subject to these.  Its DUAL gives an UPPER BOUND on meas(S7)
that is a POSITIVE/SIGNED combination of the S_r --- i.e. a Bonferroni inequality
   meas(S7) = p_0 <= sum_r y_r S_r   with the y_r the dual multipliers.
If for the EXTREMAL (consec) E this LP optimum (a clean linear comb of exact S_r)
is <= cap_k, THAT linear combination, with the proven facts S_r(E)<=S_r(consec)
(monotonicity, if it holds), is the certificate.

We:
 (A) compute S_0..S_6 exactly for E, solve the LP max p_0 (and min p_0) exactly
     via the vertices (it's a small transportation-type LP), read the dual y_r;
 (B) check the LP upper bound vs cap_k for consec_k, k=8..13;
 (C) test the MONOTONICITY hypothesis: is S_r(E) <= S_r(consec_k) for all E? (so the
     consec LP value dominates) --- the missing ingredient for a real proof.
"""
import sys, itertools, random
from math import comb, gcd
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(717)
H=F(1,14)

# ---- exact J(A,E) (forbidden sectors A; 0 in A => 0) ----
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
    """S_r = sum_{|A|=r, A subset {1..6}} J(A,E), r=0..6.  S_0=1."""
    S=[F(0)]*7
    S[0]=F(1)
    for r in range(1,7):
        s=F(0)
        for A in itertools.combinations(range(1,7), r):
            s+=measJ(A,E)
        S[r]=s
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

# ---- recover p_t EXACTLY from S_r via inverse of the C(t,r) (Stirling) matrix ----
# p_t = sum_{r>=t} (-1)^{r-t} C(r,t) S_r   (standard binomial-moment inversion).
def pvec(S):
    p=[F(0)]*7
    for t in range(7):
        s=F(0)
        for r in range(t,7):
            s += F((-1)**(r-t)*comb(r,t))*S[r]
        p[t]=s
    return p

# ---- cap_k machinery ----
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

# ============================================================
# (A) the LP: max p_0 s.t. sum_t C(t,r)p_t=S_r (r=0..R), p>=0.
#     Solve EXACTLY by enumerating basic feasible solutions for the moment LP.
#     We use the classic result: with moments S_0..S_R fixed, max/min of p_0 is
#     attained at a measure supported on <= R+1 atoms.  For a clean Bonferroni
#     bound we use ONLY S_0,S_1,S_2 (pairwise data) and the analytic dual:
#        p_0 <= 1 - S_1 + S_2        (Bonferroni, always valid)   [degree-2]
#        p_0 <= 1 - S_1 + S_2 - S_3 + ...                         (full IE = exact)
#     We report the truncated Bonferroni upper bounds B_R = sum_{r<=R}(-1)^r S_r
#     for EVEN R (these are guaranteed UPPER bounds on p_0 by Bonferroni).
# ============================================================
def bonferroni_upper(S, R):
    # B_R = sum_{r=0}^{R} (-1)^r S_r ; for EVEN R this is >= p_0 (upper bound).
    return sum(F((-1)**r)*S[r] for r in range(R+1))

print("="*92)
print("LP / Bonferroni dual certificate:  meas(S7)=p_0 <= B_R (even R) ;  exact = B_6")
print("  Columns: B_0=1, B_2=1-S1+S2, B_4, B_6(=exact meas S7), cap_k, and which B_R<=cap_k")
print("="*92)
print(f"{'k':>3}{'B2':>11}{'B4':>11}{'B6=meas':>11}{'cap_k':>11}{'B2<=cap?':>10}{'B4<=cap?':>10}")
for k in range(8,14):
    E=list(range(k)); S=Svec(E); ck=cap(k)
    B2=bonferroni_upper(S,2); B4=bonferroni_upper(S,4); B6=bonferroni_upper(S,6)
    assert B6==measS7(E), (k,B6,measS7(E))
    f2 = "YES" if B2<=ck else "no"
    f4 = "YES" if B4<=ck else "no"
    print(f"{k:>3}{float(B2):>11.5f}{float(B4):>11.5f}{float(B6):>11.5f}{float(ck):>11.5f}{f2:>10}{f4:>10}")

print("\nNOTE: B_2 = 1 - S_1 + S_2 is the degree-2 (pairwise) Bonferroni UPPER bound,")
print("valid for EVERY E by the inclusion-exclusion sign rule (even truncation).")
print("If B_2(consec_k) <= cap_k AND S_1,S_2 are controlled, that closes k by pairwise data.")
print()

# ============================================================
# (B) full vertex LP: exact max p_0 given S_0,S_1,S_2 only, p supported on {0..6}.
#     This is tighter than Bonferroni B_2 in general.  Solve by LP duality:
#     max p_0 s.t. sum p=1, sum t p_t=S_1, sum C(t,2)p_t=S_2, p>=0.
#     Dual: min y0 + y1 S_1 + y2 S_2  s.t. y0 + y1 t + y2 C(t,2) >= [t==0] for all t=0..6.
#     We grid the dual over small rational y to find the exact optimum (3 vars,
#     7 linear constraints -> optimum at a vertex; enumerate vertex triples).
# ============================================================
def lp_max_p0(S1, S2):
    """exact max of p_0 over prob vectors on t=0..6 with moments S1,S2 (or None if infeasible)."""
    import itertools as it
    ts=list(range(7))
    cons=[]  # (a0,a1,a2,rhs) meaning a0+a1 t+a2 C(t,2) ... we solve primal vertices.
    # Primal vertex: support of size <=3.  Enumerate 3-subsets, solve 3x3 for p, check >=0.
    best=None
    for sub in it.combinations(ts,3):
        # variables p_{sub}, eqns: sum p=1, sum t p =S1, sum C(t,2)p=S2
        a=[[F(1)]*3,[F(t) for t in sub],[F(comb(t,2)) for t in sub]]
        b=[F(1),S1,S2]
        # solve 3x3
        import copy
        M=[row[:]+[b[i]] for i,row in enumerate(a)]
        ok=True
        for c in range(3):
            piv=None
            for r in range(c,3):
                if M[r][c]!=0: piv=r;break
            if piv is None: ok=False;break
            M[c],M[piv]=M[piv],M[c]
            pv=M[c][c]; M[c]=[x/pv for x in M[c]]
            for r in range(3):
                if r!=c and M[r][c]!=0:
                    f=M[r][c]; M[r]=[M[r][j]-f*M[c][j] for j in range(4)]
        if not ok: continue
        p=[M[i][3] for i in range(3)]
        if any(x<0 for x in p): continue
        p0 = sum(p[i] for i in range(3) if sub[i]==0)
        if best is None or p0>best: best=p0
    # also size-1,2 supports (degenerate); covered by allowing repeated via 3-subsets? not fully.
    for sub in it.combinations(ts,2):
        a=[[F(1)]*2,[F(t) for t in sub]]
        # 2 eqns (sum=1, mean=S1); need S2 to match exactly:
        # solve p from first two, check 2nd moment
        d=a[1][0]-a[1][1]
        if d==0: continue
        p1=(S1-a[1][1])/d; p0v=F(1)-p1
        if p0v<0 or p1<0: continue
        s2=p0v*F(comb(sub[0],2))+p1*F(comb(sub[1],2))
        if s2!=S2: continue
        p0=(p0v if sub[0]==0 else F(0))+(p1 if sub[1]==0 else F(0))
        if best is None or p0>best: best=p0
    return best

print("="*92)
print("(B) exact LP max p_0 using ONLY (S_1,S_2)  [3-moment relaxation]  vs cap_k")
print("="*92)
print(f"{'k':>3}{'S_1':>11}{'S_2':>11}{'LP max p0':>12}{'meas S7':>11}{'cap_k':>11}{'LP<=cap?':>10}")
for k in range(8,14):
    E=list(range(k)); S=Svec(E); ck=cap(k)
    lp=lp_max_p0(S[1],S[2]); s7=measS7(E)
    flag="YES" if (lp is not None and lp<=ck) else "no"
    lps = f"{float(lp):.5f}" if lp is not None else "infeas"
    print(f"{k:>3}{float(S[1]):>11.5f}{float(S[2]):>11.5f}{lps:>12}{float(s7):>11.5f}{float(ck):>11.5f}{flag:>10}")

print("\nDONE.")
