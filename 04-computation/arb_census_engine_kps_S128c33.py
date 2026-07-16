#!/usr/bin/env python3
"""arb_census_engine_kps_S128c33.py -- kind-pasteur S128 cont.33.
Exhaustive n=3..7 (all tournaments via tilings): tau-vector vs odd-cycle census vs both drifts.
Q1: does (H,c5,c7) determine tau_tot / sorted tau-vector?  Q2: is tau_tot affine in census?
Q3: does tau LINEARIZE the H-drift (drift affine in (H,c5,c7,tau_tot))?
Q4: tau-drift E[Dtau_tot|T]: affine in tau_tot alone? in (tau_tot,H,c5,c7)?
Also: parity of tau_tot, transitive values, regular-stratum BEST cross-check."""
import sys, json
from math import comb, factorial
from fractions import Fraction as F
from itertools import permutations, combinations
sys.stdout.reconfigure(line_buffering=True)

def det_int(M):
    M=[r[:] for r in M]; n=len(M); sign=1; prev=1
    for k in range(n-1):
        if M[k][k]==0:
            for i in range(k+1,n):
                if M[i][k]: M[k],M[i]=M[i],M[k]; sign=-sign; break
            else: return 0
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[-1][-1]

def taus(B,n):
    L=[[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if u!=v and B[u][v]: L[v][v]+=1; L[v][u]-=1
    out=[]
    for r in range(n):
        idx=[i for i in range(n) if i!=r]
        out.append(det_int([[L[i][j] for j in idx] for i in idx]))
    return out

def ham(B,n):
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for S in range(1<<n):
        for v in range(n):
            c=dp[S][v]
            if not c or not (S>>v)&1: continue
            for u in range(n):
                if not (S>>u)&1 and B[v][u]: dp[S|1<<u][u]+=c
    return sum(dp[(1<<n)-1][v] for v in range(n))

def codd(B,n,Lc):
    if Lc>n: return 0
    tot=0
    for S in combinations(range(n),Lc):
        u=S[0]
        for perm in permutations(S[1:]):
            prev=u; ok=True
            for w in perm:
                if not B[prev][w]: ok=False; break
                prev=w
            if ok and B[prev][u]: tot+=1
    return tot

def affine_fit(rows, ncoef):
    """rows: list of (coefvec len ncoef, value). Solve exactly on independent subset, verify all.
    Returns coef list or None."""
    M=[]; piv_rows=[]
    for cv,val in rows:
        r=[F(x) for x in cv]+[F(val)]
        for pr in M:
            lead=next((i for i in range(ncoef) if pr[i]!=0),None)
            if lead is not None and r[lead]!=0:
                f=r[lead]/pr[lead]
                r=[a-f*b for a,b in zip(r,pr)]
        if any(r[i]!=0 for i in range(ncoef)):
            M.append(r)
        elif r[ncoef]!=0:
            return None   # inconsistent -> no affine law
        if len(M)==ncoef: break
    if len(M)<ncoef:
        pass  # underdetermined; still try: set free vars 0 via back-substitution on M
    # back-substitute
    sol=[F(0)]*ncoef
    for pr in reversed(M):
        lead=next(i for i in range(ncoef) if pr[i]!=0)
        s=pr[ncoef]-sum(pr[i]*sol[i] for i in range(lead+1,ncoef))
        sol[lead]=s/pr[lead]
    for cv,val in rows:
        if sum(F(c)*s for c,s in zip(cv,sol))!=val: return None
    return sol

n7drift=None
try:
    d=json.load(open(r"C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-math\f631d0eb-9f23-494b-bb86-e0501bc456e9\scratchpad\n7_drift_tuples.json"))
    n7drift={tuple(map(int,k.split(','))):F(v[0],v[1]) for k,v in d.items()}
except Exception as e:
    print("WARN: no n7 drift table:",e)

for n in (3,4,5,6,7):
    tiles=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1) if x-y>=2]
    m=len(tiles); m2=comb(n,2)
    fib_tot={}; fib_vec={}; rows_tau=[]; rows_Hdrift=[]; rows_tdrift=[]
    par={0:0,1:0}; par4={0:0,1:0,2:0,3:0}
    tdrift_by_state={}
    seenA={}
    for t in range(1<<m):
        B=[[False]*n for _ in range(n)]
        for k in range(2,n+1): B[k-1][k-2]=True
        for i,(x,y) in enumerate(tiles):
            if (t>>i)&1: B[x-1][y-1]=True
            else: B[y-1][x-1]=True
        H=ham(B,n); c5=codd(B,n,5); c7=codd(B,n,7)
        tv=taus(B,n); tt=sum(tv)
        key=(H,c5,c7)
        par[tt%2]+=1; par4[tt%4]+=1
        # Q1 fibers
        if key in fib_tot:
            if fib_tot[key]!=tt: fib_tot[key]='SPLIT'
        else: fib_tot[key]=tt
        sv=tuple(sorted(tv))
        if key in fib_vec:
            if fib_vec[key]!=sv and fib_vec[key]!='SPLIT': fib_vec[key]='SPLIT'
        else: fib_vec[key]=sv
        # tau-drift: exact sum of Dtau_tot over all flips
        s=F(0)
        for u in range(n):
            for v in range(n):
                if u!=v and B[u][v]:
                    B[u][v],B[v][u]=False,True
                    s+=sum(taus(B,n))-tt
                    B[u][v],B[v][u]=True,False
        td=F(s,m2)
        state=(tt,H,c5,c7)
        if state in tdrift_by_state:
            if tdrift_by_state[state]!=td: tdrift_by_state[state]='SPLIT'
        else: tdrift_by_state[state]=td
        # H-drift: n<=6 recompute; n=7 lookup
        if n<=6:
            sH=0
            for u in range(n):
                for v in range(n):
                    if u!=v and B[u][v]:
                        B[u][v],B[v][u]=False,True
                        sH+=ham(B,n)-H
                        B[u][v],B[v][u]=True,False
            hd=F(sH,m2)
        else:
            hd=n7drift.get(key) if n7drift else None
        if (key,tt) not in seenA:
            seenA[(key,tt)]=1
            rows_tau.append(((1,H,c5,c7),tt))
            if hd is not None: rows_Hdrift.append(((1,H,c5,c7,tt),hd))
            rows_tdrift.append(((1,tt),td))
        if n==7 and t%4096==0: print("  n=7 ... %d/32768"%t,flush=True)
    sp_tot=sum(1 for v in fib_tot.values() if v=='SPLIT')
    sp_vec=sum(1 for v in fib_vec.values() if v=='SPLIT')
    sp_td=sum(1 for v in tdrift_by_state.values() if v=='SPLIT')
    print("n=%d: census fibers %d | tau_tot splits %d | tau-vec splits %d"%(n,len(fib_tot),sp_tot,sp_vec))
    print("   tau_tot parity: even %d odd %d ; mod4 %s"%(par[0],par[1],dict(par4)))
    ft=affine_fit(rows_tau,4)
    print("   Q2 tau_tot affine in (1,H,c5,c7):", "YES "+str(ft) if ft else "NO")
    fh=affine_fit(rows_Hdrift,5) if rows_Hdrift else None
    print("   Q3 H-drift affine in (1,H,c5,c7,tau_tot):", "YES "+str(fh) if fh else "NO")
    f1=affine_fit(rows_tdrift,2)
    print("   Q4a tau-drift affine in (1,tau_tot):", "YES "+str(f1) if f1 else "NO")
    rows_td4=[(cv+ (val2,),val) for ((cv,val),( _,val2)) in zip([(r[0],r[1]) for r in rows_tdrift],[((0,),rt[1]) for rt in rows_tau])]
    # rebuild cleanly: (1,tau,H,c5,c7)
    rows_td5=[]
    seen2=set()
    for (cv,tt2) in rows_tau:
        pass
    # simpler: recollect from tdrift_by_state (only non-split states)
    rows_td5=[((1,st[0],st[1],st[2],st[3]),v) for st,v in tdrift_by_state.items() if v!='SPLIT']
    print("   Q4-pre: does (tau_tot,H,c5,c7) determine tau-drift? splits=%d of %d states"%(sp_td,len(tdrift_by_state)))
    f5=affine_fit(rows_td5,5)
    print("   Q4b tau-drift affine in (1,tau_tot,H,c5,c7):", "YES "+str(f5) if f5 else "NO")
print("DONE")
