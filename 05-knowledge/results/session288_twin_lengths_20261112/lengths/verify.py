"""Exact matching-mode and length-rigidity controls. Standard library only."""
from itertools import combinations
from math import comb, factorial
from pathlib import Path
import json, random, sys
HERE=Path(__file__).resolve().parent

def need(q,msg):
    if not q: raise RuntimeError(msg)

def table(rows):
    n=len(rows); N=1<<n; full=N-1
    start=[[0]*n for _ in range(N)]; end=[[0]*n for _ in range(N)]
    for mask in range(1,N):
        for v in range(n):
            if not mask>>v&1: continue
            rem=mask^(1<<v)
            if not rem:start[mask][v]=end[mask][v]=1
            else:
                start[mask][v]=sum(start[rem][u] for u in range(n) if rows[v]>>u&1)
                end[mask][v]=sum(end[rem][u] for u in range(n) if rows[u]>>v&1)
    B=[[0]*n for _ in range(n)]
    for mask in range(1,full):
        rem=full^mask
        for i in range(n):
            if not end[mask][i]:continue
            for j in range(n):B[i][j]+=end[mask][i]*start[rem][j]
    return sum(end[full]),B

def inflate(rows,v,m,independent=False):
    labels=[i for i in range(len(rows)) if i!=v]+[v]*m
    out=[]
    for i,a in enumerate(labels):
        out.append(sum(1<<j for j,b in enumerate(labels) if i!=j and ((a==b and not independent and i<j) or (a!=b and rows[a]>>b&1))))
    return out,labels

def matching(rows,v):
    out=[i for i in range(len(rows)) if rows[v]>>i&1]
    inside=[i for i in range(len(rows)) if i!=v and rows[i]>>v&1]
    def rec(k,used):
        if k==len(out):return 0
        return max([rec(k+1,used)]+[1+rec(k+1,used|1<<j) for j in inside if not used>>j&1 and rows[out[k]]>>j&1])
    return rec(0,0),len(out),len(inside)

def F(m,r):
    return sum((-1)**(r-j)*comb(r,j)*j**m for j in range(r+1))

def coefficients(vals):
    a=[0]
    for m,val in enumerate(vals,1):
        x=val-sum(a[r]*F(m,r) for r in range(1,m))
        need(x%factorial(m)==0,'integral skeleton');a.append(x//factorial(m))
    return a

def expcoef(a,subtract=0):
    c={j:sum((-1)**(r-j)*comb(r,j)*a[r] for r in range(j,len(a))) for j in range(1,len(a))}
    c[1]=c.get(1,0)-subtract
    return {j:x for j,x in c.items() if x}

def cutoff(c):
    if not c:return None
    M=max(c);lower=[j for j in c if j<M]
    if not lower:return 1
    q=max(lower);S=sum(abs(c[j]) for j in lower);m=1
    while abs(c[M])*M**m<=S*q**m:m+=1
    return m

def defect(H,B,sel):return H-sum(B[a][b] for a in sel for b in sel if a!=b)

def row_from_code(n,code):
    rows=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if code>>k&1:rows[i]|=1<<j
        else:rows[j]|=1<<i
    return rows

def check_core(rows,aux=False):
    n=len(rows);tested=0;seeds=0;stable=0;nonstable=0
    for v in range(n):
        nu,do,di=matching(rows,v);RH=nu+1;R=min(nu+2,1+min(do,di));old=[i for i in range(n) if i!=v]
        grids=[]
        for m in range(1,R+2):
            child,labels=inflate(rows,v,m);H,B=table(child);grids.append((H,B))
        ah=coefficients([q[0] for q in grids[:R]])
        need(max(i for i,a in enumerate(ah) if a)==RH,'matching exact HP rank')
        for mask in range(1,1<<(n-1)):
            sel=[i for i in range(n-1) if mask>>i&1]
            values=[defect(H,B,sel) for H,B in grids]
            ad=coefficients(values[:R]);c=expcoef(ad,1)
            need(sum(c.get(j,0)*j**(R+1) for j in c)==values[R]-1,'extra literal prediction')
            if aux:
                for r in range(1,R+2):
                    child,labels=inflate(rows,v,r,True);h,b=table(child)
                    ar=defect(h,b,sel)
                    need(ar==(ad[r] if r<len(ad) else 0)*factorial(r),'auxiliary independent-twin skeleton')
            if values[0]==1:
                seeds+=1
                if not c:stable+=1
                else:
                    nonstable+=1;K=cutoff(c)
                    zeros=[m for m in range(1,K) if sum(x*j**m for j,x in c.items())==0]
                    need(len(zeros)<=len(c)-1,'Chebyshev zero bound')
                    need(any(x!=1 for x in values[:R]),'failure no later than rank')
            tested+=1
    return tested,seeds,stable,nonstable

def main():
    totals=[0]*4;levels=[]
    for n in range(2,5):
        t=[0]*4
        for code in range(1<<(n*(n-1)//2)):
            a=check_core(row_from_code(n,code),aux=n<=3)
            t=[x+y for x,y in zip(t,a)]
        levels.append({'n':n,'labelled_tournaments':1<<(n*(n-1)//2),'marked_port_cases':t[0],'unit_cases':t[1],'stable_unit_ports':t[2],'unstable_unit_ports':t[3]})
    rng=random.Random(288)
    sampled=[]
    for n in (5,6,7):
        for _ in range(6):
            code=rng.randrange(1<<(n*(n-1)//2));res=check_core(row_from_code(n,code))
            sampled.append({'n':n,'code':code,'cases':res})
    # Inherited physical delayed-failure witness, now sharp for modal zero budget.
    rows=[90,28,9,16,36,79,30];v=2;O={3,6};vals=[]
    for m in range(1,7):
        child,labels=inflate(rows,v,m);h,b=table(child);sel=[i for i,a in enumerate(labels) if a in O]
        d=defect(h,b,sel);need(d==1-4*F(m,3),'physical sharp two-length witness');vals.append(d)
    need(matching(rows,v)==(2,2,4),'witness matching rank')
    # General algebraic sharpness: F_R has precisely R-1 positive integer zeros.
    sharp=[]
    for R in range(2,9):
        a=[0]*(R+1);a[R]=1;c=expcoef(a);K=cutoff(c)
        z=[m for m in range(1,K) if not sum(x*j**m for j,x in c.items())]
        need(z==list(range(1,R)),'algebraic sharpness');sharp.append({'rank':R,'integer_roots':z,'effective_cutoff':K})
    congruence=[]
    for p in (5,7,13):
        for k in range(1,7):
            m=1+(p-1)*p**(k-1);mod=p**k
            residue=(-4*(pow(3,m,mod)-3*pow(2,m,mod)+3))%mod
            need(residue==0 and m>2,'local unit agreement without literal success')
            congruence.append({'prime':p,'precision':k,'length':m,'defect_minus_one_mod_pk':residue})
    out={'padic_false_comfort':congruence,'scope':'Complete labelled tournaments orders2..4, every unselected port/nonempty exterior selected set; deterministic sampled higher cores; exact infinite statements proved separately. No all-order seed census.','levels':levels,'samples':sampled,'physical_delayed_failure':{'rows':rows,'selected':sorted(O),'port':v,'matching':2,'delta_lengths1to6':vals,'positive_unit_lengths':[1,2]},'sharp_exponential_controls':sharp}
    tag='.optimized' if not __debug__ else ''
    (HERE/('REPLAY'+tag+'.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'levels':levels,'samples':len(sampled),'delayed':vals},sort_keys=True))
if __name__=='__main__':main()
