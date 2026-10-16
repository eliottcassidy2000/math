#!/usr/bin/env python3
"""Exact odd-ray payment with original annular rows and native integer clocks.
No ROOT label is used in selecting or evaluating a controller. Run from repo root.
"""
import hashlib
import json
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from collatz_gram_slack_20261015b import close_graph, need, T

ROOT = Path(__file__).resolve().parents[2]

def w(n, residue=2):
    return 2 if n % 3 == residue else -1

def strip(n):
    a = (n & -n).bit_length() - 1
    return n >> a, a

def primitive(n):
    return n > 0 and n % 8 in (1, 3, 7)

def ray_head(n):
    j = 0
    while n % 8 == 5:
        n = (n - 1)//4
        j += 1
    return n, j

def shortcut_receipt(n, r):
    e, d = strip(n)
    for _ in range(r):
        if e == 1:
            break
        e, a = strip(3*e+1)
        d += a
    return e, (-1 if d & 1 else 1), d

def coefficient(N, b):
    """Ray telescope at primitive source b, before the last sign to U(b)."""
    if b > 4*N:
        return 0
    u, j = b, 0
    while 4*u+1 <= 4*N:
        u = 4*u+1
        j += 1
    return w(b) if u > 2*N else w(b)-w(b+j+1)

def count_residue(lo, hi, r):
    if hi <= max(lo, 0):
        return 0
    lo = max(lo, 0)
    return sum((hi-a)//24 - (lo-a)//24
               for a in range(1, 25) if primitive(a) and a % 3 == r)

def threshold(Y, p):
    return (Y-(p-1)//3)//p

def norm1_exact(N):
    """O(log N) floor count; no orbit or source enumeration."""
    total, p, j = 0, 1, 0
    while threshold(4*N, p) >= 1:
        lo, mid, hi = threshold(4*N,4*p), threshold(2*N,p), threshold(4*N,p)
        need(0 <= lo <= mid <= hi, 'nonnegative exact ray bands')
        for r in range(3):
            total += abs(w(r)-w(r+j+1))*count_residue(lo,mid,r)
            total += abs(w(r))*count_residue(mid,hi,r)
        p *= 4
        j += 1
    return total

def norm0_exact(N):
    return 4*N - 3*((2*N+3)//6)+(4*N+1)//6-(2*N+1)//6

def error_bound(N):
    return 32*(1+((4*N).bit_length()-1)//2)+1

def rows(k, cache, residue=2):
    X, b = 1 << k, 1 << (k//2)
    n = np.arange(1, 8*X+1, dtype=np.int64)
    cols = sorted({e for e,s,d in cache.values()})
    index = {e:i for i,e in enumerate(cols)}
    targets = np.array([index[cache[int(i)][0]] for i in n], dtype=np.int64)
    signed = np.array([w(int(i),residue)*cache[int(i)][1] for i in n], dtype=np.int64)
    out = []
    for a in range(X,2*X,b):
        c = np.maximum(0,np.minimum(a+b-1,n-1)-np.maximum(a,(n+3)//4)+1)
        row = np.zeros(len(cols),dtype=np.int64)
        np.add.at(row,targets,c*signed)
        out.append(row)
    return np.array(out),cols

def measurement(k, matrix, cols, exact_gram):
    X, b = 1 << k, 1 << (k//2)
    R = sum(int(s)**2 for s in np.sum(np.abs(matrix),axis=1))
    mixed = int(np.sum((np.max(matrix,axis=0)>0)&(np.min(matrix,axis=0)<0)))
    rec = dict(k=k,columns=len(cols),mixed_columns=mixed,row_envelope=[R,b],
               row_envelope_over_X3=R/(b*X**3),row_envelope_over_X2=R/(b*X**2))
    if exact_gram:
        top = int(np.max(np.abs(matrix)))
        need(len(matrix)*top*top < 2**63, 'int64 entry product bound')
        G = matrix.T @ matrix
        value = sum(int(abs(v)) for v in G.flat)
        need(value <= R, 'row envelope dominates exact Gram')
        # Independent Python-integer products at deterministic entries.
        for i in sorted({0,len(cols)//3,len(cols)//2,len(cols)-1}):
            for j in sorted({0,len(cols)//4,len(cols)//2,len(cols)-1}):
                need(int(G[i,j])==sum(int(row[i])*int(row[j]) for row in matrix), 'independent dot product')
        rec.update(gram=[value,b],gram_over_X3=value/(b*X**3),gram_over_X2=value/(b*X**2),
                   physical_clock_cancellation=R-value,
                   gram_sha256=hashlib.sha256(G.tobytes()).hexdigest())
        if mixed == 0:
            need(value == R, 'global sign cone makes row envelope exact')
    return rec

def exact_controls():
    # A direct growing annular row, independent of the floor/ray formulas.
    row=defaultdict(int)
    def add(n):
        e,s,d=shortcut_receipt(n,1)
        row[e] += w(n)*s
    for n in (2,3,4):
        add(n)
    lo,hi=None,None
    half_failures=[]
    for N in range(1,4097):
        literal=sum(abs(v) for v in row.values())
        need(literal==norm1_exact(N), 'literal source row equals logarithmic ray count')
        need(abs(F(literal)-F(16*N,7)) <= error_bound(N), 'uniform norm error')
        if 2*literal*literal > norm0_exact(N)**2:
            half_failures.append(N)
        err=7*literal-16*N
        if lo is None or err<lo[0]:lo=[err,N]
        if hi is None or err>hi[0]:hi=[err,N]
        if N<=256:
            direct=defaultdict(int)
            for n in range(N+1,4*N+1):
                e,s,d=shortcut_receipt(n,1)
                x=n
                for _ in range(d):
                    need(x!=1, 'strict first-hit native receipt')
                    x=T(x)
                need(x==e, 'integer-clock endpoint')
                direct[e]+=s*w(n)
            model=defaultdict(int)
            for b in range(1,4*N+1,2):
                if primitive(b):
                    v,a=strip(3*b+1)
                    model[v]+=((-1)**a)*coefficient(N,b)
            need({e:v for e,v in direct.items() if v}=={e:v for e,v in model.items() if v},'independent ray telescope')
        for j in (1,2,3):add(4*N+j)
    need(half_failures==[1,7,11,13,23,25,31,35,36,37,43,49], 'complete finite half-payment boundary')
    # At X=4096 the pointwise half-payment inequality is exact rational.
    X=1<<12; E=32*(1+(12+3)//2)+1
    need(2*(F(16*X,7)+E)**2 <= (F(10*X,3)-F(2,3))**2,'half payment base')
    for k in range(12,257):
        X=1<<k; E=32*(1+(k+3)//2)+1
        need(2*(F(16*X,7)+E)**2 <= (F(10*X,3)-F(2,3))**2,'dyadic half payment control')
    big=[]
    for N in (10**30+7,2**127-1,2**256+19):
        L=norm1_exact(N)
        need(abs(F(L)-F(16*N,7))<=error_bound(N),'astronomical floor control')
        big.append(dict(N=str(N),L=str(L),scaled_error=str(7*L-16*N)))
    # Honest two-odd-step clock/phase witness.
    witnesses=[]
    for N in (4,5):
        terms=[]
        for n in range(N+1,4*N+1):
            e,s,d=shortcut_receipt(n,2)
            if e==11:terms.append([n,w(n)*s,d])
        witnesses.append(dict(N=N,terms=terms,coefficient=sum(t[1] for t in terms)))
    need([x['coefficient'] for x in witnesses]==[1,-1],'second-step sign-cone failure')
    # Exact high-terminal obstruction, k=4..9, 0<=r<=k-2.
    high_checks=0
    for k in range(4,10):
        X=1<<k
        for r in range(k-1):
            inverse=defaultdict(list)
            for n in range(1,8*X+1):inverse[shortcut_receipt(n,r)[0]].append(n)
            selected=list(range(4*X+(1<<(r+1))-1,6*X+1,1<<(r+1)))
            need(len(selected)==X//(1<<r),'high branch source count')
            for n in selected:
                e,s,d=shortcut_receipt(n,r)
                need(d==r and inverse[e]==[n], 'isolated all-valuation-one endpoint')
                high_checks+=1
    # Mean-zero ternary phase alone has an unbounded cyclic-potential cost.
    ternary=[]
    for q in range(1,8):
        m=3**q;x=0;orbit=[]
        for _ in range(m):orbit.append(x);x=(4*x+1)%m
        need(x==0 and len(set(orbit))==m,'full S cycle modulo 3^q')
        half=(m-1)//2
        h=[1]*half+[0]+[-1]*half
        p=[0]
        for value in h:p.append(p[-1]+value)
        need(p[-1]==0 and max(p)-min(p)==half,'potential range')
        ternary.append(dict(q=q,period=m,min_sup_potential=[m-1,4]))
    energy_seeds={35:93862178945193578345,36:172802101944936995215}
    for k,E in energy_seeds.items():
        need(E <= 1800*(1<<k)**2, 'inherited exact energy seed fits new budget')
    need((1<<35)>36*1800*3, 'new constant budget reserves the integer seam')
    return dict(literal_N_max=4096,clock_and_ray_N_max=256,error_extrema_times7=[lo,hi],
                half_payment_failure_N=half_failures,energy_seeds=energy_seeds,
                half_payment_k_controls=[12,256],large_integer_controls=big,
                second_odd_witness=witnesses,high_terminal_checks=high_checks,ternary=ternary)

def census():
    out=[]
    for k in (8,10,12,14):
        X=1<<k
        cache={n:shortcut_receipt(n,0) for n in range(1,8*X+1)}
        for r in range(17):
            if r in (0,1,2,4,8,16):
                mat,cols=rows(k,cache)
                rec=measurement(k,mat,cols,k==8 or (k==10 and r in (1,2)))
                rec['odd_blocks']=r
                if r in (0,1):
                    need(rec['mixed_columns']==0,'one-block global sign coherence')
                    norm=norm0_exact if r==0 else norm1_exact
                    b=1<<(k//2)
                    need(rec['row_envelope'][0]==sum(sum(norm(n) for n in range(a,a+b))**2 for a in range(X,2*X,b)),'floor counts recover every coarse row')
                out.append(rec)
            ctl={1:(1,1,0)}
            for e,s,d in cache.values():
                if e not in ctl:
                    v,a=strip(3*e+1);ctl[e]=(v,(-1)**a,a)
            cache={n:(ctl[e][0],s*ctl[e][1],d+ctl[e][2]) for n,(e,s,d) in cache.items()}
    a=next(x for x in out if x['k']==8 and x['odd_blocks']==1)
    b=next(x for x in out if x['k']==8 and x['odd_blocks']==2)
    need(2*b['gram'][0]>a['gram'][0],'half payment does not automatically iterate')
    # Orthogonal mod-3 sector: halving leaves one odd shell and U is injective.
    sector=[]
    for k in (4,6,8):
        recs=[]
        for r in (0,1):
            cache={n:shortcut_receipt(n,r) for n in range(1,8*(1<<k)+1)}
            mat,cols=rows(k,cache,0);recs.append(measurement(k,mat,cols,True))
        need(recs[0]['gram']==recs[1]['gram'],'zero payment in residue-zero sector')
        sector.append(dict(k=k,gram=recs[0]['gram']))
    return out,sector

def glued_census():
    out=[]
    for k in (8,10,12,14,16):
        X=1<<k; edges={}; boundary={}
        for n in range(1,8*X+1):
            x=n
            while x%2==0:edges[x]=x//2;x//=2
            boundary[n]=x
        for r in range(1,9):
            ctl={1:1}
            for x in set(boundary.values()):
                if x==1:continue
                y=T(x);edges[x]=y
                while y%2==0:edges[y]=y//2;y//=2
                ctl[x]=y
            boundary={n:ctl[e] for n,e in boundary.items()}
            if r not in (1,2,4,8):continue
            cache,cycles=close_graph(edges,range(1,8*X+1))
            need(not cycles,'no nonROOT cycle in declared control universe')
            mat,cols=rows(k,cache)
            rec=measurement(k,mat,cols,len(cols)<=3500)
            rec.update(odd_blocks=r,observed_edges=len(edges),
                       unresolved_sources=sum(cache[n][0]!=1 for n in range(1,8*X+1)),
                       max_native_clock=max(cache[n][2] for n in range(1,8*X+1)))
            out.append(rec)
    return out

def main():
    controls=exact_controls()
    census_rows,sector=census()
    result=dict(status='PROVED one-block payment; FINITE-EXACT controls; iteration and uniform energy OPEN',
                controls=controls,fixed_depth=census_rows,orthogonal_sector=sector,glued=glued_census())
    path=ROOT/'05-knowledge/results/collatz_odd_ray_overlap_20261015d.json'
    data=json.dumps(result,indent=2,sort_keys=True)+'\n';path.write_text(data)
    print(json.dumps(dict(output=str(path),sha256=hashlib.sha256(data.encode()).hexdigest(),
                         fixed_cases=len(census_rows),glued_cases=len(result['glued'])),sort_keys=True))

if __name__=='__main__':main()
