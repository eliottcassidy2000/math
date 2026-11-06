"""Exact positive g-readout inversion by digit splitting and a fixed gcd.

The search never treats a congruence solution as a target witness: g is
evaluated exactly and the original index is reconstructed and checked.
"""
from functools import lru_cache
from math import gcd,isqrt
from pathlib import Path
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).parent/'dependencies'))
import core

def ceildiv(a,b):return -((-a)//b)

def readout(row,v):
    alpha,beta,gamma,delta=row
    return alpha*core.g(v)+beta*v+gamma*(v%2)+delta

def valid(row):
    a,b,c,d=row
    return a>=1 and a%2==1 and b>=a and -a<=c<=0 and a-b<=d<=0

def fixed_inverse(row,n,vmin=2,vmax=None,first=False):
    """Returns all indices in the requested interval, or stops at first=True.

    The original-index lower bound vmin must be at least two.
    Row order is (alpha,beta,gamma,delta), deliberately different from the
    core matrix row (alpha,gamma,beta,delta).
    """
    core.need(vmin>=2,'original index lower bound must be at least two')
    core.need(valid(row),'positive physical row domain')
    alpha,beta0,gamma0,delta0=row;d=gcd(alpha,beta0);m=alpha//d
    if n<alpha:return [],{'nodes':0,'leaves':0,'g_evaluations':0,'node_budget':1,'gcd':d}
    bound=1+(n-alpha)//beta0
    vmax=bound if vmax is None else min(bound,vmax)
    if vmin>vmax:return [],{'nodes':0,'leaves':0,'g_evaluations':0,'node_budget':1,'gcd':d}
    K=0
    while beta0*beta0*(4**K)<n*d:K+=1
    budget=2**(K+1)-1
    stats={'nodes':0,'leaves':0,'g_evaluations':0,'node_budget':budget,'gcd':d}
    stack=[(beta0,gamma0,delta0,1,0,0)];found=[]
    inverse_by_depth=[]
    if m>1:
        inverse_by_depth=[pow(2*beta0//d,-1,m)]
        inv2=pow(2,-1,m)
        for k in range(K):inverse_by_depth.append(inverse_by_depth[-1]*inv2%m)
    while stack:
        beta,gamma,delta,scale,offset,depth=stack.pop();stats['nodes']+=1
        core.need(valid((alpha,beta,gamma,delta)) and gcd(alpha,beta)==d,'preserved row and gcd')
        lo=max(1,ceildiv(n,alpha+beta),ceildiv(vmin+offset,scale))
        hi=min(1+(n-alpha)//beta,(vmax+offset)//scale)
        if lo>hi:continue
        allowed=[p for p in (0,1) if (n-delta-(beta+gamma)*p)%d==0 and ceildiv(lo-p,2)<=(hi-p)//2]
        if not allowed:continue
        if beta*beta>=n*d:
            stats['leaves']+=1;this_leaf=0
            for p in allowed:
                rhs=n-delta-(beta+gamma)*p
                if rhs%d:continue
                tlo=ceildiv(lo-p,2);thi=(hi-p)//2
                residue=0 if m==1 else (rhs//d)*inverse_by_depth[depth]%m
                t=tlo+(residue-tlo)%m
                core.need(t+m>thi,'at most one arithmetic candidate per parity')
                if t<=thi:
                    h=2*t+p;this_leaf+=1;stats['g_evaluations']+=1
                    if alpha*core.g(h)+beta*h+gamma*p+delta==n:
                        v=scale*h-offset
                        core.need(vmin<=v<=vmax and readout(row,v)==n,'exact original-index witness')
                        found.append(v)
                        if first:
                            core.need(stats['nodes']<=budget,'proved search-node budget')
                            return found,stats
            core.need(this_leaf<=2,'leaf arithmetic candidate bound')
        else:
            # h=2j-1 and h=2j. A singleton intermediate j=1 is arithmetic
            # only; every returned v>=2 is rebuilt by its canonical unit word.
            if 1 in allowed:stack.append((alpha+2*beta,-alpha,delta+gamma-beta,2*scale,scale+offset,depth+1))
            if 0 in allowed:stack.append((2*beta,0,delta,2*scale,offset,depth+1))
    core.need(stats['nodes']<=budget and len(found)==len(set(found)),'node budget and unique indices')
    return sorted(found),stats

def target_inverse(n,limit=5,first=False):
    """Complete transpose-closed central-word inverse, not a cutoff guess."""
    found=[];totals={'rows':0,'nodes':0,'g_evaluations':0}
    forms=core.form_data(limit)
    for u in range(2,isqrt(2*n)+1):
        for w,A,b in forms:
            if b*u*u>4*n:continue
            ar=core.row_for(A,u);row=(ar[0],ar[2],ar[1],ar[3])
            v,stats=fixed_inverse(row,n,u,4*n//(b*u),first)
            totals['rows']+=1
            for key in ('nodes','g_evaluations'):totals[key]+=stats[key]
            for vv in v:found.append((w,u,vv))
            if found and first:return found,totals
    return sorted(found,key=lambda z:(len(z[0]),z[0],z[1],z[2])),totals

def checks():
    controls=0;nodes=0;worst_ratio=0
    # Direct tests include valid abstract rows, not just rows already known
    # to occur as graph profiles. The numerical theorem permits that cone.
    rows=[]
    for alpha in (1,3,5,9):
        for beta in range(alpha,alpha+7):
            for gamma in (-alpha,0):
                for delta in (alpha-beta,0):rows.append((alpha,beta,gamma,delta))
    rows=sorted(set(rows))
    for row in rows:
        N=201
        expected={}
        for v in range(2,1+(N-row[0])//row[1]+1):
            n=readout(row,v)
            if n<=N:expected.setdefault(n,[]).append(v)
        for n in range(1,N+1):
            got,s=fixed_inverse(row,n)
            core.need(got==expected.get(n,[]),'complete direct fixed-row comparison')
            controls+=1;nodes+=s['nodes'];worst_ratio=max(worst_ratio,s['nodes']/s['node_budget'])
    target_checks=[]
    for n in (7,13,31,43,67,163,277,313,1877,2647,5001,17383):
        got,s=target_inverse(n)
        direct=core.inverse(n,5)
        core.need(got==direct,'complete direct/progression/digital target equality')
        target_checks.append({'target':n,'representations':len(got),**s})
    return {'status':'PASS','fixed_rows':len(rows),'fixed_row_target_controls':controls,
            'total_fixed_nodes':nodes,'maximum_node_budget_fraction':worst_ratio,
            'target_controls':target_checks,'universal_occupancy':'OPEN'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--target',type=int);p.add_argument('--all',action='store_true');a=p.parse_args()
    if a.target:
        hit,stats=target_inverse(a.target,first=not a.all)
        out={'target':a.target,'complete_search':a.all or not hit,'witnesses':hit,'statistics':stats}
    else:out=checks()
    print(json.dumps(out,indent=2,sort_keys=True))
