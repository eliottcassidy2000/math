#!/usr/bin/env python3
"""Standard-library exact controls. Run python3 -I -B verify.py; also under -O."""
import importlib.util
from pathlib import Path
import json
from fractions import Fraction as F
from collections import Counter
spec=importlib.util.spec_from_file_location('inherited',Path(__file__).with_name('inherited_verify.py'))
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
def check(ok,msg):
    if not ok: raise RuntimeError(msg)
def rank(a):
    a=[[F(x) for x in row] for row in a];r=0
    for col in range(len(a[0])):
        k=next((k for k in range(r,len(a)) if a[k][col]),None)
        if k is None: continue
        a[r],a[k]=a[k],a[r];z=a[r][col];a[r]=[x/z for x in a[r]]
        for k in range(len(a)):
            if k!=r:
                z=a[k][col];a[k]=[x-z*y for x,y in zip(a[k],a[r])]
        r+=1
    return r

def inverse(d,x):
    H,Q,s,e,a,de,T,q,b=[d[k] for k in ('H','Q','s','e','alpha','delta','T','q','b')]
    Dp,Ep,Up,dp=map(F,x);Delta=q*Q-H*b
    t=(q*Dp-H*dp)/Delta
    w=((Q-H)*dp-(b-q)*Dp)/Delta
    p=(Ep-(de-e)*t)/e
    z=w-p
    E=(Up-a*p-(2*T-a)*t-(2*a-s)*z)/s
    return [p+E,E,t+z,z]

def main():
    basis=[[0],[2,0],[2,4,1],[6,4,0]]
    X=[v.profile(a,0) for a in basis]
    check(v.determinant(X)==-4,'physical basis')
    C=[[1,0,0,0],[1,-1,1,0],[1,0,2,0],[1,-1,3,-2]]
    check(v.determinant(C)==4,'observer determinant')
    census={}; inverses=literal=perm=0;single=[]
    for n in range(1,6):
        cnt=Counter()
        for b in v.tournaments(n):
            d=v.stats(b,0);K=v.operator(b,0);r=rank(K)
            check(d==v.permutation_stats(b,0),'independent stats');perm+=1
            Delta=d['q']*d['Q']-d['H']*d['b']
            strong=v.strong(b)
            check(bool(d['s'] and d['e'])==strong,'endpoint strong equivalence')
            if n==1:
                check(r==3 and v.determinant(K)==0,'singleton boundary')
                continue
            det=d['s']*d['e']*Delta
            check(det==v.determinant(K),'determinant identity')
            check(not det or strong,'full rank requires strong')
            if strong:check(r==3+bool(Delta),'strong rank dichotomy')
            cnt[(strong,r,(Delta>0)-(Delta<0))]+=1
            if n<=4:single.append(K[0])
            for a,x in zip(basis,X):
                g,root=v.graft(a,0,b,0);xp=v.profile(g,root)
                check(xp==v.mv(K,x),'literal closure');literal+=1
                check(d['q']*xp[0]-d['H']*xp[3]==Delta*(x[2]-x[3]),'image side identity')
                check((x[2]==x[3])==(a[0]==0),'sink boundary')
                if det:
                    check(inverse(d,xp)==x,'physical inverse');inverses+=1
        census[str(n)]={'strong=%s,rank=%s,sign=%s'%k:z for k,z in sorted(cnt.items())}
    check(rank(single)==2,'single graft response rank')
    check(rank([C[0]]+single)==3,'one layer plus identity rank')
    contexts=0
    for n in range(1,6):
        for a in v.tournaments(n):
            x=v.profile(a,0);j,r=v.graft(a,0,[0]);jj,rr=v.graft(j,r,[0]);l,lr=v.graft(a,0,[2,0]);
            h=[x[0],v.profile(j,r)[0],v.profile(jj,rr)[0],v.profile(l,lr)[0]]
            check(h==v.mv(C,x),'literal context response')
            U=F(h[2]-h[0],2);E=h[0]+U-h[1];z=F(h[1]+2*U-h[3],2)
            check([h[0],E,U,z]==x,'scalar observation recovery');contexts+=1
    A=[14,28,24,16,1];B=[10,28,25,16,1]
    out=[]
    for a in (A,B):
        vals=[]
        for root in range(5):
            g,r=v.graft(a,root,v.B4);vals.append(v.profile(g,r)[0])
        out.append(vals)
    check(set(out[0]).isdisjoint(out[1]),'H9 all roots separation')
    check(out==[[423,605,357,185,261],[315,541,551,217,207]],'H9 values')
    examples={}
    for name,b in [('negative',[8,1,3,6]),('positive',[2,8,3,5]),('singular5',[18,8,3,5,14]),('triangle',[2,4,1])]:
        d=v.stats(b,0);examples[name]={'rows':b,'stats':d,'det':v.determinant(v.operator(b,0)),'rank':rank(v.operator(b,0))}
    print(json.dumps({'status':'FINITE-EXACT controls for proved statements','universe':'all labeled tournaments of orders 1..5, root 0; no filters','stats_permutation_checks':perm,'literal_graft_checks':literal,'inverse_checks':inverses,'literal_context_checks':contexts,'census':census,'examples':examples,'H9_graft_outputs':out,'observer_determinant':4},indent=2,sort_keys=True))
if __name__=='__main__':main()
