#!/usr/bin/env python3
"""Exact finite certificates WITH analytic infinite-tail cutoff, stdlib only."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
import json

def need(x,s):
    if not x: raise RuntimeError(s)

@lru_cache(None)
def w(a):
    if a<2:return Q(a)
    h=a//2
    return Q(5,2)*w(h) if a%2==0 else Q(3,2)*w(h)+w(h+1)

ROWS={1:(1,(0,)),2:(2,(0,2)),3:(4,(0,3,7,11)),4:(4,(0,4,10,16)),
      5:(8,(0,5,12,19,29,35,43,51)),6:(8,(0,6,15,23,35,43,55,65)),
      7:(8,(0,7,17,27,41,51,65,79)),8:(8,(0,8,20,32,50,62,80,98))}

def row(a,b):
    l,h=ROWS[a]
    v=w(a)*l*(b//l)+h[b%l]
    need(v.denominator==1,'row integral')
    return int(v)

@lru_cache(None)
def pats(a,maximal=True):
    ps=[]
    for p in product(range(1,a),repeat=5):
        if any(p[i]+p[(i+1)%5]>a for i in range(5)):continue
        if maximal and any(p[i]+p[(i+1)%5]<a and p[i]+p[(i-1)%5]<a for i in range(5)):continue
        ps.append(p)
    return tuple(ps)

def root_max(a,b,maximal=True):
    # In constraint-cycle order 0,2,4,1,3. Fix first port; use exact
    # max-plus path DP with closing edge. No relaxation of positive ports.
    best=(-1,None)
    if b<2:return best
    for p in pats(a,maximal):
        aa=tuple(p[i] for i in (0,2,4,1,3))
        for x in range(1,b):
            dp={x:(row(aa[0],x),(x,))}
            for j in range(1,5):
                prefix=[(-1,())]*b
                z=(-1,())
                for v in range(1,b):
                    if v in dp and dp[v][0]>z[0]:z=dp[v]
                    prefix[v]=z
                new={}
                for v in range(1,b):
                    z=prefix[b-v]
                    if z[0]>=0:new[v]=(z[0]+row(aa[j],v),z[1]+(v,))
                dp=new
            for v,z in dp.items():
                if v+x<=b and z[0]>best[0]:best=(z[0],{'kind':'c5','a':aa,'b':z[1]})
    return best

def ceil(x):return -((-x.numerator)//x.denominator)

def certificate(a):
    L,hs=ROWS[a];s=w(a)
    D=max(s*r-hs[r] for r in range(L))
    ps=pats(a);sharp=[];other=[];cut=1
    for p in ps:
        cs=[w(p[i]) for i in (0,2,4,1,3)]
        half=sum(cs)/2
        vertex=max([Q(0),half]+cs+[cs[i]+cs[(i+2)%5] for i in range(5)])
        need(vertex<=s,'inherited slope violation')
        if half==s:
            ys=[sum((-1)**j*cs[(i-j)%5] for j in range(5))/2 for i in range(5)]
            need(all(ys[i]+ys[(i-1)%5]==cs[i] for i in range(5)),'dual identity')
            need(min(ys)>0,'strict dual')
            need(all((L//2)%ROWS[x][0]==0 for x in p),'child translation period')
            S=ceil(D/min(ys))-1
            K=S+L+1
            cut=max(cut,K)
            sharp.append({'pattern':p,'dual':[str(y) for y in ys],'slack_bound':S,'cutoff':K})
        else:
            gap=s-vertex;need(gap>0,'nonsharp gap')
            K=ceil(D/gap)-1
            cut=max(cut,K)
            other.append({'pattern':p,'gap':str(gap),'cutoff':K})
    joins=[]
    for x in range(1,a):
        gap=s-w(x)-w(a-x);need(gap>0,'join gap')
        K=ceil(D/gap)-1;cut=max(cut,K)
        joins.append({'split':x,'gap':str(gap),'cutoff':K})
    finite=[];seeds=[]
    for b in range(1,cut+1):
        cr,plan=root_max(a,b)
        j=max((row(x,b)+row(a-x,b) for x in range(1,a)),default=1)
        need(cr<=row(a,b),'C5 upper certificate failed')
        need(j<=row(a,b),'join upper certificate failed')
        finite.append([b,cr,j,row(a,b)])
        if b<=L:
            opts=[(1,{'kind':'singleton'}) ,(cr,plan)]
            opts += [(row(x,b)+row(a-x,b),{'kind':'join','a':[x,a-x],'b':[b,b]}) for x in range(1,a)]
            opts += [(row(a,y)+row(a,b-y),{'kind':'union','a':[a,a],'b':[y,b-y]}) for y in range(1,b)]
            val,pl=max(opts,key=lambda z:z[0]);need(val==row(a,b),'seed missing')
            seeds.append({'a':a,'b':b,'n':val,'construction':pl})
    need(cut>=L,'all seeds checked')
    for r in range(L):
        for t in range(L):need(hs[r]+hs[t]<=s*L*((r+t)//L)+hs[(r+t)%L],'union residue certificate')
    return {'a':a,'slope':str(s),'period':L,'residue_values':hs,'max_deficit':str(D),
        'maximal_patterns':len(ps),'sharp_patterns':sharp,'other_patterns':other,
        'joins':joins,'cutoff':cut,'finite_upper':finite,'seeds':seeds,'union_checks':L*L}

def main():
    certs=[certificate(a) for a in range(5,9)]
    full=[]
    for a in range(2,7):
        for b in range(2,7):
            x=root_max(a,b,False)[0];y=root_max(a,b,True)[0]
            need(x==y,'full vs maximal root discrepancy');full.append([a,b,x])
    strict=0
    for a in range(1,129):
        for b in range(1,129):need(w(a)+w(b)<w(a+b),'strict superadditivity');strict+=1
    # Reconstructed seed stats: children are known row attainers inductively;
    # root clique/independence budget formulas are checked directly.
    seedcount=0
    for c in certs:
        for se in c['seeds']:
            a,b,n=se['a'],se['b'],se['n'];p=se['construction']
            aa,bb=p['a'],p['b'];need(sum(row(x,y) for x,y in zip(aa,bb))==n,'seed order')
            if p['kind']=='c5':
                om=max(aa[i]+aa[(i+2)%5] for i in range(5))
                al=max(bb[i]+bb[(i+1)%5] for i in range(5))
            elif p['kind']=='join':om=sum(aa);al=max(bb)
            else:om=max(aa);al=sum(bb)
            need(om<=a and al<=b,'seed feasibility');seedcount+=1
    need(row(6,1)==6 and w(6).denominator==1 and row(6,1)<w(6),'denominator hostile')
    print(json.dumps({'status':'PASS','scope':'exact all-positive-B rows5..8 via finite analytic-cutoff certificates',
        'certificates':certs,'full_pattern_controls':full,'strict_pairs':strict,
        'seed_constructions_checked':seedcount},sort_keys=True,indent=2))
if __name__=='__main__':main()
