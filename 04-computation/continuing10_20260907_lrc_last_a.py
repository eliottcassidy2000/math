#!/usr/bin/env python3
"""Remaining7200 word A: native zero-arm composite closure.

The interval/event implementation is adapted from the audited composite-wedge
engine; no external producer is imported. All universal claims are in the proof.

Standalone exact arithmetic; no mathematical producer is imported or executed.
"""
from pathlib import Path
from fractions import Fraction as F
from math import gcd, lcm
from itertools import combinations, combinations_with_replacement, product
from collections import defaultdict
from hashlib import sha256
import json
import sys
sys.stdout.reconfigure(encoding='utf-8', newline='\n')
HERE=Path(__file__).resolve()
ROOT=HERE.parent.parent if HERE.parent.name=='04-computation' else Path('C:/w/s0905')
DEST=ROOT/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
T=7200
GATES=0

def need(value,why):
    global GATES
    GATES+=1
    if not value:raise ArithmeticError(why)

def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':')).encode()

def atlas_sum(s):
    p=2
    while p*p<=s:
        k=0
        while s%p==0:s//=p;k+=1
        if k and (p%3!=2 or k>2):return False
        p+=1
    return s==1 or s%3==2

def geometry(p,q):
    """Two-pointer intersection of literal strict-danger interval lists."""
    L=14*p*q
    def cells(v,k):
        return [(max(0,(14*j-1)*k),min(L,(14*j+1)*k)) for j in range(v+1)]
    A=cells(p,q);B=cells(q,p);i=j=0;out=[]
    while i<len(A) and j<len(B):
        a,b=A[i];c,d=B[j]
        if max(a,c)<min(b,d):out.append((max(a,c),min(b,d)))
        if b<d:i+=1
        elif d<b:j+=1
        else:i+=1;j+=1
    need(out[0][0]==0 and out[-1][1]==L,'origin is an interior point of the circular intersection')
    intervals=[(out[-1][0]-L,out[0][1])]+out[1:-1]
    need(all(a<b for a,b in intervals),'strict interval lengths positive')
    return L,intervals

def count_intervals(n,L,intervals,phase):
    A,B=phase.numerator,phase.denominator
    den=B*L
    return sum(-((-(B*n*b-A*L))//den)-(B*n*a-A*L)//den-1 for a,b in intervals)

def literal_pair(n,p,q,phase,closed=False):
    A,B=phase.numerator,phase.denominator;den=B*n
    compare=(lambda r:14*min(r,den-r)<=den) if closed else (lambda r:14*min(r,den-r)<den)
    return sum(compare((p*(B*j+A))%den) and compare((q*(B*j+A))%den) for j in range(n))

PROFILE={}
def profile(n,p,q):
    p,q=sorted((p,q));key=(n,p,q)
    if key in PROFILE:return PROFILE[key]
    L,intervals=geometry(p,q);events=defaultdict(lambda:[0,0]);initial=0
    for a,b in intervals:
        initial+=-((-n*b)//L)-(n*a)//L-1
        events[n*a%L][0]+=1;events[n*b%L][1]+=1
    ordered=sorted(events.items())
    lo=hi=initial;arglo=arghi=F(0);rows=[]
    def observe(value,phase):
        nonlocal lo,hi,arglo,arghi
        need(0<=value<=n,'native count remains between zero and clock')
        if value<lo:lo=value;arglo=phase
        if value>hi:hi=value;arghi=phase
    enter0,leave0=events.get(0,(0,0));current=initial+enter0
    for index,(wall,(enter,leave)) in enumerate(ordered):
        nxt=ordered[index+1][0] if index+1<len(ordered) else ordered[0][0]+L
        if wall:
            at=current-leave;after=at+enter
        else:at=initial;after=initial+enter
        observe(at,F(wall,L));observe(after,F(wall+nxt,2*L))
        rows.append([wall,enter,leave,at,after]);current=after
    need(current-leave0==initial,'cyclic event sweep closes at the exact strict endpoint')
    need(lo==literal_pair(n,p,q,arglo),'literal grid at the exact minimum')
    need(hi==literal_pair(n,p,q,arghi),'literal grid at the exact maximum')
    for a in [F(0),F(1,2),F(1,7),arglo,arghi]:
        need(count_intervals(n,L,intervals,a)==literal_pair(n,p,q,a),'interval floor count versus independent native predicates')
    result=dict(n=n,p=p,q=q,L=L,components=len(intervals),walls=len(events),minimum=lo,maximum=hi,
        minimizer=[arglo.numerator,arglo.denominator],maximizer=[arghi.numerator,arghi.denominator],
        intervals_sha256=sha256(canonical(intervals)).hexdigest(),events_sha256=sha256(canonical(rows)).hexdigest())
    PROFILE[key]=result
    return result

def primitive(values):
    D=lcm(*(v.denominator for v in values));U=[int(D*v) for v in values];g=gcd(*U)
    return [v//g for v in U]

def connected(edges):
    seen={0}
    while True:
        more=seen|{b for a,b in edges if a in seen}|{a for a,b in edges if b in seen}
        if more==seen:return len(seen)==7
        seen=more

def main():
    path=ROOT/'05-knowledge/results/continuing8_20260906_lrc_minimum_tree_certificate.json'
    raw=path.read_bytes()
    need(sha256(raw).hexdigest()=='580a7c930103aab3bea867ad463a90b0e0208323a90ee95a685ff811a761582d','inherited full word and minimum table')
    old=next(c for c in json.loads(raw)['clocks'] if c['t']==T)
    word=[1,9,16,18,24,32,60];E=116
    need([word,E,114] in old['survivors'],'word is the actual inherited survivor')
    weights={tuple(d):r[0] for d,r in old['weights']}
    table=[[i,j,weights.get(tuple(sorted((word[i],word[j]))))] for i,j in combinations(range(7),2)]
    need(all(c is not None for i,j,c in table),'full21 possible margin pairs')
    need([(i,j,c) for i,j,c in table if i==0 and c<=E]==[(0,1,114)],'only affordable edge at margin1')
    need(min(c for i,j,c in table if c>0 and (i,j)!=(0,1))==102,'every other positive edge exceeds residual budget')
    cheap=[(i,j) for i,j,c in table if c==0 or (i,j)==(0,1)]
    need(len(cheap)==7,'all residual possibly dangerous edges')
    connected_graphs=[]
    for mask in range(1<<len(cheap)):
        edges=[e for k,e in enumerate(cheap) if mask>>k&1]
        if connected(edges) and not ((2,4) in edges and (3,4) in edges):
            connected_graphs.append(edges)
    expected=[(0,1),(1,3),(2,5),(3,4),(4,5),(4,6)]
    need(connected_graphs==[expected],'unique connected topology after proved16-24-18 exclusion')
    atlas=[(p,q) for q in range(2,356) for p in range(1,min(q,357-q)) if gcd(p,q)==1 and atlas_sum(p+q)]
    need(len(atlas)==5855,'complete strict atlas')
    arms={};low={}
    for a,b,size,low_size in [(32,24,317,14),(18,24,231,5)]:
        e=gcd(a,b);n=T//e;rows=[]
        for p,q in atlas:
            for u,v in [(p,q),(q,p)]:
                if (e*gcd(n,u),e*gcd(n,v))==(a,b):rows.append([u,v,e*profile(n,p,q)['minimum']])
        need(len(rows)==size,'entire directed arm alphabet')
        arms[a]=rows;low[a]=[r for r in rows if r[2]<=2]
        need(len(low[a])==low_size and all(c==0 for u,v,c in low[a]),'complete residual-budget arms')
    products=[]
    for (u,v,ca),(w,z,cb) in product(low[32],low[18]):
        R=primitive([F(u,v),F(1),F(w,z)])
        need(gcd(*R)==1 and len(set(R))==3,'primitive distinct actual triple')
        need([gcd(T,2*r) for r in R]==[32,24,18],'prescribed-clock normalization')
        D=gcd(R[0],R[2]);p,q=sorted((R[0]//D,R[2]//D));e=gcd(T,2*D)
        P=profile(T//e,p,q);credit=e*P['minimum']
        need(e==2 and credit==142 and credit>E,'every native endpoint closes alone')
        products.append(dict(arms=[[u,v,ca],[w,z,cb]],primitive=R,scale_gcd=2,e=e,p=p,q=q,credit=credit))
    need(len(products)==70,'complete14 by5 joint universe')
    # Actual connected, nonunit full13 control inherited from the previous exact tree bank.
    B=[827789,7450101,1105328,796194,131208,123796736,50820]
    A=[T*i for i in range(1,7)]
    need([gcd(T,u) for u in B]==word and gcd(*(A+B))==1 and len(set(A+B))==13,'nonvacuous full13 control')
    edges=[]
    for i,j in combinations(range(7),2):
        d=gcd(B[i],B[j]);p,q=sorted((B[i]//d,B[j]//d))
        if p+q<=356 and atlas_sum(p+q):edges.append((i,j))
    need(connected(edges),'actual strict graph connected on positive control')
    den=7*T
    safe=[j for j in range(T) if all(14*min(v*(7*j+1)%den,den-v*(7*j+1)%den)>=den for v in A+B)]
    need(len(safe)>0,'literal full13 safety')
    need(literal_pair(7,1,1,F(1,2))==0 and literal_pair(7,1,1,F(1,2),True)==2,'strict endpoint hostile')
    cert=dict(status='FINITE-EXACT; analytic scoped theorem in companion report',clock=T,word=word,E=E,
        inherited_pin=sha256(raw).hexdigest(),pair_table=table,cheap_edges=cheap,connected_graphs=connected_graphs,
        arms=[dict(margins=[a,24],rows=arms[a],low=low[a]) for a in [32,18]],products=products,
        profiles=[PROFILE[k] for k in sorted(PROFILE)],
        full13_control=dict(body=A,tails=B,edges=edges,phase=[1,7],safe_count=len(safe),safe_hash=sha256(canonical(safe)).hexdigest()))
    (DEST/'continuing10_20260907_lrc_last_a_certificate.json').write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Remaining A word',word,'E116: CLOSED relative to actual connected strict complement')
    print('All21 inherited minima;128 subgraphs; one residual tree',expected)
    print('Complete directed alphabets317/231; budget<=2 arms14/5; all70 native products give endpoint142')
    print('Profiles',len(PROFILE),'actual full13 safe count',len(safe),'strict endpoint hostile0 versus false2')
    print('Always-active exact gates:',GATES)

if __name__=='__main__':main()
