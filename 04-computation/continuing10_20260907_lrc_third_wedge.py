#!/usr/bin/env python3
"""Third native wedge with joint valuation exclusions and actual topology.

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

def max_tree(costs):
    par=list(range(7));edges=[]
    def find(i):
        while par[i]!=i:i=par[i]
        return i
    for cost,i,j in sorted(costs,reverse=True):
        a,b=find(i),find(j)
        if a!=b:par[a]=b;edges.append([cost,i,j])
    return sum(x[0] for x in edges),edges

def primitive(values):
    D=lcm(*(v.denominator for v in values));U=[int(D*v) for v in values];g=gcd(*U)
    return [v//g for v in U]

def projective(edges):
    adj=[[] for _ in range(7)]
    for i,j,r in edges:adj[i].append((j,r));adj[j].append((i,1/r))
    V=[None]*7;V[0]=F(1);queue=[0]
    for i in queue:
        for j,r in adj[i]:
            if V[j] is None:V[j]=V[i]*r;queue.append(j)
            else:need(V[j]==V[i]*r,'tree path consistency')
    need(len(queue)==7,'inherited selected tree connected')
    return primitive(V)

def main():
    wr=(DEST/'continuing10_20260907_lrc_composite_wedges_certificate.json').read_bytes()
    need(sha256(wr).hexdigest()=='890d00d44f0195765d62fe1d40b59ad102311f34a84a42a5fddc229d037209e9','frozen complete first-two wedges and forced24,18,30 profiles')
    P=json.loads(wr);condition=next(r for r in P['conditional_words'] if r['forced']==[24,18,30])
    need(condition['attempted']==23751 and len(condition['words'])==162 and condition['maximum']==107,'complete inherited conditional ceiling')
    raw=(ROOT/'05-knowledge/results/continuing8_20260906_lrc_minimum_tree_certificate.json').read_bytes()
    need(sha256(raw).hexdigest()=='580a7c930103aab3bea867ad463a90b0e0208323a90ee95a685ff811a761582d','complete inherited minimum-tree data')
    C=next(r for r in json.loads(raw)['clocks'] if r['t']==7200);W={tuple(d):v for d,v in C['weights']}
    atlas=[(p,q) for q in range(2,356) for p in range(1,min(q,357-q)) if gcd(p,q)==1 and atlas_sum(p+q)]
    need(len(atlas)==5855,'complete strict atlas')
    arm_bank={}
    for a in [24,30]:
        e=gcd(a,18);n=T//e;arms=[]
        for p,q in atlas:
            for u,v in [(p,q),(q,p)]:
                if (e*gcd(n,u),e*gcd(n,v))!=(a,18):continue
                R=profile(n,p,q);arms.append([u,v,e*R['minimum']])
        need(len(arms)==(231 if a==24 else 141),'complete directed arm alphabet')
        arm_bank[a]=arms
    first=0;rejected=[];accepted=[]
    for (u,v,cu),(w,z,cw) in product(arm_bank[24],arm_bank[30]):
        if cu+cw>107:first+=1;continue
        R=primitive([F(u,v),F(1),F(w,z)]);actual=[gcd(T,6*r) for r in R]
        record=dict(arms=[[u,v,cu],[w,z,cw]],primitive=R,scale_gcd=6,actual_margins=actual)
        # Any actual primitive triple has gcd(h,t)=gcd(24,18,30)=6.
        if len(set(R))<3 or actual!=[24,18,30]:
            need(len(set(R))==3,'these failures are true clock inconsistencies rather than collisions')
            need([gcd(800,6*r) for r in R]==[gcd(800,d) for d in [24,18,30]],'all failed margins differ only at prime3')
            rejected.append(record);continue
        D=gcd(R[0],R[2]);p,q=sorted((R[0]//D,R[2]//D));e=gcd(T,6*D)
        Q=profile(T//e,p,q);cost=e*Q['minimum'];forest=sum(sorted([cu,cw,cost])[-2:])
        need(forest>107,'same actual triple forest pays complete107 ceiling')
        record.update(pair=[e,p,q,cost],forest=forest);accepted.append(record)
    need((first,len(rejected),len(accepted))==(32541,26,4),'complete32571 locally compatible product partition')
    need(min(r['forest'] for r in accepted)==162,'minimum surviving native forest credit')
    need(all(r['arms'][0][2]==0 and r['arms'][1][2]==66 for r in rejected+accepted),'complete low-credit branch has exact0,66 arm costs')
    def mst(d,removed):
        edges=sorted((W[tuple(sorted((d[i],d[j])))][0],i,j) for i,j in combinations(range(7),2) if (i,j) not in removed and tuple(sorted((d[i],d[j]))) in W)
        par=list(range(7));chosen=[]
        def find(i):
            while par[i]!=i:i=par[i]
            return i
        for cost,i,j in edges:
            a,b=find(i),find(j)
            if a!=b:par[a]=b;chosen.append((cost,i,j))
        need(len(chosen)==6,'all retained finite controls have connected generous graphs')
        for c,i,j in chosen:
            seen={i};again=True
            while again:
                again=False
                for cc,u,v in chosen:
                    if (cc,u,v)==(c,i,j):continue
                    if (u in seen)!=(v in seen):seen.update([u,v]);again=True
            need(all(cc>=c for cc,u,v in edges if (u in seen)!=(v in seen)),'every minimum-tree edge passes its exact cut test')
        return sum(c for c,i,j in chosen),chosen,edges
    topology=[]
    for wi,(d,E,M) in enumerate(C['survivors']):
        policies=[(i,{18},{4,16}) for i,v in enumerate(d) if v==24]
        policies += [(i,{24},{30}) for i,v in enumerate(d) if v==18 and 24 in d and 30 in d]
        branches=[]
        for choices in product([0,1],repeat=len(policies)):
            removed=set()
            for (i,L,R),ch in zip(policies,choices):
                forbidden=L if ch==0 else R
                removed.update(tuple(sorted((i,j))) for j,v in enumerate(d) if v in forbidden)
            cost,tree,edges=mst(d,removed)
            branches.append(dict(choices=choices,removed=sorted(removed),cost=cost,tree=tree,allowed_edges=edges,closed=cost>E))
        topology.append(dict(word_index=wi,word=d,E=E,policies=[[i,sorted(L),sorted(R)] for i,L,R in policies],branches=branches,closed=all(r['closed'] for r in branches)))
    survivors=[r for r in topology if not r['closed']]
    need([r['word_index'] for r in survivors]==[0,12],'complete actual13-word deletion and exact two-word residual')
    need(sum(len(r['branches']) for r in topology)==37,'complete joint role-choice branch count')
    cert=dict(status='FINITE-EXACT complete third-wedge and topology universes; analytic candidate pending independent audit',
        first_wedges_sha256=sha256(wr).hexdigest(),minimum_tree_sha256=sha256(raw).hexdigest(),clock=T,
        arm_domains=[dict(margins=[a,18],arms=arm_bank[a]) for a in [24,30]],arm_sum_closed=first,
        valuation_rejections=rejected,accepted=accepted,profiles=[PROFILE[k] for k in sorted(PROFILE)],
        topology=topology,remaining_words=[r['word'] for r in survivors])
    data=canonical(cert)+b'\n';(DEST/(HERE.stem+'_certificate.json')).write_bytes(data)
    print('THIRD_WEDGE margins24,18,30; complete231*141=32571 native arm products')
    print('PARTITION arm-sum32541; low-credit30; prime3 inconsistencies26; jointly realizable4, minimum forest162>107')
    print('REJECTED_FIRST',rejected[0])
    print('ACCEPTED',accepted)
    print('PROFILES',len(PROFILE),'reconstructed with all strict walls and literal minimum/maximum checks')
    print('TOPOLOGY15 original words;37 complete role branches; actual connected-complement deletion13')
    for r in topology:print(r['word_index'],'E',r['E'],'branch_costs',[b['cost'] for b in r['branches']],'closed',r['closed'])
    print('REMAINING_WORDS',[r['word'] for r in survivors])
    print('SCOPE third native wedge full13 closure needs no connectivity; two-word reduction requires actual connected complement;7200 not closed')
    print('CERTIFICATE_SHA256',sha256(data).hexdigest())
    print('PASS',GATES,'always-active exact gates; raw LF')

if __name__=='__main__':main()
