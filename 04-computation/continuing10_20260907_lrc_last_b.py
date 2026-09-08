#!/usr/bin/env python3
"""Remaining B word: exact doubled zero-arm wedges and forced bridge.

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
    raw=(ROOT/'05-knowledge/results/continuing8_20260906_lrc_minimum_tree_certificate.json').read_bytes()
    need(sha256(raw).hexdigest()=='580a7c930103aab3bea867ad463a90b0e0208323a90ee95a685ff811a761582d','complete inherited atlas minimum-weight table')
    C=next(r for r in json.loads(raw)['clocks'] if r['t']==7200);W={tuple(d):v for d,v in C['weights']}
    d,E,M=C['survivors'][12]
    need((d,E,M)==([5,8,9,30,32,36,48],103,102),'exact remaining B word and budget')
    need(sum(a*((T+7*a-1)//(7*a)) for a in d)-T==103,'literal ceiling excess')
    all_edges=[(W[tuple(sorted((d[i],d[j])))][0],i,j) for i,j in combinations(range(7),2) if tuple(sorted((d[i],d[j]))) in W]
    low=sorted(e for e in all_edges if e[0]<=103)
    expected=sorted([(0,0,3),(102,2,3),(0,2,5),(0,5,6),(0,1,6),(0,1,4),(0,4,6)])
    need(low==expected,'complete generous graph below the one-edge103 budget')
    chain={(0,3),(2,3),(2,5),(5,6)};attachment={(1,6),(4,6)};trees=[]
    for candidate in combinations(low,6):
        reached={0}
        while True:
            nxt=reached|{v for c,i,j in candidate for u,v in [(i,j),(j,i)] if u in reached}
            if nxt==reached:break
            reached=nxt
        if len(reached)!=7:continue
        es={(i,j) for c,i,j in candidate}
        need(chain<=es and bool(attachment&es),'every possible spanning tree forces the full chain and one doubled wedge')
        need(sum(c for c,i,j in candidate)==102,'all three possible low trees have the102 bridge')
        trees.append(candidate)
    need(len(trees)==3,'complete seven-edge possible graph has exactly three spanning trees')
    atlas=[(p,q) for q in range(2,356) for p in range(1,min(q,357-q)) if gcd(p,q)==1 and atlas_sum(p+q)]
    need(len(atlas)==5855,'strict atlas remains complete under the changed clocks')
    arm_bank={}
    for a in [8,32,36]:
        e=gcd(a,48);n=T//e;arms=[]
        for p,q in atlas:
            for u,v in [(p,q),(q,p)]:
                if (e*gcd(n,u),e*gcd(n,v))!=(a,48):continue
                R=profile(n,p,q);arms.append([u,v,e*R['minimum']])
        need(len(arms)=={8:316,32:633,36:231}[a],'fresh complete directed arm universe at changed sheet quotient')
        need(all(c==0 or c>=e>1 for u,v,c in arms),'any positive arm credit plus the forced102 bridge exceeds103')
        zeros=[row for row in arms if row[2]==0]
        need(len(zeros)=={8:10,32:51,36:18}[a],'complete doubled zero-arm alphabet, including saturation expansion')
        arm_bank[a]=dict(margins=[a,48],e=e,n=n,arms=arms,zero_arms=zeros)
    wedges=[]
    for a,expected_count,expected_min in [(8,180,112),(32,918,140)]:
        rows=[]
        for (u,v,_),(w,z,_) in product(arm_bank[a]['zero_arms'],arm_bank[36]['zero_arms']):
            R=primitive([F(u,v),F(1),F(w,z)])
            need(len(set(R))==3 and [gcd(T,4*r) for r in R]==[a,48,36],'all joint primitive doubled paths satisfy actual margins')
            D=gcd(R[0],R[2]);p,q=sorted((R[0]//D,R[2]//D));e=gcd(T,4*D)
            need(e==4,'actual endpoint sheet multiplicity is4, not the old value2')
            Q=profile(T//e,p,q);credit=e*Q['minimum']
            need(credit>103,'fresh effective1800-clock endpoint credit closes remaining B')
            rows.append(dict(arms=[[u,v],[w,z]],primitive=R,scale_gcd=4,pair=[e,p,q,credit]))
        need(len(rows)==expected_count and min(r['pair'][-1] for r in rows)==expected_min,'complete changed-clock product count and attained lower credit')
        wedges.append(dict(margins=[a,48,36],rows=rows,minimum=expected_min))
    # Inherited old clocks are not scaled blindly: p/q578/801 is a strict test.
    old=profile(3600,578,801);new=profile(1800,578,801)
    need(2*old['minimum']==114 and 4*new['minimum']==112,'doubled margin does not preserve old credit')
    # Positive actual graph control: full row, explicit word, no physical unit.
    B=list(d);A=[T*j for j in range(1,7)];edges=[]
    for i,j in combinations(range(7),2):
        D=gcd(B[i],B[j]);p,q=sorted((B[i]//D,B[j]//D))
        if p+q<=356 and atlas_sum(p+q):edges.append([i,j])
    reached={0}
    while True:
        nxt=reached|{v for i,j in edges for u,v in [(i,j),(j,i)] if u in reached}
        if nxt==reached:break
        reached=nxt
    need(len(reached)==7 and gcd(*(A+B))==1 and len(set(A+B))==13 and min(A+B)>1,'genuine primitive distinct connected-complement full13 positive control')
    den=7*T
    safe=[j for j in range(T) if all(14*min(v*(7*j+1)%den,den-v*(7*j+1)%den)>=den for v in A+B)]
    need(bool(safe),'literal full13 weak-safe lifts at the proper-six1/7 phase')
    cert=dict(status='FINITE-EXACT remaining B closure; complete native zero-arm scope, analytic candidate pending audit',
        inherited_sha256=sha256(raw).hexdigest(),word=d,E=E,all_edges=all_edges,cheap_graph=low,all_cheap_spanning_trees=trees,
        arm_domains=[arm_bank[a] for a in [8,32,36]],wedges=wedges,
        profiles=[PROFILE[k] for k in sorted(PROFILE)],scale_control=dict(ratio=[578,801],old_credit=114,new_credit=112),
        positive_control=dict(A=A,B=B,edges=edges,alpha=[1,7],safe_count=len(safe),safe_indices_sha256=sha256(canonical(safe)).hexdigest()))
    data=canonical(cert)+b'\n';(DEST/(HERE.stem+'_certificate.json')).write_bytes(data)
    print('REMAINING_B',d,'E103; complete low graph7 edges,3 spanning trees; forced5--30--9--36--48 and one triangle attachment')
    print('FORCED_BRIDGE9--30 credit102; every positive relevant arm credit is greater than1; zero arms required under failure')
    print('FRESH_DOMAINS 8/48:316,zero10;32/48:633,zero51;36/48:231,zero18')
    print('EXACT_ZERO_WEDGES180+918, all joint margins realized; endpoint minima112 and140 on quotient clock1800, sheet4')
    print('SCALE_HOSTILE578/801: old7200-sheet2 credit114; new7200-sheet4 credit112, so no unchanged-credit scaling')
    print('PROFILES',len(PROFILE),'all strict walls; literal minimum and maximum controls')
    print('POSITIVE_FULL13 nonunit primitive row with the exact B word and connected actual graph; safe lifts',len(safe))
    print('SCOPE this complete connected-complement word closes; no claim on arbitrary doubled wedges with a larger excess')
    print('CERTIFICATE_SHA256',sha256(data).hexdigest())
    print('PASS',GATES,'always-active exact gates; raw LF')

if __name__=='__main__':main()
