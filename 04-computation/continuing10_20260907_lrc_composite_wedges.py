#!/usr/bin/env python3
"""Complete 7200 native wedge certificates and selected ratio-tree controls.

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
    raw=(ROOT/'05-knowledge/results/continuing9_20260907_lrc_ratio_tree_certificate.json').read_bytes()
    need(sha256(raw).hexdigest()=='f93d2c8b4f112027bbba23fa529742cd12ac8107e7c9bfdf57f15d6666799533','inherited complete29 orientation certificate pin')
    OLD=json.loads(raw)
    rawP=(ROOT/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
    need(sha256(rawP).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','complete inherited global profile pin')
    J=json.loads(rawP);P={int(k):{(c,tuple(w)) for c,w in x['profiles']} for k,x in J['levels'].items()}
    slots=[(I,tuple(j for j in range(7) if j not in I)) for size in range(6,0,-1) for I in combinations(range(7),size)]
    need(len(slots)==126,'all positional proper profiles retained')
    states=sorted(a for a in J['levels']['6']['gcds'] if T%a==0)
    need(len(states)==26,'complete clock-divisor alphabet')
    def valid(d):
        if gcd(*d)!=1:return False
        for I,K in slots:
            c=gcd(*(d[i] for i in I))
            if (c,tuple(sorted(gcd(c,d[j]) for j in K))) not in P[len(K)]:return False
        return True
    def excess(d):return sum(a*((T+7*a-1)//(7*a)) for a in d)-T
    words={}
    for forced,expectedN,expectedE in [((4,24,18),1033,108),((16,24,18),949,116),((24,18,30),162,107)]:
        rows=[];attempts=0
        for rest in combinations_with_replacement(states,4):
            attempts+=1;d=tuple(sorted(forced+rest))
            if valid(d):rows.append([list(d),excess(d)])
        need(attempts==23751 and len(rows)==expectedN,'complete four-free-position conditional universe')
        need(max(E for d,E in rows)==expectedE,'attained full-profile ceiling maximum')
        rows.sort();words[forced]=dict(forced=list(forced),attempted=attempts,maximum=expectedE,words=rows)
    profile_hostile=(1,9,16,18,32,60,72)
    need(valid(profile_hostile) and excess(profile_hostile)==164,'dropping middle24 loses the conditional116 bound')
    atlas=[(p,q) for q in range(2,356) for p in range(1,min(q,357-q)) if gcd(p,q)==1 and atlas_sum(p+q)]
    need(len(atlas)==5855,'complete strict pair atlas; no friendly-prime contamination')
    arm_bank={}
    for a,b in [(4,24),(16,24),(18,24)]:
        e=gcd(a,b);n=T//e;arms=[]
        for p,q in atlas:
            for u,v in [(p,q),(q,p)]:
                if (e*gcd(n,u),e*gcd(n,v))!=(a,b):continue
                R=profile(n,p,q);arms.append([u,v,e*R['minimum']])
        need(len(arms)==(231 if a==18 else 316),'complete oriented arm domain')
        arm_bank[a,b]=arms
    wedge_results=[]
    for a,E,expected_first,expected_residual,expected_min in [(4,108,72939,57,114),(16,116,72782,214,142)]:
        A=arm_bank[a,24];C=arm_bank[18,24];residual=[];first=0;zero=0
        for (u,v,cu),(w,z,cw) in product(A,C):
            if cu+cw>E:first+=1;continue
            R=primitive([F(u,v),F(1),F(w,z)])
            need(gcd(*R)==1 and len(set(R))==3,'actual primitive distinct low-credit wedge')
            # A common scale has clock gcd exactly gcd(a,24,18)=2.
            need([gcd(T,2*r) for r in R]==[a,24,18],'both actual arm labels coexist with all three margins')
            d=gcd(R[0],R[2]);p,q=sorted((R[0]//d,R[2]//d));e=gcd(T,2*d)
            CP=profile(T//e,p,q);credit=e*CP['minimum']
            costs=[cu,cw,credit];forest=sum(sorted(costs)[-2:])
            need(forest>E,'actual three-vertex forest closes every residual ratio choice')
            zero+=cu==cw==0
            residual.append(dict(arms=[[u,v,cu],[w,z,cw]],primitive=R,scale_gcd=2,pair=[e,p,q,credit],forest=forest))
        need(first==expected_first and len(residual)==expected_residual,'complete native directed wedge product partition')
        need(first+len(residual)==72996,'no pruned or forgotten arm pairs')
        need(min(r['forest'] for r in residual)==expected_min,'sharp certified lower credit on the residual bank')
        need(zero==(5 if a==4 else 30),'complete all-zero-arm wedge subtable')
        wedge_results.append(dict(margins=[a,24,18],ceiling=E,universe=72996,arm_sum_closed=first,residual=residual,minimum_residual_credit=expected_min,zero_wedges=zero))
    # Named third pattern: endpoint alone does not pay the conditioned maximum.
    third=[3476,21093,445];third_costs=[]
    for i,j in combinations(range(3),2):
        D=gcd(third[i],third[j]);p,q=sorted((third[i]//D,third[j]//D));e=gcd(T,6*D)
        R=profile(T//e,p,q);third_costs.append([i,j,e,p,q,e*R['minimum']])
    need([r[-1] for r in third_costs]==[0,96,66],'third actual path arm and composite credits')
    need(96<=107<96+66,'single-pair near miss repaired by the genuine same-triple forest')
    # Exact reconstruction of the inherited complete selected29 orientations.
    bank=[];all_count=good_count=0
    for wi,old in enumerate(OLD['fixed_tree_bank']):
        d=old['word'];choices=[]
        for i,j,e,p,q in old['tree']:
            local=[];n=T//e
            if (e*gcd(n,p),e*gcd(n,q))==(d[i],d[j]):local.append((i,j,F(q,p)))
            if (e*gcd(n,q),e*gcd(n,p))==(d[i],d[j]):local.append((i,j,F(p,q)))
            choices.append(local)
        orientations=list(product(*choices))
        need(len(orientations)==len(old['orientations']),'all allowed local orientations reconstructed')
        for oi,edges in enumerate(orientations):
            all_count+=1;U=projective(edges);prior=old['orientations'][oi]
            need(U==prior['row'],'literal path products recover frozen primitive row')
            good=len(set(U))==7 and [gcd(T,u) for u in U]==d
            need(good==prior['valid'],'distinctness and exact clipped clock margins retained')
            if not good:continue
            good_count+=1;need(valid(d),'all126 inherited profiles on selected actual word')
            pairs=[];costs=[]
            for i,j in combinations(range(7),2):
                D=gcd(U[i],U[j]);p,q=sorted((U[i]//D,U[j]//D));e=gcd(T,D)
                record=dict(i=i,j=j,e=e,p=p,q=q,evaluated=p+q<=10000)
                if record['evaluated']:
                    R=profile(T//e,p,q);record['minimum']=e*R['minimum'];record['maximum']=e*R['maximum'];costs.append((record['minimum'],i,j))
                pairs.append(record)
            M,MT=max_tree(costs)
            need(M>old['E'] and max(c for c,i,j in costs)>old['E'],'every selected valid row closes already by a short pair')
            graph={tuple(sorted((i,j))) for i,j,e,p,q in old['tree']};tags=[]
            for mid in range(7):
                if d[mid]!=24:continue
                for i,j in combinations([v for v in range(7) if tuple(sorted((v,mid))) in graph],2):
                    if sorted([d[i],d[j]]) in [[4,18],[16,18]]:tags.append(['uniform_wedge',i,mid,j])
            for indices in product(range(7),repeat=3):
                if len(set(indices))<3:continue
                h=gcd(*(U[i] for i in indices))
                if [U[i]//h for i in indices]==third and gcd(T,h)==6:tags.append(['named_third',*indices])
            need(bool(tags),'every selected realization has a proved structural family tag')
            bank.append(dict(word_index=wi,orientation_index=oi,word=d,E=old['E'],row=U,pairs=pairs,forest_credit=M,forest=MT,families=tags))
    need((all_count,good_count)==(29,20),'full inherited selected scope, not all possible ratio trees')
    # Sharp open-endpoint control and the previous common-phase zero hostile.
    endpoint=profile(7,1,1)
    need((endpoint['minimum'],endpoint['maximum'])==(0,1),'strict danger loses all overlap at the equality wall')
    need(literal_pair(7,1,1,F(1,2))==0 and literal_pair(7,1,1,F(1,2),True)==2,'closed danger would fabricate two credits')
    U=OLD['actual_row'];native=[]
    for j in range(T):
        den=2*T;native.append([14*min(u*(2*j+1)%den,den-u*(2*j+1)%den)<den for u in U])
    need(sum(not any(row) for row in native)==2020,'inherited actual worst-tail safe count')
    need(all(sum(row[i] and row[j] for row in native)==0 for i,j,*_ in OLD['actual_graph']),'all original atlas zero minima coexist')
    # Two new nonunit disconnected full13 controls, plus the inherited control.
    body=[T*i for i in range(1,7)]
    controls=[]
    for B,expectedE,expectedSafe in [([2224,264,1602,9081,9119,32416,61140],116,2541),([20856,126558,2670,9081,16208,25475,32672],107,2504),(OLD['full13_control']['tails'],None,2476)]:
        d=[gcd(T,b) for b in B]
        need(gcd(*(body+B))==1 and len(set(body+B))==13 and min(body+B)>1,'genuine primitive distinct nonunit full13 row')
        need(valid(d),'control has every projected hereditary profile')
        if expectedE is not None:need(excess(d)==expectedE,'control attains the relevant conditional maximum')
        edges=[]
        for i,j in combinations(range(7),2):
            D=gcd(B[i],B[j]);p,q=sorted((B[i]//D,B[j]//D))
            if p+q<=356 and atlas_sum(p+q):edges.append([i,j])
        need(edges==[[0,1],[1,2]],'full strict complement graph has two arms and four isolates')
        den=7*T
        safe=[j for j in range(T) if all(14*min(v*(7*j+1)%den,den-v*(7*j+1)%den)>=den for v in body+B)]
        need(len(safe)==expectedSafe,'literal full13 weak-safety preserves the six-body phase')
        controls.append(dict(body=body,tails=B,margins=d,E=excess(d),edges=edges,alpha=[1,7],safe_count=len(safe),safe_indices_sha256=sha256(canonical(safe)).hexdigest()))
    cert=dict(status='FINITE-EXACT complete stated universes; analytic candidate pending independent audit',clock=T,
        inherited_sha256=sha256(raw).hexdigest(),profile_sha256=sha256(rawP).hexdigest(),atlas=atlas,
        conditional_words=list(words.values()),arm_domains=[dict(margins=list(k),arms=v) for k,v in arm_bank.items()],
        wedges=wedge_results,third_pattern=dict(primitive=third,scale_gcd=6,costs=third_costs,forest_credit=162,ceiling=107),
        profiles=[PROFILE[k] for k in sorted(PROFILE)],selected_bank=bank,full13_controls=controls,
        dropped_middle_profile_hostile=list(profile_hostile))
    data=canonical(cert)+b'\n';(DEST/(HERE.stem+'_certificate.json')).write_bytes(data)
    print('ATLAS 5855 complete strict ratios; oriented arm domains316,316,231')
    print('PROFILE_UNIVERSES each23751; forced(4,24,18):1033/max108; (16,24,18):949/max116; (24,18,30):162/max107')
    for R in wedge_results:print('ALL_WEDGES',R['margins'],'universe72996; arm_sum',R['arm_sum_closed'],'composite',len(R['residual']),'zero',R['zero_wedges'],'minimum',R['minimum_residual_credit'])
    print('THIRD_PATTERN (3476,21093,445), scale gcd6; exact pair credits0,96,66; forest162>107')
    print('SELECTED_BANK29 orientations;20 valid on10 words; all20 short-pair closures and structural family tags; no all15-word closure')
    print('SELECTED_FOREST_CREDITS',[(r['word_index'],r['orientation_index'],r['E'],r['forest_credit']) for r in bank])
    print('NATIVE_PROFILES',len(PROFILE),'including composite ratios outside the strict atlas; full intervals/events reconstructed and semantically pinned')
    print('CONTROLS strict endpoint0 versus false closed credit2; worst-tree allsix zeros and2020 safe; nonunit disconnected full13 safe counts',[r['safe_count'] for r in controls])
    print('SCOPE conditional actual13 safety families, arbitrary other four tails, no connectivity or decoder/box premise; remaining7200 open')
    print('CERTIFICATE_SHA256',sha256(data).hexdigest())
    print('PASS',GATES,'always-active exact gates; raw LF')

if __name__=='__main__':main()
