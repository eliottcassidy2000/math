#!/usr/bin/env python3
"""Finite exact controls for GOLD_SUPPORT_THREE_HALL_BARRIER.md.
Self-contained standard-library audit. Never relies on assert or an optimizer.
Run in any directory; certificate paths resolve beside this file.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations, permutations
from collections import Counter
import json, hashlib

HERE=Path(__file__).resolve().parent
CHECKS=Counter()
def need(ok,tag):
    if not ok: raise RuntimeError(tag)
    CHECKS[tag]+=1

def close(pred):
    p=list(pred)
    for k in range(len(p)):
        for v in range(len(p)):
            if p[v]>>k&1:p[v]|=p[k]
    return None if any(p[v]>>v&1 for v in range(len(p))) else tuple(p)

def add(p,*rels):
    q=list(p)
    for a,b in rels:q[b]|=1<<a
    return close(q)

@lru_cache(None)
def count(p):
    if p is None:return 0
    n=len(p);full=(1<<n)-1
    @lru_cache(None)
    def visit(s):
        if s==full:return 1
        return sum(visit(s|1<<v)for v in range(n)if not s>>v&1 and not p[v]&~s)
    return visit(0)

def words(p):
    out=[];n=len(p)
    def visit(mask,w):
        if len(w)==n:out.append(tuple(w));return
        for v in range(n):
            if not mask>>v&1 and not p[v]&~mask:visit(mask|1<<v,w+[v])
    visit(0,[]);return out

def before(w,a,b):return w.index(a)<w.index(b)

def upsets(p):return tuple(sum(1<<v for v in range(len(p))if p[v]>>u&1)for u in range(len(p)))

def cycle(n,c):
    sig=list(range(n))
    for a,b in zip(c,c[1:]+c[:1]):sig[a]=b
    return tuple(sig)

def smallcycles(n):
    for c in combinations(range(n),2):yield c
    for c in combinations(range(n),3):
        yield c;yield (c[0],c[2],c[1])

def relabel(p,sig):
    q=[0]*len(p)
    for b in range(len(p)):
        for a in range(len(p)):
            if p[b]>>a&1:q[sig[b]]|=1<<sig[a]
    return tuple(q)

def join(p,q):return close(tuple(a|b for a,b in zip(p,q)))

def good_candidates(p):
    U=upsets(p);out=[]
    for dual,(D,V)in enumerate(((p,U),(U,p))):
        for a,b in permutations(range(len(p)),2):
            if (D[a]|V[a])>>b&1 or D[a]&~D[b]:continue
            diff=[v for v in range(len(p))if (V[b]&~V[a])>>v&1]
            if not all((D[x]>>y|D[y]>>x)&1 for x,y in combinations(diff,2)):continue
            N=count(add(p,(b,a)if dual else(a,b)))
            out.append((dual,a,b,N))
    return out

def verify_crossing_classification():
    # x=0,X=1,y=2,Y=3, sole possible outside label z=4.
    p=(0,1,0,4,0);ws=words(p);legal=set(ws);checked=0
    allowed={cycle(5,(0,4,1)),cycle(5,(4,2,3))}
    for w in ws:
        if not(before(w,0,2)and before(w,1,3)):continue
        for c in smallcycles(5):
            sig=cycle(5,c);v=tuple(sig[x]for x in w)
            if v in legal and before(v,2,0)and before(v,3,1):
                need(sig in allowed,'five-label crossing classification');checked+=1
    need(checked==2,'nonvacuous crossing patterns')
    return checked

def main():
    crossing=verify_crossing_classification()
    P=(0,0,0,6,0,19,2,18,86,151);R=add(P,(5,9))
    need(close(P)==P,'P transitive');need(R==(0,0,0,6,0,19,2,18,86,183),'R closure')
    E,b=count(P),count(R);need((E,b)==(11971,7741),'initial counts')
    need(3*E<5*b and 3*b<2*E,'strict heavy window')
    deg=[len(P)-1-(p|u).bit_count()for p,u in zip(P,upsets(P))]
    need(max(deg)==7,'exact seven-thin degree')
    gc=good_candidates(P);need(gc and all(2*z[3]>E for z in gc),'no good pair in either dual')
    table=[(a,z,count(add(R,(a,z))))for a,z in combinations(range(10),2)]
    eligible=[(a,z,c)for a,z,c in table if 2*b-E<=c<=E-b]
    need(eligible==[(0,2,3738)],'unique second pair')
    W=words(R);need(len(W)==b and W==sorted(W),'independent child enumeration and index order')
    A=[w for w in W if before(w,0,2)];B=[w for w in W if before(w,2,0)]
    need((len(A),len(B))==(3738,4003),'heavy second cells')
    need([E-b-len(A),E-b-len(B)]==[492,227],'actual Gold margins')
    x,X,y,Y=0,5,2,3;n=10
    S=add(R,(x,y),(X,Y));T=add(R,(y,x),(Y,X))
    N0=add(R,(y,x),(X,Y));need((count(S),count(N0))==(2616,1844),'Hall base counts')
    RU=upsets(R)
    Z=[z for z in range(n)if z not in(x,X,y,Y)and all(not((R[z]|RU[z])>>v&1)for v in(x,X,y,Y))]
    need(Z==[6,7],'common incomparable escape labels')
    records=[];targetsets=[]
    for z in Z:
        for c in[(x,z,X),(z,y,Y)]:
            sig=cycle(n,c);J=join(T,relabel(S,sig));d=count(J)
            records.append({'cycle':c,'predecessors':J,'count':d})
            inv=[sig.index(i)for i in range(n)]
            targets={w for w in words(J)}
            need(len(targets)==d,'escape count enumeration')
            need(all(tuple(inv[v]for v in w)in set(words(S))for w in targets),'escape inverse domains')
            targetsets.append(targets)
    need([q['count']for q in records]==[46,72,55,48],'four escape values')
    D=len(A)-count(S)+count(N0);cap=D+sum(q['count']for q in records)
    need((D,cap)==(2966,3187)and cap<2*b-E,'structural capacity obstruction')
    need(len(set().union(*targetsets))==184,'escape union overlap')
    Bid={w:j for j,w in enumerate(B)};adj=[set()for _ in A]
    for c in smallcycles(n):
        sig=cycle(n,c)
        for i,w in enumerate(A):
            j=Bid.get(tuple(sig[v]for v in w))
            if j is not None:adj[i].add(j)
    edgecount=sum(map(len,adj));need(edgecount==30897,'all support-three edges')
    need(all(adj)and len(set().union(*adj))==len(B),'no isolated vertices')
    # Independent full-word graph, sharing no cycle-generation logic.
    hamming_edges=0;four_edges=0
    for i,a in enumerate(A):
        for j,bw in enumerate(B):
            d=0
            for aa,bb in zip(a,bw):
                d+=aa!=bb
                if d>4:break
            if d<=4:four_edges+=1
            if d<=3:
                hamming_edges+=1;need(j in adj[i],'Hamming edge generated')
    need(hamming_edges==edgecount and four_edges==157835,'independent graph equality')
    certificates=[]
    for k,expected in[(3,3017),(4,3738)]:
        path=HERE/f'rigid_support{k}_certificate.json';data=json.loads(path.read_text())
        need(tuple(data['P'])==P and tuple(data['R'])==R,'certificate poset identity')
        pairs=data['matching'];L=[i for i,j in pairs];RR=[j for i,j in pairs]
        need(len(pairs)==expected and len(set(L))==expected and len(set(RR))==expected,'matching capacities')
        need(all(0<=i<len(A)and 0<=j<len(B)for i,j in pairs),'matching index range')
        need(all(sum(x!=y for x,y in zip(A[i],B[j]))<=k for i,j in pairs),'matching Hamming support')
        if k==3:
            CL=set(data['cover_left']);CR=set(data['cover_right'])
            need((len(CL),len(CR))==(1587,1430),'support-three cover sizes')
            need(all(i in CL or j in CR for i,row in enumerate(adj)for j in row),'cover every generated edge')
            need(len(CL)+len(CR)==expected,'optimal support-three matching')
        else:need(set(L)==set(range(len(A))),'support-four saturation')
        certificates.append({'support':k,'matching_size':expected,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    need(E-2*b+3738==227,'support-four Gold certificate consequence')
    return {'status':'PASS','scope':'One specified seven-thin seed and first comparison; not a poset census or universal Gold theorem','E':E,'b':b,'first':[5,9],'unique_Gold_second':eligible,'required_matching':2*b-E,'actual_heavy_margins':[492,227],'common_escape_labels':Z,'escape_counts':records,'structural_support3_cap':cap,'optimal_support3_capacity':3017,'support3_edges':edgecount,'support3_isolated_vertices':[0,0],'support4_matching':3738,'support4_edges':four_edges,'certified_support4_margin':227,'second_comparison_counts':table,'five_label_crossing_controls':crossing,'certificates':certificates,'checks':sum(CHECKS.values())}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
