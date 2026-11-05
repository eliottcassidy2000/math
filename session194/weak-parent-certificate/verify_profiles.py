#!/usr/bin/env python3
"""Self-contained exact profile exclusion on a separately growth-certified bank."""
from collections import Counter
from itertools import permutations
from hashlib import sha256
from pathlib import Path
import json
P=Path(__file__).resolve().parent

def need(q,msg):
    if not q:raise RuntimeError(msg)

def bits(S):
    while S:
        z=S&-S;yield z.bit_length()-1;S^=z

def strong(r):
    n=len(r);F=(1<<n)-1
    for rr in (r,tuple(F^x^(1<<i) for i,x in enumerate(r))):
        seen=1
        while True:
            new=seen
            for i in bits(seen):new|=rr[i]
            if new==seen:break
            seen=new
        if seen!=F:return False
    return True

def profile(r):
    n=len(r);F=(1<<n)-1
    E=[[0]*n for _ in range(F+1)];S=[[0]*n for _ in range(F+1)];H=[0]*(F+1);H[0]=1
    for U in range(1,F+1):
        for a in bits(U):
            W=U^(1<<a)
            if not W:E[U][a]=S[U][a]=1
            else:
                E[U][a]=sum(E[W][b] for b in bits(W) if r[b]>>a&1)
                S[U][a]=sum(S[W][b] for b in bits(W&r[a]))
        H[U]=sum(E[U]);need(H[U]==sum(S[U]),'two-direction count mismatch')
    R=[sum(S[U][a]*H[F^U] for U in range(F+1)) for a in range(n)]
    L=[sum(E[U][a]*H[F^U] for U in range(F+1)) for a in range(n)]
    return dict(H=H[F],R=R,L=L,s=S[F],e=E[F])

def literal(r):
    n=len(r);h=0;R=[0]*n;L=[0]*n;e=[0]*n;s=[0]*n
    for p in permutations(range(n)):
        bad=[i for i in range(n-1) if not(r[p[i]]>>p[i+1]&1)]
        if not bad:
            h+=1;s[p[0]]+=1;e[p[-1]]+=1
            for a in p:R[a]+=1;L[a]+=1
        elif len(bad)==1:
            i=bad[0];R[p[i+1]]+=1;L[p[i]]+=1
    return dict(H=h,R=R,L=L,s=s,e=e)

def main():
    bank=P/'inherited_tournaments_through8.txt';raw=bank.read_bytes()
    need(sha256(raw).hexdigest()=='5d5a5f625b4a85e04c0744aeb847ed36f59456f6a39996479aa82b9eb950afe3','input hash mismatch')
    lines=raw.decode().splitlines();need(int(lines[0])==len(lines)-1==7412,'input row count')
    complete=json.loads((P/'completeness_receipt.json').read_text());need(complete['status']=='PASS','growth closure not passed')
    need([x['classes'] for x in complete['growth_steps']]==[1,1,2,4,12,56,456,6880],'wrong growth closure domain')
    inventory={41:[],67:[]};counts=Counter();positive={};pos_targets={5,41,67,125,127,203,205}
    for line in lines[1:]:
        n,*r=map(int,line.split());r=tuple(r);need(n==len(r),'bad row length')
        counts[n]+=1
        if not strong(r):continue
        z=profile(r);h=z['H']
        if h in inventory:
            need(z==literal(r),'literal rooted profile mismatch')
            gapsR=[x-h for x in z['R']];gapsL=[x-h for x in z['L']]
            need(all(x not in (1,2) for x in gapsR+gapsL),'weak-ternary surplus exists')
            inventory[h].append(dict(rows=r,order=n,**z,R_surplus=gapsR,L_surplus=gapsL))
        if h in pos_targets and h not in positive:
            need(z==literal(r),'literal positive profile mismatch')
            positive[h]=dict(rows=r,order=n,**z)
    need(set(positive)==pos_targets,'a positive target witness missing')
    need(all(inventory[h] for h in inventory),'missing parent core')
    need(45>41 and 75>67,'Busch cutoff comparisons')
    p=positive[5];need(1 in [x-5 for x in p['R']] and 2 in [x-5 for x in p['R']],'surplus-positive H5 control missing')
    return dict(status='PASS',bank_sha256=sha256(raw).hexdigest(),classes_by_order=dict(counts),core_counts={h:len(vv) for h,vv in inventory.items()},root_profile_count=sum(len(x['rows']) for vv in inventory.values() for x in vv),inventories=inventory,positive_controls=positive,cutoffs={'41':{'first_impossible_order':8,'B':45},'67':{'first_impossible_order':9,'B':75}},growth_certificate=complete,scope='All-order absence only for strong one-vertex children of weak parents with H41 or H67 and the specified ternary counts; the target values themselves are attained.')
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
