#!/usr/bin/env python3
"""Standalone exact controls for the S194 ordered-deletion multiplier."""
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def require(q,msg):
    if not q:raise RuntimeError(msg)

def bits(S):
    while S:
        z=S&-S;yield z.bit_length()-1;S^=z

def submasks(S):
    U=S
    while True:
        yield U
        if not U:break
        U=(U-1)&S

def valid(rows):
    n=len(rows);F=(1<<n)-1
    for i,x in enumerate(rows):
        require(0<=x<=F and not x>>i&1,'invalid row')
        for j in range(i):require(((x>>j)&1)+((rows[j]>>i)&1)==1,'invalid tournament pair')

def profiles(rows):
    n=len(rows);F=(1<<n)-1
    st=[[0]*n for _ in range(F+1)];en=[[0]*n for _ in range(F+1)];h=[0]*(F+1);h[0]=1
    rev=[F^x^(1<<i) for i,x in enumerate(rows)]
    for S in range(1,F+1):
        for v in bits(S):
            U=S^(1<<v)
            if not U:st[S][v]=en[S][v]=1
            else:
                st[S][v]=sum(st[U][w] for w in bits(U&rows[v]))
                en[S][v]=sum(en[U][w] for w in bits(U&rev[v]))
        h[S]=sum(st[S]);require(h[S]==sum(en[S]),'forward/backward mismatch')
    return h,st,en

def components(rows,S):
    n=len(rows);reach=[(rows[i]&S)|(1<<i) for i in range(n)]
    for k in bits(S):
        for i in bits(S):
            if reach[i]>>k&1:reach[i]|=reach[k]
    out=[];remaining=S
    while remaining:
        i=(remaining&-remaining).bit_length()-1
        C=sum(1<<j for j in bits(remaining) if reach[i]>>j&1 and reach[j]>>i&1)
        out.append(C);remaining^=C
    snapshot=tuple(out)
    out.sort(key=lambda C:-sum(bool(rows[(C&-C).bit_length()-1]&D) for D in snapshot if D!=C))
    for i,C in enumerate(out):
        for D in out[i+1:]:require(all(rows[a]&D==D for a in bits(C)),'bad SCC order')
    return out

def q_signed(h,S):
    return sum(h[U]*h[S^U] for U in submasks(S)),sum((-1)**U.bit_count()*h[U]*h[S^U] for U in submasks(S))

def response(h,st,en,S,a):
    R=sum(st[U][a]*h[S^U] for U in submasks(S))
    L=sum(en[U][a]*h[S^U] for U in submasks(S))
    signedR=sum((-1)**(U.bit_count()-1)*st[U][a]*h[S^U] for U in submasks(S) if U)
    signedL=sum((-1)**(U.bit_count()-1)*en[U][a]*h[S^U] for U in submasks(S) if U)
    require(signedR==en[S][a] and signedL==st[S][a],'root parity failed')
    require(R>=h[S] and L>=h[S],'root split injection bound failed')
    return R,L,en[S][a],st[S][a]

def B(n):
    if n==1:return 1
    if n<3:raise ValueError(n)
    q,r=divmod(n,3)
    return (3*5**(q-1) if r==0 else 5**q if r==1 else 9*5**(q-1))

def g(n):
    if n<4:raise ValueError(n)
    if n<=8:return (2,4,6,10,16)[n-4]
    q,r=divmod(n-9,3);return (26,46,78)[r]*5**q

def literal_subsets(rows):
    n=len(rows);h=[0]*(1<<n);h[0]=1
    for k in range(1,n+1):
        for p in permutations(range(n),k):
            if all(rows[a]>>b&1 for a,b in zip(p,p[1:])):h[sum(1<<a for a in p)]+=1
    return h

def literal_selected(rows,v,cs,a,b,h):
    n=len(rows);F=(1<<n)-1;total=0;num=0
    for length in range(2,n,2):
        for P in permutations([x for x in range(n) if x!=v],length):
            if P[0]!=a or P[-1]!=b:continue
            cyc=(v,)+P
            if not all(rows[cyc[i]]>>cyc[(i+1)%len(cyc)]&1 for i in range(len(cyc))):continue
            mask=sum(1<<x for x in P)
            require(mask&cs[0] and mask&cs[-1],'selected endpoint absent')
            ranks=[next(i for i,C in enumerate(cs) if C>>x&1) for x in P]
            require(ranks==sorted(ranks),'cycle violates SCC order')
            total+=2*h[F^mask^(1<<v)];num+=1
    return total,num

def analyse(rows,h,st,en,v,literal=False):
    n=len(rows);F=(1<<n)-1;S=F^(1<<v);cs=components(rows,S)
    if len(cs)<2:return None
    aa=[a for a in bits(cs[0]) if rows[v]>>a&1]
    bb=[b for b in bits(cs[-1]) if rows[b]>>v&1]
    require(aa and bb,'strongness misses SCC incidences')
    delta=h[F]-h[S];P=1
    for C in cs[1:-1]:
        Q,signed=q_signed(h,C);require(signed==0 and Q>=2*h[C],'intermediate convolution failed');P*=Q
    counts=0;minselected=None;maxselected=None
    for a in aa:
        R,_,e,_=response(h,st,en,cs[0],a)
        for b in bb:
            _,L,_,s=response(h,st,en,cs[-1],b)
            selected=R*L*P+(e*s if len(cs)==2 else 0)
            require(selected<=delta and selected%2==0,'selected family inequality failed')
            if literal:
                z,c=literal_selected(rows,v,cs,a,b,h);require(z==selected,'literal selected root-cycle family mismatch');counts+=c
            minselected=selected if minselected is None else min(minselected,selected)
            maxselected=selected if maxselected is None else max(maxselected,selected)
    lower=h[S]+1 if len(cs)==2 else (2**(len(cs)-2))*h[S]
    require(delta>=lower and delta>h[S] and 2*delta>h[F],'strict weak-deletion ratio failed')
    if n>=4:
        require(h[F]>=B(n),'inherited Busch finite control failed')
        require(delta>g(n),'weak deletion beat or tied g')
    return dict(k=len(cs),delta=delta,H=h[F],H_deleted=h[S],lower=lower,selected_min=minselected,selected_max=maxselected,literal_cycles=counts)

def labelled(n):
    pairs=list(combinations(range(n),2))
    for z in range(1<<len(pairs)):
        r=[0]*n
        for i,(a,b) in enumerate(pairs):
            if z>>i&1:r[a]|=1<<b
            else:r[b]|=1<<a
        yield tuple(r)

def main():
    controls=Counter();bankstats=Counter();weakmins={};residualstats=Counter();qchecks=0
    for n in range(1,6):
        for rows in labelled(n):
            valid(rows);h,st,en=profiles(rows);require(literal_subsets(rows)==h,'literal Hamilton words failed');controls['labelled_tournaments']+=1
            F=(1<<n)-1
            for S in range(1,F+1):
                Q,z=q_signed(h,S);require(z==0 and Q>=2*h[S],'signed convolution failed');qchecks+=1
            if n>=3 and len(components(rows,F))==1:
                for v in range(n):
                    z=analyse(rows,h,st,en,v,True)
                    if z:
                        controls['weak_deletion_roots']+=1;controls['literal_selected_root_cycles']+=z['literal_cycles']
    bank=HERE/'inherited_tournaments_through8.txt';lines=bank.read_text().splitlines()
    require(int(lines[0])==len(lines)-1==7412,'bank size changed')
    require(sha256(bank.read_bytes()).hexdigest()=='5d5a5f625b4a85e04c0744aeb847ed36f59456f6a39996479aa82b9eb950afe3','bank hash changed')
    for line in lines[1:]:
        n,*rows=map(int,line.split());rows=tuple(rows);require(len(rows)==n,'bad bank shape');valid(rows)
        F=(1<<n)-1
        if n<3 or len(components(rows,F))!=1:continue
        h,st,en=profiles(rows)
        for v in range(n):
            bankstats[f'order_{n}_strong_roots']+=1
            z=analyse(rows,h,st,en,v)
            if z:
                bankstats[f'order_{n}_weak_deletion_roots']+=1;bankstats[f'scc_count_{z["k"]}']+=1
                if str(n) not in weakmins or z['delta']<weakmins[str(n)]['delta']:weakmins[str(n)]={**z,'rows':rows,'root':v}
            A=F^rows[v]^(1<<v);BB=rows[v];ac=components(rows,A);bc=components(rows,BB)
            extreme=any(rows[x]&ac[0] for x in bits(bc[-1]))
            if not extreme:
                residualstats[f'order_{n}_all_residual']+=1
                residualstats[f'order_{n}_weak_deletion' if z else f'order_{n}_strong_deletion']+=1
    triangle=(2,4,1);h,st,en=profiles(triangle);boundary=analyse(triangle,h,st,en,0,True)
    require(boundary['k']==2 and boundary['delta']==boundary['lower']==2,'singleton endpoint boundary failed')
    hostile=(2,4,9,19,7);h,st,en=profiles(hostile);require(h[-1]==9 and h[30]==5 and (h[7]-h[6])*h[28]==6,'merge hostile failed')
    require(len(components(hostile,7))==len(components(hostile,28))==1,'hostile factors not strong')
    for n in range(4,1001):require(2*g(n)<B(n),'B/2 versus g arithmetic failed')
    return dict(status='PASS',labelled_controls=dict(controls),signed_subset_convolution_checks=qchecks,bank_sha256=sha256(bank.read_bytes()).hexdigest(),bank_counts=dict(bankstats),weak_deletion_minima=weakmins,remaining_sector_counts=dict(residualstats),singleton_boundary=boundary,merge_hostile={'rows':hostile,'root':0,'A':7,'B':28,'H':9,'H_deleted':5,'delta':4,'false_product':6},arithmetic_control_orders=[4,1000],scope='All-order proof is in PROOF.md; finite census is a control, not a proof of unrestricted deletion minimum or spectrum completeness.')

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
