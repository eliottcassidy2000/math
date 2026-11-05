#!/usr/bin/env python3
"""Dedicated controls for STRONGER_RATIO.md; reuses the frozen local helper module."""
import json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import verify as v

def attach(rows,X):
    n=len(rows)
    return tuple(row|((not(X>>i&1))<<n) for i,row in enumerate(rows))+(X,)

def main():
    local=0;single=0
    for n in range(1,6):
        for rows in v.labelled(n):
            h,st,en=v.profiles(rows);F=(1<<n)-1
            responses=[v.response(h,st,en,F,a) for a in range(n)]
            for X in range(1,F+1):
                r=attach(rows,X);hh,_,_=v.profiles(r);d=hh[-1]-h[F]
                R=sum(responses[a][0] for a in v.bits(X));e=sum(en[F][a] for a in v.bits(X))
                v.require(d+R+e>=2*h[F],'local compensation failed');local+=1
                if X.bit_count()==1:v.require(d==R-e,'singleton local identity failed');single+=1
    stats=Counter();equal=[];checked=0
    for line in (v.HERE/'inherited_tournaments_through8.txt').read_text().splitlines()[1:]:
        n,*rows=map(int,line.split());rows=tuple(rows);F=(1<<n)-1
        if n<3 or len(v.components(rows,F))!=1:continue
        h,st,en=v.profiles(rows)
        for root in range(n):
            S=F^(1<<root);cs=v.components(rows,S)
            if len(cs)<2:continue
            D=h[F]-h[S];v.require(D>=2*h[S],'weak-deletion two-thirds failed');checked+=1
            stats[f'order_{n}_weak_roots']+=1
            if D==2*h[S]:stats[f'order_{n}_equality']+=1
            condition=False
            if len(cs)==2:
                C,E=cs;X=rows[root]&C;Y=E&~rows[root]
                dC=h[C|(1<<root)]-h[C];dE=h[E|(1<<root)]-h[E]
                R=sum(v.response(h,st,en,C,a)[0] for a in v.bits(X));e=sum(en[C][a] for a in v.bits(X))
                L=sum(v.response(h,st,en,E,a)[1] for a in v.bits(Y));s=sum(st[E][a] for a in v.bits(Y))
                v.require(D==h[E]*dC+h[C]*dE+R*L+e*s,'exact two-SCC decomposition failed')
                lower=2*h[C]*h[E]+(R-h[C])*(L-h[E])+(h[C]-e)*(h[E]-s)
                v.require(D>=lower>=2*h[S],'expanded surplus inequality failed');stats['two_scc_exact_checks']+=1
                condition=(X.bit_count()==Y.bit_count()==1 and R==h[C] and L==h[E] and (C.bit_count()==1 or E.bit_count()==1))
            v.require((D==2*h[S])==condition,'exact equality characterization failed')
            if h[S]>=3 and h[F]-3*h[S] in (2,4):
                digit=h[F]-3*h[S]
                v.require(len(cs)==2,'weak ternary birth has more than two SCCs')
                C,E=cs;X=rows[root]&C;Y=E&~rows[root]
                if C.bit_count()>1 and E.bit_count()>1:
                    condition=all(C0.bit_count()==3 for C0 in cs) and h[S]==9 and digit==4 and X.bit_count()==Y.bit_count()==1
                    stats['exceptional_9_to_31_roots']+=1
                elif E.bit_count()==1:
                    R=sum(v.response(h,st,en,C,a)[0] for a in v.bits(X))
                    condition=(X.bit_count()==1 and R==h[S]+digit//2) or (C.bit_count()==3 and X.bit_count()==2 and digit==4)
                else:
                    L=sum(v.response(h,st,en,E,a)[1] for a in v.bits(Y))
                    condition=(Y.bit_count()==1 and L==h[S]+digit//2) or (E.bit_count()==3 and Y.bit_count()==2 and digit==4)
                v.require(condition,'exact weak ternary presentation classification failed')
                stats['weak_ternary_births']+=1
    families=[]
    for m in range(1,9):
        rows=[2,0]+[0]*m
        for z in range(2,m+2):
            rows[1]|=1<<z;rows[z]|=1
            for w in range(z+1,m+2):rows[z]|=1<<w
        n=len(rows);S=tuple(x|(1<<n) for x in rows)+(0,)
        T=attach(S,1);hh,_,_=v.profiles(T)
        v.require(len(v.components(T,(1<<len(T))-1))==1,'equality family not strong')
        v.require(hh[-1]==3*(2**m+1) and hh[(1<<len(S))-1]==2**m+1,'equality family counts failed')
        families.append(dict(m=m,order=len(T),H=hh[-1],H_deleted=2**m+1,rows=T))
    return dict(status='PASS',equality_families=families,local_cut_checks=local,singleton_cut_identities=single,bank_weak_deletion_roots=checked,bank_statistics=dict(stats),hashes={name:sha256((v.HERE/name).read_bytes()).hexdigest() for name in ('PROOF.md','STRONGER_RATIO.md','verify.py','verify_stronger.py')})
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
