#!/usr/bin/env python3
"""Exact skew-double checks. Standard library; exception gates survive -O."""
from itertools import combinations, product
from collections import Counter
import json
from math import comb

def require(p,msg):
    if not p: raise RuntimeError(msg)
def tournament(n,mask):
    a=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        u,v=(i,j) if (mask>>k)&1 else (j,i)
        a[u]|=1<<v
    return a

def doubling(a):
    n=len(a); b=[0]*(2*n)
    for i in range(n):
        for j in range(n):
            if i!=j:
                if (a[i]>>j)&1:
                    b[i]|=1<<j; b[i]|=1<<(n+j); b[n+i]|=1<<j
                else: b[n+i]|=1<<(n+j)
        b[i]|=1<<(n+i)
    return b

def transitive(a,mask):
    # Independent degree-characterization: transitive iff induced outdegrees 0..k-1.
    ds=sorted((a[i]&mask).bit_count() for i in range(len(a)) if (mask>>i)&1)
    return ds==list(range(mask.bit_count()))

def profiles(a):
    n=len(a); counts=Counter(); best=0; maximizers=[]
    for s in range(1<<n):
        if transitive(a,s):
            k=s.bit_count(); counts[k]+=1
            if k>best: best=k; maximizers=[s]
            elif k==best: maximizers.append(s)
    return counts,best,maximizers

def alternating(l):
    n=2*l+1; a=[0]*n
    for i,j in combinations(range(n),2):
        if j<l: edge=(i,j)
        elif i>=l: edge=(j,i)
        else:
            uj=i+1; wi=j-l+1
            edge=(j,i) if wi<=uj else (i,j)
        a[edge[0]]|=1<<edge[1]
    return a

def copy_cnf_accept(a, prime_mask):
    for triple in combinations(range(len(a)),3):
        mask=sum(1<<i for i in triple)
        labels=[(prime_mask>>i)&1 for i in triple]
        if not transitive(a,mask):
            if sum(labels)!=2:return False
        else:
            order=sorted(triple,key=lambda i:(a[i]&mask).bit_count(),reverse=True)
            if tuple((prime_mask>>i)&1 for i in order)==(1,0,1):return False
    return True

def paired_profile(a):
    """Lossless 3-state algorithm; ignores all twin-containing subsets."""
    n=len(a); d=doubling(a); z=Counter()
    for state in product(range(3), repeat=n):
        mask=sum(1<<(i+(n if c==2 else 0)) for i,c in enumerate(state) if c)
        if transitive(d,mask): z[mask.bit_count()]+=1
    base,t,_=profiles(a)
    for k,c in base.items():
        if k: z[k+1]+=c*k*(k+1)//2
    return z

def check_tournament(a):
    n=len(a); base,t,_=profiles(a); d=doubling(a)
    twins=Counter(); total=Counter(); largest=0
    for s in range(1<<(2*n)):
        if not transitive(d,s): continue
        k=s.bit_count(); largest=max(largest,k); total[k]+=1
        common=(s&((1<<n)-1)) & (s>>n)
        require(common.bit_count()<=1,'multiple twins')
        if common:
            projected=(s&((1<<n)-1))|(s>>n)
            require(transitive(a,projected),'nontransitive twin projection')
            require(k<=t+1,'twin bound')
            twins[k]+=1
    expected={k+1:c*k*(k+1)//2 for k,c in base.items() if k}
    require(dict(twins)==expected,('twin enumerator',a,dict(twins),expected))
    require(largest<=max(n,t+1),'global bound')
    require(paired_profile(a)==total,'paired polynomial transfer')
    for p in range(1<<n):
        mask=(((1<<n)-1)^p)|(p<<n)
        require(copy_cnf_accept(a,p)==transitive(d,mask),'3CNF encoding')
    return largest

def main():
    tested=0
    for n in range(1,5):
        for m in range(1<<(n*(n-1)//2)):
            check_tournament(tournament(n,m));tested+=1
    # Every labeled 5-vertex tournament: 1,024 x 1,024 doubled subsets.
    for m in range(1<<10): check_tournament(tournament(5,m));tested+=1
    rows=[]
    for l in range(1,5):
        a=alternating(l); n=len(a); c,t,_=profiles(a)
        expected=Counter()
        for j in range(l):
            for shift,factor in enumerate((1,l+2,comb(l+2,2))):
                expected[j+shift]+=comb(l-1,j)*factor
        require(c==expected,'alternating closed polynomial')
        d=doubling(a); dc,dt,maximizers=profiles(d)
        require(t==l+1 and dt==2*l+1,'alternating maximum')
        if l>=2:
            designed=((1<<l)-1)|(((1<<(l+1))-1)<<(n+l))
            require(maximizers==[designed],'unique alternating maximum')
        rows.append({'l':l,'vertices':n,'trans_base':t,'trans_double':dt,'maximizers':len(maximizers)})
    # All 2^(2l+1) copy assignments through l=8, no doubled-subset search.
    for l in range(1,9):
        a=alternating(l);n=len(a);d=doubling(a);good=[]
        for p in range(1<<n):
            s=(((1<<n)-1)^p)|(p<<n)
            if transitive(d,s):good.append(p)
        if l>=2:require(good==[((1<<(l+1))-1)<<l],('full-support assignment',l))
        else:require(len(good)==3,'triangle assignments')
    # Hostile control: transitive bases DO attain n+1, so n bound requires nontransitivity.
    for n in range(1,6):
        require(check_tournament(tournament(n,(1<<(n*(n-1)//2))-1))==n+1,'transitive control')
    print(json.dumps({'all_labeled_tournaments_n1_to5':tested,'alternating_rows':rows,'full_support_assignment_l_max':8,'optimized_mode_safe':True,'copy_3cnf_checked':True},indent=2))
if __name__=='__main__':main()
