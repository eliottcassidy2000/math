#!/usr/bin/env python3
"""Exact certificate for Session166's shared Moser-spindle construction.
Standard library only. All checks remain active with python -O.
"""
from itertools import combinations, permutations, product
from fractions import Fraction as Q
from math import comb, factorial
from pathlib import Path
import json

BASE = Path(__file__).resolve().parent
MASKS = [18,22,26,181,81,146,145]
EDGES = [(0,1),(0,2),(1,2),(1,3),(2,3),(0,4),(0,5),(4,5),(4,6),(5,6),(3,6)]
EDGES = set(map(lambda e: tuple(sorted(e)), EDGES))

def require(b, message):
    if not b:
        raise RuntimeError(message)

def sign(p):
    return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))

def matrix(n, mask):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if mask>>k&1:a[i][j]=1
        else:a[j][i]=1
    return a

def encode(a):
    return sum(a[i][j]<<k for k,(i,j) in enumerate(combinations(range(len(a)),2)))

def relabel(a,p):
    n=len(a);b=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):b[p[i]][p[j]]=a[i][j]
    return b

def sccs(a):
    n=len(a);reach=[]
    for v in range(n):
        seen={v};todo=[v]
        while todo:
            u=todo.pop()
            for w in range(n):
                if a[u][w] and w not in seen:seen.add(w);todo.append(w)
        reach.append(seen)
    left=set(range(n));blocks=[]
    while left:
        v=min(left);c={u for u in left if u in reach[v] and v in reach[u]}
        blocks.append(c);left-=c
    # Outgoing component count gives the unique condensation order.
    blocks.sort(key=lambda c:-sum(a[next(iter(c))][next(iter(d))] for d in blocks if d!=c))
    return blocks

def submatrix(a,vs):
    return [[a[i][j]for j in vs]for i in vs]

def ordinal(blocks):
    n=sum(map(len,blocks));a=[[0]*n for _ in range(n)];s=0
    for b in blocks:
        for i in range(len(b)):
            for j in range(len(b)):a[s+i][s+j]=b[i][j]
            for j in range(s+len(b),n):a[s+i][j]=1
        s+=len(b)
    return a

def cyclic(a,b,c):
    out=ordinal([a,b,c]);na,nb,nc=len(a),len(b),len(c)
    for i in range(na):
        for j in range(na+nb,na+nb+nc):out[i][j]=0;out[j][i]=1
    return out

def profile(a):
    n=len(a);counts=[0]*(n+1);alt=0
    for p in permutations(range(n)):
        bad=sum(not a[p[i]][p[i+1]] for i in range(n-1))
        for k in range(bad+1,n+1):counts[k]+=comb(n-1-bad,k-1-bad)
        if bad==0:alt+=sign(p)
    for k in range(n+1):
        require(counts[k]%factorial(k)==0,'ordered-cover divisibility')
        counts[k]//=factorial(k)
    return counts,alt

def profile_subset(a):
    # Independent subset/first-block DP, not the permutation-and-cut method.
    n=len(a);hs=[0]*(1<<n)
    for mask in range(1,1<<n):
        vs=[i for i in range(n)if mask>>i&1]
        hs[mask]=sum(all(a[p[i]][p[i+1]]for i in range(len(p)-1))for p in permutations(vs))
    dp=[[0]*(n+1)for _ in range(1<<n)];dp[0][0]=1
    for mask in range(1,1<<n):
        first=mask&-mask;sub=mask
        while sub:
            if sub&first:
                for k in range(1,n+1):dp[mask][k]+=hs[sub]*dp[mask^sub][k-1]
            sub=(sub-1)&mask
    return dp[-1]

# Q(sqrt(3),sqrt(11)), basis [1,sqrt(3),sqrt(11),sqrt(33)].
def field(vals):return tuple(map(Q,vals))
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    z=[Q(0)]*4
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            common=i&j;coef=(3 if common&1 else 1)*(11 if common&2 else 1)
            z[i^j]+=coef*x*y
    return tuple(z)
def scalar(a,q):return tuple(x*q for x in a)
ZERO=field([0,0,0,0]);ONE=field([1,0,0,0])
A=(field([0,Q(1,2),0,0]),field([Q(1,2),0,0,0]))
B=(A[0],neg(A[1]));C=(field([0,1,0,0]),ZERO)
co=field([Q(5,6),0,0,0]);si=field([0,0,Q(1,6),0])
def rotate(p):return sub(mul(co,p[0]),mul(si,p[1])),add(mul(si,p[0]),mul(co,p[1]))
POINTS=[(ZERO,ZERO),A,B,C,rotate(A),rotate(B),rotate(C)]
def distance2(p,q):
    x,y=sub(p[0],q[0]),sub(p[1],q[1]);return add(mul(x,x),mul(y,y))

def main():
    n=6;ps=list(permutations(range(n)));pairs=list(combinations(range(n),2))
    mats=[matrix(n,t)for t in MASKS]
    require(all(len(sccs(a))==1 for a in mats),'seed strong connectivity')
    relabs=[[encode(relabel(a,p))for p in ps]for a in mats]
    dist=[[0]*7 for _ in range(7)];witnesses=[]
    for i,j in combinations(range(7),2):
        ds=[(t^MASKS[j]).bit_count()for t in relabs[i]];d=min(ds)
        dist[i][j]=dist[j][i]=d
        require(d>0,'seed isomorphism collision')
        require((d==1)==((i,j)in EDGES),'seed induced metagraph mismatch')
        if d==1:
            idx=ds.index(1);delta=relabs[i][idx]^MASKS[j]
            witnesses.append({'source':i,'target':j,'permutation':list(ps[idx]),'relabelled_mask':relabs[i][idx],'flip_pair':list(pairs[delta.bit_length()-1])})
    ds2={}
    for i,j in combinations(range(7),2):
        d=distance2(POINTS[i],POINTS[j]);require(d!=ZERO,'point collision')
        require((d==ONE)==((i,j)in EDGES),'exact geometric mismatch')
        ds2[f'{i},{j}']=[str(x)for x in d]
    # Independent exhaustive coloring check. The text also gives a short proof.
    colorings={}
    for k in (3,4):
        found=None
        for c in product(range(k),repeat=7):
            if all(c[i]!=c[j]for i,j in EDGES):found=c;break
        colorings[k]=found
    require(colorings[3]is None and colorings[4]is not None,'spindle chromatic number')
    prof=[];alts=[]
    for a in mats:
        p,alt=profile(a);require(p==profile_subset(a),'independent profile disagreement')
        require(alt%2==1,'alternating enumerator oddness');prof.append(p);alts.append(alt)
    require(prof[0]==prof[1]==[0,15,93,138,73,15,1],'blind-edge profile')
    require(alts==[-1,1,-5,3,9,-11,1],'alternating enumerator regression')
    # Check modular substitution under both cyclic and ordinal templates (order8).
    pt0=profile(cyclic([[0]],mats[0],[[0]]))[0]
    pt1=profile(cyclic([[0]],mats[1],[[0]]))[0]
    require(pt0==pt1,'cyclic module substitution')
    require(profile(ordinal([[[0]],mats[0],[[0]]]))[0]==profile(ordinal([[[0]],mats[1],[[0]]]))[0],'ordinal module substitution')
    # Literal induced Q_3 check via SCC decomposition and exact six-vertex isomorphism.
    orbit_index={t:i for i in (0,1)for t in relabs[i]}
    def family_word(a):
        bs=sccs(a)
        if len(bs)!=3 or any(len(b)!=6 for b in bs):return None
        word=[]
        for b in bs:
            t=encode(submatrix(a,sorted(b)))
            if t not in orbit_index:return None
            word.append(orbit_index[t])
        return tuple(word)
    flip_checks=0
    for word in product((0,1),repeat=3):
        a=ordinal([mats[i]for i in word]);require(family_word(a)==word,'cube word injectivity')
        neighbors=set()
        for u,v in combinations(range(18),2):
            a[u][v]^=1;a[v][u]^=1;image=family_word(a)
            a[u][v]^=1;a[v][u]^=1;flip_checks+=1
            if image is not None and image!=word:neighbors.add(image)
        expected={tuple(1-x if i==j else x for i,x in enumerate(word))for j in range(3)}
        require(neighbors==expected,'literal induced cube adjacency')
    # Score-separated wrapper guards, including every one-flip position.
    nn=6;q=2*nn+5
    guard=[[int(i!=j and 1<=(j-i)%q<=(q-1)//2)for j in range(q)]for i in range(q)]
    w=cyclic(guard,mats[0],guard);N=len(w);lo,hi=2*nn+2,3*nn+7
    require(N==5*nn+10 and len(sccs(w))==1,'wrapper order/strongness')
    require([sum(row)for row in w][:q]==[lo]*q,'wrapper low scores')
    require([sum(row)for row in w][q+nn:]==[hi]*q,'wrapper high scores')
    wrapper_checks=0
    for u,v in combinations(range(N),2):
        w[u][v]^=1;w[v][u]^=1;scores=[sum(row)for row in w]
        if not(q<=u<q+nn and q<=v<q+nn):
            require(scores.count(lo)!=q or scores.count(hi)!=q,'guard flip escaped score detection')
        w[u][v]^=1;w[v][u]^=1;wrapper_checks+=1
    # Parity action over all labelled tournaments of order5, all vertex permutations.
    n=5;pairs5=list(combinations(range(n),2));pos={e:i for i,e in enumerate(pairs5)};parity_checks=0
    for p in permutations(range(n)):
        flip=sum(p[i]>p[j]for i,j in pairs5)%2
        require((-1)**flip==sign(p),'relabelling sign formula')
        # Parity of a coordinate permutation is irrelevant; each reversed ordered pair toggles once.
        for t in range(1<<len(pairs5)):
            out=0
            for k,(i,j)in enumerate(pairs5):
                x,y=p[i],p[j];bit=(t>>k&1)^(x>y)
                out|=bit<<pos[tuple(sorted((x,y)))]
            require((-1)**out.bit_count()==sign(p)*(-1)**t.bit_count(),'cube parity equivariance')
            if out==t:require(sign(p)==1,'odd tournament automorphism')
            parity_checks+=1
    certificate={'mask_convention':'lexicographic unordered pairs, bit1 means i->j for i<j','names':['O','A','B','C','D','E','F'],'masks':MASKS,'edges':sorted(EDGES),'minimum_quotient_hamming_distances':dist,'edge_witnesses':witnesses,'proper_four_coloring':colorings[4],'profiles':prof,'alternating_hamilton_enumerators':alts,'squared_distances_basis_1_sqrt3_sqrt11_sqrt33':ds2,'cyclic_singleton_module_profile':pt0}
    (BASE/'seed_certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    report={'status':'PASS','seed_tournaments':7,'seed_pairs':21,'permutations_per_pair':720,'exact_geometric_pairs':21,'profile_algorithms':2,'literal_Q3_flip_checks':flip_checks,'wrapper_flip_checks':wrapper_checks,'parity_equivariance_checks':parity_checks,'spindle_chromatic_number':4,'external_arbitrary_girth_theorem':'Tomon arXiv:2610.03517v2 Theorem 1.5 / 5.2; not computationally established here'}
    (BASE/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
