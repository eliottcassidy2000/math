#!/usr/bin/env python3
"""Standalone exact certificate for split-interface module budgets. Standard library only."""
from itertools import combinations, permutations, product
from math import comb, factorial, gcd
from pathlib import Path
import hashlib, json
BASE=Path(__file__).resolve().parent
C3=[2,4,1]
PAIR41=[58,124,121,16,96,72,9]
SINGLE5=[2,12,9,1]
WREATH9=[58,60,57,464,480,456,135,263,71]

def need(ok,msg):
    if not ok: raise ValueError(msg)
def bits(mask):
    while mask:
        b=mask&-mask; mask-=b; yield b.bit_length()-1

def validate(T):
    n=len(T);full=(1<<n)-1
    for i,row in enumerate(T):
        need(row>=0 and row<=full and not(row>>i&1),'row bounds')
    for i,j in combinations(range(n),2):need((T[i]>>j&1)+(T[j]>>i&1)==1,'tournament pair')

def extend(T,mask):
    n=len(T);need(0<=mask<1<<n,'extension mask')
    return [row|(0 if mask>>i&1 else 1<<n)for i,row in enumerate(T)]+[mask]

def join(parts):
    n=sum(map(len,parts));T=[];offset=0
    for A in parts:
        later=((1<<n)-1)^((1<<(offset+len(A)))-1)
        T.extend((row<<offset)|later for row in A);offset+=len(A)
    return T

def substitution(Q,parts):
    labels=[(a,i)for a,A in enumerate(parts)for i in range(len(A))];T=[0]*len(labels)
    for x,y in combinations(range(len(labels)),2):
        a,i=labels[x];b,j=labels[y]
        yes=(parts[a][i]>>j&1)if a==b else(Q[a]>>b&1)
        u,v=(x,y)if yes else(y,x);T[u]|=1<<v
    return T

def path_counts(T):
    """Vertex-subset endpoint dynamic programming."""
    n=len(T);N=1<<n;E=[[0]*n for _ in range(N)]
    for v in range(n):E[1<<v][v]=1
    for mask in range(1,N):
        for v,c in enumerate(E[mask]):
            if c:
                for w in bits(T[v]&~mask):E[mask|1<<w][w]+=c
    return [1]+[sum(E[M])for M in range(1,N)]

def literal_paths(T):
    """Independent enumeration of every literal nonempty simple directed word."""
    H=[0]*(1<<len(T));H[0]=1
    def walk(mask,v):
        H[mask]+=1
        for w in bits(T[v]&~mask):walk(mask|1<<w,w)
    for v in range(len(T)):walk(1<<v,v)
    return H

def ordered_profile(T):
    """Unordered covers: isolate the path containing the least remaining name."""
    H=path_counts(T);N=len(H);n=len(T);pc=[[0]*(n+1)for _ in range(N)];pc[0][0]=1
    for mask in range(1,N):
        bit=mask&-mask;rest=mask^bit;sub=rest
        while True:
            A=sub|bit;B=mask^A
            for j in range(1,mask.bit_count()+1):pc[mask][j]+=H[A]*pc[B][j-1]
            if sub==0:break
            sub=(sub-1)&rest
    return [factorial(j)*pc[-1][j]for j in range(n+1)]

def literal_profile(T):
    """Independent permutation-and-cuts profile; only used through seven vertices."""
    n=len(T);F=[0]*(n+1)
    for p in permutations(range(n)):
        bad=sum(not(T[a]>>b&1)for a,b in zip(p,p[1:]));good=n-1-bad
        for j in range(bad+1,n+1):F[j]+=comb(good,j-1-bad)
    return F

def bins(F,d):return sum(comb(d,j)*F[j]for j in range(1,min(d,len(F)-1)+1))
def hp2(T,H=None):
    if H is None:H=path_counts(T)
    full=len(H)-1
    return H[-1],sum(H[A]*H[full^A]for A in range(full+1))
def content(F):return gcd(*F[1:])
def strong(T):
    n=len(T);full=(1<<n)-1
    for root in range(n):
        M=1<<root
        while True:
            nxt=M
            for v in bits(M):nxt|=T[v]
            if nxt==M:break
            M=nxt
        if M!=full:return False
    return True

def modules(T,blocks):
    allnames=set()
    for block in blocks:
        need(len(block)==3 and allnames.isdisjoint(block),'disjoint triples');allnames.update(block)
        need(sorted(sum(T[v]>>w&1 for w in block)for v in block)==[1,1,1],'C3 factor')
        for w in set(range(len(T)))-set(block):
            need(len({T[v]>>w&1 for v in block})==1,'old module homogeneous')

def omega(a,b):return a[0]*b[1]-a[1]*b[0]+a[2]*b[3]-a[3]*b[2]
def word_points(word):
    K=3*len(word)+1
    return [(1,0,K,0),(1,0,0,K)]+[(1,3*i,-1 if c&2 else 1,1 if c&1 else -1)for i,c in enumerate(word,1)]
def graph(points):
    T=[0]*len(points)
    for i,j in combinations(range(len(points)),2):
        w=omega(points[i],points[j]);need(w!=0,'zero pairing')
        a,b=(i,j)if w>0 else(j,i);T[a]|=1<<b
    return T

def main():
    report={'status':'PASS: exact finite certificates; arbitrary-size statements are proved in PROOF.md','counts':{}}
    need(join([C3,C3])==[58,60,57,16,32,8],'old two-block names')
    need(extend(join([C3,C3]),9)==PAIR41,'one-vertex H41 splitter')
    need(extend(C3,1)==SINGLE5,'one-vertex H5 splitter')
    seeds={}
    for name,T in [('C3',C3),('pair41',PAIR41),('single5',SINGLE5)]:
        validate(T);H=path_counts(T);need(H==literal_paths(T),'independent literal subset paths')
        F=ordered_profile(T);need(F==literal_profile(T),'literal permutation-cut seed profile')
        h,p=hp2(T,H);need((bins(F,1),bins(F,2))==(h,p),'seed profile evaluations')
        seeds[name]={'rows':T,'H':h,'P2':p,'ordered_profile':F,'content':content(F),'strong':strong(T)}
    need((seeds['C3']['H'],seeds['C3']['P2'],seeds['C3']['content'])==(3,12,3),'C3 values')
    need((seeds['pair41']['H'],seeds['pair41']['P2'],seeds['pair41']['content'])==(41,712,1),'H41 splitter values')
    need((seeds['single5']['H'],seeds['single5']['P2'],seeds['single5']['content'])==(5,32,1),'H5 splitter values')
    report['seed_certificates']=seeds
    # Exact affine rank-four geometry for all block types in the sharp family.
    geometries=[]
    for name,T,word,names in [('C3',C3,[2],[0,1,2]),('pair41',PAIR41,[1,0,3,1,0],[0,3,1,2,4,5,6]),('single5',SINGLE5,[2,2],[0,1,2,3])]:
        points=word_points(word);G=graph(points)
        permuted=[sum((T[u]>>v&1)<<j for j,v in enumerate(names))for u in names]
        need(G==permuted,'rank-four seed geometry')
        geometries.append({'name':name,'order_of_names':names,'coordinates':points})
    report['rank_four_seed_geometries']=geometries
    # Every sharp partial-budget family through four old triangle modules.
    sharp=[]
    for m in range(1,5):
        for k in range(m//2+1):
            r=m-2*k;parts=[PAIR41]*k+[C3]*r;T=join(parts);validate(T)
            old=[];offset=0;blocks=[]
            for A in parts:
                keep=range(6)if len(A)==7 else range(3)
                old.extend(offset+v for v in keep)
                blocks.extend([list(range(offset,offset+3))]+([list(range(offset+3,offset+6))]if len(A)==7 else []))
                offset+=len(A)
            oldgraph=[sum((T[u]>>v&1)<<j for j,v in enumerate(old))for u in old]
            need(oldgraph==join([C3]*m),'sharp physical old incidences')
            modules(oldgraph,[list(range(3*j,3*j+3))for j in range(m)])
            h,p=hp2(T);need((h,p)==(41**k*3**r,712**k*12**r),'sharp product H and P2')
            need(gcd(h,p)==3**r,'sharp upper gcd matches orbit lower bound')
            if len(T)<=10:need(content(ordered_profile(T))==3**r,'sharp full profile content')
            sharp.append({'old_modules':m,'births':k,'order':len(T),'H':h,'P2':p,'content':3**r})
    report['sharp_partial_budget_family']=sharp
    # Exact minimum number of births to content one for the ordered-join base.
    primitive=[]
    for m in range(1,5):
        a,r=divmod(m,2);T=join([PAIR41]*a+[SINGLE5]*r);h,p=hp2(T)
        need((h,p)==(41**a*5**r,712**a*32**r),'complete primitive repair values')
        need(gcd(h,p)==1,'complete primitive repair content')
        primitive.append({'old_modules':m,'births':a+r,'H':h,'P2':p})
    report['minimum_primitive_completion_family']=primitive
    # The three-module strong rank-four wreath. Every possible one-vertex extension.
    need(substitution(C3,[C3]*3)==WREATH9,'literal wreath graph');modules(WREATH9,[list(range(3*j,3*j+3))for j in range(3)])
    ext_hist={};full_checks=0
    masks_for_full={0,1,3,7,9,63,73,146,219,255,292,365,438,504,510,511}
    for mask in range(512):
        T=extend(WREATH9,mask);H=path_counts(T);h,p=hp2(T,H)
        need(h%3==0 and p%3==0,'all one-birth orbit-budget H and P2')
        ext_hist[h]=ext_hist.get(h,0)+1
        if mask in masks_for_full:
            F=ordered_profile(T);need(all(x%3==0 for x in F),'one-birth entire ordered profile')
            need(H==literal_paths(T),'one-birth independent literal word enumeration');full_checks+=1
    report['one_birth_wreath_census']={'attachments':512,'all_H_and_P2_divisible_by':3,'full_profile_controls':full_checks,'H_histogram':sorted(ext_hist.items())}
    # Two actual births reach content one, and one cannot by the proved bound.
    route=[];T=WREATH9
    for mask in (None,73,3):
        if mask is not None:T=extend(T,mask)
        validate(T);need(strong(T),'strong primitive-surgery route')
        H=path_counts(T);need(H==literal_paths(T),'literal route path table');F=ordered_profile(T)
        h,p=hp2(T,H);need((h,p)==(bins(F,1),bins(F,2)),'route profile evaluations')
        route.append({'new_out_mask':mask,'rows':T,'H':h,'P2':p,'content':content(F),'ordered_profile':F})
    need([(z['H'],z['P2'],z['content'])for z in route]==[(3159,48924,81),(13869,246600,9),(51703,1106588,1)],'exact sharp two-birth surgery')
    report['strong_wreath_two_birth_surgery']=route
    # Small exhaustive module quotients: full profiles on all labelled quotients through order three.
    quotient_checks=0
    for m in range(1,4):
        pairs=list(combinations(range(m),2))
        for code in range(1<<len(pairs)):
            Q=[0]*m
            for b,(i,j)in enumerate(pairs):
                a,z=(i,j)if code>>b&1 else(j,i);Q[a]|=1<<z
            T=substitution(Q,[C3]*m);F=ordered_profile(T)
            need(all(x%(3**m)==0 for x in F),'all small quotient full-profile module divisibility')
            quotient_checks+=1
    report['counts'].update({'literal_seed_profiles':3,'small_quotient_full_profiles':quotient_checks,'one_birth_wreath_H_P2':512,'one_birth_wreath_full_profiles':full_checks,'sharp_partial_budget_examples':len(sharp),'primitive_completion_examples':len(primitive),'strong_two_birth_route_profiles':3,'exact_geometric_seed_pairings':3+21+6})
    report['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=json.dumps(report,indent=2,sort_keys=True)+'\n';(BASE/'receipt.json').write_text(out);print(out,end='')
if __name__=='__main__':main()
