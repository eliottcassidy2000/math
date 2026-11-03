#!/usr/bin/env python3
"""Exact S167 divisor-lattice / aliquot-tree controls, stdlib only.
Seed and ordinal-product mechanism credited to S166 cartesian_flip/NOTE.md.
No repository scanning; this script reads only itself and writes certificate.json.
"""
from collections import deque
from hashlib import sha256
from itertools import combinations, permutations, product
from math import comb, factorial, isqrt
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def matrix(n,t):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j]=(t>>k)&1; a[j][i]=1-a[i][j]
    return a

def encode(a):
    return sum(a[i][j]<<k for k,(i,j) in enumerate(combinations(range(len(a)),2)))

def relabel(a,p):
    n=len(a); b=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n): b[p[i]][p[j]]=a[i][j]
    return b

def profile_cuts(a):
    n=len(a); hist=[0]*n
    for p in permutations(range(n)):
        hist[sum(not a[p[i]][p[i+1]] for i in range(n-1))]+=1
    return [0]+[sum(hist[b]*comb(n-1-b,k-1-b) for b in range(k))//factorial(k) for k in range(1,n+1)]

def hamilton_subsets(a):
    n=len(a);dp=[[0]*n for _ in range(1<<n)];hs=[0]*(1<<n)
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not mask>>v&1:continue
            if mask!=(1<<v):dp[mask][v]=sum(dp[mask^(1<<v)][u] for u in range(n) if mask>>u&1 and a[u][v])
        hs[mask]=sum(dp[mask])
    return hs

def profile_partitions(a):
    n=len(a);hs=hamilton_subsets(a);dp=[[0]*(n+1)for _ in range(1<<n)];dp[0][0]=1
    for mask in range(1,1<<n):
        first=mask&-mask;sub=mask
        while sub:
            if sub&first:
                for k in range(1,n+1):dp[mask][k]+=hs[sub]*dp[mask^sub][k-1]
            sub=(sub-1)&mask
    return dp[-1]

def sccs(a):
    n=len(a);reach=[]
    for v in range(n):
        seen={v};todo=[v]
        for u in todo:
            for w in range(n):
                if a[u][w] and w not in seen:seen.add(w);todo.append(w)
        reach.append(seen)
    left=set(range(n));blocks=[]
    while left:
        v=min(left);c={u for u in left if u in reach[v] and v in reach[u]};blocks.append(c);left-=c
    blocks.sort(key=lambda c:-sum(a[next(iter(c))][next(iter(d))] for d in blocks if c!=d))
    return blocks

def ordinal(blocks):
    n=sum(map(len,blocks));a=[[0]*n for _ in range(n)];s=0
    for b in blocks:
        for i in range(len(b)):
            for j in range(len(b)):a[s+i][s+j]=b[i][j]
            for j in range(s+len(b),n):a[s+i][j]=1
        s+=len(b)
    return a

def factor(n):
    out=[];p=2
    while p*p<=n:
        if n%p==0:
            a=0
            while n%p==0:n//=p;a+=1
            out.append((p,a))
        p+=1
    if n>1:out.append((n,1))
    return out

def aliquot_trial(n):
    if n<2:return 0
    total=1
    for d in range(2,isqrt(n)+1):
        if n%d==0:total+=d+(n//d if d*d!=n else 0)
    return total

def aliquot_product(n):
    if n<2:return 0
    sigma=1
    for p,a in factor(n):sigma*=sum(p**i for i in range(a+1))
    return sigma-n

def divisor_states(n):
    fac=factor(n);states={}
    for exps in product(*(range(a+1)for p,a in fac)):
        d=1;bits=[]
        for (p,a),b in zip(fac,exps):d*=p**b;bits.extend([1]*b+[0]*(a-b))
        states[d]=(exps,tuple(bits))
    return fac,states

def cycle_codes(length):
    if length==3:return [(0,),(1,),(2,)]
    require(length>=4 and length!=5,'unsupported cycle length')
    if length%2==0:
        k=length//2;word=[0]*k;out=[tuple(word)]
        for i in list(range(k))+list(range(k-1)):
            word[i]^=1;out.append(tuple(word))
        return out
    out=[(0,0,0),(0,1,0),(1,1,0),(1,1,1),(1,0,1),(2,0,1),(2,0,0)]
    while len(out)<length:
        # Replace old consecutive u,v,w by u0,u1,v1,w1,w0.
        u,v,w=out[:3]
        out=[u+(0,),u+(1,),v+(1,),w+(1,),w+(0,)]+[z+(0,)for z in out[3:]]
    return out

def main():
    masks=[18,22];mats=[matrix(6,t)for t in masks];ps=list(permutations(range(6)))
    orbits=[{encode(relabel(a,p))for p in ps}for a in mats]
    require(not(orbits[0]&orbits[1]),'seed classes collide')
    require(all(len(sccs(a))==1 for a in mats),'seed not strong')
    distance=min((t^masks[1]).bit_count()for t in orbits[0]);require(distance==1,'seed not adjacent')
    profiles=[profile_cuts(a)for a in mats]
    require(profiles[0]==profiles[1]==[0,15,93,138,73,15,1],'seed profile')
    require(all(profile_partitions(a)==profiles[i]for i,a in enumerate(mats)),'independent cover profile')
    two_block_counts={}
    for w in product((0,1),repeat=2):
        h=hamilton_subsets(ordinal([mats[x]for x in w]))[-1]
        require(h==225,'literal ordinal Hamilton count');two_block_counts[''.join(map(str,w))]=h
    # Directly examine every arc flip for each of the 8 literal Q3 representatives.
    orbit_index={t:i for i,o in enumerate(orbits)for t in o}
    def family_word(a):
        bs=sccs(a)
        if len(bs)!=3 or any(len(b)!=6 for b in bs):return None
        word=[]
        for b in bs:
            vs=sorted(b);t=encode([[a[i][j]for j in vs]for i in vs])
            if t not in orbit_index:return None
            word.append(orbit_index[t])
        return tuple(word)
    literal_flips=0
    for word in product((0,1),repeat=3):
        a=ordinal([mats[i]for i in word]);neighbors=set()
        require(family_word(a)==word,'family decode')
        for u,v in combinations(range(18),2):
            a[u][v]^=1;a[v][u]^=1;out=family_word(a);a[u][v]^=1;a[v][u]^=1;literal_flips+=1
            if out is not None and out!=word:neighbors.add(out)
        expected={tuple(1-x if i==j else x for i,x in enumerate(word))for j in range(3)}
        require(neighbors==expected,'literal cube has wrong adjacency')
    # Independently rediscovered same-session ternary seed, also in the fences packet.
    triangle=[153,177,181];triangle_mats=[matrix(6,t)for t in triangle]
    triangle_profile=[0,29,131,162,77,15,1]
    require(all(len(sccs(a))==1 for a in triangle_mats),'triangle strongness')
    require(all(profile_cuts(a)==profile_partitions(a)==triangle_profile for a in triangle_mats),'triangle profiles')
    triangle_distances={}
    for i,j in combinations(range(3),2):
        dist=min((encode(relabel(triangle_mats[i],p))^triangle[j]).bit_count()for p in ps)
        require(dist==1,'triangle quotient distance');triangle_distances[f'{triangle[i]},{triangle[j]}']=dist
    cycle_records={}
    for length in range(3,32):
        if length==5:continue
        words=cycle_codes(length)
        require(len(set(words))==length,'cycle collision')
        for i,j in combinations(range(length),2):
            h=sum(x!=y for x,y in zip(words[i],words[j]))
            require((h==1)==(j==i+1 or (i==0 and j==length-1)),'cycle has chord or missing edge')
        cycle_records[str(length)]={'dimension':len(words[0]),'words':[''.join(map(str,w))for w in words]}
    prism=list(product(range(3),range(2)));prism_five_checks=0
    for order in permutations(prism,5):
        prism_five_checks+=1
        induced=all((sum(x!=y for x,y in zip(order[i],order[j]))==1)==(j==i+1 or(i==0 and j==4))for i,j in combinations(range(5),2))
        require(not induced,'unexpected induced C5 in triangular prism')
    div_records={};div_pairs=0
    for n in [8,20,21,60,290,1729]:
        fac,states=divisor_states(n);ds=sorted(states);omega=sum(a for p,a in fac)
        require(ds==[d for d in range(1,n+1)if n%d==0],'divisor enumeration')
        for d,e in combinations(ds,2):
            hamming=sum(x!=y for x,y in zip(states[d][1],states[e][1]))
            l1=sum(abs(x-y)for x,y in zip(states[d][0],states[e][0]))
            require(hamming==l1,'grid isometry')
            require((hamming==1)==(e%d==0 and len(factor(e//d))==1 and factor(e//d)[0][1]==1),'Hasse adjacency')
            div_pairs+=1
        nontrivial=[d for d in ds if 1<d<n]
        squarefree=[d for d in nontrivial if all(b<=1 for b in states[d][0])]
        primes=[d for d in nontrivial if sum(states[d][0])==1]
        f,s,u=len(nontrivial),len(squarefree),len(primes)
        sigma=sum(ds)
        require(sigma-n==aliquot_trial(n)==aliquot_product(n),'weighted divisor evaluation')
        div_records[str(n)]={'factorization':fac,'Omega':omega,'vertices':len(ds),'F':f,'S':s,'U':u,'defect':f-s-u,'sigma':sigma,'aliquot':sigma-n,'tournament_order':6*omega,'common_H':15**omega,'divisor_addresses':{str(d):''.join(map(str,states[d][1]))for d in ds}}
    starts=[20,21,290,1729];sequences={};edges=set();vertices={0}
    for start in starts:
        seq=[start]
        while seq[-1]!=0:
            nxt=aliquot_trial(seq[-1]);require(nxt==aliquot_product(seq[-1]),'aliquot independent check')
            require(nxt not in seq,'unexpected cycle in specified starts')
            edges.add(tuple(sorted((seq[-1],nxt))));vertices.add(seq[-1]);vertices.add(nxt);seq.append(nxt)
        sequences[str(start)]=seq
    edges=sorted(edges);adj={v:set()for v in vertices}
    for i,(u,v)in enumerate(edges):adj[u].add((v,i));adj[v].add((u,i))
    codes={0:0};parent={0:None};todo=[0]
    for u in todo:
        for v,i in sorted(adj[u]):
            if v not in codes:codes[v]=codes[u]^(1<<i);parent[v]=u;todo.append(v)
    require(len(codes)==len(vertices) and len(edges)==len(vertices)-1,'aliquot union is not a tree')
    pairs_checked=0
    for u in sorted(vertices):
        dist={u:0};q=deque([u])
        while q:
            v=q.popleft()
            for w,_ in adj[v]:
                if w not in dist:dist[w]=dist[v]+1;q.append(w)
        for v in sorted(vertices):
            require((codes[u]^codes[v]).bit_count()==dist[v],'tree isometry');pairs_checked+=1
        if u:require(parent[u]==aliquot_trial(u),'aliquot direction disagrees with root')
    k=len(edges)
    report={'status':'PASS','source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'seed_masks':masks,'seed_profile':profiles[0],'seed_permutations_per_class':len(ps),'literal_two_block_H':two_block_counts,'literal_Q3_arc_flips':literal_flips,'divisor_pair_checks':div_pairs,'triangle_masks':triangle,'triangle_profile':triangle_profile,'triangle_quotient_distances':triangle_distances,'cycle_embeddings':cycle_records,'triangular_prism_C5_ordered_checks':prism_five_checks,'divisors':div_records,'aliquot':{'sequences':sequences,'vertices':len(vertices),'edges':k,'root':0,'tree_isometry_ordered_pair_checks':pairs_checked,'tournament_order':6*k,'strong_wrapper_order':30*k+10,'common_unwrapped_H':15**k,'edge_coordinates':edges,'vertex_codes':{str(v):format(codes[v],f'0{k}b')for v in sorted(vertices)}}}
    (HERE/'certificate.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in report.items()if k not in ('divisors','aliquot','cycle_embeddings')},indent=2))
    print('Aliquot tree:',len(vertices),'vertices,',k,'edges; tournament orders',6*k,'and',30*k+10)
if __name__=='__main__':main()
