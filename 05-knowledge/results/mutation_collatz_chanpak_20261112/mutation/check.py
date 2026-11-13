#!/usr/bin/env python3
"""Small exact controls for short-return SCC classes and two-root forms."""
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]

def need(q, msg):
    if not q: raise RuntimeError(msg)

def matrix(n, mask):
    B=[[0]*n for _ in range(n)]; k=0
    for i in range(n):
        for j in range(i+1,n):
            B[i][j]=1 if mask>>k&1 else -1
            B[j][i]=-B[i][j]; k+=1
    return tuple(map(tuple,B))

def induced(B, vertices):
    return tuple(tuple(B[i][j] for j in vertices) for i in vertices)

def code(B):
    return sum(int(B[i][j]>0)<<k for k,(i,j) in
               enumerate((i,j) for i in range(len(B)) for j in range(i+1,len(B))))

@lru_cache(None)
def canonical(B):
    return min(code(induced(B,p)) for p in permutations(range(len(B))))

def tournament(B):
    return all(B[i][i]==0 for i in range(len(B))) and all(
        B[i][j] in (-1,1) and B[j][i]==-B[i][j]
        for i in range(len(B)) for j in range(i+1,len(B)))

def mutate(B,k):
    return tuple(tuple(-B[i][j] if i==k or j==k else
        B[i][j]+(abs(B[i][k])*B[k][j]+B[i][k]*abs(B[k][j]))//2
        for j in range(len(B))) for i in range(len(B)))

def run(B,word):
    for k in word: B=mutate(B,k)
    return B

@lru_cache(None)
def hp(B):
    n=len(B)
    return sum(all(B[p[i]][p[i+1]]>0 for i in range(n-1))
               for p in permutations(range(n)))

def scc_word(B):
    n=len(B)
    reach=[(1<<i)|sum(1<<j for j in range(n) if B[i][j]>0) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1: reach[i]|=reach[k]
    unseen=set(range(n)); blocks=[]
    while unseen:
        i=min(unseen)
        block=tuple(j for j in unseen if reach[i]>>j&1 and reach[j]>>i&1)
        blocks.append(block); unseen.difference_update(block)
    unsorted=tuple(blocks)
    blocks.sort(key=lambda S:-sum(B[S[0]][T[0]]>0 for T in unsorted if T!=S))
    return tuple((len(S),canonical(induced(B,S))) for S in blocks)

def signature(B):
    word=scc_word(B)
    large=[i for i,s in enumerate(word) if s[0]>=4]
    if large:
        i,j=large[0],large[-1]
        return ('rigid',word[i:j+1],word[j+1:]+word[:i])
    sizes=tuple(s[0] for s in word)
    return ('necklace',min((sizes[i:]+sizes[:i] for i in range(len(sizes))),default=()))

def ordinal(A,B):
    n,m=len(A),len(B)
    return tuple(tuple(A[i][j] if i<n and j<n else B[i-n][j-n] if i>=n and j>=n
                       else 1 if i<n else -1 for j in range(n+m)) for i in range(n+m))

def append_c(R,I):
    n=len(R); B=[list(row)+[1 if I>>i&1 else -1] for i,row in enumerate(R)]
    B.append([-B[i][n] for i in range(n)]+[0]); return tuple(map(tuple,B))

def birth(B,root):
    n=len(B); C=[list(row)+[1 if i==root else -1] for i,row in enumerate(B)]
    C.append([-C[i][n] for i in range(n)]+[0]); return tuple(map(tuple,C))

def u_cover(R,I):
    n=len(R); out=0
    for bits in range(1,1<<n):
        left=[i for i in range(n) if bits>>i&1]
        right=[i for i in range(n) if not bits>>i&1]
        left_count=sum(bool(I>>p[-1]&1) and all(R[p[i]][p[i+1]]>0 for i in range(len(p)-1))
                       for p in permutations(left))
        out+=left_count*hp(induced(R,right))
    return out

def main():
    shapes=((0,), (0,1), (0,1,0), (0,1,2),
            (0,1,0,1),(0,1,0,2),(0,1,2,0),(0,1,2,1),(0,1,2,3))
    summary=[]; checks=Counter(); root_normal_forms=0
    for n in range(6):
        reps={}
        for mask in range(1<<(n*(n-1)//2)):
            B=matrix(n,mask); reps.setdefault(canonical(B),B)
            roots=[v for v in range(n) if sum(B[u][v]>0 for u in range(n))==1]
            need(len(roots)<=3,'at most three indegree-one roots')
            if len(roots)>=2:
                if len(roots)==3:
                    need(all(B[r][v]>0 for r in roots for v in range(n) if v not in roots),
                         'three roots are an initial triangle')
                else:
                    a,b=roots if B[roots[0]][roots[1]]>0 else roots[::-1]
                    c=next(v for v in range(n) if B[v][a]>0)
                    rest=[v for v in range(n) if v not in (a,b,c)]
                    need(B[b][c]>0 and all(B[a][v]>0 and B[b][v]>0 for v in rest),
                         'two-root exposed triangle')
                root_normal_forms+=1
        parent={k:k for k in reps}
        def find(k):
            while parent[k]!=k:
                parent[k]=parent[parent[k]]; k=parent[k]
            return k
        for key,B in reps.items():
            for shape in shapes:
                names=max(shape)+1
                if names>n: continue
                for assignment in permutations(range(n),names):
                    word=tuple(assignment[i] for i in shape)
                    C=run(B,word); checks['matrix_word_cases']+=1
                    if not tournament(C): continue
                    checks['tournament_returns']+=1
                    need(hp(B)==hp(C),'H neutrality control')
                    need(signature(B)==signature(C),'SCC component invariant')
                    if len(scc_word(B))==1 and n>=4:
                        need(B==C,'strong pointwise frozen')
                    target=canonical(C); parent[find(key)]=find(target)
        actual={frozenset(k for k in reps if find(k)==leader) for leader in map(find,reps)}
        predicted=defaultdict(set)
        for k,B in reps.items(): predicted[signature(B)].add(k)
        need(actual==set(map(frozenset,predicted.values())),'complete small component partition')
        summary.append({'n':n,'ordinary_classes':len(reps),'short_components':len(actual)})
    need([s['ordinary_classes'] for s in summary]==[1,1,1,2,4,12],'ordinary controls')
    need([s['short_components'] for s in summary]==[1,1,1,2,3,9],'component controls')
    two_root=[]
    for n in range(4):
        for mask in range(1<<(n*(n-1)//2)):
            R=matrix(n,mask)
            for I in range(1<<n):
                W=append_c(R,I); T=birth(birth(W,n),n+1)
                U=u_cover(R,I)
                need(hp(T)==hp(W)+2*(hp(R)+U),'two-root count formula')
                roots=[v for v in range(n+3) if sum(T[u][v]>0 for u in range(n+3))==1]
                need(len(roots)==(2 if I else 3),'root number and cut')
                two_root.append({'R_order':n,'R_mask':mask,'I':I,'H_R':hp(R),
                                 'H_W':hp(W),'U_I':U,'H_T':hp(T)})
    rows=(2,4,9,3)
    C4=tuple(tuple(0 if i==j else 1 if rows[i]>>j&1 else -1 for j in range(4)) for i in range(4))
    C3=matrix(3,5)
    original=ordinal(C4,C3); result=run(original,(4,5,6,4))
    expected=tuple(tuple(original[i][j] if i<4 and j<4 else -original[i][j]
                         for j in range(7)) for i in range(7))
    need(result==expected and tournament(result),'literal extreme triangle transfer')
    need(tuple(s[0] for s in scc_word(original))==(4,3)
         and tuple(s[0] for s in scc_word(result))==(3,4),'independent SCC order control')
    need(hp(original)==hp(result)==15,'order-seven product count')
    need(signature(original)==signature(result),'order-seven component merge')
    buffered=ordinal(matrix(1,0),C4)
    u,v,w=1,0,3
    program=(u,v,u,v,w,u,w,u,v,w,v,w,v)
    running=buffered; checkpoints=[]
    for step,pivot in enumerate(program,1):
        running=mutate(running,pivot)
        if step in (4,8,12,13):
            need(tournament(running) and hp(running)==5,'buffer tournament checkpoint')
            checkpoints.append({'step':step,'matrix':running})
    permutation=list(range(5)); permutation[u],permutation[w]=permutation[w],permutation[u]
    need(running==induced(buffered,permutation),'buffer final transposition')
    q_before=hp(induced(buffered,[i for i in range(5) if i!=u]))
    q_after=hp(induced(running,[i for i in range(5) if i!=u]))
    need((q_before,q_after)==(3,1),'fixed-name deletion response changes')
    paths=[
        '05-knowledge/results/triples_ports_sieve_20261112/ports/FOUR_RETURNS.md',
        '05-knowledge/results/triples_ports_sieve_20261112/ports/PROOF.md',
        '05-knowledge/results/adaptive_metagraph_reformulation_20261103/ANCESTRY.md',
        '05-knowledge/results/two_entry_funnels_20261103/READY.md',
    ]
    print(json.dumps({'status':'PASS; FINITE-EXACT controls',
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        'short_return_universe':summary,'checks':dict(checks),
        'normal_forms_recovered':root_normal_forms,'two_root_universe':two_root,
        'literal_triangle_transfer':{'source':original,'target':result,'word':[4,5,6,4],'H':15},
        'buffer_transport':{'source':buffered,'word':program,'checkpoints':checkpoints,
                            'fixed_name':u,'deletion_before':q_before,'deletion_after':q_after},
        'scope':'Short-return generated quotient only, not all-length mutation classes'
    },indent=2,sort_keys=True))

if __name__=='__main__': main()
