#!/usr/bin/env python3
"""Independent Session125 origin/forest audit; explicit checks survive -O.

No producer implementation is imported. Determinants use rational elimination;
rooted forests of general graphs use parent functions (not edge-mask DSU).
The path bijection is reconstructed by covering-cell geometry and root markers.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, isqrt, prod
from pathlib import Path
import hashlib
import json
import re
import sys

CHECKS = 0
def ck(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)

def square(x):
    return x >= 0 and isqrt(x)**2 == x

def edges(n):
    # Arithmetic generator independent from arrival-root enumeration.
    return tuple((u, v) for u in range(n) for v in range(u+1, n)
                 if square(u+v+2))

def det(matrix):
    a = [list(map(Fraction, row)) for row in matrix]
    answer = Fraction(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            answer = -answer
        value = a[i][i]
        answer *= value
        for j in range(i+1, len(a)):
            scale = a[j][i] / value
            a[j] = [x-scale*y for x, y in zip(a[j], a[i])]
    return answer

def lap_identity(n, es):
    a = [[int(i == j) for j in range(n)] for i in range(n)]
    for u,v in es:
        a[u][u] += 1
        a[v][v] += 1
        a[u][v] -= 1
        a[v][u] -= 1
    return a

def cofactor(a, i, j):
    return (-1)**(i+j)*det([[x for c,x in enumerate(row) if c != j]
                          for r,row in enumerate(a) if r != i])

def rooted_parent_count(n, es):
    choices = [[v] for v in range(n)]
    for u,v in es:
        choices[u].append(v)
        choices[v].append(u)
    total = 0
    for f in product(*choices):
        good = True
        for v in range(n):
            seen = set()
            while f[v] != v:
                if v in seen:
                    good = False
                    break
                seen.add(v)
                v = f[v]
            if not good:
                break
        total += good
    return total

def tilings(n):
    if n == 0:
        yield ()
    for first in (1,2):
        if first <= n:
            for tail in tilings(n-first):
                yield (first,)+tail

def decode_geometry(ts):
    # Read a complete cell covering, without producer state or placement logic.
    owner = []
    for tile in ts:
        start = len(owner)
        owner.extend([start]*tile)
    L = (len(owner)-1)//2
    digits = []
    for right in range(2*L, 0, -2):
        left = right-1
        if owner[right] == left:
            digits.append(1)  # a domino contained in this cell pair
        elif owner[left] == left-1:
            # A crossing domino: 1 if right is already covered from right,
            # 2 if right itself is a monomer.
            digits.append(1 if right+1 < len(owner) and owner[right+1] == right else 2)
        else:
            digits.append(0)
    return tuple(reversed(digits))

def tiling_blocks(ts):
    # Odd-numbered monomers mark roots, even-numbered ones mark cuts.
    blocks, length, root, monomers = [], 0, None, 0
    for tile in ts:
        if tile == 2:
            length += 1
        else:
            monomers += 1
            if monomers % 2:
                length += 1
                root = length
            else:
                blocks.append((length,root))
                length,root = 0,None
    blocks.append((length,root))
    return tuple(blocks)

def compositions(n):
    if not n:
        yield ()
    for head in range(1,n+1):
        for tail in compositions(n-head):
            yield (head,)+tail

def fib(n):
    a,b = 0,1
    for _ in range(n):
        a,b=b,a+b
    return a

def two_factors(n):
    es = edges(n)
    inc = [set(i for i,e in enumerate(es) if v in e) for v in range(n)]
    solutions=[]
    def visit(chosen, undecided):
        while True:
            changed=False
            for v in range(n):
                have=len(inc[v]&chosen)
                free=inc[v]&undecided
                need=2-have
                if need < 0 or need > len(free):
                    return
                if free and (need == 0 or need == len(free)):
                    undecided=undecided-free
                    if need:
                        chosen=chosen|free
                    changed=True
            if not changed:
                break
        if not undecided:
            solutions.append(tuple(es[i] for i in sorted(chosen)))
            return
        candidates=[]
        for v in range(n):
            free=inc[v]&undecided
            if free:
                need=2-len(inc[v]&chosen)
                candidates.append((comb(len(free),need),v,need,free))
        _,v,need,free=min(candidates,key=lambda x:(x[0],x[1]))
        for selected in combinations(sorted(free),need):
            visit(chosen|set(selected),undecided-free)
    visit(set(),set(range(len(es))))
    ck(len(set(solutions)) == len(solutions), '2-factor duplicate')
    result=[]
    for es1 in solutions:
        adj=[set() for _ in range(n)]
        for u,v in es1:
            adj[u].add(v);adj[v].add(u)
        ck(all(len(row)==2 for row in adj),'2-factor degree')
        todo=set(range(n));components=[]
        while todo:
            stack=[min(todo)];component=set()
            while stack:
                v=stack.pop()
                if v in component:continue
                component.add(v);stack.extend(adj[v]-component)
            todo-=component;components.append(len(component))
        result.append({'components':sorted(components),'edges_one_based':[[u+1,v+1] for u,v in es1]})
    return sorted(result,key=lambda row:row['components'])

# Independent origin arithmetic, boundary and exact degree controls.
for N in range(1,4001):
    actual=[v for v in range(1,N) if square(N+v)]
    predicted=[k*k-N for k in range(isqrt(N)+1,isqrt(2*N-1)+1)]
    ck(actual == predicted,'arrival neighbor set')
    ck(len(actual) == isqrt(2*N-1)-isqrt(N),'arrival degree')
for n in range(3,31):
    v=1 if n<=7 else 8 if n<=16 else 16 if n==17 else 18
    ck(sum(square(v+u) for u in range(1,n+1) if u!=v)==1,'degree-one exclusion')
for base in (31,32):
    for x,y,z in product(range(12),repeat=3):
        a,b,c=base+x,base+y,base+z
        add=lambda a,b:a+b-base
        mul=lambda a,b:base+(a-base)*(b-base)
        ck(add(a,b)==base+x+y,'transport addition')
        ck(mul(a,b)==base+x*y,'transport multiplication')
        ck(add(add(a,b),c)==add(a,add(b,c)),'addition associativity')
        ck(mul(mul(a,b),c)==mul(a,mul(b,c)),'multiplication associativity')
        ck(mul(a,add(b,c))==add(mul(a,b),mul(a,c)),'distributivity')
        ck(add(a,base)==a and mul(a,base+1)==a and mul(a,base)==base,'identities')
S=lambda n:n==0 or n>=32
ck([n for n in range(64) if S(n) and not S(n-32)] == [0]+list(range(33,64)), 'Apery set')
ck([n for n in range(1,64) if S(n) and not any(S(a) and S(n-a) for a in range(1,n))] == list(range(32,64)), 'minimal generators')
for n in range(10001):
    if not S(n):continue
    q,r=divmod(n,32)
    base=0 if not r else 32+r
    copies=q if not r else q-1
    ck(copies>=0 and n==base+32*copies,'Apery canonical form')
ck(66==32+34==33+33,'factorization nonuniqueness')

factor_rows={n:two_factors(n) for n in (31,32,33)}
ck([r['components'] for r in factor_rows[31]]==[[11,20],[14,17]],'G31 all factors')
ck([r['components'] for r in factor_rows[32]]==[[3,6,23],[3,29],[5,27],[9,23],[32]],'G32 all factors')
ck(sum(r['components']==[33] for r in factor_rows[33])==1,'G33 unique cycle')
c32=set(map(tuple,next(r['edges_one_based'] for r in factor_rows[32] if r['components']==[32])))
c33=set(map(tuple,next(r['edges_one_based'] for r in factor_rows[33] if r['components']==[33])))
ck(c32-c33 == {(5,31),(16,20)},'32-33 deleted edges')
ck(c33-c32 == {(5,20),(16,33),(31,33)},'32-33 inserted edges')

# Three independently generated finite sets and the actual geometric code.
path_rows=[]
for L in range(10):
    legal={w for w in product(range(3),repeat=L) if not re.search('21*2',''.join(map(str,w)))}
    all_ts=list(tilings(2*L+1))
    records=[(decode_geometry(t),tiling_blocks(t),t) for t in all_ts]
    expected_forests={tuple(zip(comp,roots)) for comp in compositions(L+1)
                      for roots in product(*(range(1,k+1) for k in comp))}
    ck({w for w,_,_ in records}==legal,'geometric code covers legal language')
    ck(len(records)==len(legal)==len(expected_forests),'bijection cardinality')
    ck({f for _,f,_ in records}==expected_forests,'root-marked tiling covers forests')
    hist=Counter()
    lookup={w:(f,t) for w,f,t in records}
    values=[]
    for w in sorted(legal):
        f,t=lookup[w]
        nonzero=sum(d!=0 for d in w)
        ck(nonzero==t.count(2)==L+1-len(f),'statistic is nonzero digits')
        hist[nonzero]+=1
        values.append(sum(d*fib(2*i+2) for i,d in enumerate(reversed(w))))
        w0=decode_geometry((1,1)+t)
        ck(w0==(0,)+w,'leading-zero geometry')
        ck(tiling_blocks((1,1)+t)==((1,1),)+f,'leading-zero singleton')
        ck(decode_geometry(t+(1,1))==w+(0,), 'right-zero geometry')
        ck(tiling_blocks(t+(1,1))==f+((1,1),), 'right-singleton extension')
        if w:ck(bool(w[0])==(f[0][0]!=1), 'canonical first component')
        if w and w[0] and nonzero==1:
            ck(f==((2,w[0]),)+((1,1),)*(len(w)-1), 'canonical one-edge subclass')
        value=values[-1]
        shifted=sum(d*fib(2*i+4) for i,d in enumerate(reversed(w)))
        ck(shifted==((3*value+isqrt(5*value*value))//2+1 if value else 0), 'Beatty right-singleton rank')
    ck(values==list(range(fib(2*L+2))),'fixed-length numerical rank')
    expected=[comb(2*L+1-k,k) for k in range(L+1)]
    ck([hist[k] for k in range(L+1)]==expected,'binomial coefficients')
    # Transition matrix evaluated independently by direct multiplication.
    for t in (0,1,2,3,-1):
        free,locked=1,0
        for _ in range(L):free,locked=(1+t)*free+locked,t*(free+locked)
        ck(free+locked==sum(c*t**k for k,c in enumerate(expected)),'transfer polynomial')
        a=lap_identity(L+1,[(i,i+1) for i in range(L)])
        a=[[t*(x-int(i==j))+int(i==j) for j,x in enumerate(row)] for i,row in enumerate(a)]
        ck(det(a)==free+locked,'weighted path determinant')
    path_rows.append({'L':L,'forests':len(records),'edge_coefficients':expected})
ck(decode_geometry((2,1))==(2,) and tiling_blocks((2,1))==((2,2),),'digit2 one-edge hostile')
ck(det([])==1 and fib(0)==0,'true empty versus formal zero')
def terminal_lock(word):
    locked=False
    for d in word:
        if d==0:locked=False
        elif d==2:
            if locked:raise ValueError('illegal')
            locked=True
    return locked
ck(not terminal_lock((2,0,1)) and terminal_lock((1,0,2)), '201 MSB/LSD state distinction')
ck(decode_geometry((2,1,1,1,2))==(2,0,1), '201 geometric orientation')
ck(decode_geometry((1,2,1,2,1))==(1,0,2), '102 geometric orientation')

# General matrix-forest and one-apex claims: exhaustive unfiltered labelled graphs.
counts={};graph_total=arrival_total=0
for n in range(6):
    pairs=tuple(combinations(range(n),2))
    for mask in range(1<<len(pairs)):
        es=tuple(e for i,e in enumerate(pairs) if mask>>i&1)
        value=rooted_parent_count(n,es)
        ck(det(lap_identity(n,es))==value,'matrix forest parent-function count')
        counts[n,es]=value
        graph_total+=1
for n in range(5):
    pairs=tuple(combinations(range(n),2))
    for mask in range(1<<len(pairs)):
        es=tuple(e for i,e in enumerate(pairs) if mask>>i&1)
        a=lap_identity(n,es);old=counts[n,es]
        inv=[[cofactor(a,j,i)/old for j in range(n)] for i in range(n)]
        for pm in range(1<<n):
            ports=tuple(i for i in range(n) if pm>>i&1)
            daughter=tuple(sorted(es+tuple((v,n) for v in ports)))
            actual=counts[n+1,daughter]
            m=[row[:] for row in a]
            for v in ports:m[v][v]+=1
            schur=(1+len(ports))*det(m)-sum(cofactor(m,j,i) for i in ports for j in ports)
            ck(schur==actual,'Schur adjugate response')
            h=[[inv[u][v]+int(i==j) for j,v in enumerate(ports)] for i,u in enumerate(ports)]
            green=old*(det(h)+sum(cofactor(h,j,i) for i in range(len(ports)) for j in range(len(ports))))
            ck(green==actual,'Green adjugate response')
            if len(ports)==1:
                ck(actual==2*old+cofactor(a,ports[0],ports[0]),'leaf response')
            arrival_total+=1
p3=((0,1),(1,2))
ck(counts[3,p3]==8,'P3 count')
ck(counts[4,tuple(sorted(p3+((0,3),))) ]==21,'endpoint21')
ck(counts[4,tuple(sorted(p3+((1,3),))) ]==20,'middle20')
square_counts={n:int(det(lap_identity(n,edges(n)))) for n in (12,30,31,32,33)}
ck(square_counts=={12:9240,30:4337712683055095,31:21517231556064823,32:98839029392469636,33:795810222452889290},'independent square counts')

# Nice-pair witnesses and the strict scope boundary (not all cycle transitions).
p41=(1,24,25,39,10,6,30,19,17,32,4,21,15,34,2,23,26,38,11,5,20,16,9,40,41,8,28,36,13,12,37,27,22,14,35,29,7,18,31,33,3)
p42=(1,24,25,39,10,6,3,33,31,18,7,42,22,27,37,12,13,36,28,21,15,34,30,19,17,32,4,5,11,38,26,23,2,14,35,29,20,16,9,40,41,8)
def cycle_edges(p):return {tuple(sorted(e)) for e in zip(p,p[1:]+p[:1])}
for p in (p41,p42):
    ck(set(p)==set(range(1,len(p)+1)) and len(set(p))==len(p),'nice permutation')
    ck(all(square(u+v) for u,v in cycle_edges(p)),'nice cycle')
ck({p41[i] for i in range(0,41,2)}=={p42[i] for i in range(0,42,2)},'shared phase')
ck((1,3) in cycle_edges(p41) and (1,3) not in cycle_edges(p42),'closure13 change')
nice_deletions=len(cycle_edges(p41)-cycle_edges(p42))
ck(nice_deletions>=2,'nice-pair lower-bound witness')
for v in range(1,10001):ck(not(square(v+1) and square(v+3)),'no common 1,3 neighbor bounded control')
c43=(1,15,21,43,38,26,10,6,19,30,34,2,7,18,31,5,20,29,35,14,11,25,39,42,22,3,33,16,9,27,37,12,4,32,17,8,28,36,13,23,41,40,24)
index=c43.index(5)
c44=c43[:index+1]+(44,)+c43[index+1:]
for p in (c43,c44):
    ck(set(p)==set(range(1,len(p)+1)) and len(set(p))==len(p),'insertion permutation')
    ck(all(square(u+v) for u,v in cycle_edges(p)),'insertion cycle')
ck(cycle_edges(c43)-cycle_edges(c44)=={(5,20)},'general-cycle one-deletion hostile')

here=Path(__file__).resolve().parent
inputs=[here.parent/name for name in ('PROOFS.md','forest_bridge.py','test_forest_bridge.py')]
result={'status':'PASS','checks':CHECKS,'optimization_removable_assertions':False,
        'arrival_orders':[1,4000],'apery_orders':[0,10000],
        'path_rows':path_rows,'all_labelled_graphs_orders_0_through_5':graph_total,
        'all_arrival_queries_old_orders_0_through_4':arrival_total,
        'all_two_factors':factor_rows,'square_sum_rooted_forest_counts':square_counts,
        'nice_pair_witness_old_edge_deletions':nice_deletions,
        'outside_pair_43_to_44_deletions':1,
        'input_sha256':{str(p.relative_to(here.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
output=Path(sys.argv[1] if len(sys.argv)>1 else here/'independent.normal.json')
output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':result['status'],'checks':CHECKS,'graphs':graph_total,'arrivals':arrival_total,'output':str(output)},sort_keys=True))
