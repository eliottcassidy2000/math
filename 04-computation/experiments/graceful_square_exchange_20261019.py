#!/usr/bin/env python3
"""Exact, standalone controls for graceful/square-sum path exchanges.
Only Python's standard library. No truth-bearing assert statements.
Run: python3 graceful_square_exchange_20261019.py
Optional heavier path frontier: --graceful-max 14.
"""
from collections import Counter, defaultdict
from itertools import combinations, permutations
from math import isqrt
import argparse
import hashlib
import json

def check(ok, message):
    if not ok: raise RuntimeError(message)

def canonical(p):
    p = tuple(p)
    return min(p, p[::-1])

def edges(p):
    return frozenset(tuple(sorted(e)) for e in zip(p, p[1:]))

def graceful(p):
    n = len(p)
    return set(p) == set(range(n)) and sorted(abs(x-y) for x,y in zip(p,p[1:])) == list(range(1,n))

def square_path(p):
    return set(p) == set(range(1,len(p)+1)) and all(isqrt(x+y)**2 == x+y for x,y in zip(p,p[1:]))

def graceful_paths(n):
    result = []
    def visit(p, used, differences):
        if len(p) == n:
            if p[0] < p[-1]: result.append(tuple(p))
            return
        x = p[-1]
        for y in range(n):
            if not (used >> y) & 1:
                d = abs(x-y)
                if not (differences >> d) & 1:
                    visit(p+[y], used | (1<<y), differences | (1<<d))
    for x in range(n): visit([x], 1<<x, 0)
    return result

def rotations(p, graceful_mode=False):
    result = set()
    for a in (p, p[::-1]):
        x = a[-1]
        for i in range(len(a)-2):
            v,w = a[i:i+2]
            allowed = 2*v == x+w if graceful_mode else isqrt(v+x)**2 == v+x
            if allowed: result.add(canonical(a[:i+1] + a[:i:-1]))
    return result

def splice(p):
    d = abs(p[0]-p[-1])
    i = next(i for i in range(len(p)-1) if abs(p[i]-p[i+1]) == d)
    return canonical(p[i::-1] + p[:i:-1])

class DSU:
    def __init__(self,n): self.p = list(range(n))
    def find(self,i):
        while self.p[i] != i:
            self.p[i] = self.p[self.p[i]]
            i = self.p[i]
        return i
    def join(self,i,j): self.p[self.find(i)] = self.find(j)
    def sizes(self): return sorted(Counter(self.find(i) for i in range(len(self.p))).values())

def exchange_components(paths,r):
    if not paths: return []
    n = len(paths[0]); ds = DSU(len(paths)); bank = {}
    for i,p in enumerate(paths):
        bits = [1 << (x*n+y) for x,y in sorted(edges(p))]
        mask = sum(bits)
        for removed in combinations(bits,min(r,n-1)):
            signature = mask
            for b in removed: signature ^= b
            if signature in bank: ds.join(i,bank[signature])
            else: bank[signature] = i
    return ds.sizes()

def components(paths, mode=False):
    pending = set(paths); result = []
    while pending:
        p = min(pending); pending.remove(p); queue=[p]; comp=[]
        while queue:
            p=queue.pop(); comp.append(p)
            for q in rotations(p,mode):
                if q in pending: pending.remove(q); queue.append(q)
        result.append(comp)
    return result

def walecki(n):
    return tuple(i//2 if i%2 == 0 else n-1-i//2 for i in range(n))

def square_paths(n):
    adjacency = [sum(1 << (y-1) for y in range(1,n+1)
                    if y != x and isqrt(x+y)**2 == x+y) for x in range(1,n+1)]
    full=(1<<n)-1; result=[]
    leaves=[i for i in range(n) if adjacency[i].bit_count() <= 1]
    if len(leaves) > 2: return result
    starts=leaves[:1] if leaves else range(n)
    def visit(p, remaining):
        if not remaining:
            result.append(canonical(tuple(x+1 for x in p))); return
        v=p[-1]; choices=adjacency[v]&remaining
        if not choices: return
        active=remaining|(1<<v); ones=0
        for u in range(n):
            if (remaining>>u)&1:
                degree=(adjacency[u]&active).bit_count()
                if degree == 0: return
                if degree == 1:
                    ones+=1
                    if ones > 1: return
        seen=0; frontier=1<<v
        while frontier:
            seen |= frontier; nxt=0
            while frontier:
                b=frontier&-frontier; frontier-=b
                nxt |= adjacency[b.bit_length()-1]
            frontier=nxt&active&~seen
        if seen != active: return
        while choices:
            b=choices&-choices; choices-=b
            visit(p+[b.bit_length()-1],remaining^b)
    for start in starts: visit([start],full^(1<<start))
    return sorted(set(result))

def graceful_trees(n):
    result=[]
    def visit(code,parent,d):
        if d == 0:
            result.append(tuple(reversed(code))); return
        def root(i):
            while parent[i] != i: i=parent[i]
            return i
        for a in range(n-d):
            x,y=root(a),root(a+d)
            if x != y:
                pp=parent.copy(); pp[x]=y
                visit(code+[a],pp,d-1)
    visit([],list(range(n)),n-1)
    return result

def tree_components(codes):
    index={t:i for i,t in enumerate(codes)}; ds=DSU(len(codes)); n=len(codes[0])+1
    for i,t in enumerate(codes):
        for j,a in enumerate(t):
            for b in range(a+1,n-j-1):
                q=t[:j]+(b,)+t[j+1:]
                if q in index: ds.join(i,index[q])
    return ds.sizes()

def prufer_graceful(n):
    out=set()
    from itertools import product
    for code in product(range(n),repeat=n-2):
        degree=[1]*n
        for v in code: degree[v]+=1
        ee=[]
        for v in code:
            u=next(i for i,d in enumerate(degree) if d == 1)
            ee.append((u,v)); degree[u]-=1; degree[v]-=1
        u,v=[i for i,d in enumerate(degree) if d == 1]; ee.append((u,v))
        if sorted(abs(u-v) for u,v in ee) == list(range(1,n)):
            out.add(tuple(a for d,a in sorted((abs(u-v),min(u,v)) for u,v in ee)))
    return out

def shape_summary(sizes): return sorted(Counter(sizes).items())

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--graceful-max',type=int,default=13)
    ap.add_argument('--tree-max',type=int,default=11); ap.add_argument('--square-max',type=int,default=36)
    args=ap.parse_args()
    print('Universe: paths modulo reversal; complements are distinct; trees have fixed numeric labels.')
    print('G n paths rotation_components E1_components E2_components E3_components')
    path_counts={2:1,3:2,4:2,5:4,6:12,7:16,8:20,9:60,10:148,11:324,12:664,13:1600}
    for n in range(2,args.graceful_max+1):
        ps=graceful_paths(n); lookup=set(ps)
        if n in path_counts: check(len(ps)==path_counts[n],f'graceful count {n}')
        check(all(graceful(p) for p in ps),'graceful predicate')
        if n<=8:
            independent={canonical(p) for p in permutations(range(n)) if graceful(p)}
            check(lookup==independent,'independent permutation enumeration')
        check([p for p in ps if p[0]==0] == [walecki(n)],'extreme endpoint uniqueness')
        for p in ps:
            long_edges={e for e in edges(p) if 2*abs(e[0]-e[1]) > n-1}
            for q in rotations(p,True):
                check(q in lookup,'rotation closure')
                check(long_edges <= edges(q),'long edge freezing')
            if n>=3:
                q=splice(p)
                check(q in lookup and q!=p and splice(q)==p,'splice involution')
                check(abs(q[0]-q[-1])==abs(p[0]-p[-1]),'splice endpoint gap')
        ecounts=[len(exchange_components(ps,r)) for r in (1,2,3)]
        print('G',n,len(ps),len(components(ps,True)),*ecounts,flush=True)
    for n in range(2,501):
        rr=rotations(walecki(n),True)
        check(len(rr)==int(n%3==0),'Walecki rotation residue')
        if rr:
            q=next(iter(rr)); check(edges(walecki(n))-edges(q)=={(n//3,2*n//3)},'Walecki removed edge')
            check(edges(q)-edges(walecki(n))=={(0,n//3)},'Walecki inserted edge')
    print('Walecki midpoint rotation iff 3 divides n: checked 2..500; exact proof in note.')
    for a,b,c,d in combinations(range(11),4):
        spectra=[sorted([b-a,d-c]),sorted([c-a,d-b]),sorted([d-a,c-b])]
        check(len({tuple(x) for x in spectra})==3,'matching spectrum hostile')
    p=(1,4,0,2,3);q=(1,2,4,0,3)
    check(graceful(p) and graceful(q) and (p[0],p[-1])==(q[0],q[-1]) and len(edges(p)-edges(q))==3,'five vertex trade')
    print('Minimum five-vertex fixed-endpoint graceful 3-edge trade:',p,q)
    print('T n gracefully_labeled_trees E1_components decreasing_coordinate_local_minima')
    tcounts={2:1,3:2,4:4,5:12,6:40,7:164,8:752,9:4020,10:23576,11:155632}
    for n in range(2,args.tree_max+1):
        ts=graceful_trees(n); lookup=set(ts)
        if n in tcounts: check(len(ts)==tcounts[n],f'tree count {n}')
        if n<=7: check(lookup==prufer_graceful(n),'independent Prufer census')
        minima=sum(not any(t[:j]+(b,)+t[j+1:] in lookup for j,a in enumerate(t) for b in range(a)) for t in ts)
        print('T',n,len(ts),len(tree_components(ts)),minima,flush=True)
    print('S n paths rotation_components extendable_rotation_components component_sizes')
    known={15:1,16:1,17:1,18:0,19:0,20:0,21:0,22:0,23:3,24:0,25:10,26:12,27:35,28:52,29:19,30:20,31:349,32:392,33:669,34:4041,35:17175,36:12960}
    cache={}
    for n in range(15,args.square_max+1):
        ps=square_paths(n);cache[n]=ps
        if n in known: check(len(ps)==known[n],f'square path count {n}')
        check(all(square_path(p) for p in ps),'square path predicate')
        cs=components(ps); targets={x for x in range(1,n+1) if isqrt(x+n+1)**2==x+n+1}
        hits=sum(any(p[0] in targets or p[-1] in targets for p in c) for c in cs)
        print('S',n,len(ps),len(cs),hits,shape_summary(list(map(len,cs))),flush=True)
    p=(1,24,25,11,14,22,27,9,16,20,5,4,12,13,3,6,19,17,8,28,21,15,10,26,23,2,7,18)
    q=(1,24,25,11,5,4,12,13,3,6,19,17,8,28,21,15,10,26,23,2,14,22,27,9,16,20,29,7,18)
    check(square_path(p) and square_path(q),'28 to29 witness')
    check(len(edges(p)-edges(q))==3,'three-edge insertion')
    print('28->29 removed:',sorted(edges(p)-edges(q)))
    print('28->29 added:',sorted(edges(q)-edges(p)))
    if 28 in cache and 29 in cache:
        qe=[edges(q) for q in cache[29]]
        minima=Counter(min(len(edges(p)-eq) for eq in qe) for p in cache[28])
        check(min(minima)==3,'sharp 28->29 minimum')
        check(minima==Counter({3:21,4:11,5:15,6:5}),'28->29 source histogram')
        print('28->29 minimum removed-edge histogram over all52 source paths:',sorted(minima.items()))
        frozen=hashlib.sha256(json.dumps([cache[28],cache[29]],separators=(',',':')).encode()).hexdigest()
        print('SHA256 compact JSON [Q28 paths,Q29 paths]:',frozen)
    print('ALL CHECKS PASSED. No universal conjecture is proved by these finite censuses.')

if __name__=='__main__': main()
