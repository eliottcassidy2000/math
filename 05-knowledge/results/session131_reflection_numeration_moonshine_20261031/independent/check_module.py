#!/usr/bin/env python3
"""Independent integral C2 cycle-module audit; no producer imports; stdlib only."""
from collections import deque
from itertools import combinations, combinations_with_replacement, permutations, product
import json
import random
import sys

GATES = 0

def need(condition, message):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(message)


def incidence(n, edges, z):
    ans = [0] * n
    for (u, v), k in zip(edges, z):
        ans[u] -= k
        ans[v] += k
    return ans


def forest(n, edges, order=None):
    roots = list(range(n))
    def root(v):
        while roots[v] != v:
            roots[v] = roots[roots[v]]
            v = roots[v]
        return v
    chosen = []
    for j in range(len(edges)) if order is None else order:
        u, v = edges[j]
        u, v = root(u), root(v)
        if u != v:
            roots[u] = v
            chosen.append(j)
    return chosen, len({root(v) for v in range(n)})


def fundamental(n, edges, chosen):
    adjacency = [[] for _ in range(n)]
    for j in chosen:
        u, v = edges[j]
        adjacency[u].append((v, j, 1))
        adjacency[v].append((u, j, -1))
    chosen = set(chosen)
    cycles = {}
    for j, (u, v) in enumerate(edges):
        if j in chosen:
            continue
        z = [0] * len(edges)
        z[j] = 1
        parent = {u: None}
        queue = deque([u])
        while v not in parent:
            need(bool(queue), 'forest does not span an edge')
            a = queue.popleft()
            for b, k, sign in adjacency[a]:
                if b not in parent:
                    parent[b] = (a, k, sign)
                    queue.append(b)
        b = v
        while b != u:
            a, k, sign = parent[b]
            z[k] -= sign
            b = a
        need(incidence(n, edges, z) == [0] * n, 'fundamental cycle boundary')
        cycles[j] = z
    return cycles


def determinant(matrix):
    a = [row[:] for row in matrix]
    n = len(a)
    if not n:
        return 1
    denominator, sign = 1, 1
    for k in range(n - 1):
        if not a[k][k]:
            pivot = next((j for j in range(k + 1, n) if a[j][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        value = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = value * a[i][j] - a[i][k] * a[k][j]
                need(numerator % denominator == 0, 'fraction-free division')
                a[i][j] = numerator // denominator
            a[i][k] = 0
        denominator = value
    return sign * a[-1][-1]


def examine(n, blue, black):
    qedges = blue + black
    b = len(blue)
    qforest, cq = forest(n, qedges)
    bforest, cb = forest(n, blue)
    need(qforest[:len(bforest)] == bforest, 'blue forest extension')
    bedges = list(blue)
    rperm = list(range(b))
    for uv in black:
        j = len(bedges)
        bedges.extend([uv, uv])
        rperm.extend([j + 1, j])
    lift_index = lambda j: j if j < b else b + 2 * (j - b)
    fullforest = [lift_index(j) for j in qforest]
    fundamental_cycles = fundamental(n, bedges, fullforest)
    blue_chords = [j for j in range(b) if j not in qforest]
    tree_black = [j for j in qforest if j >= b]
    chord_black = [j for j in range(b, len(qedges)) if j not in qforest]
    reflect = lambda z: [z[rperm[j]] for j in range(len(bedges))]
    basis, eigens = [], []
    for j in blue_chords:
        z = fundamental_cycles[j]
        need(reflect(z) == z, 'blue chord not invariant')
        basis.append(z)
        eigens.append(z)
    for j in tree_black:
        a = lift_index(j)
        z = fundamental_cycles[a + 1]
        need(reflect(z) == [-x for x in z], 'black tree difference not anti-invariant')
        basis.append(z)
        eigens.append(z)
    for j in chord_black:
        a = lift_index(j)
        z = fundamental_cycles[a]
        rz = reflect(z)
        need(incidence(n, bedges, rz) == [0] * n, 'reflection changes boundary')
        basis.extend([z, rz])
        eigens.extend([[x+y for x,y in zip(z,rz)], [x-y for x,y in zip(z,rz)]])
    independent_chords = sorted(fundamental_cycles)
    coordinate_matrix = [[z[j] for z in basis] for j in independent_chords]
    eigencoordinates = [[z[j] for z in eigens] for j in independent_chords]
    need(len(basis) == len(bedges)-n+cq, 'cycle rank')
    need(abs(determinant(coordinate_matrix)) == 1, 'normal form basis not unimodular')
    trivial, sign, exchange = b-n+cb, cb-cq, len(black)+cq-cb
    need((len(blue_chords),len(tree_black),len(chord_black)) == (trivial,sign,exchange), 'multiplicity formulas')
    need(abs(determinant(eigencoordinates)) == 2**exchange, 'eigenlattice index')
    need(trivial+exchange == len(qedges)-n+cq, 'invariant rank')
    need(sign+exchange == len(black), 'anti-invariant rank')
    return {'vertices':n,'blue_edges':b,'black_pairs':len(black),'c_blue':cb,'c_quotient':cq,
            'trivial':trivial,'sign':sign,'exchange':exchange,'cycle_rank':len(basis),'eigen_index':2**exchange}


def reconstruct_order4():
    pairs=list(combinations(range(4),2))
    tiles=[(0,2),(0,3),(1,3)]
    canonical=[]
    for mask in range(8):
        a=[[0]*4 for _ in range(4)]
        for u,v in pairs:
            bit=1 if v==u+1 else (mask>>tiles.index((u,v)))&1
            a[u][v]=bit;a[v][u]=1-bit
        codes=[]
        for perm in permutations(range(4)):
            code=sum(a[perm[u]][perm[v]]<<i for i,(u,v) in enumerate(pairs))
            codes.extend([code,code^63])
        canonical.append(min(codes))
    nodes={code:i for i,code in enumerate(sorted(set(canonical)))}
    blue=[];black=[];lines=[]
    for mask in range(4):
        mirror=((mask&1)<<2)|(mask&2)|((mask&4)>>2)
        endpoints=tuple(sorted((nodes[canonical[mask]],nodes[canonical[mask^7]])))
        line={'masks':[mask,mask^7],'endpoints':endpoints,'blue':mirror==mask}
        lines.append(line)
        if mirror==mask:blue.append(endpoints)
        elif mask<min(mirror,mirror^7):black.append(endpoints)
    result=examine(len(nodes),blue,black)
    need((result['trivial'],result['sign'],result['exchange'])==(1,1,0),'order4 literal reconstruction')
    return {'node_keys':sorted(nodes),'tiling_class_keys':canonical,'literal_lines':lines,'module':result}


def main():
    if '--test-failing-gate' in sys.argv:
        need(False, 'intentional failure control')
    count_by_n={}
    for n in range(4):
        slots=list(combinations_with_replacement(range(n),2))
        count=0
        for colors in product(range(4), repeat=len(slots)):
            blue=[e for e,c in zip(slots,colors) if c&1]
            black=[e for e,c in zip(slots,colors) if c&2]
            examine(n,blue,black)
            count+=1
        count_by_n[n]=count
    rng=random.Random(131)
    samples=0
    for n in range(4,13):
        for _ in range(40):
            blue=[];black=[]
            for _ in range(rng.randrange(n+5)):
                uv=(rng.randrange(n),rng.randrange(n))
                (blue if rng.randrange(2) else black).append(uv)
            examine(n,blue,black)
            samples+=1
    blue=[(0,5),(3,8),(4,6),(5,6),(5,8),(6,9),(7,9),(8,9)]
    black=[(1,4),(2,4),(2,5),(2,6),(2,7),(2,8),(4,5),(4,8),(5,7),(6,6),(6,8),(8,8)]
    actual=examine(10,blue,black)
    controls={
        'black_bridge':examine(2,[],[(0,1)]),
        'black_loop':examine(1,[],[(0,0)]),
        'blue_loop':examine(1,[(0,0)],[]),
        'blue_black_parallel':examine(2,[(0,1)],[(0,1)]),
        'parallel_blue':examine(2,[(0,1),(0,1)],[]),
        'disconnected_isolated':examine(4,[(0,1)],[(2,2)]),
    }
    need(controls['black_bridge']['sign']==1 and controls['black_bridge']['exchange']==0,'bridge sign control')
    need(controls['black_loop']['exchange']==1 and controls['black_loop']['eigen_index']==2,'loop parity hostile')
    need(controls['blue_loop']['trivial']==1,'fixed-loop trivial control')
    need(actual['trivial']==1 and actual['sign']==2 and actual['exchange']==10,'order-five normal form')
    print(json.dumps({'status':'PASS','scope':'all n<=3 graphs with zero/one blue edge and zero/one black pair per unordered pair including loops; 360 seeded larger multigraphs',
        'exhaustive_counts':count_by_n,'seeded_multigraphs':samples,'order4':reconstruct_order4(),'order5':actual,'controls':controls,'gates':GATES},indent=2,sort_keys=True))

if __name__=='__main__':main()
