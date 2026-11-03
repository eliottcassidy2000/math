#!/usr/bin/env python3
"""Exact finite controls for local/global/guard bridges. No external packages."""
from itertools import permutations, product
from collections import deque, Counter
from fractions import Fraction
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def mm(A, B):
    return tuple(tuple(sum(a*b for a,b in zip(row,col)) for col in zip(*B)) for row in A)


def mv(A, v):
    return tuple(sum(a*b for a,b in zip(row,v)) for row in A)


C=((0,1,0),(1,1,1),(1,0,0))
D=((1,0,1),(0,2,0),(0,0,0))
O=((1,1,0),(1,2,1),(3,3,2))
B=((1,-2,1),(1,2,-1),(-3,0,1))
CO=((0,1,0),(0,0,1),(1,1,1))
DO=((0,1,0),(1,3,-1),(0,3,0))
require(mm(B,O)==((2,0,0),(0,2,0),(0,0,2)), 'inverse certificate')
require(mm(O,C)==mm(CO,O), 'C conjugation')
require(mm(O,D)==mm(DO,O), 'D conjugation')
require(mm(C,((0,0,1),(1,0,0),(-1,1,-1)))==((1,0,0),(0,1,0),(0,0,1)), 'C inverse')

lattice_count=0
for h in product(range(-9,10),repeat=3):
    nums=mv(B,h)
    integral=all(x%2==0 for x in nums)
    require(integral == ((h[0]-h[2])%2==0), 'index-two criterion')
    if integral:
        require(mv(O,tuple(x//2 for x in nums))==h,'lattice round trip')
        lattice_count+=1

seed=(1,2,0)
def word_state(word):
    v=seed
    for ch in word:
        v=mv(C if ch=='C' else D,v)
    return v

v1,v2=word_state('CDC'),word_state('DDC')
require((v1,v2)==((6,9,3),(8,9,1)), 'collision inputs')
require(v1!=v2 and mv(D,v1)==mv(D,v2)==(9,18,0),'physical D collision')
require((mv(O,v1),mv(O,v2))==((15,27,51),(17,27,53)),'observer collision inputs')

# Finite controls do not establish the all-word or all-n theorems.
frontier=[seed]
word_count=0
zero_d_count=0
for depth in range(16):
    word_count+=len(frontier)
    for s,u,d in frontier:
        require(s>0 and u>0 and d>=0 and u>=s+d,'guard cone')
        h=mv(O,(s,u,d))
        require(all(x%2==1 for x in h),'odd observed counts')
        if d==0:
            zero_d_count+=1
            require(not(Fraction(1,4)<Fraction(s,u)<Fraction(3,8)), 'state gap')
            require(not(Fraction(19,11)<Fraction(h[1],h[0])<Fraction(9,5)), 'observed gap')
    frontier=[mv(M,v) for v in frontier for M in (C,D)]

# Independent endpoint DP on the literal tournament Y_n.
def path_counts(n, edges):
    dp={}
    starts=[]
    for start in range(n):
        d={(1<<start,start):1}
        for mask_size in range(1,n):
            nxt={}
            for (mask,end),num in d.items():
                for w in range(n):
                    if not(mask>>w&1) and (end,w) in edges:
                        key=(mask|1<<w,w)
                        nxt[key]=nxt.get(key,0)+num
            d=nxt
        starts.append(sum(d.values()))
    return starts

family=[]
for n in range(6,11):
    edges={(i,j) for i in range(n) for j in range(i+1,n)}
    for target in (0,2):
        edges.remove((target,n-1)); edges.add((n-1,target))
    starts=path_counts(n,edges)
    reduced={(a-1,b-1) for a,b in edges if a!=0 and b!=0}
    deleted=path_counts(n-1,reduced)
    actual=(starts[0],sum(starts[1:]),deleted[-1])
    t=2**(n-4)
    predicted=(t+1,4*t-2,0)
    require(actual==predicted,'Y_n DP')
    h=mv(O,actual)
    require(h==(5*t-1,9*t-3,15*t-3),'Y_n observations')
    require(Fraction(19,11)<Fraction(h[1],h[0])<Fraction(9,5),'Y_n forbidden interval')
    # k=1..8 strict positivity forces the last move to be C.
    v=actual
    for k in range(1,9):
        v=mv(C,v)
        require(v[2]>0,'C-tail reversal')
    family.append({'n':n,'state':actual,'observations':h})

# Permutation convention: tuple p sends i to p[i].
def compose(p,q):
    return tuple(p[q[i]] for i in range(len(p)))
def inverse(p):
    out=[0]*len(p)
    for i,j in enumerate(p):out[j]=i
    return tuple(out)
def generated_group(gens,d):
    identity=tuple(range(d)); seen={identity};todo=[identity]
    while todo:
        p=todo.pop()
        for g in gens:
            q=compose(g,p)
            if q not in seen:seen.add(q);todo.append(q)
    return seen

def cover(vertices,edges,voltages,d):
    adj={(v,i):set() for v in vertices for i in range(d)}
    identity=tuple(range(d))
    for a,b in edges:
        p=voltages.get((a,b),identity)
        for i in range(d):
            x,y=(a,i),(b,p[i]);adj[x].add(y);adj[y].add(x)
    base={v:set() for v in vertices}
    for a,b in edges:base[a].add(b);base[b].add(a)
    for (v,i),ns in adj.items():
        require(len(ns)==len(base[v]),'local cardinality')
        require({w for w,j in ns}==base[v],'local surjectivity')
    seen={next(iter(adj))};todo=list(seen)
    while todo:
        x=todo.pop()
        for y in adj[x]:
            if y not in seen:seen.add(y);todo.append(y)
    return adj,len(seen)==len(adj)

def centralizer_count(gens,d):
    return sum(all(compose(p,g)==compose(g,p) for g in gens) for p in permutations(range(d)))

def deck_root_test(adj,d):
    # All tree voltages below are identity. Thus every deck map uses a common
    # sheet permutation in all fibers. Check each candidate on every edge.
    out=[]
    for p in permutations(range(d)):
        if all({(w,p[j]) for w,j in ns}==adj[(v,p[i])] for (v,i),ns in adj.items()):out.append(p)
    return out

# Diamond: identity star tree; both non-tree edges freely receive S3 voltages.
vertices=list(range(4));edges=[(0,1),(0,2),(0,3),(1,2),(1,3)]
perms=list(permutations(range(3)))
census=Counter()
for a,b in product(perms,repeat=2):
    adj,connected=cover(vertices,edges,{(1,2):a,(1,3):b},3)
    group=generated_group([a,b],3)
    transitive={p[0] for p in group}==set(range(3))
    require(connected==transitive,'connectivity/monodromy agreement')
    decks=deck_root_test(adj,3)
    require(len(decks)==centralizer_count([a,b],3),'deck/centralizer agreement')
    census[(connected,len(decks))]+=1
require(sum(census.values())==36,'census size')
require(census[(True,1)]==18 and census[(True,3)]==8,'connected census')
a=(1,0,2);b=(1,2,0)
adj,connected=cover(vertices,edges,{(1,2):a,(1,3):b},3)
require(connected and len(deck_root_test(adj,3))==1,'sharp 12-vertex cover')
# Independent brute force: all 6^4 fiberwise permutations, not a centralizer test.
proper_color_auts=0
for pp in product(perms,repeat=4):
    if all({(w,pp[w][j]) for w,j in ns}==adj[(v,pp[v][i])] for (v,i),ns in adj.items()):proper_color_auts+=1
require(proper_color_auts==1,'brute-force colored rigidity')

degree_checks=[]
for d in range(2,8):
    a=tuple([1,0]+list(range(2,d)));b=tuple(list(range(1,d))+[0])
    adj,connected=cover(vertices,edges,{(1,2):a,(1,3):b},d)
    require(connected,'degree-d cover connectivity')
    deck_count=len(deck_root_test(adj,d))
    require(deck_count==(2 if d==2 else 1),'degree-d rigidity')
    degree_checks.append({'degree':d,'vertices':4*d,'deck_count':deck_count})

# S165's biclosed C6 x C4 bank, computed directly from denominators.
denoms=(4,6,9,5,7)
from math import lcm
closed=set()
for mask in range(1<<5):
    least=lcm(*(denoms[i] for i in range(5) if mask>>i&1))
    obs=sum(1<<i for i,b in enumerate(denoms) if least%b==0)
    if obs==mask:closed.add(mask)
bank=sorted(m for m in closed if (31^m) in closed)
bank_edges=[(a,b) for a in bank for b in bank if a<b and (a^b).bit_count()==1]
require(len(bank)==24 and len(bank_edges)==48,'S165 bank')
baseadj={v:set() for v in bank}
for a,b in bank_edges:baseadj[a].add(b);baseadj[b].add(a)
seen={bank[0]};q=deque([bank[0]]);tree=set()
while q:
    v=q.popleft()
    for w in sorted(baseadj[v]):
        if w not in seen:seen.add(w);q.append(w);tree.add(tuple(sorted((v,w))))
non_tree=[e for e in bank_edges if e not in tree]
require(len(non_tree)==25,'cycle rank')
volts={non_tree[0]:(1,0,2),non_tree[1]:(1,2,0)}
adj,connected=cover(bank,bank_edges,volts,3)
require(connected and len(deck_root_test(adj,3))==1,'rigid guard-state extension')
require(len(adj)==72 and sum(map(len,adj.values()))//2==144,'72-state cover size')

# Two-color transport: local flip constraints close iff a cycle is even.
color_cycles=[]
for length in range(3,11):
    valid=[bits for bits in product(range(2),repeat=length) if all(bits[i]!=bits[(i+1)%length] for i in range(length))]
    require(len(valid)==(2 if length%2==0 else 0),'color transport obstruction')
    color_cycles.append({'length':length,'colorings':len(valid)})

print(json.dumps({
    'status':'PASS',
    'observation_lattice_grid':{'triples':19**3,'integral_preimages':lattice_count},
    'word_controls':{'max_length':15,'words':word_count,'zero_d_words':zero_d_count},
    'physical_rank_loss':{'words':['CDC','DDC'],'before':[v1,v2],'after_D':mv(D,v1)},
    'literal_Y_family':family,
    'diamond_S3_census':[{'connected':c,'deck_count':k,'assignments':v} for (c,k),v in sorted(census.items())],
    'independent_diamond_fiber_permutations':6**4,
    'independent_diamond_color_automorphisms':proper_color_auts,
    'all_degree_controls':degree_checks,
    'guard_bank_cover':{'base_vertices':24,'base_edges':48,'cycle_rank':25,'vertices':72,'edges':144,'twisted_edges':list(volts),'connected':connected,'color_automorphisms':len(deck_root_test(adj,3))},
    'two_color_cycles':color_cycles,
},indent=2,sort_keys=True))
