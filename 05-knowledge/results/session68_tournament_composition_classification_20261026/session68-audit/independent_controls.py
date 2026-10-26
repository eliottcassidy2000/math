#!/usr/bin/env python3
"""Independent bounded controls; no producer imports, repository input, or asserts."""
from collections import Counter, deque
from itertools import permutations, product, combinations
import json

def need(ok, msg):
    if not ok: raise RuntimeError(msg)

def triangle_palette(colors, inside=(0,0,0), outer=0):
    rows = [[0]*12 for _ in range(12)]
    for u,v in combinations(range(12),2):
        b,x=divmod(u,4);c,y=divmod(v,4)
        if b!=c: wins=((c-b)%3==1)^bool(outer)
        elif x==0: wins=not colors[b]
        else: wins=((y-x)%3==1)^bool(inside[b])
        rows[u][v]=int(wins);rows[v][u]=int(not wins)
    return tuple(map(tuple,rows))

def signature(rows):
    return tuple(rows[u][v] for u,v in combinations(range(len(rows)),2))

def triangles(rows):
    return frozenset(t for t in combinations(range(len(rows)),3) if all(sum(rows[u][v] for v in t)==1 for u in t))

def image_sig(rows,p):
    out=[[0]*len(p) for _ in p]
    for u in range(len(p)):
        for v in range(len(p)):out[p[u]][p[v]]=rows[u][v]
    return signature(out)

# Full labelled node-choice space at n=12, not all 2^66 tournaments.
reals={signature(triangle_palette(a,b,c)) for a in product((0,1),repeat=3) for b in product((0,1),repeat=3) for c in (0,1)}
need(len(reals)==128,'128 distinct BCIZ node choices')
perms=[]
for p in permutations(range(3)):
    for inner in product(tuple(permutations((1,2,3))),repeat=3):
        perms.append(tuple(4*p[b]+(0 if j==0 else inner[b][j-1]) for b in range(3) for j in range(4)))
need(len(set(perms))==1296,'hypergraph group action size')
base=triangles(triangle_palette((0,0,0)))
need(len(base)==67,'67 cyclic triples')
need(all(frozenset(tuple(sorted(p[u] for u in t)) for t in base)==base for p in perms),'every candidate preserves hypergraph')
orbits=[]
for w in ((0,0,0),(0,0,1),(0,1,1),(1,1,1)):
    rows=triangle_palette(w);orb={image_sig(rows,p) for p in perms};orbits.append(orb)
need([len(o) for o in orbits]==[16,48,48,16],'four fibre orbit sizes')
need(set.union(*orbits)==reals and sum(map(len,orbits))==len(reals),'palette orbits exhaust node-choice fibre')

# W2 orbit classification from generated moves, no recurrence or all-group sum.
def rotation_map(block=None):
    if block is None:return tuple(3*((i//3+1)%3)+i%3 for i in range(9))
    return tuple(3*block+(i%3+1)%3 if i//3==block else i for i in range(9))
generators=[rotation_map()]+[rotation_map(b) for b in range(3)]
def act(x,p):return sum(((x>>i)&1)<<p[i] for i in range(9))
unseen=set(range(512));hist=Counter();orbit_sizes=[]
while unseen:
    first=min(unseen);orb={first};todo=deque([first])
    while todo:
        x=todo.popleft()
        for g in generators:
            y=act(x,g)
            if y not in orb:orb.add(y);todo.append(y)
    unseen-=orb;stab=81//len(orb);e=0
    while stab>1:need(stab%3==0,'stabilizer power 3');stab//=3;e+=1
    hist[e]+=1;orbit_sizes.append(len(orb))
need([hist[e] for e in range(5)]==[2,10,8,2,2],'P2 by move-generated BFS')

# Direct Hamilton-word DFS, independent of endpoint subset DP.
flux=[]
for colors in ((0,0,0),(1,0,0)):
    rows=triangle_palette(colors,inside=colors)
    # inside=colors makes the 1-block the literal converse Aop.
    n=12;full=(1<<n)-1;out=[sum(v<<j for j,v in enumerate(row)) for row in rows];scores=[x.bit_count() for x in out];counts=Counter()
    def visit(start,last,used):
        if used==full:counts[scores[start]-scores[last]-1]+=1;return
        nxt=out[last]&~used
        while nxt:
            bit=nxt&-nxt;nxt-=bit;visit(start,bit.bit_length()-1,used|bit)
    for s in range(n):visit(s,s,1<<s)
    flux.append(dict(sorted(counts.items())))
need(flux==[{-3:73710,-1:267624,1:121581},{-4:4293,-3:44064,-2:86535,-1:132084,0:120258,1:64017,2:11664}],'literal DFS flux histograms')

# All unordered directed path covers of C3[A,1,1] via permutations and cuts.
small=[[0]*6 for _ in range(6)]
A=triangle_palette((0,0,0))
for i,j in combinations(range(6),2):
    if j<4:win=A[i][j]
    elif i<4:win=(j==4)
    else:win=1
    small[i][j]=int(win);small[j][i]=int(not win)
other=[r[:] for r in small]
for i,j in combinations(range(4),2):other[i][j],other[j][i]=other[j][i],other[i][j]
def covers(rows):
    ans=set()
    for p in permutations(range(6)):
        for cuts in range(32):
            if any(not(cuts>>i&1) and not rows[p[i]][p[i+1]] for i in range(5)):continue
            parts=[];start=0
            for i in range(5):
                if cuts>>i&1:parts.append(p[start:i+1]);start=i+1
            parts.append(p[start:]);ans.add(tuple(sorted(parts)))
    return ans
def reverse_runs(cover):
    result=[]
    for path in cover:
        arr=[];i=0
        while i<len(path):
            j=i+1
            if path[i]<4:
                while j<len(path) and path[j]<4:j+=1
                arr.extend(reversed(path[i:j]))
            else:arr.append(path[i])
            i=j
        result.append(tuple(arr))
    return tuple(sorted(result))
c0,c1=covers(small),covers(other)
need({reverse_runs(c) for c in c0}==c1,'all run-reversal cover images')
need(all(reverse_runs(reverse_runs(c))==c for c in c0),'run-reversal involution')
print(json.dumps({'status':'PASS','full_fibre_n12':{'node_choice_realizations':128,'hypergraph_action_order':1296,'ordinary_orbit_sizes':[len(o) for o in orbits],'automorphism_orders':[1296//len(o) for o in orbits]},'W2_BFS':{'labelled_colorings':512,'ordinary_types':len(orbit_sizes),'stabilizer_exponent_histogram':dict(sorted(hist.items()))},'direct_Hamilton_DFS_flux':flux,'module_reversal_all_covers_n6':{'covers':len(c0),'profile':dict(sorted(Counter(map(len,c0)).items()))}},sort_keys=True,indent=2))
