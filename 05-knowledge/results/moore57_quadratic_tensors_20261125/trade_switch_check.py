"""Verify the Klein trades, cube signs, and overlap geometry."""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
import sys

from slot_interface_check import BRANCHES, SOURCE, leaf, read_board
from nested_core_check import graph, legal, orbit_chart

reference, raw=read_board(sys.argv[1])
m=raw[0][7]
R={(h,y):(h,*t) for h,D in enumerate(raw) if h for y,t in D.items()
   if all(x!=z for x,z in zip(t,m))}
ref={(h,y) for h,y in enumerate(reference) if h and y>=0}
coords=lambda S:{(j,z) for v in S for j,z in enumerate(R[v])}

trades=[]
for line in Path(sys.argv[2]).read_text().splitlines():
    if line.startswith('#'):continue
    v=list(map(int,line.split()))
    assert len(v)==16
    t=tuple(tuple(zip(v[j:j+8:2],v[j+1:j+8:2])) for j in (0,8))
    assert all(tuple(sorted(S))==S and len(coords(S))==16 for S in t)
    L,U=t
    assert L<U and not set(L)&set(U) and coords(L)==coords(U)
    def permutations():
        for g in range(4):
            yield tuple(next(j for j,b in enumerate(U) if R[a][g]==R[b][g]) for a in L)
    P=list(permutations())
    assert P[0]==(0,1,2,3)
    assert all(all(i!=p[i] and p[p[i]]==i for i in range(4)) for p in P[1:])
    # Boundary colors are the three coordinate directions and their xor.
    # Embed even vertices by the Klein permutations (identity and 3 nonzero directions).
    D=(1,2,4,7)
    even=[None]*4
    for p,d in zip(P,D):even[p[0]]=d^1
    odd=[u^1 for u in even]
    assert sorted(even+odd)==list(range(8))
    assert all(not x.bit_count()%2 for x in even)
    assert all(odd[P[g][i]]==even[i]^D[g] for i in range(4) for g in range(4))
    assert [sum((-1)**x.bit_count()*x**j for x in range(8)) for j in range(4)]==[0,0,0,-48]
    trades.append(t)
assert len(set(trades))==len(trades)==83
assert Counter(len((set(L)|set(U))&ref) for L,U in trades)=={0:60,1:23}

pivot={}
index={v:i for i,v in enumerate(R)}
for L,U in trades:
    w=sum(1<<index[v] for v in L+U)
    while w and w.bit_length() in pivot:w^=pivot[w.bit_length()]
    assert w
    pivot[w.bit_length()]=w
C=Counter(len(set(S)&set(T)) for a,b in combinations(trades,2) for S,T in product(a,b)
          if len(coords(set(S)|set(T)))==4*len(set(S)|set(T)))
print('Compatible side-pairs by overlap:',C)
print('Rank of trade supports over F2:',len(pivot))

# The numerical sign convention is the ordinary Thue–Morse sign (-1)^popcount.
first=trades[0]
print('First trade:',[[(BRANCHES[h],y) for h,y in S] for S in first])

assert C == {0:4856}
# Disjoint boundaries imply that every side of one move is compatible
# with every side of the other. Find the largest such trade family.
Q=[coords(L) for L,U in trades]
N=[sum(1<<j for j,b in enumerate(Q) if not a&b) for a in Q]
def unpack(P):
    while P:
        bit=P&-P
        yield bit.bit_length()-1
        P^=bit
best=()
def expand(chosen,P):
    global best
    if len(chosen)>len(best):best=tuple(chosen)
    # Greedy proper coloring of the compatibility graph: distinct
    # colors are required for every vertex of a clique.
    order=[]
    U=P
    c=0
    while U:
        T=U
        c+=1
        while T:
            b=T&-T
            i=b.bit_length()-1
            order.append((i,c))
            U^=b
            T&=~(b|N[i])
    for i,c in reversed(order):
        if len(chosen)+c<=len(best):break
        expand(chosen+[i],P&N[i])
        P^=1<<i
expand([],(1<<len(N))-1)
print('Maximum cube dimension:',len(best),'zero-based trade IDs:',best)
assert all(not Q[i]&Q[j] for i,j in combinations(best,2))
assert len(best)==6


# Realize all 64 points of the extremal cube in the actual graph.
chart,stars=orbit_chart("period_fourteen_site.txt",4,((10,45,2),(27,49,0)))
base=graph(chart,stars)
ends=set()
for choice in product((0,1),repeat=6):
    A=[R.copy() for R in base]
    v=leaf(2,7)
    assert legal(A,SOURCE,v)
    A[SOURCE].add(v)
    A[v].add(SOURCE)
    selected=tuple(sorted(q for j,b in zip(best,choice) for q in trades[j][b]))
    ends.add(selected)
    for h,y in selected:
        v=leaf(BRANCHES[h],y)
        assert legal(A,SOURCE,v)
        A[SOURCE].add(v)
        A[v].add(SOURCE)
    assert len(A[SOURCE])==29
    assert max(map(len,A))==57
assert len(ends)==64
print("All 64 patches independently embedded on 3,250 vertices.")
