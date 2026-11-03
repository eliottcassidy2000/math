"""Exact controls for guard Galois adjunction and dihedral observation modules."""
from itertools import combinations, product, permutations
from math import lcm, factorial
from collections import Counter
from fractions import Fraction
import json

def require(condition, message):
    if not condition:
        raise ValueError(message)

def observe(N, ds):
    return sum((N % d == 0) << i for i,d in enumerate(ds))

def join(mask, ds):
    return lcm(*(d for i,d in enumerate(ds) if mask >> i & 1))

def closed_states(ds):
    return {mask for mask in range(1<<len(ds)) if observe(join(mask,ds),ds)==mask}

def bichromatic(ds):
    all_bits=(1<<len(ds))-1
    C=closed_states(ds)
    return {mask for mask in C if all_bits^mask in C}

def interval_masks(m):
    return {0} | {sum(1<<i for i in range(l,r+1)) for l in range(m) for r in range(l,m)}

def cycle_labels(m):
    # Empty, increasing prefixes, full, decreasing suffixes, then empty.
    full=(1<<m)-1
    return [(1<<j)-1 for j in range(m+1)] + [full^((1<<j)-1) for j in range(1,m)]

def neighbor_mask(x,y):
    return (x^y).bit_count()==1

def closed_walk_to_complement(m, ell):
    states=cycle_labels(m); full=(1<<m)-1
    total=0
    for start in states:
        counts={start:1}
        for _ in range(ell):
            new=Counter()
            for x,c in counts.items():
                for y in states:
                    if neighbor_mask(x,y):new[y]+=c
            counts=new
        total+=counts.get(full^start,0)
    return total

def comb(n,k):
    return factorial(n)//factorial(k)//factorial(n-k)

for ds in [(4,6,9),(6,10,15),(2,4),(2,3,5),(143,85,19,17,69,29,23,62,37,31,287,43,41)]:
    m=len(ds)
    # Every subset satisfies the abstract Galois condition on a complete divisor grid.
    L=lcm(*ds)
    if L<10000:
        for S in range(1<<m):
            for N in range(1,L+1):
                require((N%join(S,ds)==0)==(S&observe(N,ds)==S),'adjunction')
        direct={observe(N,ds) for N in range(1,L+1)}
        require(direct==closed_states(ds),'arithmetic realizability')
    print('denominators',ds,'closed',len(closed_states(ds)),'biclosed',len(bichromatic(ds)))

for m in range(2,9):
    ds=[2**(m-1-i)*3**i for i in range(m)]
    C=closed_states(ds);B=bichromatic(ds);cyc=cycle_labels(m)
    require(C==interval_masks(m),'interval closure')
    require(B==set(cyc),'biclosed cycle vertices')
    require(len(cyc)==len(set(cyc))==2*m,'cycle labels distinct')
    for j,x in enumerate(cyc):
        require([k for k,y in enumerate(cyc) if neighbor_mask(x,y)]==sorted(((j-1)%(2*m),(j+1)%(2*m))),'cycle edges')
        require(x^((1<<m)-1)==cyc[(j+m)%(2*m)],'antipode')
    # Independent integer walk counts vs character trace's binomial form.
    for ell in range(2*m+2):
        rhs=2*m*sum(comb(ell,r) for r in range(ell+1) if (2*r-ell-m)%(2*m)==0)
        require(closed_walk_to_complement(m,ell)==rhs,'antipodal trace identity')
    print('rank-two module',m,'closed',len(C),'cycle',len(B),'minimum antipodal walks',closed_walk_to_complement(m,m))

# Product of disjoint-prime rank-two modules m=3,2.
ds=(4,6,9,5,7)
B=bichromatic(ds);B1=bichromatic(ds[:3]);B2=bichromatic(ds[3:])
require(B=={x|(y<<3) for x in B1 for y in B2},'prime-supported product')
full=(1<<len(ds))-1
for start in B:
    counts={start:1}
    for _ in range(5):
        new=Counter()
        for x,c in counts.items():
            for y in B:
                if neighbor_mask(x,y):new[y]+=c
        counts=new
    require(counts[full^start]==2**2*factorial(5)//(factorial(3)*factorial(2)),'product shortest paths')
print('product C6 x C4:',len(B),'states; shortest antipodal paths each = 40')

# Actual blue tilings, all n=5 states from (4,6,9).
n=5
pairs=list(combinations(range(n),2))
chords=[e for e in pairs if e[1]>e[0]+1]
reflection=lambda e:(n-1-e[1],n-1-e[0])
orbits=[];seen=set()
for e in chords:
    if e not in seen:
        orb={e,reflection(e)};orbits.append(orb);seen|=orb
assignment=(0,1,2,0)
require(len(orbits)==len(assignment),'orbit count')

def tournament(e):
    reverse={edge for i,O in enumerate(orbits) for edge in O if e>>assignment[i]&1}
    A=[[0]*n for _ in range(n)]
    for u,v in pairs:
        A[u][v]=int((u,v) not in reverse); A[v][u]=1-A[u][v]
    return A

def H(A):
    n=len(A);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for u in range(n):
            if mask>>u&1:
                for v in range(n):
                    if not mask>>v&1 and A[u][v]: dp[mask|1<<v][v]+=dp[mask][u]
    return sum(dp[-1])

def canon(A):
    n=len(A)
    return min(sum((1-A[p[i]][p[j]])<<k for k,(i,j) in enumerate(pairs)) for p in permutations(range(n)))
rows=[]
for e in cycle_labels(3):
    A=tournament(e);B=tournament(e^7)
    require(all(A[i][i+1] for i in range(n-1)),'fixed path')
    require(all(A[i][j]==A[n-1-j][n-1-i] for i,j in pairs),'actual blue fixed point')
    require(all(A[i][j]==B[i][j] if j==i+1 else A[i][j]!=B[i][j] for i,j in pairs),'actual chord complement')
    h=H(A)
    brute=sum(all(A[p[i]][p[i+1]] for i in range(n-1)) for p in permutations(range(n)))
    require(h==brute,'independent H')
    row={'guards':''.join(str(e>>i&1) for i in range(3)), 'least_integer':join(e,(4,6,9)), 'H':h,'canonical':canon(A)}
    rows.append(row)
print('blue embedding',json.dumps(rows))
require(len({canon(tournament(e)) for e in cycle_labels(3)}) < 6,'quotient collision control')
# Same observation, different successor in the stated ordered reduced program.
ds=(4,6,9)
require(observe(4,ds)==observe(16,ds)==1,'same old state')
require(observe(1,ds)==0 and observe(4,ds)==1,'different new state')
print('4 and 16 have the same guard observation but first rule 1/4 sends them to distinct observations')
print('ALL CHECKS PASS')
