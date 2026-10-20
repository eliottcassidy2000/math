#!/usr/bin/env python3
"""Separate exact implementation. Imports no producer or inherited code."""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import gcd, factorial
from pathlib import Path
from collections import Counter
import json

checks=0
def require(ok, message):
    global checks
    checks+=1
    if not ok:raise RuntimeError(message)

def farey():
    totals=[]
    for d in range(1,201):
        residues=[]
        primes=[];rem=d
        for p in range(2,d+1):
            if rem%p==0:
                primes.append(p)
                while rem%p==0:rem//=p
        for a in range(1,501):
            g=gcd(a,a+d);h=gcd(a+1,a+d+1)
            det=abs((a//g)*((a+d+1)//h)-((a+d)//g)*((a+1)//h))
            require(det*g*h==d,'reduced determinant retains content')
            require((det==1)==(a*(a+1)%d==0),'Farey adjacency divisibility iff')
            if a<=d and det==1:residues.append(a%d)
        require(len(residues)==2**len(primes),'prime-power CRT event count')
        if d in (1,6,9,15,21,27,105):totals.append({'difference':d,'residues':sorted(residues)})
    return {'rectangle':{'difference':[1,200],'a':[1,500]},'examples':totals}

def unit_graph():
    rim=[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
    pts=rim+[tuple(rim[2*i][k]+rim[2*i+1][k] for k in range(2)) for i in range(3)]+[(0,0)]
    edges=[]
    for i,j in combinations(range(10),2):
        a,b=pts[i][0]-pts[j][0],pts[i][1]-pts[j][1]
        if a*a+a*b+b*b==1:edges.append((i,j))
    require(len(edges)==18,'independent lattice norm reconstructs graph')
    ports=[(i,j) for i,j in edges if j<9]
    hist=Counter();projection=set()
    # Vertex order is six rim, then three tips. Fix wing zero to 0,1,2.
    for rest in product(range(6),repeat=6):
        colours=(0,1)+rest[:4]+(2,)+rest[4:]
        if any(colours[i]==colours[j] for i,j in ports):continue
        allowed=[c for c in range(6) if all(colours[i]!=c for i,j in edges if j==9)]
        hist[len(allowed)]+=1
        if allowed:projection.add((colours[2],colours[3],colours[7]))
    require(hist=={0:384,1:3456,2:3840,3:640,4:16},'full six-colour relation')
    require(len(projection)==100,'full two-wing projection')
    require((5**6+5)*4**3==sum(hist.values())*120,'cycle polynomial count')
    require(factorial(6)*4**3==hist[0]*120,'rainbow count')
    # All graph edges also satisfy a global three-colouring, exposing the carrier limit.
    require(all((pts[i][0]+2*pts[i][1]-pts[j][0]-2*pts[j][1])%3 for i,j in edges),
            'common triangular lattice stays three-colourable')
    return {'histogram':dict(sorted(hist.items())),'two_wing_projection':len(projection)}

def partitions(n,lower=1):
    if n==0:yield ()
    for a in range(lower,n+1):
        for tail in partitions(n-a,a):yield (a,)+tail

def partition_counts(n):
    ordinary=F(0);twisted=F(0)
    for parts in partitions(n):
        mult=Counter(parts);z=1
        for a,k in mult.items():z*=a**k*factorial(k)
        e=sum(a//2 for a in parts)+sum(gcd(a,b) for a,b in combinations(parts,2))
        if all(a%2 for a in parts):ordinary+=F(2**e,z)
        if mult[1]<=1 and all(a==1 or a%4==2 for a in parts):twisted+=F(2**e,z)
    require(ordinary.denominator==twisted.denominator==1,'Burnside integrality')
    return int(ordinary),int(twisted)

def literal_fixed(perm,twist):
    pairs=list(combinations(range(len(perm)),2));unseen=set(pairs);cycles=0
    while unseen:
        start=min(unseen);pair=start;parity=0
        while pair in unseen:
            unseen.remove(pair)
            a,b=perm[pair[0]],perm[pair[1]]
            parity ^= (a>b)^twist
            pair=tuple(sorted((a,b)))
        require(pair==start,'permutation edge orbit closes')
        if parity:return 0
        cycles+=1
    return 2**cycles

def merged_counts():
    rows=[]
    for n in range(10):
        A,S=partition_counts(n)
        if n<=6:
            for twist,expected in ((0,A),(1,S)):
                total=sum(literal_fixed(p,twist) for p in permutations(range(n)))
                require(total==factorial(n)*expected,'literal signed-edge action equals partition count')
        rows.append({'n':n,'ordinary':A,'self_converse':S,'merged':(A+S)//2})
    require([r['merged'] for r in rows[3:]]==[2,3,10,34,272,3528,97144],
            'incoming realizability count table')
    def edge_exponent(parts):
        return sum(a//2 for a in parts)+sum(gcd(a,b) for a,b in combinations(parts,2))
    for mass in range(11):
        for parts in partitions(mass):
            if not all(a%2 for a in parts):continue
            for e,c,k in product((0,1),(1,3,5),(0,1,2,3)):
                old=tuple(2*a for a in parts)+(1,)*e
                new=old+(2*c,)*k
                predicted=k*c+2*c*k*(k-1)//2+2*k*sum(gcd(c,a) for a in parts)+e*k
                require(edge_exponent(new)-edge_exponent(old)==predicted,
                        'twice-odd cycle transition including fixed-point channel')
    return rows

def path_to_exit(n):
    value=n;word=[]
    while value>=n:
        value=3*value+1;exponent=0
        while value%2==0:value//=2;exponent+=1
        word.append(exponent)
        if len(word)>1000:raise RuntimeError('first-exit cap')
    return word,value

def affine(segment):
    A=F(1);C=F(0)
    for a in sorted(segment,reverse=True):A=3*A/2**a;C=(3*C+1)/2**a
    return A,C

def block_scan(receipt):
    hard={r['source']:r for r in receipt['collatz']['hard_sources']}
    histogram=Counter();examined=0
    for n in range(3,200002,2):
        word,end=path_to_exit(n);m=len(word)
        A,C=affine(word)
        if A*n+C<n:
            require(n not in hard,'all one-block sources independently identified')
            histogram[1]+=1;continue
        require(n in hard,'producer retains every hard source')
        wanted=hard[n];maps={(i,j):affine(word[i:j]) for i in range(m) for j in range(i+1,m+1)}
        for b in range(2,wanted['minimum_blocks']+1):
            best=None
            for cuts in combinations(range(1,m),b-1):
                x=F(n);start=0
                for stop in cuts+(m,):
                    a,c=maps[start,stop];x=a*x+c;start=stop
                best=x if best is None or x<best else best;examined+=1
            if b<wanted['minimum_blocks']:require(best>=n,'every smaller block budget fails')
            else:require(best<n and best==F(wanted['upper']),'exhaustive optimum equals integer DP')
        histogram[wanted['minimum_blocks']]+=1
    require(dict(histogram)=={int(k):v for k,v in receipt['collatz']['histogram'].items()},
            'full 100000-source minimum-block histogram')
    return {'sources':100000,'histogram':dict(sorted(histogram.items())),
            'all_hard_sources':len(hard),'exhaustive_partitions':examined}

def positive_charge():
    best=F(-1);record=None;count=0
    for n in range(3,1000002,2):
        word,end=path_to_exit(n);m=len(word);K=sum(word);prefix=0
        positive=0;dividend=0
        for i,a in enumerate(word):
            suffix_power=3**(m-i-1)
            if a==1:positive+=suffix_power*(1<<(prefix+a))
            if a>=3:dividend+=((1<<a)-4)*suffix_power*(1<<prefix)
            prefix+=a
        upper=3**m*(n-1)+positive
        require(upper-dividend==(end-1)*(1<<K),'independent positive-minus-dividend identity')
        require(upper<(n-1)*(1<<K),'positive charge paid at declared finite first exits')
        ratio=F(dividend,(n-end)*(1<<K));count+=1
        if ratio>best:best=ratio;record=n
    require(best==F(215402584767891,281474976710656) and record==63,
            'exact worst dividend fraction in the declared universe')
    return {'sources':count,'worst_source':record,'worst_fraction':str(best)}

def main():
    path=Path(__file__).with_name('receipt.json')
    receipt=json.loads(path.read_text())
    result={'status':'FINITE-EXACT; separate implementation, same agent',
            'farey':farey(),'geometry':unit_graph(),'merged_counts':merged_counts(),
            'collatz':block_scan(receipt),'positive_charge':positive_charge()}
    result['checks']=checks
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
