#!/usr/bin/env python3
"""Exact controls for the projected-template cover and residue lemmas."""
from collections import Counter
from functools import lru_cache

def arc(i,j):
    return (j-i) % 5 in (1,2)

def template_counts(b):
    fixed=[j for j in range(1,5) if j!=b]
    visits=Counter()
    def recur(last, seen, a, c):
        if len(seen)==3 and a and c:
            visits[a,c]+=1
        for j in range(5):
            if last is not None and not arc(last,j):
                continue
            if j in fixed:
                if j not in seen:
                    recur(j,seen|{j},a,c)
            elif j==0:
                if a<4:
                    recur(j,seen,a+1,c)
            elif c<4:
                recur(j,seen,a,c+1)
    recur(None,set(),0,0)
    assert all(a<=3 and c<=3 for a,c in visits)
    return tuple(tuple(visits[a,c] for c in range(1,4)) for a in range(1,4))

@lru_cache(None)
def exact(lengths):
    label=[i for i,m in enumerate(lengths) for _ in range(m)]
    n=len(label)
    edges=[sum(1<<j for j in range(n)
               if (i<j if label[i]==label[j] else arc(label[i],label[j])))
           for i in range(n)]
    D=[[0]*n for _ in range(1<<n)]
    for i in range(n):
        D[1<<i][i]=1
    for S in range(1,1<<n):
        for i,q in enumerate(D[S]):
            if not q:
                continue
            T=edges[i]&~S
            while T:
                bit=T&(-T)
                T-=bit
                D[S|bit][bit.bit_length()-1]+=q
    return sum(D[-1])

def main():
    assert exact((5,1,1,1,1))==705
    assert template_counts(1)==((15,13,2),(13,14,4),(2,4,2))
    assert template_counts(2)==((15,13,2),(13,21,8),(2,8,6))
    seeds=[
        ["11222","11424","11123","12122","11114"],
        ["11111","12123","13333","11232","11133"],
        ["22222","11112","11134","11224","11242"],
        ["11113","11212","11333","11122","11324"]
    ]
    for residue4,row in enumerate(seeds):
        for residue5,seed in enumerate(row):
            sizes=tuple(map(int,seed))
            assert sum(sizes)%4==residue4
            assert exact(sizes)%5==residue5
    for a in range(1,9):
        assert exact((a,1,1,1,1))==2*3**a+7*2**a-5
    print("Exact template matrices, H=705 witness, 20 fixed-order residue seeds, and eight one-variable controls passed.")

from itertools import permutations
from collections import defaultdict, Counter

def degree_four(n):
    forest = {}
    def add(paths):
        arcs = [(a,b) for p in paths for a,b in zip(p,p[1:])]
        E = tuple(sorted(tuple(sorted(e)) for e in arcs))
        c = (-1)**sum(a>b for a,b in arcs) * 2**len(paths)
        if E in forest:
            assert forest[E] == c
        forest[E] = c
    for p in permutations(range(n),5):
        add((p,))
    for p in permutations(range(n),6):
        add((p[:3],p[3:]))
    tri = {(0,1):1,(1,2):1,(0,2):-1}
    terms = defaultdict(int)
    for E,c in forest.items():
        internal = [e for e in E if e in tri]
        if len(internal) == 1:
            terms[tuple(e for e in E if e not in tri)] += 2*c*tri[internal[0]]
    return Counter(abs(c) for E,c in terms.items() if c
        and {v for e in E for v in e if v>=3} == set(range(3,n)))

want = [{},{},{4:12,8:6},{4:36,8:18},{}]
for k,h in enumerate(want):
    assert degree_four(3+k) == h
print("Exact degree-four numerator: 576*binom(n-3,2)+1728*binom(n-3,3)")

def hp(out):
    n=len(out)
    D=[[0]*n for _ in range(1<<n)]
    for v in range(n):
        D[1<<v][v]=1
    for S in range(1,1<<n):
        for i,c in enumerate(D[S]):
            if not c: continue
            T=out[i]&~S
            while T:
                b=T&-T
                T-=b
                D[S|b][b.bit_length()-1]+=c
    return sum(D[-1])

def diagonal():
    from math import isqrt
    out=[0]*14
    for u in range(14):
        for v in range(u+1,14):
            d=(v%7-u%7)%7
            forward=(d in (1,2,4) if u//7==v//7 else d in (0,1))
            if forward: out[u]|=1<<v
            else: out[v]|=1<<u
    q=2163467
    assert hp(out)==7*q and all(q%d for d in range(2,isqrt(q)+1))

def dyadic():
    from itertools import product
    expected=[(1,0,0,0),(1,0,0,0),(1,1,0,0),(1,2,5,0),
              (19,14,23,8),(284,168,372,200),(9136,8382,8216,7034)]
    A=[{1}]
    for n,want in enumerate(expected,1):
        free=[(i,j) for i in range(n) for j in range(i+2,n)]
        freq=[0]*4
        for bits in product((0,1),repeat=len(free)):
            out=[(1<<(i+1) if i+1<n else 0) for i in range(n)]
            for (i,j),b in zip(free,bits):
                if b: out[i]|=1<<j
                else: out[j]|=1<<i
            freq[(hp(out)//2)%4]+=1
        assert tuple(freq)==want
        A.append({2*i+1 for i,z in enumerate(freq) if z%2})
    for mod in (4,8):
        a=[set() for _ in A]
        for n,row in enumerate(A):
            for r in row: a[n]^={r%mod}
        want=({1},{1},set(),{3},{1},{1,3}) if mod==4 else (
            {5},{5},set(),{3},{5},{3,5})
        S=[set()]
        for n in range(1,31):
            cur=set(a[n]) if n<len(a) else set()
            for j in range(1,min(n,len(a))):
                for r in a[j]:
                    for t in S[n-j]:
                        cur^={r*t%mod}
            S.append(cur)
            assert cur==({1} if n==1 else want[n%6])
    print("Direct fixed-path cubes through n=7 and both strong residue cycles passed.")

if __name__ == "__main__":
    main()
    diagonal()
    dyadic()
