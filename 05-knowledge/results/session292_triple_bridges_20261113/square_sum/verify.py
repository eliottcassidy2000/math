#!/usr/bin/env python3
"""Exact guards and independent literal graph controls; no third-party packages."""
from itertools import combinations, product
from math import isqrt
import json
from pathlib import Path


def require(ok, message):
    if not ok: raise RuntimeError(message)


def direct(word, n, cycle=False):
    if cycle and len(word)<3: return []
    out=[]
    for x in range(1,n+1):
        vs=[x]
        for a in word: vs.append(a-vs[-1])
        tested=vs[:-1] if cycle else vs
        if all(1<=v<=n for v in vs) and len(set(tested))==len(tested):
            if not cycle or vs[-1]==x: out.append(x)
    return out


def compiled(word,n,cycle=False):
    if cycle and len(word)<3: return []
    cs=[0]
    for a in word: cs.append(a-cs[-1])
    lo,hi=1,n
    for i,c in enumerate(cs):
        low,high=(c-n,c-1) if i%2 else (1-c,n-c)
        lo,hi=max(lo,low),min(hi,high)
    end=len(word)
    if cycle:
        if end%2:
            if cs[-1]%2: return []
            lo,hi=max(lo,cs[-1]//2),min(hi,cs[-1]//2)
        elif cs[-1]:return []
    blocked=set()
    for i,j in combinations(range(end if cycle else end+1),2):
        if i%2==j%2:
            if cs[i]==cs[j]:return []
        else:
            # ((-1)^i-(-1)^j)x=c_j-c_i
            d=(1 if i%2==0 else -1)-(1 if j%2==0 else -1)
            if (cs[j]-cs[i])%d==0: blocked.add((cs[j]-cs[i])//d)
    return [x for x in range(lo,hi+1) if x not in blocked]


def hex_count(a,b,c,n):
    p,q=b-a,c-b
    parity=(a+q)%2
    radius=min(a+b-c-2,2*n-b-c+a)
    if radius<0:return 0
    return (radius+parity)//2-sum(d<=radius and d%2==parity for d in {p,q,p+q})


def canonical(vs):
    return min(tuple(v[i:]+v[:i]) for v in [vs,list(reversed(vs))] for i in range(len(v)))


def hex_cycles(a,b,c,n):
    cycles=set()
    for x in direct([a,b,c,a,b,c],n,True):
        vs=[x]
        for s in [a,b,c,a,b]:vs.append(s-vs[-1])
        cycles.add(canonical(vs))
    return sorted(cycles)


def graph(n):
    e=[(x,y) for x in range(1,n+1) for y in range(x+1,n+1) if isqrt(x+y)**2==x+y]
    adj={x:set() for x in range(1,n+1)}
    for x,y in e:adj[x].add(y);adj[y].add(x)
    seen=set();components=[]
    for x in adj:
        if x in seen:continue
        stack=[x];seen.add(x);comp=[]
        while stack:
            y=stack.pop();comp.append(y)
            for z in adj[y]-seen:seen.add(z);stack.append(z)
        components.append(sorted(comp))
    squares=[s*s for s in range(2,isqrt(2*n-1)+1)]
    hexes=[{"colours":[a,b,c],"cycles":hex_cycles(a,b,c,n)} for a,b,c in combinations(squares,3) if hex_count(a,b,c,n)]
    triangles=[list(v) for v in combinations(range(1,n+1),3) if all(y in adj[x] for x,y in combinations(v,2))]
    return {"n":n,"edges":len(e),"components":components,"cycle_rank":len(e)-n+len(components),"leaves":[x for x in adj if len(adj[x])==1],"hexagons_ABCABC":hexes,"triangles":triangles}


def main():
    for word,n in [([4],2),([4,4],3)]:
        require(direct(word,n,True)==compiled(word,n,True)==[], "short-cycle hostile")
    word_tests=0
    for length in range(1,7):
        for word in product([4,9,16],repeat=length):
            for n in [1,3,7,12,16,23,25]:
                for cycle in [False,True]:
                    require(compiled(word,n,cycle)==direct(word,n,cycle),(word,n,cycle))
                    word_tests+=1
    hex_tests=0
    for a,b,c in combinations(range(1,26),3):
        for n in range(1,31):
            require(hex_count(a,b,c,n)==len(hex_cycles(a,b,c,n)),(a,b,c,n))
            hex_tests+=1
    births=[]
    for t in range(4,101):
        a,b,c=t*t,(t+1)**2,(t+2)**2
        n=(t*t+6*t+5)//2+1
        require(hex_count(a,b,c,n-1)==0 and hex_count(a,b,c,n)==1,t)
        if t<=15:
            require(len(hex_cycles(a,b,c,n))==1,t)
            births.append({"t":t,"n":n,"cycle":hex_cycles(a,b,c,n)[0]})
    atlas=[graph(n) for n in range(1,26)]
    require(not any(row['hexagons_ABCABC'] for row in atlas[:22]),'hex early')
    require(atlas[22]['hexagons_ABCABC'][0]['cycles']==[(2,14,22,3,13,23)],'first hex canonical')
    require(graph(29)['triangles']==[] and graph(30)['triangles']==[[6,19,30]],'triangle boundary')
    out={"status":"PASS", "word_tests":word_tests,"hex_tests":hex_tests,"consecutive_square_birth_controls":97,"first25":atlas,"consecutive_square_examples":births,"scope":"All-word and all-N proofs in PROOF.md; finite controls independently compare direct traversal. No Hamiltonicity or Chan-Pak closure."}
    Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['first25','consecutive_square_examples']},indent=2))

if __name__=='__main__':main()
