#!/usr/bin/env python3
"""No repository imports. Literal reachability and mirror tests through n=7."""
from pathlib import Path

def formulas(limit):
    total=[1]+[1<<((n-1)*(n-2)//2) for n in range(1,limit+1)]
    blue=[1]+[1<<((n-1)**2//4) for n in range(1,limit+1)]
    strong=[0]*(limit+1)
    strongblue=[0]*(limit+1)
    for n in range(1,limit+1):
        strong[n]=total[n]-sum(strong[k]*total[n-k] for k in range(1,n))
        strongblue[n]=blue[n]-sum(strong[k]*blue[n-2*k] for k in range(1,n//2+1))
    return total,blue,strong,strongblue

def literal(n):
    tiles=[(x,y) for x in range(n) for y in range(x-1)]
    index={edge:i for i,edge in enumerate(tiles)}
    mirror=[index[n-1-y,n-1-x] for x,y in tiles]
    full=(1<<n)-1
    def reaches_all(adj):
        seen=1
        while True:
            new=seen
            for v in range(n):
                if seen>>v&1: new|=adj[v]
            if new==seen: return seen==full
            seen=new
    strengths=[]
    blues=[]
    for mask in range(1<<len(tiles)):
        adj=[0]*n
        rev=[0]*n
        for x in range(1,n):
            adj[x]|=1<<(x-1)
            rev[x-1]|=1<<x
        for i,(x,y) in enumerate(tiles):
            u,v=(x,y) if mask>>i&1 else (y,x)
            adj[u]|=1<<v
            rev[v]|=1<<u
        strengths.append(reaches_all(adj) and reaches_all(rev))
        blues.append(all(((mask>>i)^(mask>>j))&1==0 for i,j in enumerate(mirror)))
    maximum=len(strengths)-1
    cut=core=bluecut=bluecore=periphery=0
    for t in range((maximum+1)//2):
        u=t^maximum
        assert blues[t]==blues[u]
        if strengths[t] != strengths[u]:
            cut+=1; bluecut+=blues[t]
        elif strengths[t]:
            core+=1; bluecore+=blues[t]
        else: periphery+=1
    assert periphery==0
    return (sum(strengths),sum(s and b for s,b in zip(strengths,blues)),cut,bluecut,core,bluecore)

if __name__=='__main__':
    total,blue,strong,strongblue=formulas(12)
    print('n total blue strong strongblue cut bluecut core bluecore blackcut blackcore')
    for n in range(3,13):
        s,b=strong[n],strongblue[n]
        cut=total[n]-s; bluecut=blue[n]-b
        core=s-total[n]//2; bluecore=b-blue[n]//2
        expected=(s,b,cut,bluecut,core,bluecore)
        if n<=7: assert literal(n)==expected,(n,literal(n),expected)
        assert (cut-bluecut)%2==0 and (core-bluecore)%2==0
        print(n,total[n],blue[n],*expected,cut-bluecut,core-bluecore)
    print('PASS: independent bitset reachability and coordinate-mirror enumeration for every tiling n3..7.')
    print('n8..12 are recurrence-derived values, not direct enumeration in this script.')
