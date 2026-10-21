#!/usr/bin/env python3
"""Independent signed-equation and literal marked-tournament audit.

No producer imports. Actual permutation equations and generator orbit search
are compared with a separately coded multiplicity-vector channel sum.
"""
import argparse
import hashlib
import itertools
import json
import math
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

CHECKS = 0

def require(condition, description):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(description)

@lru_cache(None)
def channel(size, a, b):
    def visit(length, left, cycles):
        if left == 0:
            counts = Counter(cycles)
            denominator = math.prod(l**c * math.factorial(c) for l,c in counts.items())
            exponent = sum((l-1)//2 for l in cycles)
            exponent += sum(math.gcd(x,y) for i,x in enumerate(cycles) for y in cycles[:i])
            return Fraction(2**(a*exponent+b*len(cycles)), denominator)
        if length > left:
            return Fraction(0)
        return sum((visit(length+2,left-c*length,cycles+(length,)*c)
                    for c in range(left//length+1)), Fraction(0))
    return visit(1,size,())

def expected(n,k):
    ordinary = 2**((k-1)*(k-2)//2)*channel(n-k,1,k)
    r,delta = divmod(k,2)
    s,epsilon = divmod(n-k,2)
    converse = Fraction(0) if delta*epsilon else 2**(r*(r-1+delta+epsilon))*channel(s,2,k+epsilon)
    require(ordinary.denominator == converse.denominator == 1, f"integrality n={n},k={k}")
    return int(ordinary),int(converse)

def fixed_by_equations(n,k,permutation,anti):
    """Solve x_target XOR x_source = endpoint_flip XOR anti, then path pins."""
    edges=list(itertools.combinations(range(n),2))
    index={edge:i for i,edge in enumerate(edges)}
    anchor=len(edges)
    parent=list(range(anchor+1))
    parity=[0]*(anchor+1)
    def root(v):
        value=0
        while parent[v]!=v:
            value ^= parity[v]
            v=parent[v]
        return v,value
    def join(u,v,change):
        a,x=root(u);b,y=root(v)
        if a==b:
            return x^y==change
        parent[a]=b
        parity[a]=x^y^change
        return True
    for j,(u,v) in enumerate(edges):
        x,y=permutation[u],permutation[v]
        target=index[min(x,y),max(x,y)]
        if not join(j,target,int(x>y)^anti):
            return 0
    for v in range(k-1):
        if not join(index[v,v+1],anchor,1):
            return 0
    free=len({root(v)[0] for v in range(anchor+1)})-1
    return 2**free

def direct_burnside(bound):
    rows=[]
    for n in range(1,bound+1):
        for k in range(1,n+1):
            total=[0,0]
            calls=0
            for tail in itertools.permutations(range(k,n)):
                calls+=1
                for anti in (0,1):
                    head=tuple(range(k-1,-1,-1)) if anti else tuple(range(k))
                    total[anti]+=fixed_by_equations(n,k,head+tail,anti)
            got=tuple(Fraction(v,math.factorial(n-k)) for v in total)
            target=expected(n,k)
            require(got==target,f"signed Burnside mismatch {(n,k,got,target)}")
            rows.append({"n":n,"k":k,"ordinary":int(got[0]),"converse":int(got[1]),"permutations":calls})
    return rows

def involution_transform(n,p):
    """Pull back all ordered orientation queries under an involution p."""
    edges=list(itertools.combinations(range(n),2))
    index={edge:i for i,edge in enumerate(edges)}
    sources=[]
    bias=0
    for target,(u,v) in enumerate(edges):
        x,y=p[u],p[v]
        source=index[min(x,y),max(x,y)]
        sources.append((source,1<<target))
        if x>y:
            bias|=1<<target
    # Independent byte lookup construction from pullback query sources.
    tables=[]
    for offset in range(0,len(edges),8):
        table=[]
        for value in range(256):
            out=0
            for source,targetbit in sources:
                if offset<=source<offset+8 and (value>>(source-offset))&1:
                    out|=targetbit
            table.append(out)
        tables.append((offset,table))
    def apply(value):
        out=bias
        for offset,table in tables:
            out^=table[(value>>offset)&255]
        return out
    return apply

def literal_census(n,k):
    edges=list(itertools.combinations(range(n),2))
    pinned=sum(1<<j for j,(u,v) in enumerate(edges) if v==u+1 and v<k)
    full=(1<<len(edges))-1
    swaps=[]
    for v in range(k,n-1):
        p=list(range(n));p[v],p[v+1]=p[v+1],p[v]
        swaps.append(involution_transform(n,p))
    reversal=involution_transform(n,tuple(range(k-1,-1,-1))+tuple(range(k,n)))
    visited=bytearray(full+1)
    ordinary=converse=0
    histogram=Counter()
    for mask in range(full+1):
        if mask&pinned!=pinned or visited[mask]:
            continue
        orbit={mask}
        stack=[mask]
        visited[mask]=1
        while stack:
            value=stack.pop()
            for action in swaps:
                image=action(value)
                if image not in orbit:
                    orbit.add(image)
                    stack.append(image)
                    visited[image]=1
        require(all(value&pinned==pinned for value in orbit),"literal pin preservation")
        require(math.factorial(n-k)%len(orbit)==0,"orbit-stabilizer divisibility")
        ordinary+=1
        converse+=reversal(full^mask) in orbit
        histogram[len(orbit)]+=1
    target=expected(n,k)
    require((ordinary,converse)==target,f"literal mismatch {(n,k,ordinary,converse,target)}")
    return {"n":n,"k":k,"ordinary":ordinary,"converse":converse,"orbit_sizes":dict(sorted(histogram.items()))}

def correction_controls():
    # n=8, lambda=(5,3): four unordered pair orbits, no fixed vertices.
    true_layer=Fraction(8*2**4,5*3)
    printed_formula=Fraction(8*2**4,5*3*math.factorial(2))
    require(true_layer==Fraction(128,15),"mixed-cycle correction")
    require(printed_formula==Fraction(64,15),"printed formula hostile")
    rooted=[int(channel(n-1,1,1)) for n in range(1,9)]
    unrooted=[int(channel(n,1,0)) for n in range(1,9)]
    deficit=[n*unrooted[n-1]-rooted[n-1] for n in range(1,9)]
    require(rooted==[1,2,4,12,48,296,3040,54256],"THM260 rooted prefix")
    require(deficit==[0,0,2,4,12,40,152,784],"THM260 deficit prefix")
    # An odd self-converse tournament need not have one distinguished vertex.
    # The directed 3-cycle has three rotationally equivalent possible centers.
    n=3;k=1
    require(expected(n,k)[1]==2,"root center count")
    return {"mixed_cycle_true":str(true_layer),"printed_formula":str(printed_formula),
            "rooted":rooted,"deficit":deficit,"unique_vertex_hostile":"Directed 3-cycle: all three vertices are possible anti-involution centers in one Aut orbit."}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--receipt",default="audit-normal.json")
    parser.add_argument("--literal-bound",type=int,default=7)
    parser.add_argument("--burnside-bound",type=int,default=9)
    args=parser.parse_args()
    rows=direct_burnside(args.burnside_bound)
    print(json.dumps({"phase":"signed-equation Burnside","cells":len(rows),"permutations":sum(x['permutations'] for x in rows)}),flush=True)
    literal=[]
    for n in range(1,args.literal_bound+1):
        for k in range(1,n+1):
            row=literal_census(n,k)
            literal.append(row)
            print(json.dumps({"phase":"literal","n":n,"k":k,"ordinary":row['ordinary'],"converse":row['converse']}),flush=True)
    correction=correction_controls()
    receipt={"status":"PASS","checks":CHECKS,"signed_burnside":rows,"literal_orbits":literal,
             "correction":correction,"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(args.receipt).write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"status":"PASS","checks":CHECKS,"receipt":args.receipt}),flush=True)

if __name__=="__main__":
    main()
