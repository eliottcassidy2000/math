#!/usr/bin/env python3
"""Independent finite audit of the proposed symmetric 7p threshold.
Matrix tournaments + cached recursive path counting, no producer imports.
"""
from functools import lru_cache
from itertools import product
import json

def need(c,m):
    if not c: raise RuntimeError(m)

def frommask(n,m):
    a=[[0]*n for _ in range(n)];k=0
    for i in range(n):
        for j in range(i+1,n):a[i][j]=(m>>k)&1;a[j][i]=1-a[i][j];k+=1
    return a

def fromrows(rs):return [[(r>>j)&1 for j in range(len(rs))] for r in rs]
def valid(a):return all(a[i][i]==0 and all(a[i][j] in (0,1) and a[i][j]+a[j][i]==1 for j in range(i)) for i in range(len(a)))
def strong(a):
    for start in range(len(a)):
        seen={start};front={start}
        while front:
            front={j for i in front for j in range(len(a)) if a[i][j]}-seen
            seen|=front
        if len(seen)<len(a):return False
    return True

def H(a):
    n=len(a)
    @lru_cache(None)
    def walk(last,remaining):
        if not remaining:return 1
        return sum(walk(v,remaining^(1<<v)) for v in range(n) if remaining>>v&1 and a[last][v])
    return sum(walk(v,((1<<n)-1)^(1<<v)) for v in range(n))
def circulant(n,choice):
    ds={i if b else n-i for i,b in enumerate(choice,1)}
    return [[int((j-i)%n in ds) for j in range(n)] for i in range(n)]
def substitute(outer,inner):
    n=len(inner);k=len(outer);a=[[0]*(n+k-1) for _ in range(n+k-1)]
    for i in range(len(a)):
        for j in range(len(a)):
            a[i][j]=inner[i][j] if i<n and j<n else outer[0 if i<n else i-n+1][0 if j<n else j-n+1]
    need(valid(a),'substituted tournament validity')
    return a

def automorphism(a,p):return all(a[i][j]==a[p[i]][p[j]] for i in range(len(a)) for j in range(len(a)))
def main():
    triangle=[[0,1,0],[0,0,1],[1,0,0]]
    orbit={};c7=[]
    for n in (5,7):
        pairs=set()
        for bits in product((0,1),repeat=(n-1)//2):
            a=circulant(n,bits);need(valid(a) and strong(a),'circulant validity/strongness')
            pairs.add((H(a),H(substitute(triangle,a))))
            if n==7:c7.append(a)
        orbit[n]=sorted(pairs)
    need(orbit[5]==[(15,135)],'five-orbit figures')
    need(orbit[7]==[(175,2163),(189,2205)],'seven-orbit figures')
    outs=[frommask(4,m) for m in range(64) if strong(frommask(4,m))]
    need(len(outs)==24,'strong4 coverage')
    values=[]
    for out in outs:
        for inn in c7:values.append(H(substitute(out,inn)))
    need(len(values)==192 and min(values)==4151,'rooted substitution minimum')
    rows=[142,156,184,240,225,195,135,512,255,383]
    a=fromrows(rows);p=list(range(1,7))+[0,7,8,9]
    need(valid(a) and strong(a) and H(a)==4151,'threshold witness')
    need(automorphism(a,p),'order-seven symmetry')
    need(4151==7*593 and all(593%d for d in range(2,25)),'prime quotient')
    bad=[row[:] for row in a];bad[0][1],bad[1][0]=bad[1][0],bad[0][1]
    need(valid(bad) and not automorphism(bad,p),'mutated symmetry hostile')
    print(json.dumps({'status':'PASS','circulant_and_triangle_values':orbit,'rooted_strong4_substitutions':192,
      'substitution_values':sorted(set(values)),'minimum':4151,'prime_quotient':593,'threshold_witness_rows':rows,
      'cyclic_automorphism':p,'hostile_one_arc_breaks_symmetry':True},indent=2,sort_keys=True))
if __name__=='__main__':main()
