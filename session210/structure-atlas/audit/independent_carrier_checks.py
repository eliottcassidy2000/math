#!/usr/bin/env python3
"""Independent exact checks: no imports from session 210 worker implementations."""
from itertools import combinations, permutations, product
from collections import Counter
import json
from pathlib import Path

ROOT=Path(__file__).parent

def require(x, context):
    if not x: raise RuntimeError(context)

def tournament(n, bits):
    a=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        u,v=(i,j) if bits>>k&1 else (j,i)
        a[u]|=1<<v
    return a

def hp(a, subset=None):
    n=len(a); full=(1<<n)-1 if subset is None else subset
    if not full:return 1
    dp=[{} for _ in range(1<<n)]
    for v in range(n):
        if full>>v&1:dp[1<<v][v]=1
    for m in range(1,full+1):
        if m&~full:continue
        for v,c in dp[m].items():
            avail=a[v]&full&~m
            while avail:
                bit=avail&-avail;avail-=bit;w=bit.bit_length()-1
                dp[m|bit][w]=dp[m|bit].get(w,0)+c
    return sum(dp[full].values())

def cycles(a, odd=False):
    n=len(a);out=[]
    for k in range(3,n+1):
        if odd and k%2==0:continue
        for sub in combinations(range(n),k):
            r=sub[0]
            for p in permutations(sub[1:]):
                q=(r,)+p
                if all(a[q[i]]>>q[(i+1)%k]&1 for i in range(k)):
                    out.append((q,sum(1<<v for v in sub)))
    return out

def packings(odd):
    out=[()]
    def rec(start,mask,ids):
        for i in range(start,len(odd)):
            s=odd[i][1]
            if not mask&s:
                p=ids+(i,);out.append(p);rec(i+1,mask|s,p)
    rec(0,0,())
    return out

def signed_cells(packs):
    cells=[]
    for p in packs:
        if p:
            for s in product((-1,1), repeat=len(p)-1):
                cells.append(tuple(zip(p,(1,)+s)))
    return cells

def canonical(face):
    if not face:return ()
    return tuple((v,s*face[0][1]) for v,s in face)

def mod2_rank(columns):
    basis={}
    for c in columns:
        while c:
            p=c.bit_length()-1
            if p not in basis:basis[p]=c;break
            c^=basis[p]
    return len(basis)

def topology(cells):
    if not cells:return {'f':[],'betti_mod2':[],'euler':0}
    d=max(map(len,cells))-1
    levels=[[c for c in cells if len(c)==k+1] for k in range(d+1)]
    ranks=[0]
    for k in range(1,d+1):
        pos={c:i for i,c in enumerate(levels[k-1])};cols=[]
        for c in levels[k]:
            v=0
            for j in range(len(c)):
                face=canonical(c[:j]+c[j+1:]);v^=1<<pos[face]
            cols.append(v)
        ranks.append(mod2_rank(cols))
    ranks.append(0)
    b=[len(levels[k])-ranks[k]-ranks[k+1] for k in range(d+1)]
    return {'f':list(map(len,levels)),'betti_mod2':b,'euler':sum((-1)**k*len(v) for k,v in enumerate(levels))}

def gate(a):
    n=len(a);g=list(a)+[1<<(n+1),(1<<n)-1]
    for v in range(n):g[v]|=1<<n
    return g

def gates(a):
    n=len(a);ind=[sum(row>>v&1 for row in a) for v in range(n)]
    return [(u,v) for u in range(n) if a[u].bit_count()==1 for v in range(n) if a[u]>>v&1 and ind[v]==1]

def join(a,b):
    k=len(a);return [row|(((1<<len(b))-1)<<k) for row in a]+[row<<k for row in b]

def main():
    out={'scope':'All labelled tournaments through order 5; all 32768 at order 6 for OCF and projective cells. Gate and full cycle-polynomial checks through core order 5.','orders':[]}
    for n in range(1,7):
        spectrum=set();cases=1<<(n*(n-1)//2)
        for bits in range(cases):
            a=tournament(n,bits);h=hp(a);odd=cycles(a,odd=True);p=packings(odd);counts=Counter(map(len,p));cells=signed_cells(p)
            require(h==sum(c*2**r for r,c in counts.items()),('ocf',n,bits))
            require(len(cells)==(h-1)//2,('projective count',n,bits))
            require(all(canonical(c)==c for c in cells),('sign normalization',n,bits))
            # Every cell closure has distinct quotient vertices and all faces exist.
            cellset=set(cells)
            for c in cells:
                require(len(set(v for v,s in c))==len(c),('injectivity',n,bits,c))
                for j in range(len(c)):
                    if len(c)>1:require(canonical(c[:j]+c[j+1:]) in cellset,('face',n,bits,c))
            spectrum.add(h)
            if n<=5:
                g=gate(a);allc=cycles(g);cg=Counter(len(q) for q,m in allc);orig=Counter(len(q) for q,m in cycles(a))
                for m in range(1,1<<n):orig[m.bit_count()+2]+=hp(a,m)
                require(cg==orig,('gate cycle polynomial',n,bits,cg,orig))
                require(cg[n+2]==h,('gate HC',n,bits))
                require(len(gates(g))==(3 if n==1 else 1),('unique gate',n,bits,gates(g)))
        out['orders'].append({'n':n,'tournaments':cases,'spectrum':sorted(spectrum),'projective_indices':sorted((h-1)//2 for h in spectrum)})
        print('passed n',n,flush=True)
    controls=[];a=[];c3=[2,4,1]
    for r in range(1,5):
        a=join(a,c3);p=packings(cycles(a,odd=True));cs=signed_cells(p);top=topology(cs)
        require(top['betti_mod2']==[1]*r,('RP topology',r,top))
        require(len(cs)==(3**r-1)//2,('RP cells',r))
        controls.append({'triangle_factors':r,'cells':len(cs),**top})
    out['projective_controls']=controls
    # Abstract carriers can have forbidden indices: r isolated cycle vertices.
    out['unrestricted_complex_hostiles']=[{'isolated_vertices':r,'total_projective_cells':len(signed_cells([()]+[(j,) for j in range(r)]))} for r in (3,10)]
    out['status']='PASS'
    (ROOT/'independent_carrier_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
