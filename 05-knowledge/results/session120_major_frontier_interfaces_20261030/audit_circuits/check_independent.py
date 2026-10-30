#!/usr/bin/env python3
"""Independent audit: standard library only, no producer imports/asserts."""
import itertools as it
import json

def require(ok, message):
    if not ok: raise RuntimeError(message)

def span(rows):
    s={0}
    for r in rows: s |= {x^r for x in tuple(s)}
    return s

def forest(n, es):
    # Component growth, rather than the producer's union/find test.
    seen=set(); components=0
    for v in range(n):
        if v in seen: continue
        components+=1; todo=[v]; seen.add(v)
        while todo:
            a=todo.pop()
            for x,y in es:
                b=y if x==a else x if y==a else None
                if b is not None and b not in seen: seen.add(b); todo.append(b)
    return len(es)==n-components

def matching(es):
    for k in range(len(es),-1,-1):
        if any(len({v for edge in chosen for v in edge})==2*k for chosen in it.combinations(es,k)):
            return k

def truth_linear(row,n):
    return sum(((row&a).bit_count()%2)<<a for a in range(1<<n))

def essentials(table,n):
    return sum(any(((table>>a)^(table>>(a^(1<<j))))&1 for a in range(1<<n)) for j in range(n))

def build(n,es,delta):
    # Highest leaf, opposite producer choice. Literal gate DAG and shear records.
    remain=set(es); gates=[]; products=[]; shears=[]; pivots=[]
    def gate(op,a,b):
        require(a<n+delta+len(gates) and b<n+delta+len(gates),'acyclic predecessor')
        gates.append((op,a,b)); return n+delta+len(gates)-1
    while remain:
        adj={i:set() for i in range(n)}
        for a,b in remain: adj[a].add(b); adj[b].add(a)
        x=max(i for i in adj if len(adj[i])==1); y=next(iter(adj[x])); w=x
        for z in sorted(adj[y]-{x},reverse=True):
            w=gate('xor',w,z); shears.append((x,z))
        products.append(gate('and',y,w)); pivots.append((x,y))
        remain={edge for edge in remain if x not in edge and y not in edge}
    terms=products+([n] if delta else []); out=terms.pop()
    while terms: out=gate('xor',out,terms.pop())
    return gates,out,shears,pivots

def main():
    graphs=assignments=first_pivots=0
    for n in range(2,7):
        universe=list(it.combinations(range(n),2))
        for e in range(1,n):
            for es in it.combinations(universe,e):
                if len({v for edge in es for v in edge})!=n or not forest(n,es): continue
                graphs+=1; nu=matching(es)
                adj=[sum(1<<b for a,b in es if a==i)+sum(1<<a for a,b in es if b==i) for i in range(n)]
                require(len(span(adj))==1<<(2*nu),'polar rank')
                # Check every legal first leaf choice, not just the constructed trace.
                for x in range(n):
                    neighbors=[y for y in range(n) if adj[x]>>y&1]
                    if len(neighbors)==1:
                        y=neighbors[0]; residual=[edge for edge in es if x not in edge and y not in edge]
                        require(matching(residual)==nu-1,'arbitrary leaf matching exchange'); first_pivots+=1
                for delta in (0,1):
                    N=n+delta; gates,out,shears,pivots=build(n,es,delta)
                    require(len(gates)==e+nu+delta-1,'fixed inventory')
                    require(len(shears)==e-nu and len(pivots)==nu,'affine inventory')
                    inputs=[truth_linear(1<<i,N) for i in range(N)]; expected=0
                    for a,b in es: expected^=inputs[a]&inputs[b]
                    if delta: expected^=inputs[n]
                    vals=inputs[:]
                    for op,a,b in gates: vals.append(vals[a]^vals[b] if op=='xor' else vals[a]&vals[b])
                    require(vals[out]==expected,'all-assignment fixed semantics')
                    coords=inputs[:]; matrix=[1<<i for i in range(N)]
                    for x,z in shears: coords[x]^=coords[z]; matrix[x]^=matrix[z]
                    require(len(span(matrix))==1<<N,'invertible shear matrix')
                    normal=0
                    for x,y in pivots: normal^=coords[x]&coords[y]
                    if delta: normal^=coords[n]
                    require(normal==expected,'all-assignment affine semantics')
                    for x,z in reversed(shears): coords[x]^=coords[z]
                    require(coords==inputs,'reverse shear inverse')
                    canonical=0
                    for x,y in pivots: canonical^=inputs[x]&inputs[y]
                    if delta: canonical^=inputs[n]
                    require(essentials(canonical,N)==2*nu+delta,'endpoint I count')
                    require(essentials(expected,N)==N,'original essential inputs')
                    require((2*e+delta-1)-(len(gates)-len(shears))==2*(e-nu),'total charge saving')
                    assignments+=1<<N
    # Exhaust all GL(4,2) charts for four structurally different representatives.
    cases=[([(0,1),(1,2),(2,3)],0), ([(0,1),(0,2),(0,3)],0),
           ([(0,1),(2,3)],0), ([(0,1),(1,2)],1)]
    minima=[4]*len(cases); charts=0
    for rows in it.product(range(1,16),repeat=4):
        if len(span(rows))!=16: continue
        charts+=1; vals=[truth_linear(r,4) for r in rows]
        for j,(es,delta) in enumerate(cases):
            t=0
            for a,b in es: t^=vals[a]&vals[b]
            if delta: t^=vals[3]
            minima[j]=min(minima[j],essentials(t,4))
    require(charts==20160 and minima==[4,2,4,3],'exhaustive affine chart essentials')
    # Empty forest terminal exceptions are outside theorem: constant/literal cost zero.
    require(essentials(0,1)==0 and essentials(2,1)==1,'terminal exceptions')
    print(json.dumps({'status':'PASS','forests':graphs,'assignments':assignments,
        'all_first_leaf_choices':first_pivots,'GL4_charts_per_case':charts,
        'chart_min_essential_inputs':minima,'max_vertices':6},sort_keys=True))
if __name__=='__main__': main()
