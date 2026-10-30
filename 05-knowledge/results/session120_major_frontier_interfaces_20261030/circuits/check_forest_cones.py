#!/usr/bin/env python3
"""Independent exhaustive finite checks; no producer imports, no assert guards."""
from itertools import combinations
import json

def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)

def forest(n, edges):
    p=list(range(n))
    def root(a):
        while p[a]!=a: a=p[a]
        return a
    for a,b in edges:
        a,b=root(a),root(b)
        if a==b: return False
        p[a]=b
    return True

def matching(edges):
    if not edges: return 0
    a,b=edges[0]
    return max(matching(edges[1:]),1+matching([e for e in edges[1:] if a not in e and b not in e]))

def rank(rows):
    rows=rows[:]; out=0
    for bit in range(max((r.bit_length() for r in rows),default=0)):
        k=next((k for k in range(out,len(rows)) if rows[k]>>bit&1),None)
        if k is None: continue
        rows[k],rows[out]=rows[out],rows[k]
        for j in range(len(rows)):
            if j!=out and rows[j]>>bit&1: rows[j]^=rows[out]
        out+=1
    return out

def build(n,edges,delta):
    remaining=set(edges); steps=[]; gates=[]; products=[]
    def gate(op,a,b):
        gates.append((op,a,b)); return n+delta+len(gates)-1
    while remaining:
        neighbors={i:[] for i in range(n)}
        for a,b in remaining:
            neighbors[a].append(b); neighbors[b].append(a)
        x=next(i for i in range(n) if len(neighbors[i])==1)
        y=neighbors[x][0]; others=sorted(z for z in neighbors[y] if z!=x)
        w=x
        for z in others: w=gate('xor',w,z)
        products.append(gate('and',w,y)); steps.append((x,y,others))
        remaining={e for e in remaining if x not in e and y not in e}
    terms=products+([n] if delta else [])
    out=terms[0]
    for t in terms[1:]: out=gate('xor',out,t)
    return steps,gates,out

def cycle_controls():
    # Odd cycle: rank alone does not account for the radical-linear bit.
    need(rank([0b110,0b101,0b011])==2,'triangle rank')
    images=set()
    for bits in range(8):
        x,y,z=(bits>>i&1 for i in range(3))
        q=(x&y)^(y&z)^(z&x)
        xp,yp,zp=x^z,y^z,z
        images.add((xp,yp,zp))
        need(q==((xp&yp)^zp),'triangle affine identity')
        need((xp^zp,yp^zp,zp)==(x,y,z),'triangle inverse')
        opposite=((x^1)&(y^1))^((y^1)&(z^1))^((z^1)&(x^1))
        need((q^opposite)==1,'triangle radical derivative')
    need(len(images)==8,'triangle bijection')
    need(matching([(0,1),(1,2),(0,2)])==1,'triangle matching')
    # Even cycle: matching size no longer determines polar rank.
    need(rank([0b1010,0b0101,0b1010,0b0101])==2,'C4 rank')
    images=set()
    for bits in range(16):
        a,b,c,d=(bits>>i&1 for i in range(4))
        q=(a&b)^(b&c)^(c&d)^(d&a)
        p,r=a^c,b^d
        images.add((p,r,c,d))
        need(q==(p&r),'C4 affine identity')
        need((p^c,r^d,c,d)==(a,b,c,d),'C4 inverse')
    need(len(images)==16,'C4 bijection')
    need(matching([(0,1),(1,2),(2,3),(0,3)])==2,'C4 matching')
    return 24

def run():
    graphs=assignments=star_graphs=0
    for n in range(2,7):
        possible=list(combinations(range(n),2))
        for mask in range(1,1<<len(possible)):
            edges=[e for j,e in enumerate(possible) if mask>>j&1]
            if len(edges)>=n or len(set(sum(([a,b] for a,b in edges),[])))!=n: continue
            if not forest(n,edges): continue
            graphs+=1; e=len(edges); nu=matching(edges)
            rows=[sum(1<<b for a,b in edges if a==i)+sum(1<<a for a,b in edges if b==i) for i in range(n)]
            need(rank(rows)==2*nu,'polar rank')
            # A forest with nu=number of components consists of stars.
            is_star=(nu==n-e)
            star_graphs+=is_star
            for delta in (0,1):
                steps,gates,out=build(n,edges,delta)
                need(len(steps)==nu,'matching recursion')
                need(sum(len(s[2]) for s in steps)==e-nu,'transform inventory')
                need(len(gates)==e+nu+delta-1,'fixed gate inventory')
                need((2*e+delta-1)-(2*nu+delta-1)==2*(e-nu),'affine savings')
                if is_star: need(len(gates)==n+delta-1,'star sharpness')
                for bits in range(1<<(n+delta)):
                    original=[bits>>i&1 for i in range(n+delta)]
                    expected=sum(original[a]*original[b] for a,b in edges)%2
                    if delta: expected^=original[n]
                    vals=original[:]
                    for op,a,b in gates:
                        need(a<len(vals) and b<len(vals),'DAG topological order')
                        vals.append(vals[a]^vals[b] if op=='xor' else vals[a]&vals[b])
                    need(vals[out]==expected,'fixed-coordinate semantics')
                    coordinates=original[:]
                    for x,y,others in steps:
                        for z in others: coordinates[x]^=coordinates[z]
                    canonical=sum(coordinates[x]*coordinates[y] for x,y,_ in steps)%2
                    if delta: canonical^=coordinates[n]
                    need(canonical==expected,'canonical semantics')
                    recovered=coordinates[:]
                    for x,y,others in reversed(steps):
                        for z in others: recovered[x]^=recovered[z]
                    need(recovered==original,'invertible embedding')
                    assignments+=1
    cycle_assignments=cycle_controls()
    print(json.dumps(dict(status='PASS',forests=graphs,star_forests=star_graphs,assignments=assignments,cycle_control_assignments=cycle_assignments,max_vertices=6),sort_keys=True))

if __name__=='__main__': run()
