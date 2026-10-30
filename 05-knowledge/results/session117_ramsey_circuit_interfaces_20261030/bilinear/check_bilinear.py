#!/usr/bin/env python3
"""Exact, dependency-free checks; proof supplies the all-m conclusions."""
import json

def require(ok, label):
    if not ok:
        raise RuntimeError(label)

def dot(a,b): return (a&b).bit_count()&1

def span(values):
    s={0}
    for v in values: s |= {u^v for u in tuple(s)}
    return frozenset(s)

def subsets(n):
    return [frozenset(i for i in range(n) if mask>>i&1) for mask in range(1,1<<n)]

def subspaces(d):
    level={frozenset({0})}; out=set(level)
    for _ in range(d):
        nxt=set()
        for u in level:
            for v in range(1<<d):
                if v not in u: nxt.add(frozenset(u | {x^v for x in u}))
        out |= nxt; level=nxt
    return out

def run(m):
    N=1<<m; vals=range(N); full=frozenset(vals)
    ip=lambda z: dot(z&(N-1),z>>m)
    pol=lambda z,w: ip(z^w)^ip(z)^ip(w)
    H=[[(-1)**dot(x,y) for y in vals] for x in vals]
    for x in vals:
        for y in vals:
            require(sum(H[x][z]*H[y][z] for z in vals)==(N if x==y else 0),'orthogonality')
    for a in vals:
        for b in vals:
            require(sum((-1)**(dot(x,y)^dot(a,x)^dot(b,y)) for x in vals for y in vals)==N*((-1)**dot(a,b)),'Walsh')
    sets=subsets(N); mono=[0,0]; maxima=[0,0]; eq=[0,0]
    for A in sets:
        for B in sets:
            a=min(A); b=min(B); c=dot(a,b)
            if any(dot(x,y)!=c for x in A for y in B): continue
            mono[c]+=1; area=len(A)*len(B); bound=N if c==0 else N//2
            require(area<=bound,'rectangle bound'); maxima[c]=max(maxima[c],area)
            if area==bound:
                eq[c]+=1
                if c==0:
                    U=span(A)
                    require(A==U and B==frozenset(y for y in vals if all(dot(x,y)==0 for x in U)),'zero equality')
                else:
                    U=span(x^a for x in A); V=span(y^b for y in B)
                    require(A==frozenset(a^x for x in U) and B==frozenset(b^y for y in V),'one affine equality')
                    require(V==frozenset(y for y in vals if not dot(a,y) and all(not dot(x,y) for x in U)),'one equality parameterization')
    require(maxima==[N,N//2],'sharp rectangles')
    spaces=subspaces(2*m); maxdim=[-1,-1]; nummaxzero=0; lag=0; flats=0
    for U in spaces:
        dim=len(U).bit_length()-1
        isotropic=all(not pol(u,v) for u in U for v in U)
        remaining=set(range(N*N)); constant=[]
        while remaining:
            a=min(remaining); C={a^u for u in U}; remaining-=C; flats+=1
            colors={ip(z) for z in C}
            if len(colors)==1:
                c=next(iter(colors)); constant.append(c); maxdim[c]=max(maxdim[c],dim)
                require(isotropic and dim<=m,'constant isotropy')
                if dim==m:
                    require(c==0,'maximum flat color'); nummaxzero+=1
        if isotropic and dim==m:
            lag+=1; require(constant==[0],'unique constant coset')
    expected=1
    for j in range(1,m+1): expected*=2**j+1
    require(maxdim==[m,m-1] and lag==expected and nummaxzero==expected,'flat classification')
    # Literal gate construction and sensitivity witnesses, not a circuit search.
    gates=[('AND',i,m+i) for i in vals if i<m]
    out=2*m
    for i in range(1,m):
        gates.append(('XOR',out,2*m+i)); out=2*m+len(gates)-1
    require(len(gates)==2*m-1,'gate size')
    for z in range(N*N):
        wires=[(z>>i)&1 for i in range(2*m)]
        for op,a,b in gates: wires.append((wires[a]&wires[b]) if op=='AND' else wires[a]^wires[b])
        require(wires[out]==ip(z),'circuit evaluation')
    for i in range(2*m):
        witness=1<<((i+m)%(2*m))
        require(ip(witness)!=ip(witness^(1<<i)),'essential variable')
    return {'m':m,'rectangle_pairs':len(sets)**2,'monochromatic_rectangles':mono,'maximum_rectangles':eq,'linear_subspaces':len(spaces),'affine_flats':flats,'maximum_constant_dimensions':maxdim,'maximum_zero_flats':nummaxzero,'gate_count':len(gates)}

if __name__=='__main__':
    print(json.dumps({'status':'PASS','range':[1,2,3],'results':[run(m) for m in range(1,4)]},indent=2,sort_keys=True))
