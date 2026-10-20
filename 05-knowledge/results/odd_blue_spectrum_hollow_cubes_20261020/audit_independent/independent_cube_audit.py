#!/usr/bin/env python3
"""Exhaustive independent free-cube audit of all hollow frames through n=9."""
from independent_pairing_audit import *


def audit(m):
    ps=list(it.permutations(range(m)));frames=packets=0;dims=co.Counter();Bcounts=[co.Counter(),co.Counter()]
    for e in range(1<<(m*(m-1))):
        A,B=make_e(m,e);codes=[ecode(A,B,p) for p in ps]
        if e!=min(codes):continue
        frames+=1;H=[p for p,c in zip(ps,codes) if c==e]
        remaining=set(range(m));orbits=[]
        while remaining:
            O={p[min(remaining)] for p in H};remaining-=O;orbits.append(O)
        assert all(len(O)%2 for O in orbits)
        r=len(orbits);assert r%2==m%2
        flips=[sum(1<<i for i in O) for O in orbits]
        translations=[0]
        for f in flips:translations += [b^f for b in translations]
        canon=lambda b:min(perm_b(b,p) for p in H)
        G0=graph(A,B,0);base=c3(G0)
        coefficients=[2*m-1-2*sum(G0[i]) for i in range(m)]
        h=len(sym_paths(G0));assert h%len(H)==0 and h%2
        b0=h//len(H);mass=0;shore_mass=[0,0];unseen={canon(b) for b in range(1<<m)}
        while unseen:
            b=min(unseen);cube={canon(b^t) for t in translations}
            assert len(cube)==1<<r and cube<=unseen
            unseen-=cube;dims[r]+=1;packets+=1
            stabs=[tuple(p for p in H if perm_b(v,p)==v) for v in cube]
            assert all(s==stabs[0] for s in stabs)
            parity_vectors={tuple(sum((v>>i)&1 for i in O)%2 for O in orbits) for v in cube}
            assert len(parity_vectors)==1<<r
            blue=h//len(stabs[0]);parities=[]
            for v in cube:
                triangles=c3(graph(A,B,v))
                assert triangles==base+sum(c*((v>>i)&1) for i,c in enumerate(coefficients))
                parity=triangles%2;parities.append(parity)
                assert parity==(m-v.bit_count())%2
                Bcounts[parity][blue]+=1;mass+=blue;shore_mass[parity]+=blue
            assert sum(parities)==1<<(r-1)
        assert mass==b0*(1<<m) and shore_mass==[b0*(1<<(m-1))]*2
    assert Bcounts[0]==Bcounts[1]
    if m%2==0:assert all(v%2==0 for v in Bcounts[0].values())
    return {'n':2*m+1,'hollow_frame_classes':frames,'cube_packets':packets,'cube_dimensions':dict(sorted(dims.items())),'each_shore_spectrum':dict(sorted(Bcounts[0].items()))}

if __name__=='__main__':
    rows=[audit(m) for m in range(1,5)]
    for r in rows:print(json.dumps(r))
    (OUT/'independent_cube_audit.json').write_text(json.dumps(rows,indent=2)+'\n')
