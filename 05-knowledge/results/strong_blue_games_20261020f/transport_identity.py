"""Exhaustive small controls for the paired-triangle relocation identity."""
from itertools import combinations
from pathlib import Path
import hashlib
import json
from boundary_transport_probe import G

HERE=Path(__file__).resolve().parent


def c3(t):
    return sum(all(sum(bool(t[u]>>v&1) for v in vs if v!=u)==1 for u in vs)
               for vs in combinations(range(len(t)),3))


def main():
    counts=[]
    for n in (5,6):
        fixed=[(i,i+1) for i in range(n-1)]
        free=[(u,v) for u,v in combinations(range(n),2) if v!=u+1]
        tested=0;orientations=0
        for mask in range(1<<len(free)):
            t=[0]*n
            for u,v in fixed:t[u]|=1<<v
            for bit,(u,v) in enumerate(free):
                if mask>>bit&1:t[u]|=1<<v
                else:t[v]|=1<<u
            t=tuple(t);P=tuple(range(n));S=G.complement(t,P)
            orientations+=1
            for x in range(1,n-3):
                for z in range(x+2,n-1):
                    Q=P[:x]+P[x+1:z+1]+(x,)+P[z+1:]
                    if not G.is_path(t,Q):continue
                    a,b,c=x-1,x+1,z+1
                    triangles=((a,x,b),(x,z,c))
                    changed=set();replay=list(S)
                    for tri in triangles:
                        for u,v in zip(tri,tri[1:]+tri[:1]):
                            G.need(bool(S[u]>>v&1),'directed triangle')
                            e=tuple(sorted((u,v)));G.need(e not in changed,'disjoint edge supports');changed.add(e)
                            replay[u]^=1<<v;replay[v]^=1<<u
                    target=G.complement(t,Q)
                    G.need(tuple(replay)==target and len(changed)==6,'literal paired triangle reversal')
                    G.need([r.bit_count() for r in S]==[r.bit_count() for r in target],'labelled score preservation')
                    G.need(c3(S)==c3(target),'direct triangle count preservation')
                    tested+=1
        G.need(tested>0,'positive control at every order')
        counts.append(dict(n=n,orientations=orientations,legal_relocations=tested))
    out=dict(status='FINITE-EXACT controls for an all-order identity',
             source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             universe='all standard-path tournaments n=5,6; all legal interior relocations with |B|>=2',counts=counts)
    (HERE/'transport_identity.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
