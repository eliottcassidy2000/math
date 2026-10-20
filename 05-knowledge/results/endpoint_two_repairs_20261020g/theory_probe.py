"""Bounded universal-block ledger search for the second-position port."""
from importlib.util import spec_from_file_location,module_from_spec
from itertools import permutations
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
s=spec_from_file_location('endpoint_old',HERE.parent/'strong_blue_games_20261020f/detours_probe.py')
D=module_from_spec(s);s.loader.exec_module(D);G=D.G


def endpoint_two(middle,u):
    a=list(G.join([(0,),middle,(0,)]));n=len(a);v=n-1-u
    for x,y in ((0,n-1),(0,u),(v,n-1)):
        a[x]^=1<<y;a[y]^=1<<x
    return tuple(a)


def graph_path(edges,n):
    if len(edges)!=n-1:return None
    adj=[[] for _ in range(n)]
    for u,v in edges:adj[u].append(v);adj[v].append(u)
    if any(not vs or len(vs)>2 for vs in adj):return None
    ends=[u for u,vs in enumerate(adj) if len(vs)==1]
    if len(ends)!=2:return None
    p=[ends[0]]
    while len(p)<n:
        nxt=[v for v in adj[p[-1]] if len(p)<2 or v!=p[-2]]
        if not nxt:return None
        p.append(nxt[0])
    return tuple(p) if len(set(p))==n else None


def main():
    # P=(a,b,x,U,y,c,z), U has two actual vertices.
    n=8;tokens=((0,),(1,),(2,),(3,4),(5,),(6,),(7,));P=tuple(range(n))
    sources=[endpoint_two(m,2) for m in G.blue_tilings(n-2)
             if m[-1]&1 and m[0]>>(n-4)&1]
    words=[sum(parts,()) for parts in permutations(tokens)]
    first=[q for q in words if all(G.is_path(a,q) for a in sources)]
    certs=[];valid=0
    for q in first:
        after=[G.complement(a,q) for a in sources]
        if any(G.blue(a,q) for a in sources):continue
        for r in words:
            if not all(G.is_path(a,r) for a in after):continue
            if any(G.blue(a,r) for a in after):continue
            sp=graph_path(G.edge_set(P)^G.edge_set(q)^G.edge_set(r),n)
            if sp is None:continue
            after2=[G.complement(a,r) for a in after]
            if not G.is_path(after2[0],sp):sp=sp[::-1]
            if not all(G.is_path(a,sp) for a in after2):continue
            valid+=1
            if any(G.blue(a,sp) for a in after2):continue
            if not all(G.complement(a,sp)==G.complement(t,P) for a,t in zip(after2,sources)):
                raise RuntimeError('parity failure')
            certs.append((q,r,sp))
    out=dict(n=n,sources=len(sources),first_words=len(first),valid=valid,certificates=certs)
    (HERE/'theory_probe.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
