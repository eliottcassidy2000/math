"""Exact full automorphism orders by injective adjacency-preserving backtracking."""
from pathlib import Path
import json
P=Path(__file__).resolve().parent

def automorphisms(T):
    n=len(T);deg=[r.bit_count() for r in T];maps=[];visits=0
    def rec(mapping,used):
        nonlocal visits
        visits+=1
        if len(mapping)==n:
            maps.append(tuple(mapping[i] for i in range(n)));return
        best=None
        for v in range(n):
            if v in mapping:continue
            cand=[w for w in range(n) if not(used>>w&1) and deg[v]==deg[w] and all((T[v]>>a&1)==(T[w]>>b&1) for a,b in mapping.items())]
            if not cand:return
            if best is None or len(cand)<len(best[1]):best=(v,cand)
        v,cand=best
        for w in cand:
            mapping[v]=w;rec(mapping,used|1<<w);del mapping[v]
    rec({},0)
    return maps,visits

def main():
    zs=json.loads((P/'input_audit.json').read_text())['hostile']['rows'];out=[]
    for T,expected in zip(zs,[81,27]):
        autos,visits=automorphisms(T)
        if len(autos)!=expected:raise RuntimeError('full order mismatch')
        if len(set(autos))!=len(autos):raise RuntimeError('duplicate permutation')
        out.append({'full_automorphism_order':len(autos),'search_nodes':visits,'automorphisms':autos})
    (P/'automorphism_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print([(x['full_automorphism_order'],x['search_nodes']) for x in out])
if __name__=='__main__':main()
