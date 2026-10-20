"""One-port hostile controls on the already authenticated sixteen-point gate."""
from pathlib import Path
from collections import Counter
import json

def need(ok,message):
    if not ok:raise RuntimeError(message)

def main():
    old=Path(__file__).resolve().parents[1]/'first_exit_complement_games_20261020/palette-glue.json'
    data=json.loads(old.read_text()); edges=[(a,b) for a,b in data['edges'] if a<16 and b<16]
    adj=[set() for _ in range(16)]
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    core_edges={(a,b) for a,b in edges if 1<=a<=6 and 1<=b<=6}
    triangles=[{1,2,3},{4,5,6}]
    need(core_edges=={(1,2),(1,3),(2,3),(4,5),(4,6),(5,6),(3,6)},'two triangles plus bridge')
    base={v:(v-7)//3 for v in range(7,16)}
    def counts(boundary):
        if any(boundary[a]==boundary[b] for a,b in edges if a in boundary and b in boundary):return None
        cs=boundary.copy()
        def go():
            if len(cs)==16:return 1
            choices=[(set(range(6))-{cs[j] for j in adj[i] if j in cs},i) for i in range(16) if i not in cs]
            allowed,i=min(choices,key=lambda p:(len(p[0]),-len(adj[p[1]])))
            total=0
            for c in allowed:cs[i]=c;total+=go()
            cs.pop(i,None);return total
        out=[]
        for c in range(6):cs[0]=c;out.append(go());cs.pop(0)
        return out
    need(counts(base)==[24,24,24,0,0,0],'inherited gate control')
    changes=[];hist=Counter()
    for v in range(7,16):
        for c in range(6):
            if c==base[v]:continue
            b=base.copy();b[v]=c;out=counts(b)
            if out is None:hist['improper_boundary']+=1
            else:
                escape=sum(out[3:]);hist['repairs_disjoint_join' if escape else 'still_blocked']+=1
                changes.append({'port':v,'old':base[v],'new':c,'output_counts':out,
                                'extensions_with_other_palette_345':24*escape})
    # A universal attempted synchronizer in a four-colour carrier must survive this.
    c4=data['four_colouring'][:16];need(all(c4[a]!=c4[b] for a,b in edges),'inherited 4-colouring')
    recolour_controls=0
    for a in range(7,16):
        for b in range(a+1,16):
            c=c4[:];c[a]=4;c[b]=5
            need(all(c[u]!=c[v] for u,v in edges),'two fresh-colour ports')
            recolour_controls+=1
    print(json.dumps({'status':'FINITE-EXACT, externally conditioned ports',
        'internal_24_decoder':'two labelled triangles joined by one edge: 3*2*2*2',
        'one_port_changes':changes,'histogram':dict(hist),
        'fresh_two_port_controls':recolour_controls},indent=2))
if __name__=='__main__':main()
