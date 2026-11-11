"""Exact marked-seed deletion DAG through seven ordinary base vertices."""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import json
import port_transport as pt

HERE=Path(__file__).resolve().parent


def delete(rows,mask,v):
    labels=[i for i in range(len(rows)) if i!=v]
    out=tuple(sum(1<<j for j,b in enumerate(labels) if rows[a]>>b&1) for a in labels)
    return out,sum(1<<j for j,a in enumerate(labels) if mask>>a&1)


@lru_cache(None)
def state(rows,mask):
    if not mask:return None
    h,u,w,_,_,_=pt.ib.tc.path_tables(rows)
    selected=[v for v in range(len(rows)) if mask>>v&1]
    if h-sum(w[a][b] for a,b in combinations(selected,2))!=1:return None
    q=sum(u[v] for v in selected)
    c=q
    while c%2==0:c//=2
    return h,q,c


@lru_cache(None)
def reachable(rows,mask):
    if state(rows,mask) is None:return False
    if len(rows)==1:return True
    return any(reachable(*delete(rows,mask,v)) for v in range(len(rows)))


def main():
    bank=[[]];levels=[];atoms=[];edges=[]
    for n in range(1,8):
        codes={pt.ib.cc.canonical(pt.ib.cc.child(r,m)) for r in bank for m in range(1<<(n-1))}
        bank=[pt.ib.cc.decode(n,c) for c in sorted(codes)]
        seeds,irreducible,accessible=0,0,0
        for code,rr in zip(sorted(codes),bank):
            rows=tuple(rr)
            for mask in range(1,1<<n):
                profile=state(rows,mask)
                if profile is None:continue
                seeds+=1
                predecessors=[]
                for v in range(n):
                    child=delete(rows,mask,v)
                    smaller=state(*child)
                    if smaller:
                        predecessors.append(dict(deleted=v,selected=bool(mask>>v&1),H=smaller[0],U=smaller[1],c=smaller[2]))
                ok=reachable(rows,mask)
                accessible+=ok
                if not predecessors:
                    irreducible+=1
                    atoms.append(dict(n=n,code=code,rows=rr,selected=mask,H=profile[0],U=profile[1],c=profile[2]))
                edges.append(dict(n=n,code=code,selected=mask,H=profile[0],U=profile[1],c=profile[2],
                                  reachable_from_singleton=ok,predecessors=predecessors))
        levels.append(dict(n=n,seeds=seeds,deletion_atoms=irreducible,reachable_from_singleton=accessible))
    report=dict(levels=levels,atoms=atoms,states=edges)
    (HERE/'birth_atoms.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'levels':levels,'first_nontrivial_atoms':[a for a in atoms if a['n']>1][:5]},indent=2))


if __name__=='__main__':main()
