"""Independent edge-constraint certificates for Q24 and literal Q25 repair."""
from pathlib import Path
from math import isqrt
from collections import Counter
import json,gzip,sys,importlib.util
HERE=Path(__file__).resolve().parent

def need(p,m):
    if not p:raise RuntimeError(m)

def graph(n):
    edges=[(a,b) for a in range(1,n+1) for b in range(a+1,n+1) if isqrt(a+b)**2==a+b]
    inc=[0]*(n+1)
    for i,(a,b) in enumerate(edges):inc[a]|=1<<i;inc[b]|=1<<i
    return edges,inc

def solve(n,ends):
    edges,inc=graph(n);full=(1<<len(edges))-1
    target=[0]+[1 if v in ends else 2 for v in range(1,n+1)]
    reasons=Counter();leaves=[]
    def visit(yes,no,choices):
        while True:
            old=(yes,no)
            for v in range(1,n+1):
                have=(yes&inc[v]).bit_count();free=inc[v]&~(yes|no)
                if have>target[v] or have+free.bit_count()<target[v]:
                    reasons['degree']+=1;leaves.append({'choices':choices,'contradiction':'degree','vertex':v});return []
                if have==target[v]:no|=free
                elif have+free.bit_count()==target[v]:yes|=free
            # No Hamilton path contains a selected cycle, even a spanning one.
            parent=list(range(n+1))
            def root(v):
                while parent[v]!=v:v=parent[v]
                return v
            for i,(a,b) in enumerate(edges):
                if yes>>i&1:
                    a,b=root(a),root(b)
                    if a==b:
                        reasons['cycle']+=1;leaves.append({'choices':choices,'contradiction':'cycle'});return []
                    parent[a]=b
            if (yes,no)==old:break
        todo=full&~(yes|no)
        if not todo:
            need(yes.bit_count()==n-1,'degree completion is a spanning tree')
            return [yes]
        # Branch at the most constrained remaining vertex.
        v=min((v for v in range(1,n+1) if inc[v]&todo),key=lambda v:(inc[v]&todo).bit_count())
        bit=inc[v]&todo;bit&=-bit;i=bit.bit_length()-1
        return visit(yes|bit,no,choices+[[i,1]])+visit(yes,no|bit,choices+[[i,0]])
    solutions=visit(0,0,[])
    return {'endpoints':ends,'solutions':solutions,'leaves':leaves,'reasons':dict(reasons)}

def run():
    edges,inc=graph(24);leaves=[v for v in range(1,25) if inc[v].bit_count()==1]
    need(leaves==[18],'18 is the compulsory endpoint')
    records=[solve(24,[18,v]) for v in range(1,25) if v!=18]
    need(all(not r['solutions'] for r in records),'every possible second endpoint is excluded')
    path=HERE.parent/'joint_incidence_clocked_frames_20261102/path_kernel.py'
    spec=importlib.util.spec_from_file_location('old_square_paths',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    endpoint_rows=[];repairs=[]
    for n in [23,24,25]:
        _,paths=m.graph_paths(n,True)
        pairs=Counter(tuple(sorted((p[0]+1,p[-1]+1))) for p in paths)
        endpoint_rows.append({'n':n,'directed_paths':len(paths),'endpoint_pairs':[[*p,c//2] for p,c in sorted(pairs.items())]})
        if n==25:
            for p in paths:
                if p>p[::-1]:continue
                i=p.index(24);need(0<i<len(p)-1,'25 is always internal')
                need({p[i-1]+1,p[i+1]+1}=={11,24},'the two actual entry ports')
                parts=[[v+1 for v in p[:i]],[v+1 for v in p[i+1:]]]
                repairs.append(parts)
    need([r['directed_paths'] for r in endpoint_rows]==[6,0,20],'independent first-birth path counts')
    out={'scope':'All23 possible second endpoints of Q24, with18 forced by degree1. Every binary search leaf has a degree or cycle contradiction. Independent literal path enumeration checks Q23..25; deleting25 exhibits all10 unoriented repairs through ports11,24.','edges':edges,'endpoint_exclusions':records,'path_counts':endpoint_rows,'Q25_repairs':repairs}
    tag='' if __debug__ else '.optimized';(HERE/('square_endpoint_certificate'+tag+'.json.gz')).write_bytes(gzip.compress((json.dumps(out,separators=(',',':'))+'\n').encode(),mtime=0))
    print('Q24 contradiction leaves',sum(len(r['leaves']) for r in records),'max branch depth',max(len(l['choices']) for r in records for l in r['leaves']),'endpoints',endpoint_rows,flush=True)
if __name__=='__main__':run()
