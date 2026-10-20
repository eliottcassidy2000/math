"""Exact ABSTRACT MACROGRAPH control; edges are proved black paths, not edges
of the induced tournament metagraph. Every golden SCC word through order18
is generated without sampling, then converse-merged by literal word reversal.
"""
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
WEIGHTS={'a':1,'b':3,'c':4}


def need(x,msg):
    if not x:
        raise RuntimeError(msg)


@lru_cache(None)
def words(n):
    if n==0:
        return ('',)
    return tuple(letter+tail for letter,size in WEIGHTS.items() if n>=size
                 for tail in words(n-size))


def active(w):
    return 'c' in w or w!=w[::-1]


def canonical(w):
    return min(w,w[::-1])


def mass(w):
    return sum(WEIGHTS[x] for x in w),w.count('b')+w.count('c')


def main():
    rows=[]
    boundary={'interior':0,'left_end':0,'right_end':0,'whole_word':0}
    special={}
    all_node_hash=sha256()
    all_edge_hash=sha256()
    for n in range(19):
        ordinary=words(n)
        by_f={}
        for w in ordinary:
            if active(w):
                by_f.setdefault(mass(w)[1],set()).add(canonical(w))
        for f,node_set in sorted(by_f.items()):
            nodes=sorted(node_set)
            parent={w:w for w in nodes}
            def find(w):
                while parent[w]!=w:
                    parent[w]=parent[parent[w]]
                    w=parent[w]
                return w
            edges=set()
            for w in nodes:
                all_node_hash.update(f'{n}:{f}:{w}\n'.encode())
                candidates=[]
                for i,x in enumerate(w):
                    if x=='c':
                        for replacement in ('ab','ba'):
                            candidates.append((w[:i]+replacement+w[i+1:],i,1))
                for i in range(len(w)-1):
                    if w[i:i+2] in ('ab','ba'):
                        candidates.append((w[:i]+w[i:i+2][::-1]+w[i+2:],i,2))
                for v,i,width in candidates:
                    if not active(v):
                        continue
                    need(mass(v)==(n,f),'macro changed endpoint inventory')
                    v=canonical(v)
                    need(v in node_set,'macro target is missing from universe')
                    key=tuple(sorted((w,v)))
                    if key[0]!=key[1]:
                        edges.add(key)
                        uroot,vroot=find(w),find(v)
                        parent[uroot]=vroot
                    left,right=i>0,i+width<len(w)
                    location='interior' if left and right else 'right_end' if left else 'left_end' if right else 'whole_word'
                    boundary[location]+=1
            components=len({find(w) for w in nodes})
            need(components==1,'fixed-inventory macrograph disconnected')
            target=canonical('a'*(n-3*f)+'b'*f)
            need(target in node_set,'normal target absent or pure')
            for u,v in sorted(edges):
                all_edge_hash.update(f'{n}:{f}:{u}:{v}\n'.encode())
            row=dict(n=n,f=f,ordinary_golden_words=len(ordinary),
                     active_merged_nodes=len(nodes),simple_macro_edges=len(edges),
                     components=components,normal_target=target)
            rows.append(row)
            if n in (4,5):
                special[str(n)]=dict(nodes=nodes,edges=sorted(edges),target=target)
    need(special['4']['nodes']==['ab','c'],'n4 boundary universe changed')
    need(special['5']['nodes']==['aab','ac'],'n5 boundary universe changed')
    need(all(boundary.values()),'an endpoint category was not tested')
    report=dict(scope='complete abstract macrograph on active golden words through vertex order18, modulo reversal',
                alphabet=dict(a='singleton',b='C3',c='R4'),
                total_ordinary_words=sum(len(words(n)) for n in range(19)),
                total_active_merged_nodes=sum(r['active_merged_nodes'] for r in rows),
                fixed_inventory_strata=len(rows),strata=rows,
                boundary_macro_presentations=boundary,n4_n5_controls=special,
                node_universe_sha256=all_node_hash.hexdigest(),
                simple_edge_universe_sha256=all_edge_hash.hexdigest(),
                source_sha256_lf=sha256(Path(__file__).read_text().encode()).hexdigest())
    output=json.dumps(report,indent=2)+'\n'
    (HERE/'arithmetic_wordgraph.json').write_text(output,encoding='utf-8')
    print(output,end='')


if __name__=='__main__':
    main()
