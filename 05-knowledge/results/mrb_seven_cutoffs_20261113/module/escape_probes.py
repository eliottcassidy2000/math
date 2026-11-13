"""Complete small escape menus past the fixed uniform-module obstruction."""
import json,itertools
from pathlib import Path
import module_certificates as z
HERE=Path(__file__).resolve().parent

def run():
    bank=json.loads((HERE/'verification.json').read_text())['search']['candidates']
    sources=[r for r in bank if r['target']==881]
    z.need(len(sources)==2,'both classified881 sources')
    selective=[];extensions=[];predecessor=[];tests=0;ready=0
    for q in sources:
        n=q['order']-1;core=z.f.a.unrank(n,q['rank'])
        for k in range(1,5):
            parent,ix=z.inflate(core,0,k)
            for t in sorted(set([0,k-1])):
                p=ix[0]+t;child,new=z.f.a.birth(parent,p)
                predecessor.append(dict(words=q['words'],k=k,mark=t,H=z.count(child),defect=z.count(parent,p)))
        for v in range(1,n):
            parent,ix=z.inflate(core,v,2)
            for bits in itertools.product((0,1),repeat=2):
                rows=list(parent);p=ix[0]
                for t,bit in enumerate(bits):
                    w=ix[v]+t
                    if bit:rows[p]|=1<<w;rows[w]&=~(1<<p)
                    else:rows[w]|=1<<p;rows[p]&=~(1<<w)
                child,new=z.f.a.birth(tuple(rows),p)
                selective.append(dict(words=q['words'],module=v,port_bits=bits,H=z.count(child),defect=z.count(tuple(rows),p)))
        for mask in range(1<<n):
            parent=tuple(row|((1<<n) if mask>>i&1 else 0) for i,row in enumerate(core))+(sum(1<<i for i in range(n) if not mask>>i&1),)
            delta=z.count(parent,0);tests+=1
            if delta!=1:continue
            ready+=1;child,new=z.f.a.birth(parent,0);h=z.count(child)
            if h==1691:extensions.append(dict(words=q['words'],mask=mask,H=h,defect=delta))
    z.need(len(selective)==68 and not any(r['H']==1691 for r in selective),'complete one-selective-predecessor menu')
    z.need(tests==1536 and not extensions,'complete one-new-vertex menu')
    z.need(any(r['defect']==1 for r in selective),'positive ready controls')
    z.need(any(r['defect']>1 for r in predecessor),'hostile rootclone controls')
    z.need(predecessor[3]['H']!=1691,'target obstruction control')
    result=dict(status='FINITE-EXACT',selective_scope='Both complete881 strong ready sources; every nonport vertex cloned once; all four comparison words against sole predecessor; other comparisons uniform; transitive clone pair',selective_tests=len(selective),selective_ready=sum(r['defect']==1 for r in selective),selective_rows=selective,extension_scope='Every one-new-vertex induced extension of both classified881 strong ready sources; same physical indegree-one root and predecessor; every orientation against all other old vertices; filter exact defect1',extension_tests=tests,extension_ready=ready,extension_target_hits=extensions,predecessor_clone_scope='TT_k substitutes the old predecessor in the parent; marked first or last clone becomes the new root sole predecessor; k1..4',predecessor_clone_rows=predecessor)
    name='escapes.json' if __debug__ else 'escapes.optimized.json'
    (HERE/name).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS',len(selective),'selective comparisons;',tests,'induced extensions;',ready,'ready extensions; no1691')

if __name__=='__main__':run()
