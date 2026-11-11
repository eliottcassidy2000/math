"""Independent physical graph counts and hostile controls for the new packet.

Uses the earlier forward Hamilton DP/ear-curvature verifier, not this
packet's quotient-word endpoint computation. Same-agent algorithmic audit.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
from math import comb
import hashlib,json,sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'endpoint_coverage_compiler_20261111'))
import verify as independent
import construct as old
sys.path.insert(0,str(HERE))
import realize
import cones
import birth_atoms as ba

require=independent.require


@lru_cache(None)
def marked_state(rows,mask):
    if mask==0:return None
    selected=[i for i in range(len(rows)) if mask>>i&1]
    h,w,u=independent.endpoint_state(list(rows),selected)
    if h-w!=1:return None
    c=u
    while c%2==0:c//=2
    return h,u,c


@lru_cache(None)
def reachable(rows,mask):
    if marked_state(rows,mask) is None:return False
    return len(rows)==1 or any(reachable(*ba.delete(rows,mask,v)) for v in range(len(rows)))


def main():
    completion=json.loads((HERE/'head_completion.json').read_text())
    source=HERE.parents[2]/completion['inventory_source']
    require(hashlib.sha256(source.read_bytes().replace(b'\r\n',b'\n')).hexdigest()==completion['inventory_sha256_lf'],'inventory pin')
    kernels=prefix=0
    for name,head in completion['completed_heads'].items():
        w=head['witness'];values=[]
        # At most n visits to one inflated port in a single spanning path.
        for a in range(1,len(w['rows'])+1):
            _,rows=old.ib.inflate(w['rows'],w['port'],a)
            values.append(independent.count(rows)[0])
        coefficients=independent.invert(values)
        require(coefficients==head['kernel']+[0]*(len(coefficients)-len(head['kernel'])),'independent full kernel '+name)
        kernels+=1
        for item in head['prefix']:
            require(independent.count(item['rows'])[0]==item['H'],'literal finite prefix')
            prefix+=1
    constructions=backtracks=0
    for family in ['three','twentyone','mixed11']:
        for r in range(1,5):
            for s in range(1,7):
                if family=='mixed11' and s==2:continue
                item=realize.realize(family,r,s)
                h,_=independent.count(item['rows'])
                require(h==item['expected_H'],'actual completed-line realization')
                if len(item['rows'])<=10:
                    require(independent.recursive_count(item['rows'])==h,'recursive path count')
                    backtracks+=1
                constructions+=1
    print('Independent Hamilton kernels and constructors passed.',flush=True)
    simultaneous=0
    for item in json.loads((HERE/'cones.json').read_text())['records']:
        sizes=[2]*len(item['ports'])
        selected=[v for v in range(item['n']) if item['selected']>>v&1]
        rows,es,_=cones.expand(item['rows'],selected,item['ports'],sizes)
        h,w,u=independent.endpoint_state(rows,es)
        polynomial={tuple(term[:-1]):term[-1] for term in item['slope_terms']}
        require(h-w==1 and cones.eval_kernel(polynomial,sizes)==u,'independent simultaneous seed')
        simultaneous+=1
    transport=json.loads((HERE/'port_transport.json').read_text())['decisions']
    for item in transport:
        if item['stable_sizes']=='all':
            require(item['defect_tensor']==[[1,1,1]],'joint stability tensor')
        else:
            certificate=item['certificate'];poly=dict(certificate['exponential_coefficients'])
            top=max(poly);a=item['no_roots_from']
            require(abs(poly[top])*top**a>sum(abs(c)*b**a for b,c in poly.items() if b!=top),'dominance inequality')
            def column(j,x):
                return sum(d*sum((-1)**(i-b)*comb(i,b)*b**x for b in range(1,i+1))
                           for i,k,d in item['defect_tensor'] if k==j)
            bound=max(k for _,k,_ in item['defect_tensor'])
            roots=[x for x in range(1,a) if all(column(j,x)==0 for j in range(2,bound+1))]
            require(roots==item['stable_sizes'],'all pre-dominance roots')
    print('Simultaneous seeds and all-size transport certificates passed.',flush=True)
    births=json.loads((HERE/'birth_atoms.json').read_text());states=0
    for item in births['states']:
        rows=tuple(old.ib.cc.decode(item['n'],item['code']));mask=item['selected']
        require(marked_state(rows,mask)==(item['H'],item['U'],item['c']),'independent seed state')
        predecessors=[]
        for v in range(len(rows)):
            small=marked_state(*ba.delete(rows,mask,v))
            if small:predecessors.append(dict(deleted=v,selected=bool(mask>>v&1),H=small[0],U=small[1],c=small[2]))
        require(predecessors==item['predecessors'],'all immediate deletion choices')
        require(reachable(rows,mask)==item['reachable_from_singleton'],'full ancestral deletion DAG')
        states+=1
    # A Hamilton doubling port at the first parameter does not transport.
    hostile_rows=[8,33,3,6,15,29]
    actual=[]
    for a,b in [(1,1),(1,2),(2,1),(2,2)]:
        rows,_,_=cones.expand(hostile_rows,[],[2,0],[a,b])
        actual.append(independent.count(rows)[0])
    require(actual==[23,45,59,129] and actual[3]!=2*actual[2]-1,'doubling transport hostile')
    for family,r,s in [('three',0,1),('twentyone',1,0),('mixed11',1,2)]:
        try:realize.realize(family,r,s)
        except ValueError:pass
        else:raise RuntimeError('invalid or unclosed constructor input accepted')
    report=dict(independent_full_Hamilton_kernels=kernels,finite_prefixes=prefix,
                literal_line_constructions=constructions,recursive_backtracks=backtracks,
                simultaneous_seeds=simultaneous,all_size_transport_certificates=len(transport),
                independently_checked_marked_states=states,doubling_hostile_counts=actual,
                separate_person_audit=False)
    (HERE/'verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
